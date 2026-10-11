"""ClinVar data ingestion module.

Ingests ClinVar variant records from tab-delimited exports while strictly
preserving source files unchanged, calculating SHA-256 hashes, distinguishing
aggregate VCV records, preserving conflicting clinical assertions verbatim,
and validating reference residues against canonical UniProt P13569.
"""

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from sqlalchemy.orm import Session

from backend.app.db.models import ClinVarVariant, DataSource
from backend.app.ingestion.models import ImportSummary
from backend.app.scientific.normalization import normalize_protein_change
from backend.app.scientific.reference import validate_wt_residue


def compute_file_sha256(file_path: Path) -> str:
    """Compute SHA-256 hash of a file without modifying it."""
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def ingest_clinvar(
    file_path: Path,
    db_session: Optional[Session] = None,
    dry_run: bool = False,
) -> ImportSummary:
    """Ingest ClinVar dataset with strict provenance, normalization, and validation.
    
    Args:
        file_path: Path to ClinVar TSV file (e.g. data/raw/clinvar_download_*.tsv).
        db_session: Optional SQLAlchemy Session for database storage.
        dry_run: If True, executes parsing and validation without writing to DB.
        
    Returns:
        ImportSummary with row counts, validation metrics, and exclusion reasons.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"ClinVar source file not found at: {path}")

    file_sha256 = compute_file_sha256(path)

    summary = ImportSummary(
        source_name="ClinVar",
        source_file=path.name,
        file_sha256=file_sha256,
        source_version=path.stem.replace("clinvar_download_", ""),
    )

    seen_variation_ids = set()
    db_records: List[ClinVarVariant] = []

    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            summary.rows_read += 1

            variation_id_raw = row.get("VariationID", "").strip()
            name = row.get("Name", "").strip()
            vcv_accession = row.get("Accession", "").strip()

            if not variation_id_raw:
                summary.skipped_rows += 1
                summary.add_exclusion_reason("MISSING_VARIATION_ID")
                continue

            try:
                variation_id = int(variation_id_raw)
            except ValueError:
                summary.skipped_rows += 1
                summary.add_exclusion_reason("INVALID_VARIATION_ID")
                continue

            if variation_id in seen_variation_ids:
                summary.skipped_rows += 1
                summary.duplicate_identifiers += 1
                summary.add_exclusion_reason("DUPLICATE_VARIATION_ID")
                continue
            seen_variation_ids.add(variation_id)

            protein_change_raw = row.get("Protein change", "").strip()
            norm = normalize_protein_change(protein_change_raw, fallback_name=name)

            is_valid_wt = None
            val_flag = "UNVALIDATED"
            val_message = None

            if norm.success and norm.pos and norm.ref_aa:
                val_res = validate_wt_residue(norm.pos, norm.ref_aa)
                is_valid_wt = val_res.is_valid
                val_flag = val_res.flag
                val_message = val_res.message
                if val_res.flag == "MATCH":
                    summary.wt_matches += 1
                elif val_res.flag == "WT_MISMATCH":
                    summary.wt_mismatches += 1
                elif val_res.flag == "OUT_OF_BOUNDS":
                    summary.wt_out_of_bounds += 1
            elif norm.success:
                summary.unvalidated_rows += 1
                val_flag = "NO_POSITION_OR_REF"
            else:
                summary.normalization_failures += 1
                summary.add_exclusion_reason(norm.failure_reason or "NORMALIZATION_FAILURE")
                val_flag = "NORMALIZATION_FAILED"

            summary.valid_rows += 1

            if not dry_run and db_session:
                record = ClinVarVariant(
                    variation_id=variation_id,
                    vcv_accession=vcv_accession or None,
                    name=name,
                    protein_change_raw=protein_change_raw or None,
                    canonical_protein=norm.canonical_protein,
                    uniprot_pos=norm.pos,
                    ref_aa=norm.ref_aa,
                    alt_aa=norm.alt_aa,
                    variant_type=row.get("Variant type", "unknown"),
                    molecular_consequence=row.get("Molecular consequence") or None,
                    is_valid_wt=is_valid_wt,
                    validation_flag=val_flag,
                    validation_message=val_message,
                    germline_classification=row.get("Germline classification") or None,
                    review_status=row.get("Germline review status") or None,
                    date_last_evaluated=row.get("Germline date last evaluated") or None,
                    conditions=row.get("Condition(s)") or None,
                    grch38_location=row.get("GRCh38Location") or None,
                    grch37_location=row.get("GRCh37Location") or None,
                    dbsnp_id=row.get("dbSNP ID") or None,
                    allele_id=row.get("AlleleID(s)") or None,
                    raw_data_json=json.dumps(row),
                )
                db_records.append(record)

    if not dry_run and db_session:
        try:
            existing_source = (
                db_session.query(DataSource)
                .filter_by(source_name="ClinVar", file_sha256=file_sha256)
                .first()
            )
            if existing_source:
                data_source = existing_source
                data_source.source_version = summary.source_version
                data_source.record_count = len(db_records)
                data_source.summary_json = json.dumps(summary.to_dict())
                data_source.imported_at = datetime.now(timezone.utc)
            else:
                data_source = DataSource(
                    source_name="ClinVar",
                    filename=path.name,
                    file_sha256=file_sha256,
                    source_version=summary.source_version,
                    record_count=len(db_records),
                    summary_json=json.dumps(summary.to_dict()),
                )
                db_session.add(data_source)
            db_session.flush()

            existing_records = {
                r.variation_id: r
                for r in db_session.query(ClinVarVariant).all()
            }

            for rec in db_records:
                rec.source_id = data_source.id
                if rec.variation_id in existing_records:
                    target = existing_records[rec.variation_id]
                    target.vcv_accession = rec.vcv_accession
                    target.name = rec.name
                    target.protein_change_raw = rec.protein_change_raw
                    target.canonical_protein = rec.canonical_protein
                    target.uniprot_pos = rec.uniprot_pos
                    target.ref_aa = rec.ref_aa
                    target.alt_aa = rec.alt_aa
                    target.variant_type = rec.variant_type
                    target.molecular_consequence = rec.molecular_consequence
                    target.is_valid_wt = rec.is_valid_wt
                    target.validation_flag = rec.validation_flag
                    target.validation_message = rec.validation_message
                    target.germline_classification = rec.germline_classification
                    target.review_status = rec.review_status
                    target.date_last_evaluated = rec.date_last_evaluated
                    target.conditions = rec.conditions
                    target.grch38_location = rec.grch38_location
                    target.grch37_location = rec.grch37_location
                    target.dbsnp_id = rec.dbsnp_id
                    target.allele_id = rec.allele_id
                    target.source_id = data_source.id
                    target.raw_data_json = rec.raw_data_json
                else:
                    db_session.add(rec)
            db_session.commit()
        except Exception:
            db_session.rollback()
            raise

    return summary
