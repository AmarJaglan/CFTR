"""CFTR2 data ingestion module.

Ingests CFTR2 variant records from Excel exports while strictly preserving
source files unchanged, calculating SHA-256 hashes, validating reference
residues against canonical UniProt P13569, and preserving original clinical
interpretations verbatim.
"""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import xml.etree.ElementTree as ET
import zipfile

from sqlalchemy.orm import Session

from backend.app.db.models import CFTR2Variant, DataSource
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


def _read_xlsx_rows(file_path: Path) -> Tuple[List[Dict[str, str]], Optional[str]]:
    """Parse rows from the primary CFTR2 worksheet using standard library zipfile/xml.
    
    Returns:
        Tuple of (data_rows, detected_version_string).
    """
    with zipfile.ZipFile(file_path, "r") as z:
        # Load shared strings table
        ss_root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        ns = {"ns": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        shared_strings = []
        for si in ss_root.findall(".//ns:si", ns):
            t_elem = si.find(".//ns:t", ns)
            if t_elem is not None and t_elem.text is not None:
                shared_strings.append(t_elem.text)
            else:
                texts = [t.text for t in si.findall(".//ns:t", ns) if t.text]
                shared_strings.append("".join(texts))

        # Read sheet1
        sheet_root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        header_map: Dict[str, str] = {}
        header_row_idx: Optional[int] = None
        data_rows: List[Dict[str, str]] = []
        detected_version: Optional[str] = None

        for r in sheet_root.findall(".//ns:row", ns):
            r_idx = int(r.attrib.get("r", "0"))
            row_dict: Dict[str, str] = {}
            for c in r.findall("./ns:c", ns):
                coord = c.attrib.get("r", "")
                col = "".join(filter(str.isalpha, coord))
                t_attr = c.attrib.get("t")
                v = c.find("./ns:v", ns)
                val = v.text if v is not None and v.text is not None else ""
                if t_attr == "s" and val:
                    val = shared_strings[int(val)]
                row_dict[col] = val.strip()

            # Attempt to extract version from early metadata rows
            if r_idx < 10:
                for val in row_dict.values():
                    if "determination" in val.lower() or "version" in val.lower() or "date:" in val.lower():
                        if not detected_version:
                            detected_version = val

            # Locate header row containing 'Variant legacy name'
            if not header_map:
                for col, val in row_dict.items():
                    if "Variant legacy name" in val:
                        header_map = row_dict
                        header_row_idx = r_idx
                        # Extract version from header column if available
                        for h_val in row_dict.values():
                            if "current version" in h_val.lower():
                                detected_version = h_val.replace("Variant final determination", "").replace("(current version)", "").strip()
                        break
            elif header_row_idx is not None and r_idx > header_row_idx:
                if any(row_dict.values()):
                    # Map row by header names
                    named_row: Dict[str, str] = {}
                    for col, val in row_dict.items():
                        h_name = header_map.get(col, col)
                        named_row[h_name] = val
                    data_rows.append(named_row)

        return data_rows, detected_version


def ingest_cftr2(
    file_path: Path,
    db_session: Optional[Session] = None,
    dry_run: bool = False,
) -> ImportSummary:
    """Ingest CFTR2 dataset with strict provenance, normalization, and validation.
    
    Args:
        file_path: Path to CFTR2 Excel file (e.g. data/raw/CFTR2_30January2026.xlsx).
        db_session: Optional SQLAlchemy Session for database storage.
        dry_run: If True, executes parsing and validation without writing to DB.
        
    Returns:
        ImportSummary with row counts, validation metrics, and exclusion reasons.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"CFTR2 source file not found at: {path}")

    file_sha256 = compute_file_sha256(path)
    raw_rows, detected_version = _read_xlsx_rows(path)

    summary = ImportSummary(
        source_name="CFTR2",
        source_file=path.name,
        file_sha256=file_sha256,
        source_version=detected_version,
    )

    seen_identifiers = set()
    db_records: List[CFTR2Variant] = []

    for row in raw_rows:
        summary.rows_read += 1

        # Extract fields matching expected CFTR2 headers
        legacy_name = ""
        protein_name = ""
        cdna_name = ""
        alternative_names = ""
        alleles_reported_raw = ""
        allele_freq_raw = ""
        determination_prev = ""
        determination_curr = ""
        change_from_prev = ""

        for h, v in row.items():
            h_lower = h.lower()
            if "legacy name" in h_lower:
                legacy_name = v
            elif "protein name" in h_lower:
                protein_name = v
            elif "cdna name" in h_lower:
                cdna_name = v
            elif "alternative names" in h_lower:
                alternative_names = v
            elif "alleles reported" in h_lower:
                alleles_reported_raw = v
            elif "allele frequency" in h_lower:
                allele_freq_raw = v
            elif "previous version" in h_lower:
                determination_prev = v
            elif "current version" in h_lower or "final determination" in h_lower:
                determination_curr = v
            elif "change from previous" in h_lower:
                change_from_prev = v

        if not legacy_name:
            summary.skipped_rows += 1
            summary.add_exclusion_reason("MISSING_LEGACY_NAME")
            continue

        # Skip non-variant footer / disclaimer rows
        legacy_lower = legacy_name.lower()
        if (
            legacy_lower.startswith(("*represents", "permitted use", "please use", "©copyright", "copyright", "note:"))
            or (not cdna_name and not determination_curr and not protein_name)
        ):
            summary.skipped_rows += 1
            summary.add_exclusion_reason("FOOTER_DISCLAIMER_ROW")
            continue

        if legacy_name in seen_identifiers:
            summary.skipped_rows += 1
            summary.duplicate_identifiers += 1
            summary.add_exclusion_reason("DUPLICATE_LEGACY_NAME")
            continue
        seen_identifiers.add(legacy_name)

        # Parse numeric metrics safely
        alleles_reported = None
        if alleles_reported_raw:
            try:
                alleles_reported = int(alleles_reported_raw.replace(",", ""))
            except ValueError:
                pass

        allele_freq = None
        if allele_freq_raw:
            try:
                allele_freq = float(allele_freq_raw)
            except ValueError:
                pass

        # Normalize protein change
        norm = normalize_protein_change(protein_name, fallback_name=legacy_name)

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
            record = CFTR2Variant(
                legacy_name=legacy_name,
                protein_name=protein_name or None,
                cdna_name=cdna_name or None,
                canonical_protein=norm.canonical_protein,
                uniprot_pos=norm.pos,
                ref_aa=norm.ref_aa,
                alt_aa=norm.alt_aa,
                variant_type=norm.variant_type,
                is_valid_wt=is_valid_wt,
                validation_flag=val_flag,
                validation_message=val_message,
                alleles_reported=alleles_reported,
                allele_frequency=allele_freq,
                determination_current=determination_curr or None,
                determination_previous=determination_prev or None,
                change_from_previous=change_from_prev or None,
                alternative_names=alternative_names or None,
                raw_data_json=json.dumps(row),
            )
            db_records.append(record)

    if not dry_run and db_session:
        try:
            # Check or create DataSource provenance entry
            existing_source = (
                db_session.query(DataSource)
                .filter_by(source_name="CFTR2", file_sha256=file_sha256)
                .first()
            )
            if existing_source:
                data_source = existing_source
                data_source.source_version = detected_version
                data_source.record_count = len(db_records)
                data_source.summary_json = json.dumps(summary.to_dict())
                data_source.imported_at = datetime.now(timezone.utc)
            else:
                data_source = DataSource(
                    source_name="CFTR2",
                    filename=path.name,
                    file_sha256=file_sha256,
                    source_version=detected_version,
                    record_count=len(db_records),
                    summary_json=json.dumps(summary.to_dict()),
                )
                db_session.add(data_source)
            db_session.flush()

            # Upsert records by legacy_name to prevent duplicate rows on re-import
            existing_records = {
                r.legacy_name: r
                for r in db_session.query(CFTR2Variant).all()
            }

            for rec in db_records:
                rec.source_id = data_source.id
                if rec.legacy_name in existing_records:
                    target = existing_records[rec.legacy_name]
                    target.protein_name = rec.protein_name
                    target.cdna_name = rec.cdna_name
                    target.canonical_protein = rec.canonical_protein
                    target.uniprot_pos = rec.uniprot_pos
                    target.ref_aa = rec.ref_aa
                    target.alt_aa = rec.alt_aa
                    target.variant_type = rec.variant_type
                    target.is_valid_wt = rec.is_valid_wt
                    target.validation_flag = rec.validation_flag
                    target.validation_message = rec.validation_message
                    target.alleles_reported = rec.alleles_reported
                    target.allele_frequency = rec.allele_frequency
                    target.determination_current = rec.determination_current
                    target.determination_previous = rec.determination_previous
                    target.change_from_previous = rec.change_from_previous
                    target.alternative_names = rec.alternative_names
                    target.source_id = data_source.id
                    target.raw_data_json = rec.raw_data_json
                else:
                    db_session.add(rec)
            db_session.commit()
        except Exception:
            db_session.rollback()
            raise

    return summary
