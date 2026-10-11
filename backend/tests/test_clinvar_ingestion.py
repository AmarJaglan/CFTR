"""Automated tests for ClinVar ingestion using synthetic test fixtures."""

from pathlib import Path
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.app.db.models import Base, ClinVarVariant, DataSource
from backend.app.ingestion.clinvar import compute_file_sha256, ingest_clinvar


@pytest.fixture
def in_memory_session():
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(bind=engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        yield session
    finally:
        session.close()


def create_synthetic_clinvar_tsv(dest_path: Path, rows: list) -> None:
    """Generate a TSV file with synthetic ClinVar records."""
    headers = [
        "Name",
        "Gene(s)",
        "Protein change",
        "Condition(s)",
        "Accession",
        "GRCh37Chromosome",
        "GRCh37Location",
        "GRCh38Chromosome",
        "GRCh38Location",
        "VariationID",
        "AlleleID(s)",
        "dbSNP ID",
        "Canonical SPDI",
        "Variant type",
        "Molecular consequence",
        "Germline classification",
        "Germline date last evaluated",
        "Germline review status",
        "Somatic clinical impact",
        "Somatic clinical impact date last evaluated",
        "Somatic clinical impact review status",
        "Oncogenicity classification",
        "Oncogenicity date last evaluated",
        "Oncogenicity review status",
    ]
    with open(dest_path, "w", encoding="utf-8") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows:
            f.write("\t".join(str(val) for val in r) + "\n")


def test_synthetic_clinvar_ingestion(tmp_path, in_memory_session):
    synthetic_file = tmp_path / "synthetic_clinvar.tsv"

    # Synthetic rows:
    # 1. Matching missense: M1V (Met1 matches canonical WT)
    # 2. Matching with conflicting classifications: G551D (must preserve "Conflicting classifications of pathogenicity")
    # 3. Mismatch: P508del (reported P at 508, canonical is F) -> must flag WT_MISMATCH
    # 4. Duplicate VariationID: 990001 (should be skipped)
    # 5. Missing VariationID: empty (should be skipped)
    # 6. Non-coding variant: empty protein change (should record unvalidated)
    test_rows = [
        [
            "NM_000492.4(CFTR):c.1A>G (p.Met1Val)", "CFTR", "M1V", "Cystic fibrosis",
            "VCV000990001", "7", "117120149", "7", "117480095", "990001", "880001", "rs99901",
            "NC_000007.14:117480094:A:G", "single nucleotide variant", "missense variant",
            "Pathogenic", "Mar 17, 2020", "reviewed by expert panel", "", "", "", "", "", "",
        ],
        [
            "NM_000492.4(CFTR):c.1652G>A (p.Gly551Asp)", "CFTR", "G551D", "Cystic fibrosis",
            "VCV000990002", "7", "117120150", "7", "117480096", "990002", "880002", "rs99902",
            "NC_000007.14:117480095:G:A", "single nucleotide variant", "missense variant",
            "Conflicting classifications of pathogenicity", "Aug 7, 2024",
            "criteria provided, conflicting classifications", "", "", "", "", "", "",
        ],
        [
            "NM_000492.4(CFTR):c.1522_1524del (p.Pro508del)", "CFTR", "P508del", "Cystic fibrosis",
            "VCV000990003", "7", "117120151", "7", "117480097", "990003", "880003", "rs99903",
            "NC_000007.14:117480096:T:C", "deletion", "inframe deletion",
            "Uncertain significance", "Jul 24, 2023", "criteria provided, single submitter",
            "", "", "", "", "", "",
        ],
        [
            "NM_000492.4(CFTR):c.1A>G (p.Met1Val)", "CFTR", "M1V", "Cystic fibrosis",
            "VCV000990001", "7", "117120149", "7", "117480095", "990001", "880001", "rs99901",
            "NC_000007.14:117480094:A:G", "single nucleotide variant", "missense variant",
            "Pathogenic", "Mar 17, 2020", "reviewed by expert panel", "", "", "", "", "", "",
        ],
        [
            "NM_000492.4(CFTR):c.50G>A (p.Arg17His)", "CFTR", "R17H", "Cystic fibrosis",
            "VCV000990005", "7", "117120152", "7", "117480098", "", "880005", "rs99905",
            "", "single nucleotide variant", "missense variant",
            "Likely pathogenic", "Jan 1, 2021", "criteria provided, single submitter",
            "", "", "", "", "", "",
        ],
        [
            "NM_000492.4(CFTR):c.-887C>T", "CFTR", "", "Cystic fibrosis",
            "VCV000990006", "7", "117119262", "7", "117479208", "990006", "880006", "rs99906",
            "", "single nucleotide variant", "non-coding transcript variant",
            "Benign", "Jan 1, 2022", "criteria provided, single submitter",
            "", "", "", "", "", "",
        ],
    ]

    create_synthetic_clinvar_tsv(synthetic_file, test_rows)
    initial_sha = compute_file_sha256(synthetic_file)

    summary = ingest_clinvar(synthetic_file, db_session=in_memory_session, dry_run=False)

    # 1. Verify file immutability: SHA-256 did not change
    assert compute_file_sha256(synthetic_file) == initial_sha
    assert summary.file_sha256 == initial_sha

    # 2. Verify row accounting
    assert summary.rows_read == 6
    assert summary.valid_rows == 4
    assert summary.skipped_rows == 2  # 1 missing VariationID, 1 duplicate
    assert summary.duplicate_identifiers == 1
    assert "MISSING_VARIATION_ID" in summary.exclusion_reasons
    assert "DUPLICATE_VARIATION_ID" in summary.exclusion_reasons

    # 3. Verify scientific validation
    assert summary.wt_matches == 2  # M1V and G551D
    assert summary.wt_mismatches == 1  # P508del (reported P, canonical F)
    assert summary.wt_out_of_bounds == 0

    # 4. Verify database records and verbatim clinical assertions
    records = in_memory_session.query(ClinVarVariant).all()
    assert len(records) == 4

    rec_m1v = in_memory_session.query(ClinVarVariant).filter_by(variation_id=990001).first()
    assert rec_m1v is not None
    assert rec_m1v.canonical_protein == "M1V"
    assert rec_m1v.is_valid_wt is True
    assert rec_m1v.validation_flag == "MATCH"
    assert rec_m1v.germline_classification == "Pathogenic"

    # Conflicting assertions preserved verbatim without artificial consensus
    rec_conflict = in_memory_session.query(ClinVarVariant).filter_by(variation_id=990002).first()
    assert rec_conflict is not None
    assert rec_conflict.canonical_protein == "G551D"
    assert rec_conflict.is_valid_wt is True
    assert rec_conflict.germline_classification == "Conflicting classifications of pathogenicity"
    assert rec_conflict.review_status == "criteria provided, conflicting classifications"

    # WT mismatch flagged, never silently corrected
    rec_mismatch = in_memory_session.query(ClinVarVariant).filter_by(variation_id=990003).first()
    assert rec_mismatch is not None
    assert rec_mismatch.is_valid_wt is False
    assert rec_mismatch.validation_flag == "WT_MISMATCH"
    assert "Reference mismatch at position 508" in rec_mismatch.validation_message

    # DataSource provenance
    data_source = in_memory_session.query(DataSource).first()
    assert data_source is not None
    assert data_source.source_name == "ClinVar"
    assert data_source.record_count == 4
    assert data_source.file_sha256 == initial_sha
