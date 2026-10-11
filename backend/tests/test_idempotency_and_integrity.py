"""Automated tests for database integrity, import idempotency, and transaction rollback."""

from pathlib import Path
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.app.db.models import Base, CFTR2Variant, ClinVarVariant, DataSource
from backend.app.ingestion.cftr2 import ingest_cftr2
from backend.app.ingestion.clinvar import ingest_clinvar
from backend.tests.test_cftr2_ingestion import create_synthetic_cftr2_xlsx
from backend.tests.test_clinvar_ingestion import create_synthetic_clinvar_tsv


@pytest.fixture
def db_session():
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(bind=engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        yield session
    finally:
        session.close()


def test_cftr2_idempotent_reimport(tmp_path, db_session):
    """Verify that re-importing CFTR2 does not create duplicate rows."""
    xlsx_file = tmp_path / "cftr2_test.xlsx"
    rows = [
        ["SYN_1", "p.Gly551Asp", "c.1652G>A", "", "10", "0.01", "CF-causing", "CF-causing", "No"],
        ["SYN_2", "p.Phe508del", "c.1521_1523del", "", "100", "0.5", "CF-causing", "CF-causing", "No"],
    ]
    create_synthetic_cftr2_xlsx(xlsx_file, rows)

    # First import
    sum1 = ingest_cftr2(xlsx_file, db_session=db_session, dry_run=False)
    assert sum1.valid_rows == 2
    assert db_session.query(CFTR2Variant).count() == 2
    assert db_session.query(DataSource).count() == 1

    # Second import of identical file
    sum2 = ingest_cftr2(xlsx_file, db_session=db_session, dry_run=False)
    assert sum2.valid_rows == 2
    # Verify NO duplicate rows created
    assert db_session.query(CFTR2Variant).count() == 2
    assert db_session.query(DataSource).count() == 1


def test_clinvar_idempotent_reimport(tmp_path, db_session):
    """Verify that re-importing ClinVar does not create duplicate rows or throw IntegrityError."""
    tsv_file = tmp_path / "clinvar_test.tsv"
    rows = [
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
    ]
    create_synthetic_clinvar_tsv(tsv_file, rows)

    # First import
    sum1 = ingest_clinvar(tsv_file, db_session=db_session, dry_run=False)
    assert sum1.valid_rows == 2
    assert db_session.query(ClinVarVariant).count() == 2
    assert db_session.query(DataSource).count() == 1

    # Second import
    sum2 = ingest_clinvar(tsv_file, db_session=db_session, dry_run=False)
    assert sum2.valid_rows == 2
    # Verify row count is unchanged
    assert db_session.query(ClinVarVariant).count() == 2
    assert db_session.query(DataSource).count() == 1


def test_multiple_variants_with_same_protein_change_remain_distinct(tmp_path, db_session):
    """Verify that different source records with identical protein change are NOT collapsed."""
    tsv_file = tmp_path / "clinvar_same_protein.tsv"
    # Two distinct cDNA variants (c.3G>A and c.3G>T) causing the same protein change (M1I)
    rows = [
        [
            "NM_000492.4(CFTR):c.3G>A (p.Met1Ile)", "CFTR", "M1I", "Cystic fibrosis",
            "VCV000053869", "7", "117120151", "7", "117480097", "53869", "68289", "rs397508476",
            "NC_000007.14:117480095:G:A", "single nucleotide variant", "missense variant",
            "Pathogenic", "Jan 1, 2020", "criteria provided, single submitter", "", "", "", "", "", "",
        ],
        [
            "NM_000492.4(CFTR):c.3G>T (p.Met1Ile)", "CFTR", "M1I", "Cystic fibrosis",
            "VCV000053870", "7", "117120151", "7", "117480097", "53870", "68290", "rs397508476",
            "NC_000007.14:117480095:G:T", "single nucleotide variant", "missense variant",
            "Pathogenic", "Jan 1, 2020", "criteria provided, single submitter", "", "", "", "", "", "",
        ],
    ]
    create_synthetic_clinvar_tsv(tsv_file, rows)

    summary = ingest_clinvar(tsv_file, db_session=db_session, dry_run=False)
    assert summary.valid_rows == 2

    records = db_session.query(ClinVarVariant).all()
    assert len(records) == 2
    # Confirm both records share canonical_protein 'M1I' but have distinct VariationIDs
    assert records[0].canonical_protein == "M1I"
    assert records[1].canonical_protein == "M1I"
    assert {r.variation_id for r in records} == {53869, 53870}


def test_cftr2_footer_rows_skipped(tmp_path, db_session):
    """Verify that spreadsheet footer disclaimer rows are skipped and counted as exclusions."""
    xlsx_file = tmp_path / "cftr2_footers.xlsx"
    rows = [
        ["SYN_REAL", "p.Gly551Asp", "c.1652G>A", "", "10", "0.01", "CF-causing", "CF-causing", "No"],
        ["*Represents allele frequency note", "", "", "", "", "", "", "", ""],
        ["Permitted use disclaimer", "", "", "", "", "", "", "", ""],
        ["©Copyright 2011 US CF Foundation", "", "", "", "", "", "", "", ""],
    ]
    create_synthetic_cftr2_xlsx(xlsx_file, rows)

    summary = ingest_cftr2(xlsx_file, db_session=db_session, dry_run=False)
    assert summary.rows_read == 4
    assert summary.valid_rows == 1
    assert summary.skipped_rows == 3
    assert summary.exclusion_reasons.get("FOOTER_DISCLAIMER_ROW") == 3
    assert db_session.query(CFTR2Variant).count() == 1
