"""Automated tests for CFTR2 ingestion using synthetic test fixtures."""

import io
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.app.db.models import Base, CFTR2Variant, DataSource
from backend.app.ingestion.cftr2 import compute_file_sha256, ingest_cftr2


def create_synthetic_cftr2_xlsx(dest_path: Path, rows: list) -> None:
    """Generate a valid minimal XLSX file with synthetic CFTR2 records."""
    headers = [
        "Variant legacy name",
        "Variant protein name",
        "Variant cDNA name",
        "Variant alternative names",
        "# alleles reported in CFTR2",
        "Allele frequency in CFTR2*",
        "Variant final determination 25 September 2024 (previous version)",
        "Variant final determination 30 January 2026 (current version)",
        "Change from previous?",
    ]

    all_rows = [headers] + rows

    # Collect shared strings
    shared_strings = []
    string_map = {}
    for r in all_rows:
        for val in r:
            s_val = str(val)
            if s_val not in string_map:
                string_map[s_val] = len(shared_strings)
                shared_strings.append(s_val)

    # Build XML
    with zipfile.ZipFile(dest_path, "w") as z:
        z.writestr(
            "[Content_Types].xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>
</Types>""",
        )
        z.writestr(
            "xl/workbook.xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<sheets><sheet name="Sheet1" sheetId="1" r:id="rId1" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"/></sheets>
</workbook>""",
        )

        ss_entries = "".join(f"<si><t>{s}</t></si>" for s in shared_strings)
        z.writestr(
            "xl/sharedStrings.xml",
            f'<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="{len(shared_strings)}">{ss_entries}</sst>',
        )

        cols = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
        sheet_rows = []
        for r_idx, r_data in enumerate(all_rows, start=12):
            cells = []
            for c_idx, val in enumerate(r_data):
                col_letter = cols[c_idx]
                s_idx = string_map[str(val)]
                cells.append(f'<c r="{col_letter}{r_idx}" t="s"><v>{s_idx}</v></c>')
            sheet_rows.append(f'<row r="{r_idx}">{"".join(cells)}</row>')

        z.writestr(
            "xl/worksheets/sheet1.xml",
            f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<sheetData>{"".join(sheet_rows)}</sheetData>
</worksheet>""",
        )


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


def test_synthetic_cftr2_ingestion(tmp_path, in_memory_session):
    synthetic_file = tmp_path / "synthetic_cftr2.xlsx"

    # Synthetic rows designed specifically for test verification:
    # 1. Matching missense: G551D (Gly551 matches canonical WT)
    # 2. Matching deletion: F508del (Phe508 matches canonical WT)
    # 3. Mismatch: Ala508Val (reported Ala at 508, canonical is Phe) -> must flag WT_MISMATCH
    # 4. Duplicate identifier: SYN_G551D (should be skipped as duplicate)
    # 5. Missing identifier: empty legacy name (should be skipped)
    # 6. Non-coding / splice: p.? with legacy SYN_SPLICE -> preserved determination
    test_rows = [
        ["SYN_G551D", "p.Gly551Asp", "c.1652G>A", "", "3500", "0.018", "CF-causing", "CF-causing", "No"],
        ["SYN_F508del", "p.Phe508del", "c.1521_1523del", "", "130000", "0.65", "CF-causing", "CF-causing", "No"],
        ["SYN_MISMATCH", "p.Ala508Val", "c.1523C>T", "", "10", "0.0001", "Varying clinical consequence", "Varying clinical consequence", "No"],
        ["SYN_G551D", "p.Gly551Asp", "c.1652G>A", "", "3500", "0.018", "CF-causing", "CF-causing", "No"],
        ["", "p.Arg117His", "c.350G>A", "", "50", "0.0005", "Varying clinical consequence", "Varying clinical consequence", "No"],
        ["SYN_SPLICE", "p.?", "c.3718-2477C>T", "", "200", "0.001", "CF-causing", "Non CF-causing", "Yes"],
    ]

    create_synthetic_cftr2_xlsx(synthetic_file, test_rows)
    initial_sha = compute_file_sha256(synthetic_file)

    summary = ingest_cftr2(synthetic_file, db_session=in_memory_session, dry_run=False)

    # 1. Verify file immutability: SHA-256 did not change
    assert compute_file_sha256(synthetic_file) == initial_sha
    assert summary.file_sha256 == initial_sha

    # 2. Verify row accounting
    assert summary.rows_read == 6
    assert summary.valid_rows == 4
    assert summary.skipped_rows == 2  # 1 missing legacy name, 1 duplicate
    assert summary.duplicate_identifiers == 1
    assert "MISSING_LEGACY_NAME" in summary.exclusion_reasons
    assert "DUPLICATE_LEGACY_NAME" in summary.exclusion_reasons

    # 3. Verify scientific validation
    assert summary.wt_matches == 2  # SYN_G551D and SYN_F508del
    assert summary.wt_mismatches == 1  # SYN_MISMATCH (Ala508 instead of Phe508)
    assert summary.wt_out_of_bounds == 0

    # 4. Verify database persistence and preserved determinations
    records = in_memory_session.query(CFTR2Variant).all()
    assert len(records) == 4

    rec_g551d = in_memory_session.query(CFTR2Variant).filter_by(legacy_name="SYN_G551D").first()
    assert rec_g551d is not None
    assert rec_g551d.canonical_protein == "G551D"
    assert rec_g551d.is_valid_wt is True
    assert rec_g551d.validation_flag == "MATCH"
    assert rec_g551d.determination_current == "CF-causing"
    assert rec_g551d.alleles_reported == 3500

    rec_mismatch = in_memory_session.query(CFTR2Variant).filter_by(legacy_name="SYN_MISMATCH").first()
    assert rec_mismatch is not None
    assert rec_mismatch.is_valid_wt is False
    assert rec_mismatch.validation_flag == "WT_MISMATCH"
    assert "Reference mismatch at position 508" in rec_mismatch.validation_message
    assert rec_mismatch.determination_current == "Varying clinical consequence"

    rec_splice = in_memory_session.query(CFTR2Variant).filter_by(legacy_name="SYN_SPLICE").first()
    assert rec_splice is not None
    assert rec_splice.determination_current == "Non CF-causing"
    assert rec_splice.is_valid_wt is None

    # 5. Verify DataSource record
    data_source = in_memory_session.query(DataSource).first()
    assert data_source is not None
    assert data_source.source_name == "CFTR2"
    assert data_source.record_count == 4
    assert data_source.file_sha256 == initial_sha
