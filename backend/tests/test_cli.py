"""Automated tests for CFTR-StructMap CLI commands."""

from pathlib import Path
import subprocess
import sys

import pytest

from backend.app.cli import main
from backend.tests.test_cftr2_ingestion import create_synthetic_cftr2_xlsx
from backend.tests.test_clinvar_ingestion import create_synthetic_clinvar_tsv


def test_cli_help():
    result = subprocess.run(
        [sys.executable, "-m", "backend.app.cli", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "cftr-structmap-cli" in result.stdout
    assert "import-cftr2" in result.stdout
    assert "import-clinvar" in result.stdout
    assert "map-structures" in result.stdout
    assert "extract-features" in result.stdout
    assert "status" in result.stdout


def test_cli_dry_run_synthetic(tmp_path):
    cftr2_file = tmp_path / "cftr2.xlsx"
    clinvar_file = tmp_path / "clinvar.tsv"

    create_synthetic_cftr2_xlsx(
        cftr2_file,
        [["SYN_1", "p.Gly551Asp", "c.1652G>A", "", "10", "0.01", "CF-causing", "CF-causing", "No"]],
    )
    create_synthetic_clinvar_tsv(
        clinvar_file,
        [[
            "NM_000492.4(CFTR):c.1A>G (p.Met1Val)", "CFTR", "M1V", "Cystic fibrosis",
            "VCV000990001", "7", "117120149", "7", "117480095", "990001", "880001", "rs99901",
            "NC_000007.14:117480094:A:G", "single nucleotide variant", "missense variant",
            "Pathogenic", "Mar 17, 2020", "reviewed by expert panel", "", "", "", "", "", "",
        ]],
    )

    db_file = tmp_path / "test.db"

    # Test import-cftr2 dry-run
    res_cftr2 = subprocess.run(
        [sys.executable, "-m", "backend.app.cli", "--db", str(db_file), "import-cftr2", "--file", str(cftr2_file), "--dry-run"],
        capture_output=True,
        text=True,
    )
    assert res_cftr2.returncode == 0
    assert "CFTR2 Ingestion Summary" in res_cftr2.stdout

    # Test import-clinvar dry-run
    res_clinvar = subprocess.run(
        [sys.executable, "-m", "backend.app.cli", "--db", str(db_file), "import-clinvar", "--file", str(clinvar_file), "--dry-run"],
        capture_output=True,
        text=True,
    )
    assert res_clinvar.returncode == 0
    assert "ClinVar Ingestion Summary" in res_clinvar.stdout
