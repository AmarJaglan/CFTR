"""Automated tests for Phase 3 CFTR Structural Mapping Engine.

Verifies mmCIF parsing with Gemmi, canonical sequence correspondence,
unresolved residue detection, author vs label numbering separation,
reference mismatches (E1371Q), Delta-F508 deletion in 8EIQ, negative control
filtering (1AV5, 7JQ4), database idempotency, and audit generation.
"""

from pathlib import Path
import tempfile
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.app.db.models import Base, ResidueMapping, StructuralChain, Structure
from backend.app.scientific.mapping import (
    STATUS_MAPPED_OBSERVED,
    STATUS_MAPPED_UNRESOLVED,
    STATUS_NOT_COVERED,
    STATUS_REFERENCE_MISMATCH,
    TYPE_DOMAIN,
    TYPE_FULL_LENGTH,
    TYPE_NON_CFTR,
    export_audit_csvs,
    generate_audit_markdown,
    map_all_structures,
    map_structure,
    persist_mapping_results,
)
from backend.app.scientific.mmcif import (
    ChainMetadata,
    PolymerResidueInfo,
    StructureMetadata,
    parse_mmcif_file,
)
from backend.app.scientific.reference import CFTR_CANONICAL_SEQUENCE


RAW_DATA_DIR = Path("data/raw")
EXPECTED_CIF_FILES = [
    RAW_DATA_DIR / "1AV5.cif",
    RAW_DATA_DIR / "5UAK(1).cif",
    RAW_DATA_DIR / "6MSM(1).cif",
    RAW_DATA_DIR / "7jq4.cif",
    RAW_DATA_DIR / "8EIQ(1).cif",
]


@pytest.fixture
def available_cif_files():
    """Return available CIF files, resolving any symlink or naming differences."""
    resolved = []
    for p in EXPECTED_CIF_FILES:
        if p.exists():
            resolved.append(p)
        else:
            alt = p.parent / p.name.replace("(1)", "")
            if alt.exists():
                resolved.append(alt)
    return resolved


@pytest.fixture
def in_memory_session():
    """Provide an isolated in-memory SQLite database session."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(bind=engine)
    session_factory = sessionmaker(bind=engine)
    session = session_factory()
    try:
        yield session
    finally:
        session.close()


def test_parse_real_structures_metadata(available_cif_files):
    """Verify all input structures parse with Gemmi and extract metadata."""
    assert len(available_cif_files) == 5, f"Expected 5 CIF files, found {len(available_cif_files)}"

    parsed = [parse_mmcif_file(p) for p in available_cif_files]
    pdb_ids = {m.structure_id for m in parsed}

    assert "5UAK" in pdb_ids
    assert "6MSM" in pdb_ids
    assert "8EIQ" in pdb_ids
    assert "1AV5" in pdb_ids
    assert "7JQ4" in pdb_ids

    for m in parsed:
        assert m.file_sha256 is not None
        assert len(m.file_sha256) == 64
        assert m.experimental_method is not None
        assert m.chains, f"Structure {m.structure_id} has no chains"


def test_missing_file_handling():
    """Verify FileNotFoundError on non-existent CIF file."""
    with pytest.raises(FileNotFoundError):
        parse_mmcif_file(Path("data/raw/non_existent_structure.cif"))


def test_malformed_cif_handling():
    """Verify ValueError when parsing a corrupt CIF file."""
    with tempfile.NamedTemporaryFile(suffix=".cif", mode="w", delete=False) as tf:
        tf.write("data_corrupt\n_invalid_syntax { [[[ not a cif")
        temp_path = Path(tf.name)

    try:
        with pytest.raises(ValueError, match="Failed to parse CIF syntax"):
            parse_mmcif_file(temp_path)
    finally:
        temp_path.unlink()


def test_chain_and_entity_classification(available_cif_files):
    """Verify CFTR vs non-CFTR chains and negative controls."""
    results = map_all_structures(available_cif_files)
    res_by_id = {r.structure_id: r for r in results}

    # 5UAK, 6MSM, 8EIQ must have CFTR = True, Type = FULL_LENGTH
    assert res_by_id["5UAK"].has_cftr is True
    assert res_by_id["5UAK"].structure_type == TYPE_FULL_LENGTH
    assert res_by_id["6MSM"].has_cftr is True
    assert res_by_id["6MSM"].structure_type == TYPE_FULL_LENGTH
    assert res_by_id["8EIQ"].has_cftr is True
    assert res_by_id["8EIQ"].structure_type == TYPE_FULL_LENGTH

    # 1AV5 and 7JQ4 must have CFTR = False, Type = NON_CFTR
    assert res_by_id["1AV5"].has_cftr is False
    assert res_by_id["1AV5"].structure_type == TYPE_NON_CFTR
    assert res_by_id["7JQ4"].has_cftr is False
    assert res_by_id["7JQ4"].structure_type == TYPE_NON_CFTR

    # Negative controls must not inject rows into CFTR residue table
    assert len(res_by_id["1AV5"].residue_mappings) == 0
    assert len(res_by_id["7JQ4"].residue_mappings) == 0

    # Chain B/R in 5UAK, 6MSM, 8EIQ must be non-CFTR
    for sid in ["5UAK", "6MSM", "8EIQ"]:
        r = res_by_id[sid]
        chain_a = next(cs for cs in r.chain_summaries if cs.label_asym_id == "A")
        assert chain_a.is_cftr is True
        chain_b = next(cs for cs in r.chain_summaries if cs.label_asym_id == "B")
        assert chain_b.is_cftr is False


def test_canonical_sequence_alignment_and_numbering(available_cif_files):
    """Verify canonical 1..1480 residue mappings for CFTR structures."""
    results = map_all_structures(available_cif_files)
    res_5uak = next(r for r in results if r.structure_id == "5UAK")

    assert len(res_5uak.residue_mappings) == 1480

    # Check position order
    positions = [rm.canonical_position for rm in res_5uak.residue_mappings]
    assert positions == list(range(1, 1481))

    # Check canonical amino acids match UniProt P13569
    for rm in res_5uak.residue_mappings:
        expected_aa = CFTR_CANONICAL_SEQUENCE[rm.canonical_position - 1]
        assert rm.canonical_amino_acid == expected_aa


def test_unresolved_residues_detection(available_cif_files):
    """Verify that unmodeled residues have coordinates flagged as False and status MAPPED_UNRESOLVED."""
    results = map_all_structures(available_cif_files)
    res_5uak = next(r for r in results if r.structure_id == "5UAK")
    map_by_pos = {rm.canonical_position: rm for rm in res_5uak.residue_mappings}

    # Residues 1..4 are unresolved N-terminus in 5UAK
    for pos in [1, 2, 3, 4]:
        rec = map_by_pos[pos]
        assert rec.has_observed_residue is False
        assert rec.has_ca_atom is False
        assert rec.mapping_status == STATUS_MAPPED_UNRESOLVED
        assert rec.reference_identity_match is True

    # Residues in ordered regions (e.g. pos 5..10) are observed
    for pos in [5, 6, 7, 8, 9, 10]:
        rec = map_by_pos[pos]
        assert rec.has_observed_residue is True
        assert rec.has_ca_atom is True
        assert rec.mapping_status == STATUS_MAPPED_OBSERVED

    # Ensure no unobserved residue is reported as MAPPED_OBSERVED
    for rm in res_5uak.residue_mappings:
        if not rm.has_observed_residue:
            assert rm.mapping_status != STATUS_MAPPED_OBSERVED


def test_author_vs_label_numbering_and_shifts(available_cif_files):
    """Verify distinction between canonical, label_seq_id, and auth_seq_id."""
    results = map_all_structures(available_cif_files)
    res_8eiq = next(r for r in results if r.structure_id == "8EIQ")
    map_8eiq = {rm.canonical_position: rm for rm in res_8eiq.residue_mappings}

    # Position 509 (Glycine) in 8EIQ:
    # construct label_seq_id is 508 (due to Delta-F508 deletion)
    # auth_seq_id is '509' (author numbering maintained standard UniProt)
    rec_509 = map_8eiq[509]
    assert rec_509.canonical_amino_acid == "G"
    assert rec_509.entity_sequence_position == 508
    assert rec_509.label_seq_id == 508
    assert rec_509.auth_seq_id == "509"
    assert rec_509.label_seq_id != int(rec_509.auth_seq_id)

    # Position 1371 in 8EIQ:
    # construct label_seq_id is 1370
    # auth_seq_id is '1371'
    rec_1371 = map_8eiq[1371]
    assert rec_1371.label_seq_id == 1370
    assert rec_1371.auth_seq_id == "1371"


def test_reference_mismatch_e1371q(available_cif_files):
    """Verify that catalytic mutant E1371Q in 6MSM and 8EIQ is detected as REFERENCE_MISMATCH."""
    results = map_all_structures(available_cif_files)

    for sid in ["6MSM", "8EIQ"]:
        r = next(res for res in results if res.structure_id == sid)
        map_by_pos = {rm.canonical_position: rm for rm in r.residue_mappings}
        rec = map_by_pos[1371]

        assert rec.canonical_amino_acid == "E"
        assert rec.structure_amino_acid == "Q"
        assert rec.reference_identity_match is False
        assert rec.mapping_status == STATUS_REFERENCE_MISMATCH
        assert rec.has_observed_residue is True
        assert rec.has_ca_atom is True
        assert "E1371Q" in rec.mapping_warning


def test_delta_f508_deletion_in_8eiq(available_cif_files):
    """Verify authentic Delta-F508 patient deletion in 8EIQ is marked NOT_COVERED."""
    results = map_all_structures(available_cif_files)
    res_8eiq = next(r for r in results if r.structure_id == "8EIQ")
    map_by_pos = {rm.canonical_position: rm for rm in res_8eiq.residue_mappings}

    rec_508 = map_by_pos[508]
    assert rec_508.canonical_amino_acid == "F"
    assert rec_508.structure_amino_acid is None
    assert rec_508.entity_sequence_position is None
    assert rec_508.label_seq_id is None
    assert rec_508.auth_seq_id is None
    assert rec_508.has_observed_residue is False
    assert rec_508.has_ca_atom is False
    assert rec_508.mapping_status == STATUS_NOT_COVERED
    assert "Delta-F508" in rec_508.mapping_warning

    # Check chain summary detected deletions
    chain_a = next(cs for cs in res_8eiq.chain_summaries if cs.label_asym_id == "A")
    assert "Delta-F508" in chain_a.detected_deletions


def test_synthetic_cftr_domain_classification():
    """Verify that an isolated CFTR domain (NBD1) is classified as DOMAIN, not FULL_LENGTH."""
    # NBD1 sequence: residues 389..678 (290 aa)
    nbd1_seq = CFTR_CANONICAL_SEQUENCE[388:678]

    residues = [
        PolymerResidueInfo(
            seq_id=idx + 1,
            mon_id="ALA",
            one_letter_code=aa,
            auth_seq_id=str(389 + idx),
            auth_comp_id="ALA",
            auth_asym_id="A",
            label_asym_id="A",
            is_observed=True,
            has_ca=True,
        )
        for idx, aa in enumerate(nbd1_seq)
    ]

    chain = ChainMetadata(
        label_asym_id="A",
        auth_asym_id="A",
        entity_id="1",
        entity_description="CFTR NBD1 Domain",
        organism="Homo sapiens",
        polymer_sequence=nbd1_seq,
        polymer_length=len(nbd1_seq),
        residues=residues,
        observed_count=len(residues),
        ca_count=len(residues),
    )

    structure_meta = StructureMetadata(
        structure_id="TEST_NBD1",
        source_filename="test_nbd1.cif",
        file_sha256="0" * 64,
        title="CFTR NBD1 Domain Construct",
        chains={"A": chain},
        polymer_entities={"1": {"description": "CFTR NBD1 Domain"}},
        uniprot_accessions={"1": "P13569"},
    )

    res = map_structure(structure_meta)
    assert res.has_cftr is True
    assert res.structure_type == TYPE_DOMAIN  # Explicitly domain-level!
    assert res.chain_summaries[0].structure_type == TYPE_DOMAIN
    assert res.chain_summaries[0].canonical_coverage < 0.50
    assert len(res.residue_mappings) == 1480

    # Outside NBD1: NOT_COVERED
    map_by_pos = {rm.canonical_position: rm for rm in res.residue_mappings}
    assert map_by_pos[100].mapping_status == STATUS_NOT_COVERED
    # Inside NBD1: MAPPED_OBSERVED
    assert map_by_pos[500].mapping_status == STATUS_MAPPED_OBSERVED


def test_idempotent_database_persistence(available_cif_files, in_memory_session):
    """Verify that repeatedly persisting mapping results produces zero duplicates."""
    results = map_all_structures(available_cif_files)

    # First persistence run
    persist_mapping_results(results, in_memory_session)

    struct_count_1 = in_memory_session.query(Structure).count()
    chain_count_1 = in_memory_session.query(StructuralChain).count()
    res_count_1 = in_memory_session.query(ResidueMapping).count()

    assert struct_count_1 == 5
    assert chain_count_1 == 9
    assert res_count_1 == 4440

    # Second persistence run (idempotent)
    persist_mapping_results(results, in_memory_session)

    struct_count_2 = in_memory_session.query(Structure).count()
    chain_count_2 = in_memory_session.query(StructuralChain).count()
    res_count_2 = in_memory_session.query(ResidueMapping).count()

    assert struct_count_2 == struct_count_1
    assert chain_count_2 == chain_count_1
    assert res_count_2 == res_count_1


def test_audit_csv_and_markdown_generation(available_cif_files, tmp_path):
    """Verify audit CSVs and markdown summary report generation."""
    results = map_all_structures(available_cif_files)

    csv_dict = export_audit_csvs(results, tmp_path)
    md_file = generate_audit_markdown(results, tmp_path / "structural_mapping_summary.md")

    # Check files exist
    assert csv_dict["structure_inventory"].exists()
    assert csv_dict["chain_summary"].exists()
    assert csv_dict["residue_mapping"].exists()
    assert csv_dict["warnings"].exists()
    assert md_file.exists()

    # Check non-empty content
    assert csv_dict["structure_inventory"].stat().st_size > 500
    assert csv_dict["chain_summary"].stat().st_size > 500
    assert csv_dict["residue_mapping"].stat().st_size > 100_000  # 4440 rows
    assert csv_dict["warnings"].stat().st_size > 500
    assert md_file.stat().st_size > 3000

    # Check markdown contents
    md_text = md_file.read_text(encoding="utf-8")
    assert "Delta-F508" in md_text
    assert "E1371Q" in md_text
    assert "1AV5" in md_text
    assert "7JQ4" in md_text
    assert "PreScission" in md_text or "SNSLEVLFQ" in md_text
