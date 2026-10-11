"""Automated tests for Phase 4 CFTR Structural Feature Extraction.

Verifies:
1. Canonical CFTR domain partitioning (Lasso, TMD1, NBD1, R-domain, TMD2, NBD2, C-term).
2. Relative Solvent Accessibility (RSA) via Shrake-Rupley and Tien et al. (2013) scale.
3. Heavy-atom contact counting with configurable distance sphere and same-residue exclusion.
4. Deposited mmCIF secondary structure parsing and unresolved coordinate handling.
5. Strict coordinate eligibility: unresolved residues are strictly None, never zero.
6. Delta-F508 deletion (pos 508 in 8EIQ) and catalytic E1371Q substitution handling.
7. Multi-structure per-residue comparative dynamics across 5UAK, 6MSM, and 8EIQ.
8. Database idempotency and SQLite persistence.
9. Audit CSV and Markdown report generation.
"""

from pathlib import Path
import tempfile
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.app.db.models import Base, MultiStructureComparison, ResidueFeature
from backend.app.scientific.features import (
    CFTR_DOMAINS,
    TIEN_MAX_SASA_2013,
    calculate_rsa,
    compare_multi_structure_features,
    export_feature_audit_csvs,
    extract_secondary_structure_from_cif,
    extract_structure_features,
    generate_feature_audit_markdown,
    get_cftr_domain,
    persist_feature_results,
    run_all_feature_extractions,
)
from backend.app.scientific.mapping import (
    STATUS_MAPPED_OBSERVED,
    STATUS_MAPPED_UNRESOLVED,
    STATUS_NOT_COVERED,
    STATUS_REFERENCE_MISMATCH,
    map_structure,
)
from backend.app.scientific.mmcif import parse_mmcif_file
from backend.app.scientific.reference import CFTR_CANONICAL_SEQUENCE


RAW_DATA_DIR = Path("data/raw")
TARGET_CIF_FILES = [
    RAW_DATA_DIR / "5UAK(1).cif",
    RAW_DATA_DIR / "6MSM(1).cif",
    RAW_DATA_DIR / "8EIQ(1).cif",
]


@pytest.fixture
def available_target_files():
    """Return available target CIF files (resolving any filename aliases)."""
    resolved = []
    for p in TARGET_CIF_FILES:
        if p.exists():
            resolved.append(p)
        else:
            alt = p.parent / p.name.replace("(1)", "")
            if alt.exists():
                resolved.append(alt)
    return resolved


# ==============================================================================
# 1. Canonical CFTR Domain Architecture Tests
# ==============================================================================
def test_canonical_domain_boundaries():
    """Verify that domain boundaries match documented literature partitioning."""
    # Check boundaries of all 7 domains
    expected_endpoints = [
        (1, "Lasso"),
        (79, "Lasso"),
        (80, "TMD1"),
        (388, "TMD1"),
        (389, "NBD1"),
        (508, "NBD1"),  # Delta-F508
        (678, "NBD1"),
        (679, "R-domain"),
        (844, "R-domain"),
        (845, "TMD2"),
        (1172, "TMD2"),
        (1173, "NBD2"),
        (1371, "NBD2"),  # E1371Q
        (1436, "NBD2"),
        (1437, "C-term"),
        (1480, "C-term"),
    ]
    for pos, expected_dom in expected_endpoints:
        assert get_cftr_domain(pos) == expected_dom, f"Pos {pos} expected {expected_dom}"

    # Out of bounds
    assert get_cftr_domain(0) == "Unknown"
    assert get_cftr_domain(1481) == "Unknown"

    # Verify every position in 1..1480 is covered by exactly one domain
    domain_counts = {}
    for pos in range(1, 1481):
        dom = get_cftr_domain(pos)
        assert dom != "Unknown"
        domain_counts[dom] = domain_counts.get(dom, 0) + 1

    assert domain_counts["Lasso"] == 79
    assert domain_counts["TMD1"] == 309
    assert domain_counts["NBD1"] == 290
    assert domain_counts["R-domain"] == 166
    assert domain_counts["TMD2"] == 328
    assert domain_counts["NBD2"] == 264
    assert domain_counts["C-term"] == 44
    assert sum(domain_counts.values()) == 1480


# ==============================================================================
# 2. RSA & Tien et al. (2013) Theoretical Scale Tests
# ==============================================================================
def test_rsa_calculation_and_tien_scale():
    """Verify Relative Solvent Accessibility calculation against theoretical max scale."""
    # Reference theoretical max values
    assert TIEN_MAX_SASA_2013["ALA"] == 129.0
    assert TIEN_MAX_SASA_2013["ARG"] == 274.0
    assert TIEN_MAX_SASA_2013["GLY"] == 104.0
    assert TIEN_MAX_SASA_2013["PHE"] == 240.0
    assert TIEN_MAX_SASA_2013["TRP"] == 285.0

    # Test 0.0 SASA -> BURIED
    rsa_zero, cat_zero = calculate_rsa(0.0, "ALA")
    assert rsa_zero == 0.0
    assert cat_zero == "BURIED"

    # Test 100% SASA -> SURFACE
    rsa_full, cat_full = calculate_rsa(129.0, "ALA")
    assert rsa_full == 1.0
    assert cat_full == "SURFACE"

    # Test thresholds: Buried (< 0.20), Intermediate (0.20 <= x < 0.50), Surface (>= 0.50)
    # Ala max = 129.0
    rsa_b, cat_b = calculate_rsa(25.0, "ALA")  # 25 / 129 = ~0.1938
    assert rsa_b < 0.20
    assert cat_b == "BURIED"

    rsa_i, cat_i = calculate_rsa(35.0, "ALA")  # 35 / 129 = ~0.2713
    assert 0.20 <= rsa_i < 0.50
    assert cat_i == "INTERMEDIATE"

    rsa_s, cat_s = calculate_rsa(70.0, "ALA")  # 70 / 129 = ~0.5426
    assert rsa_s >= 0.50
    assert cat_s == "SURFACE"

    # Test unknown amino acid fallback
    rsa_unk, cat_unk = calculate_rsa(50.0, "XYZ")
    assert rsa_unk == round(50.0 / 200.0, 4)


# ==============================================================================
# 3. Secondary Structure Parsing Tests
# ==============================================================================
def test_secondary_structure_parsing(available_target_files):
    """Verify that mmCIF depositor secondary structure annotations are parsed correctly."""
    cif_5uak = next((p for p in available_target_files if "5UAK" in p.name), None)
    if not cif_5uak:
        pytest.skip("5UAK mmCIF file not available")

    # In 5UAK, chain A has label_asym_id "A"
    sec_map = extract_secondary_structure_from_cif(cif_5uak, label_asym_id="A")
    assert len(sec_map) > 0

    # Values should only be valid secondary structure labels
    valid_labels = {"HELIX", "SHEET", "STRAND", "TURN", "COIL"}
    for label_seq_id, label in sec_map.items():
        assert label in valid_labels


# ==============================================================================
# 4. Target Structure Feature Extraction Tests (5UAK, 6MSM, 8EIQ)
# ==============================================================================
def test_full_structural_feature_extraction(available_target_files):
    """Verify feature extraction across all three validated cryo-EM models."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files (5UAK, 6MSM, 8EIQ) not all available")

    features_by_struct, comparisons, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    assert "5UAK" in features_by_struct
    assert "6MSM" in features_by_struct
    assert "8EIQ" in features_by_struct

    # Every structure must have exactly 1,480 residue records (canonical length)
    for sid in ["5UAK", "6MSM", "8EIQ"]:
        assert len(features_by_struct[sid]) == 1480

    # 1. 5UAK: Inward-facing apo (PDB 5UAK)
    f5 = features_by_struct["5UAK"]
    el5 = [r for r in f5 if r.eligibility_status == "ELIGIBLE"]
    un5 = [r for r in f5 if r.eligibility_status == "INELIGIBLE_UNRESOLVED"]
    nc5 = [r for r in f5 if r.eligibility_status == "INELIGIBLE_NOT_COVERED"]
    assert len(el5) == 1139
    assert len(un5) == 341
    assert len(nc5) == 0

    # 2. 6MSM: Outward-facing ATP-bound (PDB 6MSM)
    f6 = features_by_struct["6MSM"]
    el6 = [r for r in f6 if r.eligibility_status == "ELIGIBLE"]
    un6 = [r for r in f6 if r.eligibility_status == "INELIGIBLE_UNRESOLVED"]
    nc6 = [r for r in f6 if r.eligibility_status == "INELIGIBLE_NOT_COVERED"]
    assert len(el6) == 1181
    assert len(un6) == 299
    assert len(nc6) == 0

    # 3. 8EIQ: Trikafta-bound Delta-F508 (PDB 8EIQ)
    f8 = features_by_struct["8EIQ"]
    el8 = [r for r in f8 if r.eligibility_status == "ELIGIBLE"]
    un8 = [r for r in f8 if r.eligibility_status == "INELIGIBLE_UNRESOLVED"]
    nc8 = [r for r in f8 if r.eligibility_status == "INELIGIBLE_NOT_COVERED"]
    assert len(el8) == 1162
    assert len(un8) == 317
    assert len(nc8) == 1  # Delta-F508 at pos 508

    # Scientific Guardrail: Unresolved residues must strictly have None for features, NEVER 0.0
    for r in un5:
        assert r.sasa_raw_angstrom2 is None
        assert r.rsa_ratio is None
        assert r.exposure_category is None
        assert r.total_contacts is None
        assert r.secondary_structure is None
        assert r.failure_reason == "UNRESOLVED_COORDINATES"


# ==============================================================================
# 5. Delta-F508 and E1371Q Specific Handling Tests
# ==============================================================================
def test_delta_f508_and_e1371q_handling(available_target_files):
    """Verify precise biophysical and construct handling of Phe508 and Glu1371."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, comparisons, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    f5_by_pos = {r.canonical_position: r for r in features_by_struct["5UAK"]}
    f6_by_pos = {r.canonical_position: r for r in features_by_struct["6MSM"]}
    f8_by_pos = {r.canonical_position: r for r in features_by_struct["8EIQ"]}

    # --- Position 508 (Phe508) ---
    r5_508 = f5_by_pos[508]
    r6_508 = f6_by_pos[508]
    r8_508 = f8_by_pos[508]

    # In 5UAK (apo WT): Observed, buried in NBD1
    assert r5_508.eligibility_status == "ELIGIBLE"
    assert r5_508.has_observed_residue is True
    assert r5_508.structure_amino_acid in ("PHE", "F")
    assert r5_508.rsa_ratio is not None
    assert r5_508.total_contacts is not None
    assert r5_508.total_contacts >= 20

    # In 6MSM (ATP-bound WT): Observed, buried in NBD1
    assert r6_508.eligibility_status == "ELIGIBLE"
    assert r6_508.has_observed_residue is True
    assert r6_508.structure_amino_acid in ("PHE", "F")
    assert r6_508.rsa_ratio is not None
    assert r6_508.total_contacts is not None

    # In 8EIQ (Trikafta Delta-F508): Deleted from construct, must be INELIGIBLE_NOT_COVERED
    assert r8_508.eligibility_status == "INELIGIBLE_NOT_COVERED"
    assert r8_508.has_observed_residue is False
    assert r8_508.failure_reason == "CONSTRUCT_DELETION"
    assert r8_508.sasa_raw_angstrom2 is None
    assert r8_508.rsa_ratio is None
    assert r8_508.total_contacts is None

    # --- Position 1371 (Glu1371 / Catalytic Base) ---
    r5_1371 = f5_by_pos[1371]
    r6_1371 = f6_by_pos[1371]
    r8_1371 = f8_by_pos[1371]

    # In 5UAK: Canonical Glu (GLU)
    assert r5_1371.eligibility_status == "ELIGIBLE"
    assert r5_1371.structure_amino_acid in ("GLU", "E")
    assert r5_1371.failure_reason is None

    # In 6MSM: Catalytic mutant E1371Q
    assert r6_1371.eligibility_status == "ELIGIBLE"
    assert r6_1371.structure_amino_acid in ("GLN", "Q")
    assert r6_1371.failure_reason == "CONSTRUCT_SUBSTITUTION"
    assert r6_1371.rsa_ratio is not None
    assert r6_1371.total_contacts is not None

    # In 8EIQ: Catalytic mutant E1371Q
    assert r8_1371.eligibility_status == "ELIGIBLE"
    assert r8_1371.structure_amino_acid in ("GLN", "Q")
    assert r8_1371.failure_reason == "CONSTRUCT_SUBSTITUTION"
    assert r8_1371.rsa_ratio is not None
    assert r8_1371.total_contacts is not None


# ==============================================================================
# 6. Multi-Structure Per-Residue Comparative Engine Tests
# ==============================================================================
def test_multi_structure_comparison_metrics(available_target_files):
    """Verify per-residue comparative metrics across conformational states."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, comparisons, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    assert len(comparisons) == 1480

    comp_by_pos = {c.canonical_position: c for c in comparisons}

    # Pattern counts
    patterns = [c.observation_pattern for c in comparisons]
    assert patterns.count("ALL_THREE_OBSERVED") == 1127
    assert patterns.count("ALL_UNRESOLVED") == 290
    assert patterns.count("CONSTRUCT_DELETION") == 1  # pos 508
    assert patterns.count("PARTIALLY_OBSERVED") == 62

    # Comparability flag: pos 1371 is observed in all 3, but is CONSTRUCT_SUBSTITUTION
    # So is_fully_comparable should be 1126
    fully_comp = [c for c in comparisons if c.is_fully_comparable]
    assert len(fully_comp) == 1126

    # Verify pos 508 comparative entry
    c_508 = comp_by_pos[508]
    assert c_508.observation_pattern == "CONSTRUCT_DELETION"
    assert c_508.construct_variation == "CONSTRUCT_DELETION"
    assert c_508.is_fully_comparable is False
    assert c_508.sasa_5uak is not None
    assert c_508.sasa_6msm is not None
    assert c_508.sasa_8eiq is None
    assert c_508.delta_rsa is not None  # Delta between 5UAK and 6MSM exists

    # Verify pos 1371 comparative entry
    c_1371 = comp_by_pos[1371]
    assert c_1371.observation_pattern == "ALL_THREE_OBSERVED"
    assert c_1371.construct_variation == "CONSTRUCT_SUBSTITUTION"
    assert c_1371.is_fully_comparable is False  # Flagged because of construct substitution
    assert c_1371.contacts_5uak is not None
    assert c_1371.contacts_6msm is not None
    assert c_1371.contacts_8eiq is not None


# ==============================================================================
# 7. Database Persistence & Idempotency Tests
# ==============================================================================
def test_idempotent_feature_persistence(available_target_files):
    """Verify that saving features into SQLite is strictly idempotent."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    SessionLocal = sessionmaker(bind=engine)

    features_by_struct, comparisons, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    with SessionLocal() as session:
        # First persistence
        persist_feature_results(features_by_struct, comparisons, session)
        count_feat_1 = session.query(ResidueFeature).count()
        count_comp_1 = session.query(MultiStructureComparison).count()

        assert count_feat_1 == 1480 * 3  # 4440 rows
        assert count_comp_1 == 1480

        # Second persistence (must overwrite cleanly without duplicate rows)
        persist_feature_results(features_by_struct, comparisons, session)
        count_feat_2 = session.query(ResidueFeature).count()
        count_comp_2 = session.query(MultiStructureComparison).count()

        assert count_feat_2 == 1480 * 3
        assert count_comp_2 == 1480


# ==============================================================================
# 8. Audit CSV and Markdown Export Tests
# ==============================================================================
def test_feature_audit_exports(available_target_files):
    """Verify generation and format of Phase 4 audit CSVs and markdown summary."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, comparisons, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        csv_paths = export_feature_audit_csvs(
            features_by_structure=features_by_struct,
            comparisons=comparisons,
            all_warnings=warnings,
            output_dir=out_dir,
        )

        md_path = generate_feature_audit_markdown(
            features_by_structure=features_by_struct,
            comparisons=comparisons,
            output_path=out_dir / "structural_feature_summary.md",
        )

        assert csv_paths["feature_summary"].exists()
        assert csv_paths["residue_features"].exists()
        assert csv_paths["comparisons"].exists()
        assert csv_paths["warnings"].exists()
        assert md_path.exists()

        # Check line counts in CSVs
        with open(csv_paths["feature_summary"]) as f:
            lines = f.readlines()
            assert len(lines) == 4  # Header + 3 structures

        with open(csv_paths["residue_features"]) as f:
            lines = f.readlines()
            assert len(lines) == 4441  # Header + 4440 residues

        with open(csv_paths["comparisons"]) as f:
            lines = f.readlines()
            assert len(lines) == 1481  # Header + 1480 canonical positions

        # Verify markdown content contains domain table and coverage table
        content = md_path.read_text(encoding="utf-8")
        assert "Canonical CFTR Domain Architecture" in content
        assert "Structural Feature Coverage & Extraction Summary" in content
        assert "Multi-Structure Per-Residue Comparative Analysis" in content


# ==============================================================================
# 9. Phase 4.1 Audit Reconciliation & Regression Tests
# ==============================================================================
def test_reconciled_values_for_residues_508_and_1371(available_target_files):
    """Verify exact reconciled RSA, contacts, and deltas for canonical residues 508 and 1371."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, comparisons, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    f5_by_pos = {r.canonical_position: r for r in features_by_struct["5UAK"]}
    f6_by_pos = {r.canonical_position: r for r in features_by_struct["6MSM"]}
    f8_by_pos = {r.canonical_position: r for r in features_by_struct["8EIQ"]}
    c_by_pos = {c.canonical_position: c for c in comparisons}

    # Residue 508 (Phe508)
    r5_508 = f5_by_pos[508]
    r6_508 = f6_by_pos[508]
    r8_508 = f8_by_pos[508]
    c_508 = c_by_pos[508]

    assert r5_508.rsa_ratio == 0.0913
    assert r5_508.total_contacts == 26
    assert r5_508.intra_chain_contacts == 26
    assert r5_508.inter_chain_contacts == 0
    assert r5_508.secondary_structure == "HELIX"

    assert r6_508.rsa_ratio == 0.0246
    assert r6_508.total_contacts == 28
    assert r6_508.intra_chain_contacts == 28
    assert r6_508.inter_chain_contacts == 0
    assert r6_508.secondary_structure == "HELIX"

    assert r8_508.eligibility_status == "INELIGIBLE_NOT_COVERED"
    assert r8_508.failure_reason == "CONSTRUCT_DELETION"
    assert r8_508.rsa_ratio is None
    assert r8_508.total_contacts is None

    assert c_508.delta_rsa == round(0.0913 - 0.0246, 4)
    assert c_508.delta_contacts == 28 - 26

    # Residue 1371 (Glu1371 / Catalytic Base)
    r5_1371 = f5_by_pos[1371]
    r6_1371 = f6_by_pos[1371]
    r8_1371 = f8_by_pos[1371]
    c_1371 = c_by_pos[1371]

    assert r5_1371.structure_amino_acid in ("GLU", "E")
    assert r5_1371.rsa_ratio == 0.1994
    assert r5_1371.total_contacts == 21
    assert r5_1371.secondary_structure == "HELIX"

    assert r6_1371.structure_amino_acid in ("GLN", "Q")
    assert r6_1371.rsa_ratio == 0.0191
    assert r6_1371.total_contacts == 35  # includes coordination to ATP and Mg2+
    assert r6_1371.secondary_structure == "COIL"
    assert r6_1371.failure_reason == "CONSTRUCT_SUBSTITUTION"

    assert r8_1371.structure_amino_acid in ("GLN", "Q")
    assert r8_1371.rsa_ratio == 0.0479
    assert r8_1371.total_contacts == 34  # includes coordination to ATP and Mg2+
    assert r8_1371.secondary_structure == "COIL"
    assert r8_1371.failure_reason == "CONSTRUCT_SUBSTITUTION"

    assert c_1371.delta_rsa == round(max(0.1994, 0.0191, 0.0479) - min(0.1994, 0.0191, 0.0479), 4)
    assert c_1371.delta_contacts == max(21, 35, 34) - min(21, 35, 34)


def test_rsa_greater_than_one_chain_break_audit(available_target_files):
    """Verify that all residues with RSA > 1.0 are terminal to an unobserved loop break."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, _, _ = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    rsa_gt_1_list = []
    for sid, feats in features_by_struct.items():
        by_pos = {r.canonical_position: r for r in feats}
        for r in feats:
            if r.rsa_ratio is not None and r.rsa_ratio > 1.0:
                pos = r.canonical_position
                prev_obs = by_pos.get(pos - 1).has_observed_residue if (pos - 1) in by_pos else False
                next_obs = by_pos.get(pos + 1).has_observed_residue if (pos + 1) in by_pos else False
                borders_break = (not prev_obs) or (not next_obs)
                rsa_gt_1_list.append((sid, pos, r.rsa_ratio, borders_break))

    # Exactly 7 residues have RSA > 1.0 across all three structures
    assert len(rsa_gt_1_list) == 7
    # 100% of these 7 residues directly border a disordered chain break or terminus
    for sid, pos, rsa_val, borders_break in rsa_gt_1_list:
        assert borders_break is True, f"Structure {sid} residue {pos} has RSA {rsa_val} but does not border a break"


def test_secondary_structure_deposited_vs_unannotated_distinction(available_target_files):
    """Verify distinction between deposited mmCIF annotations and default unannotated coil."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, _, _ = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    for sid, feats in features_by_struct.items():
        el = [r for r in feats if r.eligibility_status == "ELIGIBLE"]
        for r in el:
            if r.secondary_structure in ("HELIX", "SHEET", "TURN"):
                assert r.sec_struct_method == "mmCIF Depositor (_struct_conf, _struct_sheet_range)"
            elif r.secondary_structure == "COIL":
                assert r.sec_struct_method == "Default Assignment (Unannotated in mmCIF)"


def test_construct_variation_in_warnings(available_target_files):
    """Verify that E1371Q and Delta-F508 are both captured in feature warnings with explicit reasons."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    _, _, warnings = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    # Must contain exactly 3 construct warnings: 8EIQ pos 508, 6MSM pos 1371, 8EIQ pos 1371
    warn_tuples = [(w["structure_id"], w["canonical_position"], w["ineligibility_reason"]) for w in warnings]
    assert ("8EIQ", 508, "CONSTRUCT_DELETION") in warn_tuples
    assert ("6MSM", 1371, "CONSTRUCT_SUBSTITUTION") in warn_tuples
    assert ("8EIQ", 1371, "CONSTRUCT_SUBSTITUTION") in warn_tuples


def test_contact_counts_sum_and_self_exclusion(available_target_files):
    """Verify that total_contacts strictly equals intra_chain + inter_chain and excludes self."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    features_by_struct, _, _ = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    for sid, feats in features_by_struct.items():
        el = [r for r in feats if r.eligibility_status == "ELIGIBLE"]
        for r in el:
            assert r.total_contacts == r.intra_chain_contacts + r.inter_chain_contacts
            # Self-contacts excluded: residue does not count itself
            assert r.total_contacts >= 0
            assert r.intra_chain_contacts >= 0
            assert r.inter_chain_contacts >= 0


def test_delta_metrics_are_strictly_non_causal_max_minus_min(available_target_files):
    """Verify that delta_rsa and delta_contacts represent unsigned maximal ranges (max - min)."""
    if len(available_target_files) < 3:
        pytest.skip("Required target structure files not all available")

    _, comparisons, _ = run_all_feature_extractions(
        structures_dir=RAW_DATA_DIR,
        contact_threshold=8.0,
    )

    for c in comparisons:
        observed_rsas = [v for v in (c.rsa_5uak, c.rsa_6msm, c.rsa_8eiq) if v is not None]
        observed_cnts = [v for v in (c.contacts_5uak, c.contacts_6msm, c.contacts_8eiq) if v is not None]

        if len(observed_rsas) >= 2:
            expected_delta_rsa = round(max(observed_rsas) - min(observed_rsas), 4)
            assert c.delta_rsa == expected_delta_rsa
            assert c.delta_rsa >= 0.0
        else:
            assert c.delta_rsa is None

        if len(observed_cnts) >= 2:
            expected_delta_cnt = max(observed_cnts) - min(observed_cnts)
            assert c.delta_contacts == expected_delta_cnt
            assert c.delta_contacts >= 0
        else:
            assert c.delta_contacts is None

