"""Automated tests for Phase 5: Variant-Level Structural Evidence Integration.

Validates:
1. Missense normalization and canonical reference validation.
2. Synonymous and unsupported variant class handling.
3. Source-record identity preservation (no merging of distinct records; ClinVar VCV aggregate documented).
4. Duplicate normalized protein changes across distinct source records (e.g. F311L).
5. CFTR2 and ClinVar verbatim classifications preserved without consensus synthesis.
6. M470V historical canonical-reference mismatch handling.
7. F508del authentic patient deletion absent in 8EIQ (CONSTRUCT_DELETION, null features).
8. E1371Q catalytic construct substitution handling in 6MSM and 8EIQ.
9. R-domain unmodeled/unresolved coordinates (evidence availability issue, not negative impact).
10. Domain context and per-structure feature joins.
11. Cross-structure comparability metrics and dynamics (non-causal observational ranges).
12. Idempotent SQLite persistence.
13. Audit CSV and Markdown summary export reproducibility.
14. Read-only FastAPI endpoints (/variants and /variants/{identifier}).
"""

import csv
from pathlib import Path
import tempfile
import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from backend.app.db.models import (
    Base,
    CFTR2Variant,
    ClinVarVariant,
    MultiStructureComparison,
    ResidueFeature,
    ResidueMapping,
    VariantStructuralEvidence,
    VariantStructuralSummary,
)
from backend.app.db.session import get_db, get_engine, get_session_factory
from backend.app.main import app
from backend.app.scientific.variant_evidence import (
    ELIGIBILITY_ELIGIBLE_MISSENSE,
    ELIGIBILITY_IN_FRAME_DELETION,
    ELIGIBILITY_REFERENCE_MISMATCH,
    ELIGIBILITY_SYNONYMOUS,
    ELIGIBILITY_UNRESOLVED_CHANGE,
    ELIGIBILITY_UNSUPPORTED_CLASS,
    build_variant_evidence,
    classify_cftr2_variant,
    classify_clinvar_variant,
    export_variant_audit_artifacts,
    persist_variant_evidence,
)


@pytest.fixture(scope="module")
def db_session():
    """Provide a database session connected to the real processed SQLite database."""
    engine = get_engine()
    session_factory = get_session_factory(engine)
    session = session_factory()
    yield session
    session.close()


@pytest.fixture(scope="module")
def api_client():
    """FastAPI TestClient fixture."""
    return TestClient(app)


# ==============================================================================
# 1. Variant Eligibility Classification Tests
# ==============================================================================


def test_classify_cftr2_missense():
    """Test CFTR2 missense variant with valid WT match is ELIGIBLE_MISSENSE."""
    var = CFTR2Variant(
        legacy_name="G551D",
        variant_type="missense",
        uniprot_pos=551,
        ref_aa="G",
        alt_aa="D",
        is_valid_wt=True,
        validation_flag="MATCH",
    )
    status, reason = classify_cftr2_variant(var)
    assert status == ELIGIBILITY_ELIGIBLE_MISSENSE
    assert reason is None


def test_classify_cftr2_synonymous():
    """Test CFTR2 synonymous variant is classified as SYNONYMOUS."""
    var = CFTR2Variant(
        legacy_name="T854T",
        variant_type="synonymous",
        uniprot_pos=854,
        ref_aa="T",
        alt_aa="T",
        is_valid_wt=True,
        validation_flag="MATCH",
    )
    status, reason = classify_cftr2_variant(var)
    assert status == ELIGIBILITY_SYNONYMOUS
    assert "Synonymous" in reason


def test_classify_cftr2_inframe_deletion():
    """Test CFTR2 in-frame deletion is classified as IN_FRAME_DELETION."""
    var = CFTR2Variant(
        legacy_name="F508del",
        variant_type="inframe_deletion",
        uniprot_pos=508,
        ref_aa="F",
        alt_aa=None,
        is_valid_wt=True,
        validation_flag="MATCH",
    )
    status, reason = classify_cftr2_variant(var)
    assert status == ELIGIBILITY_IN_FRAME_DELETION
    assert "deletion" in reason.lower()


def test_classify_cftr2_reference_mismatch():
    """Test CFTR2 M470V is classified as REFERENCE_MISMATCH."""
    var = CFTR2Variant(
        legacy_name="M470V",
        variant_type="missense",
        uniprot_pos=470,
        ref_aa="M",
        alt_aa="V",
        is_valid_wt=False,
        validation_flag="WT_MISMATCH",
    )
    status, reason = classify_cftr2_variant(var)
    assert status == ELIGIBILITY_REFERENCE_MISMATCH
    assert "mismatch" in reason.lower()


def test_classify_cftr2_unsupported_classes():
    """Test CFTR2 nonsense, frameshift, and unparsed are classified appropriately."""
    nonsense_var = CFTR2Variant(variant_type="nonsense", validation_flag="MATCH", uniprot_pos=542)
    assert classify_cftr2_variant(nonsense_var)[0] == ELIGIBILITY_UNSUPPORTED_CLASS

    frameshift_var = CFTR2Variant(variant_type="frameshift", validation_flag="MATCH", uniprot_pos=684)
    assert classify_cftr2_variant(frameshift_var)[0] == ELIGIBILITY_UNSUPPORTED_CLASS

    unparsed_var = CFTR2Variant(variant_type="unknown", validation_flag="NORMALIZATION_FAILED")
    assert classify_cftr2_variant(unparsed_var)[0] == ELIGIBILITY_UNRESOLVED_CHANGE


def test_classify_clinvar_snv():
    """Test ClinVar SNV missense variant with valid WT match is ELIGIBLE_MISSENSE."""
    var = ClinVarVariant(
        variation_id=7120,
        variant_type="single nucleotide variant",
        molecular_consequence="missense variant",
        uniprot_pos=551,
        ref_aa="G",
        alt_aa="D",
        is_valid_wt=True,
        validation_flag="MATCH",
    )
    status, reason = classify_clinvar_variant(var)
    assert status == ELIGIBILITY_ELIGIBLE_MISSENSE
    assert reason is None


def test_classify_clinvar_haplotype_and_indel():
    """Test ClinVar Haplotype and Indel records are marked as UNSUPPORTED_VARIANT_CLASS."""
    haplo = ClinVarVariant(
        variation_id=209047,
        variant_type="Haplotype",
        validation_flag="MATCH",
        uniprot_pos=117,
        ref_aa="R",
        alt_aa="H",
    )
    assert classify_clinvar_variant(haplo)[0] == ELIGIBILITY_UNSUPPORTED_CLASS

    indel = ClinVarVariant(
        variation_id=53596,
        variant_type="Indel",
        validation_flag="MATCH",
        uniprot_pos=973,
        ref_aa="L",
        alt_aa="F",
    )
    assert classify_clinvar_variant(indel)[0] == ELIGIBILITY_UNSUPPORTED_CLASS


# ==============================================================================
# 2. Source Identity & Verbatim Preservation Tests
# ==============================================================================


def test_source_identity_preservation_and_counts(db_session):
    """Verify all source records are retained without merging or row loss."""
    cftr2_count = db_session.query(CFTR2Variant).count()
    clinvar_count = db_session.query(ClinVarVariant).count()
    total_source = cftr2_count + clinvar_count

    ev_count = db_session.query(VariantStructuralEvidence).count()
    assert ev_count == total_source, f"Expected {total_source} records, got {ev_count}"

    cftr2_ev = db_session.query(VariantStructuralEvidence).filter_by(source_db="CFTR2").count()
    clinvar_ev = db_session.query(VariantStructuralEvidence).filter_by(source_db="ClinVar").count()
    assert cftr2_ev == cftr2_count
    assert clinvar_ev == clinvar_count


def test_duplicate_protein_changes_preserved_separately(db_session):
    """Verify identical protein changes across distinct ClinVar VariationIDs remain separate rows."""
    f311l_records = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="ClinVar", canonical_protein="F311L")
        .all()
    )
    assert len(f311l_records) == 3, "Expected 3 distinct ClinVar records for F311L"

    variation_ids = {r.source_record_id for r in f311l_records}
    assert len(variation_ids) == 3, "Each record must have a distinct VariationID"
    assert "7153" in variation_ids
    assert "1474955" in variation_ids
    assert "2871602" in variation_ids

    # Check that each preserves its authentic review status
    expert_panel = [r for r in f311l_records if r.source_record_id == "7153"][0]
    assert expert_panel.clinvar_review_status == "reviewed by expert panel"
    assert expert_panel.source_determination == "Pathogenic"


def test_verbatim_classifications_no_consensus(db_session):
    """Verify CFTR2 and ClinVar classifications are preserved verbatim without consensus synthesis."""
    g551d_cftr2 = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="CFTR2", canonical_protein="G551D")
        .first()
    )
    g551d_clinvar = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="ClinVar", canonical_protein="G551D")
        .first()
    )

    assert g551d_cftr2.source_determination == "CF-causing"
    assert g551d_clinvar.source_determination == "Pathogenic"

    # Overlapping variant with divergent wording
    r117h_cftr2 = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="CFTR2", canonical_protein="R117H")
        .first()
    )
    r117h_clinvar = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="ClinVar", source_record_id="7109")
        .first()
    )
    assert r117h_cftr2.source_determination == "Varying clinical consequence"
    assert r117h_clinvar.source_determination == "Pathogenic"


# ==============================================================================
# 3. Structural Evidence & Special Cases Tests
# ==============================================================================


def test_f508del_handling_and_8eiq_absence(db_session):
    """Verify F508del has wild-type features in 5UAK/6MSM and CONSTRUCT_DELETION (null) in 8EIQ."""
    f508del = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="CFTR2", source_record_id="F508del")
        .first()
    )
    assert f508del is not None
    assert f508del.canonical_position == 508
    assert f508del.eligibility_status == ELIGIBILITY_IN_FRAME_DELETION

    # 5UAK (apo full-length WT)
    assert f508del.mapping_status_5uak == "MAPPED_OBSERVED"
    assert f508del.feature_status_5uak == "ELIGIBLE"
    assert f508del.rsa_5uak == pytest.approx(0.0913, abs=0.001)
    assert f508del.contacts_5uak == 26

    # 6MSM (ATP-bound WT at pos 508)
    assert f508del.mapping_status_6msm == "MAPPED_OBSERVED"
    assert f508del.feature_status_6msm == "ELIGIBLE"
    assert f508del.rsa_6msm == pytest.approx(0.0246, abs=0.001)
    assert f508del.contacts_6msm == 28

    # 8EIQ (authentic Delta-F508 patient deletion construct)
    assert f508del.mapping_status_8eiq == "NOT_COVERED"
    assert f508del.feature_status_8eiq == "INELIGIBLE_NOT_COVERED"
    assert f508del.failure_reason_8eiq == "CONSTRUCT_DELETION"
    assert f508del.rsa_8eiq is None
    assert f508del.sasa_8eiq is None
    assert f508del.contacts_8eiq is None
    assert f508del.sec_struct_8eiq is None
    assert "Delta-F508" in f508del.construct_note_8eiq

    # Cross-structure comparative metrics
    assert f508del.observation_pattern == "CONSTRUCT_DELETION"
    assert f508del.construct_variation == "CONSTRUCT_DELETION"
    assert f508del.is_fully_comparable is False


def test_e1371q_catalytic_construct_substitution(db_session):
    """Verify position 1371 correctly flags catalytic E1371Q substitution in 6MSM and 8EIQ."""
    ev = (
        db_session.query(VariantStructuralEvidence)
        .filter(VariantStructuralEvidence.canonical_position == 1371)
        .first()
    )
    assert ev is not None
    assert ev.domain_name == "NBD2"

    # 5UAK is wild-type Glu
    assert ev.mapping_status_5uak == "MAPPED_OBSERVED"
    assert ev.structure_aa_5uak == "E"
    assert ev.construct_note_5uak is None

    # 6MSM has engineered Gln (E1371Q)
    assert ev.mapping_status_6msm == "REFERENCE_MISMATCH"
    assert ev.structure_aa_6msm == "Q"
    assert ev.failure_reason_6msm == "CONSTRUCT_SUBSTITUTION"
    assert "E1371Q" in ev.construct_note_6msm

    # 8EIQ has engineered Gln (E1371Q)
    assert ev.mapping_status_8eiq == "REFERENCE_MISMATCH"
    assert ev.structure_aa_8eiq == "Q"
    assert ev.failure_reason_8eiq == "CONSTRUCT_SUBSTITUTION"
    assert "E1371Q" in ev.construct_note_8eiq

    assert ev.construct_variation == "CONSTRUCT_SUBSTITUTION"
    assert ev.is_fully_comparable is False


def test_m470v_reference_mismatch(db_session):
    """Verify CFTR2 M470V is preserved verbatim and flagged as REFERENCE_MISMATCH with null features."""
    m470v = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="CFTR2", source_record_id="M470V")
        .first()
    )
    assert m470v is not None
    assert m470v.canonical_position == 470
    assert m470v.ref_aa == "M"
    assert m470v.alt_aa == "V"
    assert m470v.is_valid_wt is False
    assert m470v.validation_flag == "WT_MISMATCH"
    assert m470v.eligibility_status == ELIGIBILITY_REFERENCE_MISMATCH
    assert "mismatch" in m470v.ineligibility_reason.lower()

    # Features must be strictly None (not mapped to wild-type Valine coordinates)
    assert m470v.rsa_5uak is None
    assert m470v.rsa_6msm is None
    assert m470v.rsa_8eiq is None
    assert m470v.feature_status_5uak == "INELIGIBLE_REFERENCE_MISMATCH"

    # Contrast with ClinVar V470M (VariationID 7130) which correctly matches canonical Val
    v470m = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(source_db="ClinVar", source_record_id="7130")
        .first()
    )
    assert v470m is not None
    assert v470m.validation_flag == "MATCH"
    assert v470m.eligibility_status == ELIGIBILITY_ELIGIBLE_MISSENSE
    assert v470m.rsa_5uak is not None


def test_r_domain_unresolved_coordinates(db_session):
    """Verify variants in intrinsically disordered R-domain have null features with explicit reason."""
    r_domain_variants = (
        db_session.query(VariantStructuralEvidence)
        .filter(VariantStructuralEvidence.domain_name == "R-domain")
        .all()
    )
    assert len(r_domain_variants) > 0

    for r in r_domain_variants[:20]:
        # Features should be null across all structures
        assert r.rsa_5uak is None
        assert r.rsa_6msm is None
        assert r.rsa_8eiq is None
        assert r.contacts_5uak is None
        assert r.observation_pattern == "ALL_UNRESOLVED"
        assert r.feature_status_5uak in ("INELIGIBLE_UNRESOLVED", None)


def test_cross_structure_comparability(db_session):
    """Verify cross-structure deltas and comparability flags are populated correctly."""
    g551d = (
        db_session.query(VariantStructuralEvidence)
        .filter_by(canonical_protein="G551D")
        .first()
    )
    assert g551d is not None
    assert g551d.domain_name == "NBD1"
    assert g551d.observation_pattern == "ALL_THREE_OBSERVED"
    assert g551d.construct_variation == "CANONICAL_WT"
    assert g551d.is_fully_comparable is True
    assert g551d.delta_rsa is not None
    assert g551d.delta_contacts is not None
    assert g551d.sec_struct_consistency is not None


# ==============================================================================
# 4. Database Idempotency & Audit Outputs Tests
# ==============================================================================


def test_idempotent_persistence(db_session):
    """Verify persist_variant_evidence runs cleanly twice without duplicate rows."""
    count_before_ev = db_session.query(VariantStructuralEvidence).count()
    count_before_sum = db_session.query(VariantStructuralSummary).count()

    evidence_list, summary_list, _, _ = build_variant_evidence(db_session)
    n_ev, n_sum = persist_variant_evidence(evidence_list, summary_list, db_session)

    assert n_ev == count_before_ev
    assert n_sum == count_before_sum
    assert db_session.query(VariantStructuralEvidence).count() == count_before_ev
    assert db_session.query(VariantStructuralSummary).count() == count_before_sum


def test_reproducible_audit_export(db_session):
    """Verify export_variant_audit_artifacts generates valid CSVs and Markdown report."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        evidence_list, summary_list, comp_list, warn_list = build_variant_evidence(db_session)

        paths = export_variant_audit_artifacts(
            evidence_list=evidence_list,
            summary_list=summary_list,
            comparisons_list=comp_list,
            warnings_list=warn_list,
            output_dir=tmp_path,
        )

        for name, p in paths.items():
            assert p.exists(), f"Expected artifact {name} at {p}"
            assert p.stat().st_size > 0, f"Artifact {name} is empty"

        # Check evidence CSV row count (header + 4756 rows)
        with open(paths["evidence"], newline="", encoding="utf-8") as f:
            reader = list(csv.reader(f))
            assert len(reader) == len(evidence_list) + 1

        # Check summary CSV row count
        with open(paths["summary"], newline="", encoding="utf-8") as f:
            reader = list(csv.reader(f))
            assert len(reader) == len(summary_list) + 1

        # Check Markdown content
        md_text = paths["markdown"].read_text(encoding="utf-8")
        assert "CFTR-StructMap — Phase 5" in md_text
        assert "4,756" in md_text
        assert "F508del" in md_text
        assert "E1371Q" in md_text
        assert "M470V" in md_text


# ==============================================================================
# 5. Read-Only FastAPI Endpoints Tests
# ==============================================================================


def test_api_list_variants(api_client):
    """Test GET /variants pagination and filtering."""
    r = api_client.get("/variants?limit=10&offset=0")
    assert r.status_code == 200
    data = r.json()
    assert data["total"] == 4756
    assert len(data["results"]) == 10

    # Filter by position
    r_pos = api_client.get("/variants?position=508")
    assert r_pos.status_code == 200
    assert any(x["source_record_id"] == "F508del" for x in r_pos.json()["results"])

    # Filter by domain
    r_dom = api_client.get("/variants?domain=NBD1&limit=5")
    assert r_dom.status_code == 200
    for res in r_dom.json()["results"]:
        assert res["domain_name"] == "NBD1"


def test_api_get_variant_by_protein(api_client):
    """Test GET /variants/G551D returns complete structural and classification payload."""
    r = api_client.get("/variants/G551D")
    assert r.status_code == 200
    data = r.json()
    assert data["canonical_protein"] == "G551D"
    assert data["canonical_position"] == 551
    assert data["domain_name"] == "NBD1"
    assert data["source_associations"]["cftr2_count"] >= 1
    assert data["source_associations"]["clinvar_count"] >= 1

    # Check structure payloads
    assert "5UAK" in data["structures"]
    assert "6MSM" in data["structures"]
    assert "8EIQ" in data["structures"]
    assert data["structures"]["5UAK"]["rsa_ratio"] == pytest.approx(0.0769, abs=0.001)

    # Check guardrail notice
    assert "does not predict clinical pathogenicity" in data["guardrail_notice"].lower()


def test_api_get_variant_f508del(api_client):
    """Test GET /variants/F508del surfaces CONSTRUCT_DELETION in 8EIQ."""
    r = api_client.get("/variants/F508del")
    assert r.status_code == 200
    data = r.json()
    assert data["canonical_position"] == 508
    assert data["structures"]["8EIQ"]["feature_status"] == "INELIGIBLE_NOT_COVERED"
    assert data["structures"]["8EIQ"]["failure_reason"] == "CONSTRUCT_DELETION"
    assert data["structures"]["8EIQ"]["rsa_ratio"] is None
    assert any(w["type"] == "CONSTRUCT_DELETION" for w in data["warnings"])


def test_api_get_variant_m470v(api_client):
    """Test GET /variants/M470V surfaces REFERENCE_MISMATCH warning."""
    r = api_client.get("/variants/M470V")
    assert r.status_code == 200
    data = r.json()
    assert any(w["type"] == "REFERENCE_MISMATCH" for w in data["warnings"])
    assert data["structures"]["5UAK"]["feature_status"] == "INELIGIBLE_REFERENCE_MISMATCH"


def test_api_get_variant_position_1371(api_client):
    """Test GET /variants/1371 surfaces CONSTRUCT_SUBSTITUTION warning."""
    r = api_client.get("/variants/1371")
    assert r.status_code == 200
    data = r.json()
    assert any(w["type"] == "CONSTRUCT_SUBSTITUTION" for w in data["warnings"])
    assert data["structures"]["6MSM"]["mapping_status"] == "REFERENCE_MISMATCH"


def test_api_get_variant_not_found(api_client):
    """Test GET /variants/NONEXISTENT returns 404."""
    r = api_client.get("/variants/NONEXISTENT_XYZ")
    assert r.status_code == 404
