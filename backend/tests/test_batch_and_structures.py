"""Comprehensive test suite for Phase 7 batch analysis and 3D structure endpoints."""

import pytest
from fastapi.testclient import TestClient
from backend.app.main import app


@pytest.fixture
def client():
    return TestClient(app)


# ==============================================================================
# Batch Analysis Tests
# ==============================================================================

def test_batch_variants_valid_multiple(client):
    """Test batch analysis with multiple valid identifiers."""
    payload = {"identifiers": ["G551D", "F508del", "M470V"]}
    resp = client.post("/variants/batch", json=payload)
    assert resp.status_code == 200
    data = resp.json()

    assert data["total_input_rows"] == 3
    assert data["unique_identifiers"] == 3
    assert data["resolved_count"] == 3
    assert data["unresolved_count"] == 0
    assert data["evidence_available_count"] == 3
    assert len(data["results"]) == 3

    assert data["results"][0]["row_number"] == 1
    assert data["results"][0]["canonical_protein"] == "G551D"
    assert data["results"][0]["status"] == "RESOLVED"
    assert data["results"][0]["has_structural_evidence"] is True

    assert data["results"][1]["row_number"] == 2
    assert data["results"][1]["canonical_protein"] == "F508del"
    assert data["results"][1]["status"] == "RESOLVED"

    assert data["results"][2]["row_number"] == 3
    assert data["results"][2]["canonical_protein"] == "M470V"


def test_batch_variants_mixed_statuses_and_duplicates(client):
    """Test mixed valid, invalid, duplicate, and blank rows while preserving order and row numbers."""
    payload = {
        "identifiers": [
            "G551D",
            "",
            "p.Gly551Asp",  # HGVS normalization to G551D
            "INVALID_NOTATION_999",
            "G551D",  # Duplicate row
            "E1371Q",  # Engineered construct mutation
        ]
    }
    resp = client.post("/variants/batch", json=payload)
    assert resp.status_code == 200
    data = resp.json()

    assert data["total_input_rows"] == 6
    assert data["unique_identifiers"] == 4  # G551D, p.Gly551Asp, INVALID_NOTATION_999, E1371Q
    assert data["resolved_count"] == 4  # G551D, p.Gly551Asp, G551D (dup), E1371Q
    assert data["unresolved_count"] == 2  # empty row, invalid notation

    results = data["results"]
    assert len(results) == 6

    # Row 1: G551D
    assert results[0]["row_number"] == 1
    assert results[0]["original_input"] == "G551D"
    assert results[0]["status"] == "RESOLVED"

    # Row 2: Blank
    assert results[1]["row_number"] == 2
    assert results[1]["status"] == "EMPTY"
    assert results[1]["failure_reason"] == "Blank input row."

    # Row 3: p.Gly551Asp
    assert results[2]["row_number"] == 3
    assert results[2]["original_input"] == "p.Gly551Asp"
    assert results[2]["normalized_identifier"] == "G551D"
    assert results[2]["status"] == "RESOLVED"

    # Row 4: INVALID_NOTATION_999
    assert results[3]["row_number"] == 4
    assert results[3]["status"] == "UNSUPPORTED_FORMAT"
    assert "Unable to resolve" in results[3]["failure_reason"]

    # Row 5: Duplicate G551D
    assert results[4]["row_number"] == 5
    assert results[4]["original_input"] == "G551D"
    assert results[4]["status"] == "RESOLVED"

    # Row 6: E1371Q
    assert results[5]["row_number"] == 6
    assert results[5]["original_input"] == "E1371Q"
    assert results[5]["canonical_position"] == 1371
    assert results[5]["status"] == "RESOLVED"


def test_batch_variants_f508del_null_integrity(client):
    """Verify that batch analysis preserves strictly null structural values for deleted F508 in 8EIQ."""
    payload = {"identifiers": ["F508del"]}
    resp = client.post("/variants/batch", json=payload)
    assert resp.status_code == 200
    res = resp.json()["results"][0]

    structs = res["structures"]
    assert structs["5UAK"]["sasa_angstrom2"] is not None
    assert structs["6MSM"]["sasa_angstrom2"] is not None

    # In 8EIQ, F508 is deleted: SASA, RSA, contacts must be null, never 0.0
    assert structs["8EIQ"]["sasa_angstrom2"] is None
    assert structs["8EIQ"]["rsa_ratio"] is None
    assert structs["8EIQ"]["total_contacts"] is None
    assert structs["8EIQ"]["mapping_status"] == "NOT_COVERED"
    assert structs["8EIQ"]["failure_reason"] == "CONSTRUCT_DELETION"


def test_batch_variants_g551d_structural_evidence(client):
    """Verify that batch analysis returns genuine computed zero SASA and ligand warnings for G551D."""
    payload = {"identifiers": ["G551D"]}
    resp = client.post("/variants/batch", json=payload)
    assert resp.status_code == 200
    res = resp.json()["results"][0]

    structs = res["structures"]
    # 5UAK apo SASA ~ 8.0
    assert structs["5UAK"]["sasa_angstrom2"] == 8.0
    assert structs["5UAK"]["total_contacts"] == 19

    # 6MSM and 8EIQ computed zero SASA
    assert structs["6MSM"]["sasa_angstrom2"] == 0.0
    assert structs["6MSM"]["rsa_ratio"] == 0.0
    assert structs["6MSM"]["total_contacts"] == 30
    assert structs["8EIQ"]["sasa_angstrom2"] == 0.0

    # Ligand warnings present
    assert any(w["type"] == "DIMER_INTERFACE_LIGANDS" for w in res["warnings"])


def test_batch_variants_empty_payload_rejected(client):
    """Verify that empty batch list is rejected with 400 Bad Request."""
    resp = client.post("/variants/batch", json={"identifiers": []})
    assert resp.status_code == 400
    assert "cannot be empty" in resp.json()["detail"]


def test_batch_variants_size_limit_exceeded(client):
    """Verify that batch size exceeding 500 items is rejected with 400 Bad Request."""
    oversized = ["G551D"] * 501
    resp = client.post("/variants/batch", json={"identifiers": oversized})
    assert resp.status_code == 400
    assert "exceeds maximum limit of 500" in resp.json()["detail"]


# ==============================================================================
# 3D Structure Coordinate Delivery Tests
# ==============================================================================

def test_structure_coordinates_valid(client):
    """Verify that supported structure coordinates are served with correct MIME type and caching headers."""
    for sid in ["5UAK", "6MSM", "8EIQ", "1AV5", "7JQ4"]:
        resp = client.get(f"/structures/{sid}/coordinates")
        assert resp.status_code == 200
        assert "chemical/x-mmcif" in resp.headers.get("content-type", "")
        assert "max-age=86400" in resp.headers.get("cache-control", "")
        assert len(resp.content) > 100000


def test_structure_coordinates_invalid_and_traversal(client):
    """Verify that unsupported structure IDs and path traversal attempts are safely rejected with 404."""
    # Unsupported structure
    resp_unknown = client.get("/structures/UNKNOWN_PDB/coordinates")
    assert resp_unknown.status_code == 404
    assert "allowlist" in resp_unknown.json()["detail"]

    # Traversal attempt
    resp_traversal = client.get("/structures/..%2F..%2Fetc%2Fpasswd/coordinates")
    assert resp_traversal.status_code == 404


def test_structure_residue_mapping(client):
    """Verify canonical residue mapping to structure author chain and numbering."""
    # 5UAK pos 551 -> Chain A, Res 551, Glycine, MAPPED_OBSERVED
    resp = client.get("/structures/5UAK/residues/551")
    assert resp.status_code == 200
    data = resp.json()
    assert data["structure_id"] == "5UAK"
    assert data["canonical_position"] == 551
    assert data["auth_asym_id"] == "A"
    assert data["auth_seq_id"] == "551"
    assert data["structure_amino_acid"] == "G"
    assert data["mapping_status"] == "MAPPED_OBSERVED"
    assert data["has_observed_residue"] is True

    # 8EIQ pos 508 -> Deleted construct, NOT_COVERED
    resp_8eiq = client.get("/structures/8EIQ/residues/508")
    assert resp_8eiq.status_code == 200
    d_8eiq = resp_8eiq.json()
    assert d_8eiq["mapping_status"] == "NOT_COVERED"
    assert d_8eiq["has_observed_residue"] is False
    assert "Delta-F508" in d_8eiq["mapping_warning"]


def test_structure_metadata(client):
    """Verify structure metadata and polymer chains endpoint."""
    resp = client.get("/structures/5UAK/metadata")
    assert resp.status_code == 200
    data = resp.json()
    assert data["id"] == "5UAK"
    assert data["has_cftr"] is True
    assert len(data["chains"]) > 0
    cftr_chains = [c for c in data["chains"] if c["is_cftr"]]
    assert len(cftr_chains) >= 1
    assert cftr_chains[0]["auth_asym_id"] == "A"


# ==============================================================================
# E1371Q and E1371X Separation & Structural Evidence Integrity Tests
# ==============================================================================

def test_e1371q_engineered_construct_single_endpoint(client):
    """Verify explicit E1371Q single lookup preserves identity, features from DB, and zero patient records."""
    resp = client.get("/variants/E1371Q")
    assert resp.status_code == 200
    data = resp.json()

    assert data["canonical_protein"] == "E1371Q"
    assert data["canonical_position"] == 1371
    assert data["ref_aa"] == "E"
    assert data["alt_aa"] == "Q"
    assert data["domain_name"] == "NBD2"

    # Zero patient records attached
    assert data["source_associations"]["cftr2_count"] == 0
    assert data["source_associations"]["cftr2_records"] == []
    assert data["source_associations"]["clinvar_count"] == 0
    assert data["source_associations"]["clinvar_records"] == []

    # Real structural features from ResidueFeature
    s5 = data["structures"]["5UAK"]
    assert s5["structure_aa"] == "E"
    assert s5["mapping_status"] == "MAPPED_OBSERVED"
    assert s5["sasa_angstrom2"] == pytest.approx(44.47, abs=0.01)
    assert s5["rsa_ratio"] == pytest.approx(0.1994, abs=0.001)
    assert s5["total_contacts"] == 21
    assert s5["secondary_structure"] == "HELIX"

    s6 = data["structures"]["6MSM"]
    assert s6["structure_aa"] == "Q"
    assert s6["mapping_status"] == "REFERENCE_MISMATCH"
    assert s6["failure_reason"] == "CONSTRUCT_SUBSTITUTION"
    assert s6["sasa_angstrom2"] == pytest.approx(4.29, abs=0.01)
    assert s6["rsa_ratio"] == pytest.approx(0.0191, abs=0.001)
    assert s6["total_contacts"] == 35
    assert s6["secondary_structure"] == "COIL"
    assert "E1371Q" in s6["construct_note"]

    s8 = data["structures"]["8EIQ"]
    assert s8["structure_aa"] == "Q"
    assert s8["mapping_status"] == "REFERENCE_MISMATCH"
    assert s8["failure_reason"] == "CONSTRUCT_SUBSTITUTION"
    assert s8["sasa_angstrom2"] == pytest.approx(10.78, abs=0.01)
    assert s8["rsa_ratio"] == pytest.approx(0.0479, abs=0.001)
    assert s8["total_contacts"] == 34
    assert s8["secondary_structure"] == "COIL"
    assert "E1371Q" in s8["construct_note"]

    # Warnings and comparative metrics
    assert any(w["type"] == "CONSTRUCT_SUBSTITUTION" for w in data["warnings"])
    comp = data["comparative_metrics"]
    assert comp["construct_variation"] == "CONSTRUCT_SUBSTITUTION"
    assert comp["is_fully_comparable"] is False
    assert comp["delta_rsa"] == pytest.approx(0.1803, abs=0.001)
    assert comp["delta_contacts"] == 14
    assert comp["sec_struct_consistency"] == "DIVERGENT"


def test_e1371q_hgvs_synonym_resolution(client):
    """Verify HGVS synonym p.Glu1371Gln resolves identically to E1371Q."""
    resp = client.get("/variants/p.Glu1371Gln")
    assert resp.status_code == 200
    data = resp.json()
    assert data["canonical_protein"] == "E1371Q"
    assert data["canonical_position"] == 1371
    assert data["ref_aa"] == "E"
    assert data["alt_aa"] == "Q"
    assert data["source_associations"]["cftr2_count"] == 0


def test_e1371x_patient_nonsense_integrity(client):
    """Verify E1371X resolves exclusively to patient nonsense record with null structural metrics."""
    resp = client.get("/variants/E1371X")
    assert resp.status_code == 200
    data = resp.json()

    assert data["canonical_protein"] == "E1371X"
    assert data["canonical_position"] == 1371
    assert data["ref_aa"] == "E"
    assert data["alt_aa"] == "X"

    # Patient record present
    assert data["source_associations"]["cftr2_count"] == 1
    assert data["source_associations"]["cftr2_records"][0]["determination_verbatim"] == "CF-causing"

    # Structural metrics null due to nonsense class
    assert data["structures"]["5UAK"]["sasa_angstrom2"] is None
    assert data["structures"]["6MSM"]["sasa_angstrom2"] is None
    assert data["structures"]["8EIQ"]["sasa_angstrom2"] is None


def test_batch_single_consistency_and_separation(client):
    """Verify batch and single endpoint consistency across E1371Q, E1371X, G551D, and F508del."""
    batch_resp = client.post(
        "/variants/batch",
        json={"identifiers": ["E1371Q", "E1371X", "G551D", "F508del"]},
    )
    assert batch_resp.status_code == 200
    results = batch_resp.json()["results"]

    # 1. E1371Q
    r_e1371q = results[0]
    single_e1371q = client.get("/variants/E1371Q").json()
    assert r_e1371q["status"] == "RESOLVED"
    assert r_e1371q["canonical_protein"] == "E1371Q"
    assert r_e1371q["ref_aa"] == "E"
    assert r_e1371q["alt_aa"] == "Q"
    assert r_e1371q["source_record_count"] == 0
    assert r_e1371q["source_determinations"] == []
    assert r_e1371q["structures"]["5UAK"]["sasa_angstrom2"] == single_e1371q["structures"]["5UAK"]["sasa_angstrom2"]
    assert r_e1371q["structures"]["6MSM"]["sasa_angstrom2"] == single_e1371q["structures"]["6MSM"]["sasa_angstrom2"]
    assert r_e1371q["structures"]["8EIQ"]["sasa_angstrom2"] == single_e1371q["structures"]["8EIQ"]["sasa_angstrom2"]

    # 2. E1371X
    r_e1371x = results[1]
    single_e1371x = client.get("/variants/E1371X").json()
    assert r_e1371x["status"] == "RESOLVED"
    assert r_e1371x["canonical_protein"] == "E1371X"
    assert r_e1371x["ref_aa"] == "E"
    assert r_e1371x["alt_aa"] == "X"
    assert r_e1371x["source_record_count"] == 1
    assert r_e1371x["source_determinations"] == ["CF-causing"]
    assert r_e1371x["structures"]["5UAK"]["sasa_angstrom2"] is None

    # 3. G551D
    r_g551d = results[2]
    single_g551d = client.get("/variants/G551D").json()
    assert r_g551d["status"] == "RESOLVED"
    assert r_g551d["canonical_protein"] == "G551D"
    assert r_g551d["structures"]["5UAK"]["sasa_angstrom2"] == single_g551d["structures"]["5UAK"]["sasa_angstrom2"]
    assert r_g551d["structures"]["6MSM"]["sasa_angstrom2"] == 0.0
    assert r_g551d["structures"]["8EIQ"]["sasa_angstrom2"] == 0.0

    # 4. F508del
    r_f508del = results[3]
    single_f508del = client.get("/variants/F508del").json()
    assert r_f508del["status"] == "RESOLVED"
    assert r_f508del["canonical_protein"] == "F508del"
    assert r_f508del["structures"]["8EIQ"]["sasa_angstrom2"] is None
    assert r_f508del["structures"]["8EIQ"]["failure_reason"] == "CONSTRUCT_DELETION"
