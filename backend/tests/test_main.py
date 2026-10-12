import pytest
from fastapi.testclient import TestClient

from backend.app.main import app
from backend.app.config import settings


@pytest.fixture
def client():
    return TestClient(app)


def test_root_endpoint(client):
    response = client.get("/")
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == settings.PROJECT_NAME
    assert data["version"] == settings.VERSION
    assert data["status"] == "online"
    assert data["docs_url"] == "/docs"


def test_health_endpoint(client):
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "healthy"
    assert data["service"] == settings.PROJECT_NAME
    assert data["version"] == settings.VERSION


def test_docs_accessible(client):
    response = client.get("/docs")
    assert response.status_code == 200
    assert "swagger-ui" in response.text.lower() or "html" in response.headers.get("content-type", "")


def test_openapi_schema(client):
    response = client.get("/openapi.json")
    assert response.status_code == 200
    schema = response.json()
    assert schema["info"]["title"] == settings.PROJECT_NAME
    assert schema["info"]["version"] == settings.VERSION
    assert "/health" in schema["paths"]
    assert "/" in schema["paths"]


def test_cors_allowed_origin(client):
    response = client.get(
        "/health",
        headers={"Origin": "http://localhost:3000"},
    )
    assert response.status_code == 200
    assert response.headers.get("access-control-allow-origin") == "http://localhost:3000"


def test_cors_disallowed_origin(client):
    response = client.get(
        "/health",
        headers={"Origin": "http://unauthorized-origin.com"},
    )
    assert response.status_code == 200
    assert "access-control-allow-origin" not in response.headers


def test_stats_endpoint(client):
    """Test GET /stats returns verified live aggregate database metrics."""
    response = client.get("/stats")
    assert response.status_code == 200
    data = response.json()
    assert data["database_connected"] is True
    assert data["evidence_loaded"] is True
    assert data["total_evidence_records"] == 4756
    assert data["grouped_protein_changes"] == 3353
    assert data["cftr2_records"] == 2092
    assert data["clinvar_records"] == 2664
    assert data["cftr_structures_count"] == 3
    assert data["control_structures_count"] == 2
    assert data["comparable_positions"] > 0


def test_variant_lookup_hgvs_normalization(client):
    """Test lookup by HGVS protein notation (e.g. p.Gly551Asp -> G551D)."""
    response = client.get("/variants/p.Gly551Asp")
    assert response.status_code == 200
    data = response.json()
    assert data["canonical_protein"] == "G551D"
    assert data["canonical_position"] == 551

    # Comparative metrics decomposition
    metrics = data["comparative_metrics"]
    assert metrics["all_structures_observed"] is True
    assert metrics["canonical_identity_match"] is True
    assert metrics["construct_context_caveat"] is True
    assert metrics["environment_comparability"] == "QUALIFIED_DIFFERENT_CONFORMATIONAL_STATES"
    assert len(metrics["comparison_caveats"]) >= 2

    # Verify contact definition
    assert "Heavy-atom" in data["structures"]["5UAK"]["contact_definition"]


def test_variant_lookup_e1371q_construct_context(client):
    """Test lookup by engineered construct mutation E1371Q resolves to position 1371 with E-to-Q identity."""
    response = client.get("/variants/E1371Q")
    assert response.status_code == 200
    data = response.json()
    assert data["canonical_protein"] == "E1371Q"
    assert data["canonical_position"] == 1371
    assert data["ref_aa"] == "E"
    assert data["alt_aa"] == "Q"
    assert data["source_associations"]["cftr2_count"] == 0
    assert data["source_associations"]["clinvar_count"] == 0
    assert data["structures"]["5UAK"]["structure_aa"] == "E"
    assert data["structures"]["6MSM"]["structure_aa"] == "Q"
    assert data["structures"]["8EIQ"]["structure_aa"] == "Q"
    assert data["structures"]["6MSM"]["failure_reason"] == "CONSTRUCT_SUBSTITUTION"
    assert data["structures"]["8EIQ"]["failure_reason"] == "CONSTRUCT_SUBSTITUTION"
    assert any(w["type"] == "CONSTRUCT_SUBSTITUTION" for w in data["warnings"])
    assert data["comparative_metrics"]["environment_comparability"] == "INCOMPARABLE_CONSTRUCT_SUBSTITUTION"


def test_zero_versus_null_distinction_api(client):
    """Verify that computed zero (G551) is distinct from absent/unmodeled null (F508 in 8EIQ)."""
    # G551 has genuine calculated 0.0 SASA and 0.0 RSA in 6MSM and 8EIQ
    res_g551 = client.get("/variants/G551D")
    assert res_g551.status_code == 200
    d_g551 = res_g551.json()
    assert d_g551["structures"]["6MSM"]["sasa_angstrom2"] == 0.0
    assert d_g551["structures"]["6MSM"]["rsa_ratio"] == 0.0
    assert d_g551["structures"]["6MSM"]["total_contacts"] == 30
    assert d_g551["structures"]["8EIQ"]["sasa_angstrom2"] == 0.0
    assert d_g551["structures"]["8EIQ"]["rsa_ratio"] == 0.0
    assert d_g551["structures"]["8EIQ"]["total_contacts"] == 30

    # F508 is deleted in 8EIQ: values MUST be null, not 0.0
    res_f508 = client.get("/variants/F508del")
    assert res_f508.status_code == 200
    d_f508 = res_f508.json()
    assert d_f508["structures"]["8EIQ"]["sasa_angstrom2"] is None
    assert d_f508["structures"]["8EIQ"]["rsa_ratio"] is None
    assert d_f508["structures"]["8EIQ"]["total_contacts"] is None
    assert d_f508["structures"]["8EIQ"]["secondary_structure"] is None
    assert d_f508["structures"]["8EIQ"]["failure_reason"] == "CONSTRUCT_DELETION"


def test_g551_dimer_ligand_warning(client):
    """Verify that G551 surfaces DIMER_INTERFACE_LIGANDS warning regarding bound ATP/Mg."""
    response = client.get("/variants/G551D")
    assert response.status_code == 200
    data = response.json()
    warnings = data["warnings"]
    assert any(w["type"] == "DIMER_INTERFACE_LIGANDS" for w in warnings)
    assert any("ATP1504" in w["message"] for w in warnings)


