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
