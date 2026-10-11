# CFTR-StructMap

CFTR-StructMap is a scientific bioinformatics application for reproducible structural-context analysis of CFTR (Cystic Fibrosis Transmembrane Conductance Regulator) variants. It integrates variant annotations with experimentally determined 3D protein structures to generate interpretable structural-context evidence and prioritization features (not clinical pathogenicity predictions).

## Project Setup

### Prerequisites

- macOS (Apple Silicon / arm64) or Linux
- Python 3.13+
- Virtual environment (`.venv`)

### Environment Initialization

Activate the virtual environment:

```bash
source .venv/bin/activate
```

Dependencies are managed in `backend/requirements.txt`:

```bash
pip install -r backend/requirements.txt
```

## Running the Backend

Start the FastAPI development server:

```bash
uvicorn backend.app.main:app --reload
```

The server binds to `http://127.0.0.1:8000` by default.

### Interactive API Documentation

- Swagger UI: [http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs)
- ReDoc: [http://127.0.0.1:8000/redoc](http://127.0.0.1:8000/redoc)
- Health check: [http://127.0.0.1:8000/health](http://127.0.0.1:8000/health)

### Configuration

CORS origins can be customized using the `CORS_ORIGINS` environment variable:

```bash
export CORS_ORIGINS="http://localhost:3000,http://127.0.0.1:3000"
```

## Running Tests

Execute automated tests with pytest:

```bash
pytest -v
```

Or directly via the virtual environment binary:

```bash
./.venv/bin/pytest -v
```

## Scientific Data Pipeline CLI

CFTR-StructMap includes reproducible CLI commands for data ingestion, mmCIF parsing, and residue mapping:

```bash
# Ingest CFTR2 and ClinVar clinical variant records
python -m backend.app.cli import-all

# Parse mmCIF structures with Gemmi and map to canonical CFTR (P13569)
python -m backend.app.cli map-structures

# Check database counts and provenance status
python -m backend.app.cli status
```

## Repository Structure

```text
CFTR-StructMap/
├── backend/
│   ├── app/
│   │   ├── __init__.py
│   │   ├── config.py
│   │   └── main.py
│   ├── tests/
│   │   ├── __init__.py
│   │   └── test_main.py
│   ├── requirements.txt
│   └── __init__.py
├── data/
│   ├── raw/           # Raw downloaded files (ClinVar, CFTR2, mmCIF) - excluded from git
│   └── processed/     # Derived tabular/SQLite artifacts - excluded from git
├── docs/              # Scientific documentation & architectural specs
├── pytest.ini         # Pytest configuration
├── .gitignore         # Ignores caches, .venv, data, and DB artifacts
└── README.md
```
