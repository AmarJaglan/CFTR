# CFTR-StructMap frontend

React + TypeScript + Vite workspace for the existing FastAPI backend. The frontend is intentionally additive: Python files, datasets, SQLite data, audit artifacts, and scientific logic remain in the repository's existing locations.

## Prerequisites

- Node.js 20+
- Python 3.11+ with the dependencies in `backend/requirements.txt`
- The backend database populated according to the existing project workflow

## Development

From the repository root, start FastAPI in one terminal:

```bash
uvicorn backend.app.main:app --reload --port 8000
```

Then start the frontend:

```bash
cd frontend
npm install
npm run dev
```

Vite proxies `/api/*` to `http://localhost:8000/*`, so browser code never hardcodes a development server URL. The current backend exposes `/health`, `/variants`, and `/variants/{identifier}`; the UI uses those contracts directly.

## Configuration

Set `VITE_API_BASE_URL` for a deployed API, for example:

```bash
VITE_API_BASE_URL=https://api.example.org npm run build
```

When it is unset, the frontend uses `/api` and therefore expects a reverse proxy in production. If the API is hosted on a different origin, configure FastAPI `CORS_ORIGINS` to include the frontend origin and expose the API through a secure deployment boundary.

## Tests and build

```bash
npm run test
npm run build
```

## Current limitations

- The backend has no batch-analysis endpoint, so Batch Analysis is an explicit unavailable state.
- The backend has no coordinate or mesh endpoint, so the Structural Laboratory does not pretend to provide a 3D viewer.
- Overview counts that are not exposed by an endpoint are labelled `API` rather than fabricated.
- Provenance and methods screens summarize verified repository documentation; they do not claim browser access to local audit files.
