# CFTR-StructMap

**CFTR-StructMap** is a research-oriented bioinformatics application for structural-context analysis of CFTR (*Cystic Fibrosis Transmembrane Conductance Regulator*) variants. It integrates curated registry and assertion records from **CFTR2** and **ClinVar** with residue-level structural evidence from high-resolution experimental Cryo-EM models (5UAK, 6MSM, 8EIQ).

> **Scientific Boundary:** This platform provides mechanistic structural evidence observations only. It does **not** predict clinical pathogenicity, compute automated damage scores, or generate clinical consensus classifications.

---

## 1. System Requirements & Prerequisites

- **OS:** macOS (Apple Silicon arm64 / Intel) or Linux (x86_64)
- **Python:** 3.13+ with standard virtual environment (`.venv`)
- **Node.js:** 18+ and `npm` 9+
- **Biophysical Tooling:** Gemmi 0.7.5+, Biopython 1.88+, SQLAlchemy 2.0+

---

## 2. Quickstart & Local Execution

### Backend Setup (FastAPI)
```bash
# 1. Activate virtual environment
source .venv/bin/activate

# 2. Install backend dependencies
pip install -r backend/requirements.txt

# 3. Start FastAPI server (binds to http://127.0.0.1:8000)
uvicorn backend.app.main:app --host 127.0.0.1 --port 8000 --reload
```
- Interactive Swagger UI: [http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs)
- Health Check: [http://127.0.0.1:8000/health](http://127.0.0.1:8000/health)
- Database Statistics: [http://127.0.0.1:8000/stats](http://127.0.0.1:8000/stats)

### Frontend Setup (React / Vite)
```bash
# 1. Install frontend dependencies
npm install --prefix frontend

# 2. Start Vite development server
npm run dev --prefix frontend
# Alternatively from repository root:
npm run dev
```
- Local Web Interface: [http://localhost:3000](http://localhost:3000) (or `http://localhost:3001` if port 3000 is occupied).
- The development server automatically proxies `/api/*` requests to `http://127.0.0.1:8000/*`.

---

## 3. Configuration & Environment Variables

| Variable | Default Value | Description |
| :--- | :--- | :--- |
| `DATABASE_URL` | `sqlite:///data/processed/cftr_structmap.db` | SQLAlchemy database connection string. |
| `CORS_ORIGINS` | `http://localhost:3000,http://127.0.0.1:3000,http://localhost:3001,http://127.0.0.1:3001` | Allowed origins for direct API access. |
| `VITE_API_BASE_URL` | `/api` | Base path for frontend API calls (proxied locally). |

---

## 4. Scientific Data Pipeline CLI

When rebuilding the SQLite database from raw inputs, execute commands from the repository root in this order:

```bash
source .venv/bin/activate

# Step 1: Ingest and validate CFTR2 and ClinVar clinical records
python -m backend.app.cli import-all

# Step 2: Parse mmCIF structures with Gemmi and map to canonical CFTR (P13569)
python -m backend.app.cli map-structures

# Step 3: Extract residue-level features (SASA, RSA, contacts, DSSP)
python -m backend.app.cli extract-features

# Step 4: Integrate variant-level structural evidence and multi-structure metrics
python -m backend.app.cli integrate-evidence

# Step 5: Check database counts and provenance status
python -m backend.app.cli status
```

---

## 5. Automated Testing & Verification

### Running Backend Tests (Pytest)
```bash
source .venv/bin/activate
pytest backend/tests/ -v
```
*Current test suite: **92 passed** covering clinical ingestion, sequence normalization, mmCIF mapping, biophysical features, zero vs. null integrity, batch analysis (`POST /variants/batch`), and coordinate streaming endpoints.*

### Running Frontend Tests (Vitest)
```bash
npm test --prefix frontend
```
*Current test suite: **9 passed** covering typed API calls, statistics, query serialization, batch payload submission, coordinate URL generation, and the `BatchAnalysis` React UI component.*

### Production Build Verification
```bash
npm run build --prefix frontend
```
*Compiles TypeScript and bundles client assets via Vite.*

---

## 6. Functional Capabilities (Phase 7)

### 1. High-Throughput Batch Variant Analysis
- **Endpoint:** `POST /variants/batch`
- **Supported Formats:**
  - Canonical protein changes: `G551D`, `F508del`, `M470V`
  - HGVS protein notations: `p.Gly551Asp`, `p.Phe508del`
  - Engineered construct variants: `E1371Q`
  - Source IDs / Legacy names: `F508del`, `7120`
- **Input Channels:** Multiline textarea (one per line) or CSV file upload (auto-detects `variant`, `identifier`, or 1st column).
- **Batch Limits:** Up to 500 rows per request, max 2 MB CSV file upload.
- **Duplicate Policy:** Database lookups are deduplicated in memory, while each input occurrence is preserved in output with its original row number and text.
- **Resolution Statuses:**
  - `RESOLVED`: Found in CFTR2/ClinVar database with structural evidence.
  - `NOT_FOUND`: Valid protein syntax but absent from ingested database.
  - `UNSUPPORTED_FORMAT`: Unparseable input string.
  - `EMPTY`: Blank input row.
- **Export Schema:** Client-side CSV download containing 29 columns with RFC 4180 escaping. Null/unmodeled coordinates remain blank (never converted to 0.0).

### 2. Interactive 3D Molecular Structure Viewer
- **Viewer Technology:** Integrated **Mol\*** (version 5.13.1) WebGL molecular canvas inside `Structural Laboratory`.
- **Coordinate Delivery:** Served securely via `GET /structures/{structure_id}/coordinates` (`chemical/x-mmcif`) with directory traversal protection and HTTP cache headers.
- **Supported Models:**
  - `5UAK`: Apo, dephosphorylated monomer (3.87 Å)
  - `6MSM`: ATP-bound dimer with E1371Q catalytic trap (3.20 Å)
  - `8EIQ`: ATP-bound dimer carrying authentic ΔF508 patient deletion and E1371Q (3.40 Å)
- **Residue-Level Mapping:** `GET /structures/{structure_id}/residues/{canonical_position}` resolves canonical positions 1–1480 to structure-specific author chain (`auth_asym_id`) and sequence numbers (`auth_seq_id`), accounting for coordinate numbering shifts.
- **Camera & Interaction Controls:** Reset camera, zoom, pan, rotate, and 3D sequence strip navigation.

---

## 7. Scientific Interpretation Guidelines

### 1. Zero vs. Null Distinction
- `0.0 Å²` / `0.0%`: A **genuine calculated physical value** (e.g. residue G551 in 6MSM and 8EIQ is completely buried/occluded by NBD2 and bound ATP/Mg at the closed catalytic interface).
- `null` / `—`: **Missing, unmodeled, or inapplicable data** (e.g. F508 in 8EIQ is deleted from the construct and thus strictly `null`, never zero; R-domain residues 679–844 lack coordinates due to disorder).

### 2. Contact Neighborhood Semantics
- **Distance Criterion:** 8.0 Å heavy-atom sphere radius excluding the residue itself.
- **Entity Inclusions:** In monomeric apo state (5UAK), contacts comprise surrounding protein residues (19 contacts for G551). In ATP-bound dimer states (6MSM and 8EIQ), the 30 contacts include 28 protein residues plus catalytic ligands (`ATP1504` and `MG1502`).

### 3. Multi-Structure Comparability & Construct Caveats
- `ALL_THREE_OBSERVED`: Coordinates resolved in 5UAK, 6MSM, and 8EIQ.
- `CANONICAL_WT`: Residue identity matches canonical UniProt P13569 in all models.
- **G551D Context:** Cryo-EM models capture reference wild-type Glycine 551; mutant Aspartate sidechains are not modeled in silico.
- **E1371Q Trap:** Catalytic substitution introduced to disable ATP hydrolysis and enable Cryo-EM imaging.
- **ΔF508 Deletion:** Authentic patient deletion present in 8EIQ; residue 508 coordinates are absent.

### 4. Canonical Domain Architecture (P13569)
- **Lasso:** 1–79
- **TMD1:** 80–388
- **NBD1:** 389–678
- **R-domain:** 679–844 (intrinsically disordered in Cryo-EM maps)
- **TMD2:** 845–1172
- **NBD2:** 1173–1436
- **C-terminal:** 1437–1480

---

## 8. Capability Boundaries
- **Research Platform Only:** CFTR-StructMap provides biophysical and structural evidence observations. It does not predict clinical pathogenicity or generate automated consensus classifications.
- **Non-CFTR Controls:** Isolated crystal structures (`1AV5`) and non-CFTR ABC transporters (`7JQ4`) are parsed as negative controls but excluded from full-length human CFTR mapping.
