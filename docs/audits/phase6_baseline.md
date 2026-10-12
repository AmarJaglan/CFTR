# CFTR-StructMap — Phase 6 Baseline Integrity Audit

**Date:** 2026-10-12
**Auditor:** Senior Bioinformatics Software Engineer & QA
**Repository Working Directory:** `/Users/amar/CFTR-StructMap`

---

## 1. Actual Git Branch and Commit

- **Branch:** `main` (tracked to `origin/main`)
- **Current HEAD Commit:** `484ea0f` (*Merge pull request #3 from AmarJaglan/v0/fix-vercel-output-directory*)
- **Recent Commit Log:**
  - `484ea0f` Merge pull request #3 from AmarJaglan/v0/fix-vercel-output-directory
  - `644c24e` Configure Vercel output directory for Vite
  - `c88d5a4` Merge pull request #2 from AmarJaglan/v0/fix-vercel-build
  - `2b04575` Fix pnpm workspace build configuration
  - `8b1c19f` Merge pull request #1 from AmarJaglan/v0/cftr-frontend

---

## 2. Initial Working-Tree Status

Working tree contains pre-existing uncommitted changes implementing the integration of frontend and backend:
- `modified: backend/app/main.py`: Added search query parameter `search: str | None` supporting text and position queries on `GET /variants`.
- `modified: frontend/src/lib/api.ts`: Typed data contracts (`VariantRecord`, `VariantList`, `StructureEvidence`, `CFTR2Record`, `ClinVarRecord`, `ComparativeMetrics`, `VariantDetail`).
- `modified: frontend/src/main.tsx`: Complete single-page application shell with all 7 navigation panels, light theme, dark navy sidebar, teal accents.
- `modified: frontend/vite.config.ts` & `frontend/vite.config.js`: Host `127.0.0.1`, port `3000`, `/api` reverse proxy to `http://127.0.0.1:8000`.
- `modified: package.json`: Standardized npm scripts to `--prefix frontend`.
- `untracked: frontend/node_modules/`, `node_modules/`, build info caches.

All existing work is preserved intact without resetting or discarding changes.

---

## 3. Discovered Components

- **Backend Framework:** FastAPI application in `backend/app/main.py`, configuration in `backend/app/config.py`, SQLite session in `backend/app/db/session.py`.
- **Database Layer:** SQLAlchemy models in `backend/app/db/models.py`.
- **Scientific Pipelines:**
  - `backend/app/scientific/reference.py`: Canonical CFTR UniProt P13569 sequence (1,480 aa) and 7 structural domains.
  - `backend/app/scientific/normalization.py`: HGVS cDNA and protein normalization.
  - `backend/app/scientific/mapping.py` / `structural_mapping.py`: Gemmi-based mmCIF residue mapping.
  - `backend/app/scientific/features.py`: Shrake-Rupley SASA, Tien (2013) RSA, 8.0 Å heavy-atom contacts, mmCIF secondary structure.
- **Frontend Stack:** React 19, TypeScript, Vite 8.3, Tailwind CSS v4, Lucide icons, Vitest test runner.
- **Auditing Suite:** 19 baseline audit CSVs and markdown reports in `docs/audits/`.

---

## 4. Scientific Inputs Availability

| Resource | Status | Path & Details |
| :--- | :--- | :--- |
| **5UAK Structure** | Available | `data/raw/5UAK(1).cif` -> `5UAK.cif` (977 KB, Cryo-EM apo) |
| **6MSM Structure** | Available | `data/raw/6MSM(1).cif` -> `6MSM.cif` (1.13 MB, Cryo-EM ATP dimer, E1371Q) |
| **8EIQ Structure** | Available | `data/raw/8EIQ(1).cif` -> `8EIQ.cif` (1.12 MB, Cryo-EM ΔF508 ATP dimer, E1371Q) |
| **1AV5 (Control)** | Available | `data/raw/1AV5.cif` (248 KB, non-CFTR PKC-interacting protein) |
| **7JQ4 (Control)** | Available | `data/raw/7jq4.cif` (552 KB, non-CFTR 3C-like protease) |
| **CFTR2 Dataset** | Available | `data/raw/CFTR2_30January2026.xlsx` (353 KB) |
| **ClinVar Dataset** | Available | `data/raw/clinvar_download_20261009_143804.tsv` (831 KB) |
| **Processed SQLite DB** | Available | `data/processed/cftr_structmap.db` (13.6 MB) |

---

## 5. Database Status and Verified Record Counts

The SQLite database at `data/processed/cftr_structmap.db` is verified and populated:

- `cftr2_variants`: **2,092** rows
- `clinvar_variants`: **2,664** rows
- `structures`: **5** rows (3 CFTR models, 2 non-CFTR controls)
- `structural_chains`: **9** rows
- `residue_mappings`: **4,440** rows (1,480 positions × 3 CFTR structures)
- `residue_features`: **4,440** rows
- `multi_structure_comparisons`: **1,480** rows (1,480 canonical residues)
- `variant_structural_evidence`: **4,756** rows
- `variant_structural_summaries`: **3,353** rows (3,293 unique protein changes + 60 unresolved)

---

## 6. Current Test and Build Status

- **Frontend Vitest (`npm test --prefix frontend`):** **2 passed** (100% pass rate in 0.63s).
- **Frontend Production Build (`npm run build --prefix frontend`):** **Passed** in 356ms (`dist/index.html`, `dist/assets/index-*.css`, `dist/assets/index-*.js`).
- **Backend Pytest (`pytest backend/tests/`):** **77 passed**, 1 warning in 57.29s.

---

## 7. Discovered Problems and Severity

1. **[Medium] Contact Count Semantics (Protein vs. Ligand):**
   - The current contact counting in `compute_structure_biopython_features` uses `model.get_atoms()`, which includes hetero-ligands (`ATP1504`, `MG1502`). For residue 551 in 6MSM and 8EIQ, the 30 contacts include 28 protein residues and 2 ligand entities. The UI and documentation should explicitly state this or partition protein vs ligand neighbors.
2. **[Medium] Comparability Concept Overloading:**
   - `is_fully_comparable` is currently a boolean checking only local residue identity match and coordinate presence. It ignores whether the surrounding structure has an engineered catalytic mutation (e.g. E1371Q, which directly contacts G551) or a disease deletion ($\Delta$F508 in 8EIQ). A multi-faceted comparability model (`all_structures_observed`, `canonical_identity_match`, `construct_context_caveat`, `comparison_caveats`) is needed.
3. **[Low] Missing Aggregate Endpoint (`GET /stats`):**
   - The Overview panel presently displays live total evidence records from `GET /variants` (4,756), but must label grouped counts as unavailable because no `GET /stats` route exists on the backend. Adding a small read-only `/stats` endpoint backed by actual DB queries will provide verified counts for grouped protein changes without hardcoding.
4. **[Low] UI Typography and Special Character Formatting:**
   - Verify proper rendering of mathematical and chemical entities ($\Delta\text{RSA}$, $\Delta\text{Contacts}$, $\text{\AA}^2$, $\text{Mg}^{2+}$, $\Delta\text{F508}$, $\neq$) without LaTeX syntax bleed or broken symbols.
