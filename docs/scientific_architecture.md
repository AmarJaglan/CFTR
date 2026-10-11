# CFTR-StructMap: Scientific Architecture Specification

## 1. System Mission & Boundary Conditions

CFTR-StructMap is a specialized bioinformatics platform designed to provide reproducible, structural-context analysis for variants in the Human Cystic Fibrosis Transmembrane Conductance Regulator (*CFTR* / *ABCC7*, UniProt: [P13569](https://www.uniprot.org/uniprotkb/P13569/entry)).

### Core Principles & Boundaries
1. **Structural-Context Evidence, Not Pathogenicity Prediction**: The platform quantifies objective biophysical, topological, and microenvironmental context (e.g. solvent exposure, secondary structure, proximity to pore or nucleotide-binding pockets, contact disruption). It **does not** emit a clinical pathogenicity classification or pseudo-diagnostic score.
2. **Provenance & Source Traceability**: Raw data files are preserved verbatim in `data/raw/` with cryptographic SHA-256 hashes. Ingested records retain original accessions, submitter evidence, and version stamps.
3. **Strict Coordinate Truth**: Sequence positions are explicitly distinguished from experimentally observed atomic coordinates. If a residue resides in an unresolved region (e.g. flexible loops, disordered parts of the R-domain), calculations are not fabricated or filled with zeros; they are flagged as uncomputable due to missing coordinates.
4. **Context-Specific Structure Treatment**: Experimental structures represent specific biological states. For example, **PDB 8EIQ** represents CFTR in a ΔF508 mutant, triple-combination modulator-bound context (ivacaftor/tezacaftor/elexacaftor or equivalent therapeutic modulators), not an ordinary wild-type reference.

---

## 2. Modular Architecture Overview

```text
backend/app/
├── core/                       # Core configuration, logging, exceptions
│   ├── config.py
│   └── logging.py
├── domain/                     # Domain schemas & value objects (Pydantic)
│   ├── variants.py             # HGVS, sequence positions, variant classes
│   ├── structures.py           # PDB/mmCIF representations, chains, conformations
│   ├── features.py             # Structural feature representations & metrics
│   └── provenance.py           # Source and run metadata
├── ingestion/                  # External data ingestion & parsing (read-only raw)
│   ├── clinvar.py              # ClinVar VCV vs SCV distinction & parsing
│   ├── cftr2.py                # CFTR2 clinical & functional ingestion
│   └── uniprot.py              # Canonical P13569 sequence & domain definitions
├── scientific/                 # Bio-computational algorithms & parsers
│   ├── normalization.py        # Variant normalization & WT residue validation
│   ├── mapping.py              # UniProt <-> PDB residue mapping & SIFTS alignment
│   ├── mmcif.py                # Gemmi-based mmCIF parsing & model inspection
│   └── features/               # Structural-context calculators
│       ├── environment.py      # Spatial neighbor search, contact networks
│       ├── accessibility.py    # Surface accessibility & burial
│       └── functional_sites.py # NBD/Walker motifs, channel pore, modulator pockets
├── db/                         # Persistence layer (SQLAlchemy + SQLite)
│   ├── session.py              # Engine and session lifecycle
│   ├── models.py               # Relational tables with foreign-key provenance
│   └── repository.py           # Data access objects
└── api/                        # FastAPI presentation layer
    └── v1/
        ├── router.py           # Main v1 router assembly
        ├── endpoints/
        │   ├── variants.py     # Variant search and retrieval
        │   ├── structures.py   # Structure metadata & state inspection
        │   ├── mapping.py      # Sequence-to-structure residue mapping
        │   ├── features.py     # Structural feature queries
        │   └── provenance.py   # Data versioning & execution lineage
```

---

## 3. Detailed Component Specifications

### 3.1. Variant Normalization and Validation (`domain/variants.py`, `scientific/normalization.py`)
- **HGVS Parsing**: Support HGVS cDNA (`c.`) and protein (`p.`) nomenclature, as well as legacy CFTR nomenclature (e.g. `ΔF508` / `F508del`, `G551D`, `W1282X`, `N1303K`).
- **Canonical Sequence Reference**: UniProt P13569 (1480 amino acids).
- **Validation**:
  - Enforce that the reported wild-type amino acid matches UniProt P13569 at the specified 1-based index (e.g. index 508 must be Phenylalanine `F`).
  - Flag non-canonical or ambiguous coordinates immediately prior to any calculation.
  - Categorize variant type: Missense, In-frame deletion, In-frame insertion, Frameshift, Nonsense, Synonymous.

### 3.2. ClinVar and CFTR2 Source Ingestion (`ingestion/clinvar.py`, `ingestion/cftr2.py`)
- **ClinVar Model**:
  - Distinguish **Aggregate Variation records (VCV)** from **Individual Submissions (SCV)**.
  - Track submitter name, assertion criteria, review status (0 to 4 stars), clinical significance, and evaluation date.
  - Do not summarize conflicting assertions into a single synthetic certainty.
- **CFTR2 Model**:
  - Ingest clinical consequence categories (*CF-causing*, *Varying clinical consequence*, *Non CF-causing*, *Unknown significance*).
  - Ingest aggregate phenotypic metrics: sweat chloride distributions, pancreatic insufficiency percentages, and functional channel activity assays.
- **Raw File Safeguards**:
  - Raw source files in `data/raw/` are treated as immutable.
  - The ingestion engine records SHA-256 hashes, source URLs, and download timestamps in SQLite.

### 3.3. UniProt and PDB Residue Mapping (`scientific/mapping.py`)
- **Coordinate Systems**:
  - Canonical UniProt position (`1..1480`).
  - mmCIF Author residue number (`auth_seq_id`) and chain (`auth_asym_id`).
  - mmCIF Canonical residue number (`label_seq_id`) and chain (`label_asym_id`).
- **Engineered Mutations & Discrepancies**:
  - Experimental structures frequently harbor mutations (e.g., E1371Q in ATP-hydrolysis-deficient cryo-EM models, thermostabilizing mutations, purification tags).
  - Residue mapping explicitly compares the experimental residue (`auth_comp_id`) to the UniProt reference and sets `is_engineered_mutation = True` if they diverge.
- **Observability / Coordinate Presence**:
  - Flag whether coordinates exist in the model for that residue (`is_resolved = True / False`).

### 3.4. mmCIF Structure Parsing with Gemmi (`scientific/mmcif.py`)
- **Gemmi Library**:
  - Utilizes Gemmi's mmCIF parser to access `_atom_site`, `_chem_comp`, `_struct_conn`, and header metadata directly without lossy PDB-format conversions.
- **Structural Conformation Annotation**:
  - `5UAK`: Inward-facing, dephosphorylated, ATP-free (wild-type reference).
  - `6MS4`: Outward-facing, phosphorylated, ATP-bound (wild-type reference with catalytic mutation E1371Q).
  - `8EIQ`: ΔF508 mutant in complex with corrector/potentiator modulators (therapeutic-stabilized state; **explicitly flagged as mutant/modulator context**).
  - Header extraction: Resolution, cryo-EM map details, deposition date, ligand entities (ATP, Mg²⁺, modulators).

### 3.5. Structural Feature Calculations (`scientific/features/`)
- **Guard Rail**: A feature is calculated **if and only if** the corresponding residue has resolved coordinates in the selected structure model.
- **Calculated Features**:
  1. *Local Microenvironment*: Number of contacting residues within 4.0 Å, 6.0 Å, and 8.0 Å spheres; local atom density.
  2. *Interaction Network*: Salt bridges (e.g. Arg/Lys to Asp/Glu within 4.0 Å), hydrogen bonds, aromatic interactions.
  3. *Solvent Exposure*: Relative solvent accessibility (buried vs surface-exposed).
  4. *Functional Domain & Site Proximity*:
     - Domain: Lasso, MSD1, NBD1, Regulatory (R) domain, MSD2, NBD2.
     - Proximity (Å) to Walker A, Walker B, ABC signature motif, Cl⁻ permeation pore, and drug-binding cavities.
  5. *Status Reporting*: If unresolved, output `feature_status: "UNRESOLVED_COORDINATES"`, with `null` metrics and diagnostic explanation.

### 3.6. Evidence Provenance and Data-Version Tracking (`domain/provenance.py`, `db/models.py`)
- Every database record and API response includes a provenance payload:
  - Database schema version.
  - Ingestion timestamp and raw source hash.
  - Calculation engine version and parameter set.
  - Git commit SHA of the running software.

### 3.7. SQLite Persistence (`db/models.py`, `db/repository.py`)
- Schema design:
  - `sources`: Source identifier, release tag, SHA-256 hash, fetch date.
  - `variants`: Normalized HGVS, UniProt position, WT/ALT amino acids.
  - `clinvar_assertions`: SCV level records, submitter, review stars, clinical classification.
  - `cftr2_records`: Phenotypic and functional annotations.
  - `structures`: PDB ID, resolution, state, mutant/modulator context flags.
  - `residue_mappings`: Exact per-structure residue coordinates and resolved flags.
  - `structural_features`: Precomputed or cached structural feature metrics.
  - `provenance_runs`: Execution logs and algorithm versions.

### 3.8. API Routes (`api/v1/`)
- `GET /api/v1/variants`: Search by HGVS or legacy name (e.g. `F508del`).
- `GET /api/v1/variants/{id}`: Detailed variant record with ClinVar & CFTR2 annotations.
- `GET /api/v1/structures`: Available PDB models with conformational and mutant context annotations.
- `GET /api/v1/structures/{pdb_id}/residues/{pos}`: Residue-level coordinate and resolution status.
- `GET /api/v1/features/{variant_id}`: Comprehensive structural context across all structures.
- `GET /api/v1/provenance`: Data source releases, software commit, and audit trails.

### 3.9. Automated Scientific and Software Tests (`backend/tests/`)
- Unit and integration tests:
  - Validating WT residue matching (e.g. F508 == Phe; G551 == Gly).
  - Checking that missing coordinates produce explicit unavailable responses rather than default zeros.
  - Verifying 8EIQ is properly tagged as ΔF508/modulator context and not wild-type.
  - Ensuring ClinVar SCV submitter assertions are not collapsed into fabricated consensuses.
