# CFTR-StructMap: Data Dictionary & Provenance Policy

## 1. Source Provenance Policy

### 1.1. Raw Source Immutability
- Raw input files located in `data/raw/` are treated as immutable, read-only scientific records.
- Ingestion pipelines open all raw data files in read-only mode (`rb` or `r`).
- Cryptographic SHA-256 hashes are calculated upon ingestion and persisted in the `data_sources` table. Any modification to a raw file triggers an audit mismatch.
- Original files are never modified, overwritten, trimmed, or re-formatted in place.

### 1.2. Tracking Identifiers & Versions
- **CFTR2**: Tracked by release version (e.g., `30 January 2026`) and legacy variant identifier (`Variant legacy name`).
- **ClinVar**: Tracked by stable NCBI `VariationID` and aggregate `Accession` (`VCV...`).
- Ingestion runs record:
  - Source filename and full path.
  - Ingestion timestamp (UTC ISO-8601).
  - Software version and execution summary metrics.

---

## 2. Scientific Data Dictionaries

### 2.1. CFTR2 Schema (`cftr2_variants`)

| Field Name | Type | Description | Scientific Notes |
| :--- | :--- | :--- | :--- |
| `id` | Integer (PK) | Internal surrogate primary key. | Autoincrementing integer. |
| `legacy_name` | String (Unique) | Historical CFTR variant designation (e.g. `F508del`, `G551D`, `M470V`). | Primary unique key in CFTR2 patient registry. |
| `protein_name` | String | Reported HGVS protein description (e.g. `p.Gly551Asp`). | Preserved verbatim from source export. |
| `cdna_name` | String | Reported HGVS cDNA description (e.g. `c.1652G>A`). | Preserved verbatim from source export. |
| `canonical_protein` | String | Standardized compact notation (e.g. `G551D`, `F508del`, `H139=`). | Normalized via `backend.app.scientific.normalization`. `NULL` for unparseable. |
| `uniprot_pos` | Integer | 1-based residue coordinate in canonical UniProt P13569. | Bounded `1..1480` for coding variants; `NULL` for non-coding/splice. |
| `ref_aa` | String (1-aa) | Reported reference amino acid. | Single-letter amino acid code. |
| `alt_aa` | String (1-aa) | Reported alternate amino acid. | Single-letter code, `NULL` for deletions, `X` for nonsense. |
| `variant_type` | String | Variant classification category. | `missense`, `inframe_deletion`, `nonsense`, `frameshift`, `synonymous`, `unknown`. |
| `is_valid_wt` | Boolean | Agreement with UniProt P13569 sequence at `uniprot_pos`. | `True` if matches, `False` if mismatch, `NULL` if unvalidated. |
| `validation_flag` | String | Machine-readable validation status. | `MATCH`, `WT_MISMATCH`, `OUT_OF_BOUNDS`, `NORMALIZATION_FAILED`. |
| `validation_message`| Text | Diagnostic explanation of validation finding. | Explicit rationale (e.g. flags known polymorphisms such as `M470V`). |
| `alleles_reported` | Integer | Total patient alleles observed in CFTR2 registry. | Cohort-level clinical abundance metric. |
| `allele_frequency` | Float | Allele frequency in CFTR2 cohort. | Expressed as fraction of total identified alleles. |
| `determination_current` | String | Current clinical consequence determination. | **Preserved verbatim**: `CF-causing`, `Varying clinical consequence`, `Non CF-causing`, `No interpretation available`. |
| `determination_previous`| String | Historical determination in previous release. | Preserved verbatim for longitudinal consistency. |
| `change_from_previous` | String | Longitudinal status change flag (`Yes` / `No`). | Direct indication of reclassification. |
| `alternative_names`| Text | Historical aliases and alternate legacy names. | Raw comma-delimited aliases. |
| `raw_data_json` | Text (JSON) | Complete verbatim row record from source file. | Audit preservation. |

---

### 2.2. ClinVar Schema (`clinvar_variants`)

| Field Name | Type | Description | Scientific Notes |
| :--- | :--- | :--- | :--- |
| `variation_id` | Integer (PK) | Stable NCBI ClinVar VariationID (e.g. `53423`). | Primary stable identifier across ClinVar releases. |
| `vcv_accession` | String | ClinVar Variation record accession (e.g. `VCV000053423`). | Aggregate VCV identifier. |
| `name` | String | Full ClinVar variant title descriptor. | e.g. `NM_000492.4(CFTR):c.1652G>A (p.Gly551Asp)`. |
| `protein_change_raw` | String | Raw protein notation from ClinVar export (e.g. `G551D`). | Preserved verbatim. |
| `canonical_protein` | String | Standardized compact notation (e.g. `G551D`). | Normalized representation. |
| `uniprot_pos` | Integer | 1-based residue position on UniProt P13569. | Bounded `1..1480`. |
| `ref_aa` | String (1-aa) | Reported reference amino acid. | Single-letter code. |
| `alt_aa` | String (1-aa) | Reported alternate amino acid. | Single-letter code. |
| `variant_type` | String | Raw ClinVar variant type. | `single nucleotide variant`, `Indel`, `Haplotype`. |
| `molecular_consequence` | String | Molecular effect descriptor. | e.g. `missense variant`, `initiator_codon_variant`. |
| `is_valid_wt` | Boolean | Strict agreement with UniProt P13569. | Validated against canonical sequence. |
| `validation_flag` | String | Machine-readable validation status. | `MATCH`, `WT_MISMATCH`, `OUT_OF_BOUNDS`. |
| `validation_message`| Text | Diagnostic explanation of validation finding. | Mismatches are never silently corrected. |
| `germline_classification` | String | Aggregate germline clinical classification. | **Preserved verbatim**: `Pathogenic`, `Likely pathogenic`, `Conflicting classifications of pathogenicity`, `Uncertain significance`, `Benign`. |
| `review_status` | String | ClinVar review confidence status. | `reviewed by expert panel`, `criteria provided, multiple submitters, no conflicts`, `criteria provided, conflicting classifications`, etc. |
| `date_last_evaluated` | String | Date of last clinical curation evaluation. | Human-readable date string from ClinVar. |
| `conditions` | Text | Associated clinical conditions / phenotypes. | Pipe-delimited list of indicated disorders. |
| `grch38_location` | String | Genomic coordinate on GRCh38 assembly. | Chromosome and position. |
| `grch37_location` | String | Genomic coordinate on GRCh37 assembly. | Historical genomic alignment coordinate. |
| `dbsnp_id` | String | dbSNP accession (e.g. `rs397508328`). | Cross-reference identifier. |
| `allele_id` | String | NCBI AlleleID. | Cross-reference identifier. |
| `raw_data_json` | Text (JSON) | Complete verbatim row record from source file. | Audit preservation. |

---

## 3. Protein-Change Normalization Rules

The normalization engine (`backend.app.scientific.normalization.normalize_protein_change`) applies strict, ordered pattern-matching without guessing or heuristic reduction:

1. **Nonsense / Stop-gain**:
   - Syntax: `p.Gly542X`, `G542X`, `p.Trp1282*`, `p.Trp1282Ter`
   - Normalization: `[Ref1][Pos]X` (e.g., `G542X`, `W1282X`)
   - Type: `nonsense` (`is_missense=False`, `alt_aa='X'`)
2. **Three-letter HGVS Missense**:
   - Syntax: `p.Gly551Asp`, `Gly551Asp`, `(p.Met1Val)`
   - Normalization: `[Ref1][Pos][Alt1]` (e.g., `G551D`, `M1V`)
   - Type: `missense` (or `synonymous` if `Ref1 == Alt1`)
3. **One-letter Missense**:
   - Syntax: `G551D`, `M1V`, `N1303K`
   - Normalization: `[Ref1][Pos][Alt1]`
   - Type: `missense`
4. **Three-letter and One-letter In-frame Deletion**:
   - Syntax: `p.Phe508del`, `F508del`
   - Normalization: `[Ref1][Pos]del` (e.g., `F508del`)
   - Type: `inframe_deletion` (`alt_aa=None`, `is_missense=False`)
5. **Frameshift Indels**:
   - Syntax: `p.Leu167fs`, `p.Phe430SerfsX12`
   - Normalization: `[Ref1][Pos]fs` (e.g., `L167fs`)
   - Type: `frameshift` (`is_missense=False`)
6. **Synonymous Coding Changes**:
   - Syntax: `p.His139=`, `p.Thr854=`, `p.H139=`
   - Normalization: `[Ref1][Pos]=` (e.g., `H139=`, `T854=`)
   - Type: `synonymous` (`ref_aa == alt_aa`, `is_missense=False`)
7. **Complex Alleles, Multi-residue Deletions, and Indels**:
   - Patterns: `S912L, G1244V`, `p.[Arg117His;7T]`, `p.Ile1023_Val1024del`, `p.Gln359_Thr360delinsLysLys`
   - Normalization: **Not normalized to single residue**; emitted with `canonical_protein = NULL`, `variant_type = "complex_unparsed"`. Preserves biological complexity without data loss.

---

## 4. Reference Sequence Validation & The M470 Polymorphism

- **Canonical Sequence Reference**: UniProt **P13569** (`CFTR_HUMAN`, 1480 amino acids, SHA-256 `a86f7dbeb9c5f48401a1a452b440d782c35067d2da152180a09ce26f76d11639`).
- **Coordinate System**: All positions use 1-based biological numbering externally (`1..1480`).
- **Validation Rule**: At position `pos`, the reported wild-type residue `ref_aa` is evaluated against `P13569[pos - 1]`.
  - If identical: `validation_flag = "MATCH"`, `is_valid_wt = True`.
  - If different: `validation_flag = "WT_MISMATCH"`, `is_valid_wt = False`. Mismatches are **never silently corrected**.
  - If outside 1..1480: `validation_flag = "OUT_OF_BOUNDS"`, `is_valid_wt = False`.
- **M470V Reference Discrepancy Findings**:
  - The canonical sequence encodes **Valine (`V`)** at position 470 (`rs1800089`, codon `GTG` / c.1408G).
  - The original 1989 cDNA clone (GenBank `M28668.1`) had `1540A` (c.1408A) encoding Methionine (`M`).
  - CFTR2 contains records under both legacy naming (`M470V`, `p.Met470Val`, `c.1408A>G`) and canonical naming (`V470M`, `p.Val470Met`, `c.1408G>A`).
  - In CFTR2 Sheet 2, the curators explicitly state: *"old nomenclature; reference sequence now has V as reference so genomic coordinates are not provided"*.
  - The system correctly flags `M470V` as `WT_MISMATCH` with an explicit diagnostic note, while `V470M` passes validation as a `MATCH`.

---

## 5. Source-Level vs. Variant-Level Accounting

| Concept | Definition | Rule |
| :--- | :--- | :--- |
| **Source Record** | An individual row or accession in a source database file. | Every source row has distinct metadata (phenotype, frequency, submitter). Must never be collapsed upon import. |
| **Unique Normalized Variant** | A distinct amino acid change (e.g., `G551D`). | Multiple source records can map to the same protein change (e.g. different cDNA mutations). |
| **Cross-Dataset Counts** | Summing records across CFTR2 and ClinVar. | **Prohibited**: CFTR2 records and ClinVar records overlap substantially. They must be reported separately. |

---

## 6. ClinVar VCV vs. SCV Granularity

- **Available Export**: `clinvar_download_20261009_143804.tsv` is an aggregate **Variation-level export (VCV)**.
- **Accessions**: Records possess `VCV...` accessions and unique `VariationID`s.
- **Review Status**: Reflects aggregate consensus tiers (`reviewed by expert panel`, `criteria provided, multiple submitters, no conflicts`, `criteria provided, conflicting classifications`).
- **Individual Submissions (SCV)**: SCV records (submitter names, individual submission dates, specific assertion methods) are **not present** in this aggregate TSV file. Submitter-level history must not be claimed unless a separate SCV export (ClinVar XML or VCF) is ingested.
- **Classification Preservation**: Classification terms (including `Conflicting classifications of pathogenicity`) are preserved verbatim without artificial consensus.

---

## 7. Known Limitations & Audit Findings

1. **Non-Coding and Splicing Records**:
   - 393 CFTR2 records represent intronic or splice-site mutations (e.g., `3849+10kbC->T`, `621+1G->T`) with protein consequence reported as `p.?`.
   - These are retained in the database with their full clinical determinations and `canonical_protein = NULL`. They cannot be mapped to 3D atomic coordinates.
2. **Complex Alleles in cis**:
   - 87 CFTR2 records and 11 ClinVar records represent complex alleles or haplotypes (e.g. `R117H;7T`, `G576A;R668C`, `F508del;I1027T`).
   - These are preserved intact without artificial reduction to single-residue substitutions.
3. **Spreadsheet Footer Rows**:
   - 5 rows in the raw CFTR2 spreadsheet represent legal disclaimers, copyright notices, and footnotes. These are explicitly excluded from variant ingestion, yielding exactly 2,092 genuine biological records.
