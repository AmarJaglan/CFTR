# CFTR-StructMap: Normalization Audit & Ingestion Integrity Summary

**Phase**: 2.1 Data Quality Audit & Ingestion Integrity  
**Date**: October 2026  
**Auditor**: Lead Software Engineer  
**Scope**: Ingestion modules, scientific reference validation, protein normalization, SQLite database integrity, idempotency, and source provenance.

---

## 1. Executive Summary

A comprehensive, rigorous audit of the Phase 2 data ingestion pipeline was performed across both primary datasets:
- **CFTR2**: `data/raw/CFTR2_30January2026.xlsx` (SHA-256: `7b3c13b20fd4243ea707413c43653099698e672421dbd15ee9c1f5c56f851ad9`)
- **ClinVar**: `data/raw/clinvar_download_20261009_143804.tsv` (SHA-256: `a5291e5d9bab3e957dc9f51be035d3d0305893525ad6d7de9d49c27874353132`)

### Key Findings:
1. **M470V Reference Discrepancy Resolved**:
   - The reported wild-type mismatch at position 470 (`M470V`, `p.Met470Val`, `c.1408A>G`) is biologically and historically authentic.
   - Canonical human UniProt [P13569](https://www.uniprot.org/uniprotkb/P13569/entry) and NCBI RefSeq NM_000492 / NP_000483 standardize on **Valine (`V`)** at residue 470 (`rs1800089`).
   - The original 1989 cloned CFTR cDNA (GenBank `M28668.1`) contained `1540A` (c.1408A) encoding Methionine (`Met`). In CFTR2 Sheet 2, the curators explicitly annotate `M470V`: *"old nomenclature; reference sequence now has V as reference so genomic coordinates are not provided"*.
   - The validator correctly flagged this record with `WT_MISMATCH`. The mismatch was preserved and contextualized, never silently masked.
2. **CFTR2 Spreadsheet Footer Anomaly Uncovered**:
   - The raw CFTR2 sheet contains 2,097 rows below the column headers. The last 5 rows are spreadsheet footnotes, copyright disclaimers, and citation notices (`*Represents...`, `Permitted use...`, `Please use...`, `©Copyright...`, `NOTE:...`).
   - The initial pipeline ingested these 5 rows as variant records with unparseable protein changes.
   - Excluding these non-variant rows brings the valid record count to **exactly 2,092**, aligning perfectly with the CFTR2 document's stated header: *"Number of variants reported in CFTR2: 2,092"*.
3. **Synonymous Notation Normalization (`p.Xxx123=`)**:
   - 40 CFTR2 records used standard HGVS synonymous notation (e.g., `p.His139=`, `p.Glu528=`, `p.Thr854=`).
   - All 40 were verified against UniProt P13569 with 100% agreement. Support was implemented for standard HGVS synonymous notation (`H139=`), without altering the original source text.
4. **Authentic Normalization Failures**:
   - **CFTR2**: 560 real records cannot and should not be collapsed into single-residue substitutions (393 splicing/intronic `p.?`, 87 complex alleles/haplotypes, 40 duplications, 13 multi-residue deletions, 12 complex delins, 6 unassigned `p.?`, 5 complex frameshifts, 3 UTR, 1 insertion).
   - **ClinVar**: Exactly 13 records fail normalization: 11 are multi-variant haplotypes (`Type: Haplotype`) and 2 are multi-residue deletion-insertions (`p.Gln359_Thr360delinsLysLys`, `p.Met1101_Arg1102delinsIleGly`).
5. **Database Idempotency and Transactional Integrity**:
   - Both ingestion modules were upgraded to support idempotent upserts by stable identifier (`legacy_name` for CFTR2, `variation_id` for ClinVar). Re-running imports does not duplicate rows or trigger constraint crashes. All writes are guarded by transaction rollbacks on failure.

---

## 2. Source-Record Counts vs. Unique Normalized Protein Changes

> [!IMPORTANT]
> **Source-record counts** must never be confused with **unique protein-variant counts**. Multiple distinct nucleotide alterations in a gene can encode the identical amino acid substitution. Furthermore, records from CFTR2 and ClinVar must never be summed directly as "unique variants".

| Metric | CFTR2 (`30 January 2026`) | ClinVar (`20261009_143804`) | Combined Interpretation |
| :--- | :--- | :--- | :--- |
| **Source Records Read** | 2,097 (includes 5 footer rows) | 2,664 | Raw input lines parsed |
| **Valid Source Records Persisted** | **2,092** | **2,664** | Genuine variant records in database |
| **Unique Stable Identifiers** | 2,092 (`legacy_name`) | 2,664 (`variation_id`) | 100% unique primary identifiers |
| **Parseable Protein Changes** | 1,532 | 2,651 | Normalized to compact protein format |
| **Unique Normalized Protein Changes** | **1,437** | **2,571** | Distinct protein alterations |
| **Duplicate Protein Change Labels** | 82 distinct labels | 75 distinct labels | Multiple nucleotide variants yielding same protein change |
| **Wild-type Reference Mismatches** | 1 (`M470V`) | 0 | Flagged against canonical UniProt P13569 |
| **Out-of-Range Positions** | 0 | 0 | All positions within 1..1480 |
| **Normalization Failures (Persisted with `canonical_protein = NULL`)** | **560** | **13** | Genuine non-coding, complex, or multi-residue records |

---

## 3. Detailed Audit: The M470V Discrepancy

### 3.1. Evidence Synthesis
- **CFTR2 Record**:
  - Legacy name: `M470V`
  - Protein notation: `p.Met470Val`
  - cDNA notation: `c.1408A>G`
  - Alternative names: `1540A>G`
  - Determination: `Non CF-causing`
- **Canonical Reference**:
  - UniProt: `P13569` (Entry `CFTR_HUMAN`, 1480 aa)
  - RefSeq: `NM_000492.4` / `NP_000483.3`
  - Sequence at position 470: **Valine (`V`)**, encoded by codon `GTG` (c.1408G).
  - dbSNP: `rs1800089` (GRCh38 chr7:117559479 G>A).
- **ClinVar Record**:
  - VariationID: `7130` (`VCV000007130`)
  - Title: `NM_000492.4(CFTR):c.1408G>A (p.Val470Met)`
  - Classification: `Benign/Likely benign`
- **CFTR2 Sheet 2 Curator Note**:
  - `c.1408A>G` (`M470V`): *"old nomenclature; reference sequence now has V as reference so genomic coordinates are not provided"*
  - `c.1408G>A` (`V470M`): `grch38_pos: 117559479`, `ref: G`, `alt: A`.

### 3.2. Evaluation Conclusion
When CFTR was originally cloned in 1989 (Riordan et al., *Science* 1989), the sequenced cDNA (GenBank `M28668`) had an adenine at position 1540 (`1540A`, corresponding to `c.1408A`), specifying Methionine. Subsequent population genetics demonstrated that Valine (`c.1408G`) is the major/ancestral allele. Modern reference genomes, RefSeq, and UniProt standardized on Val470.

The validator's behavior is scientifically sound: `M470V` reports reference Met, whereas the canonical reference is Val. The system correctly records `validation_flag = "WT_MISMATCH"` and contextualizes the finding without suppressing the warning or altering the source data.

---

## 4. Categorization of Normalization Failures

### 4.1. CFTR2 Normalization Failures (560 Records)

Full inventory recorded in [`docs/audits/cftr2_normalization_audit.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/cftr2_normalization_audit.csv).

| Category | Count | Representative Examples | Description & Recommended Action |
| :--- | :--- | :--- | :--- |
| **Splicing / Intronic (`p.?`)** | 393 | `3849+10kbC->T`, `621+1G->T`, `2789+5G->A`, `1717-1G->A` | Non-coding or splice-junction variants with no direct coding residue change. Retained with `canonical_protein = NULL` and full clinical determination. |
| **Complex Allele / Haplotype** | 87 | `R117H;7T`, `G576A;R668C`, `F508del;I1027T`, `M470V;T854T;Q1463Q` | Co-occurring mutations reported in *cis*. Must not be collapsed into single amino acid changes; downstream structural mapping can decompose individual elements. |
| **In-frame Duplication** | 40 | `p.Leu138dup` (`c.413_415dup`), `p.Val456dup`, `p.Phe337_Ile340dup` | Duplication of wild-type residues. Cannot be represented as single-residue substitutions. Retained as in-frame duplication indels. |
| **Multi-residue In-frame Deletion** | 13 | `p.Ile1023_Val1024del` (`3199del6`), `p.Met607_Gln634del` (`1949del84`) | Deletions spanning multiple residues. Retained as multi-residue deletions with start/end residue bounds. |
| **Complex Delins** | 12 | `p.Gln359_Thr360delinsLysLys`, `p.His939_Leu967delinsGln` | Simultaneous deletion and insertion across multiple codons. Retained as complex delins. |
| **Missing / Unresolved Protein Change** | 6 | `M1V` (`c.1A>G`, `p.?`), `M1T` (`c.2T>C`, `p.?`), `TG10` (`polyT/TG`) | Coding cDNA reported in CFTR2 but protein consequence unassigned in export. Retained without guessing. |
| **Complex Frameshift Allele** | 5 | `p.[Tyr577Phe;Leu581PhefsX8]`, `p.[Leu365TrpfsX16;Ser1235Arg]` | Combinations of frameshifts with missense variants in cis. Retained as complex frameshift alleles. |
| **UTR / Regulatory Variant** | 3 | `c.-887C>T` (`5UTR-816C->T`), `polyT/TG` | Untranslated region or repeat variants. Retained as non-coding regulatory records. |
| **In-frame Insertion** | 1 | `p.Val1272_Ser1273insTrp` (`c.3815_3816insTTG`) | Insertion of novel amino acids between residues 1272 and 1273. Retained as in-frame insertion. |

#### 4.1.1. Cross-Tabulation & Reconciliation: Source Protein Notation (`p.?`) vs. Mutation Mechanism

> [!NOTE]
> **Reconciliation of 482 `p.?` records with failure categories**:
> In the raw CFTR2 source file, exactly **482 records** report unknown or missing protein consequence (`protein_name` in `p.?`, `?`, or empty). The remaining **78 unparseable records** provide explicit HGVS protein notation that cannot be represented as a single-residue missense substitution.
>
> The failure categories classify variants by their **biological and molecular mutation mechanism** (determined from cDNA coordinates and syntax). Consequently, the 482 `p.?` records are partitioned across six mechanisms, with zero category overlap:

| Biological / Molecular Category | Explicit Protein String (Non-`p.?`) | Missing / Unknown Notation (`p.?`) | Total Unparseable Records | cDNA / Structural Definition |
| :--- | :---: | :---: | :---: | :--- |
| **Splicing / Intronic** | 0 | **393** | **393** | cDNA coordinates with `+` or `-` (intronic splice junctions) |
| **Complex Allele / Haplotype** | 48 | **39** | **87** | Co-occurring mutations in *cis* (`;` in cDNA or protein notation) |
| **In-frame Duplication** | 6 | **34** | **40** | Duplications of one or more residues (`dup` in cDNA or protein) |
| **Multi-residue In-frame Deletion** | 13 | **0** | **13** | In-frame deletions spanning $>1$ codon (`_...del`) |
| **Complex Delins** | 5 | **7** | **12** | Deletion-insertions spanning $>1$ codon (`delins` in cDNA or protein) |
| **Missing / Unresolved Protein Change** | 0 | **6** | **6** | Single-nucleotide coding cDNAs (`c.1A>G`, `c.2T>C`) with unassigned protein |
| **Complex Frameshift Allele** | 5 | **0** | **5** | Multi-variant alleles combining frameshifts with missense changes |
| **UTR / Regulatory Variant** | 0 | **3** | **3** | Non-coding 5' UTR or repeat alterations (`c.-...`) |
| **In-frame Insertion** | 1 | **0** | **1** | Multi-residue insertions (`ins` in cDNA or protein) |
| **Total** | **78** | **482** | **560** | **100% of CFTR2 normalization failures accounted for** |

- **Reconciliation Proof**: `393 (Splicing) + 39 (Complex) + 34 (Duplication) + 7 (Delins) + 6 (Missing/Unresolved) + 3 (UTR) = 482 (p.?)`.
- **Complete Total**: `482 (p.?) + 78 (Non-p.?) = 560 unparseable records`.
- **Database Alignment**: Exactly matches `SELECT count(*) FROM cftr2_variants WHERE canonical_protein IS NULL;` (560).

### 4.2. ClinVar Normalization Failures (13 Records)

Full inventory recorded in [`docs/audits/clinvar_normalization_audit.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/clinvar_normalization_audit.csv).

| VariationID | Accession | Variant Name | Protein Change | Variant Type | Molecular Consequence | Explanation & Action |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **7239** | `VCV000007239` | `c.[2735C>T;3731G>T]` | `S912L, G1244V` | Haplotype | missense variant | Multi-site haplotype in *cis*. Preserved without collapsing. |
| **40200** | `VCV000040200` | `c.[1075C>A;1079C>A]` | `Q359K, T360K` | Haplotype | missense variant | Multi-site haplotype in *cis*. Preserved without collapsing. |
| **53180** | `VCV000053180` | `c.1075_1079delinsAAAAA (p.Gln359_Thr360delinsLysLys)` | `None` | Indel | missense variant | Two-codon deletion-insertion. Preserved as delins. |
| **634835** | `VCV000634835` | `c.[1397C>G;3209G>A]` | `S466*, R1070Q` | Haplotype | nonsense\|missense | Multi-site nonsense/missense haplotype. Preserved. |
| **634837** | `VCV000634837` | `c.[1521_1523delCTT;3080T>C]` | `None` | Haplotype | inframe_deletion\|missense | Deletion combined with missense variant in *cis*. Preserved. |
| **870319** | `VCV000870319` | `c.[1327G>T;1727G>C;2002C>T]` | `G576A, R668C, D443Y` | Haplotype | missense variant | Triple missense haplotype in *cis*. Preserved. |
| **916697** | `VCV000916697` | `c.[1727G>C;2002C>T]` | `G576A, R668C` | Haplotype | missense variant | Double missense haplotype in *cis*. Preserved. |
| **977735** | `VCV000977735` | `c.[220C>T;3808G>A]` | `D1270N, R74W` | Haplotype | missense variant | Double missense haplotype in *cis*. Preserved. |
| **996223** | `VCV000996223` | `c.[2816A>G;2846A>T]` | `H939R, H949L` | Haplotype | missense variant | Double missense haplotype in *cis*. Preserved. |
| **1300167** | `VCV001300167` | `c.[1523T>G;3752G>A]` | `F508C, S1251N` | Haplotype | missense variant | Double missense haplotype in *cis*. Preserved. |
| **1300168** | `VCV001300168` | `c.[3846G>A;3848G>T]` | `W1282*, R1283M` | Haplotype | nonsense\|missense | Nonsense plus missense haplotype in *cis*. Preserved. |
| **1320029** | `VCV001320029` | `c.[4C>T;7A>T]` | `Q2*, R3W` | Haplotype | nonsense\|missense | Nonsense plus missense haplotype in *cis*. Preserved. |
| **2032842** | `VCV002032842` | `c.3303_3304delinsAG (p.Met1101_Arg1102delinsIleGly)` | `None` | Indel | missense variant | Two-codon deletion-insertion. Preserved as delins. |

---

## 5. Verification of Database Integrity and Provenance

1. **Idempotent Upsert Verification**:
   - Re-running `import-all` against an existing database was verified to execute with zero `IntegrityError` exceptions and zero unintended row additions.
   - The database contains **2,092** CFTR2 records and **2,664** ClinVar records across repeated executions.
2. **Stable Identifier Preservation**:
   - In ClinVar, primary keys are strictly anchored to NCBI `VariationID`. Multiple variants encoding the same protein change (e.g. `VariationID 53869` and `VariationID 53870`, both `p.Met1Ile`) remain distinct rows.
   - In CFTR2, records are uniquely keyed by `legacy_name`.
3. **Cryptographic Provenance**:
   - `data_sources` preserves SHA-256 fingerprints calculated from raw binary reads of `data/raw/` files:
     - `CFTR2_30January2026.xlsx`: `7b3c13b20fd4243ea707413c43653099698e672421dbd15ee9c1f5c56f851ad9`
     - `clinvar_download_20261009_143804.tsv`: `a5291e5d9bab3e957dc9f51be035d3d0305893525ad6d7de9d49c27874353132`
4. **Verbatim Clinical Assertions**:
   - ClinVar germline classifications (including `Conflicting classifications of pathogenicity`) and review tiers (e.g. `reviewed by expert panel`) are preserved verbatim.
   - No majority-voting or manufactured consensus is performed.
