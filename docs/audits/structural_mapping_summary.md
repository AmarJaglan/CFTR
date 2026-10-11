# CFTR-StructMap: Structural Mapping Engine Summary

**Phase**: Phase 3 Structural Mapping Engine  
**Canonical Reference**: Human CFTR UniProt [P13569](https://www.uniprot.org/uniprotkb/P13569/entry) (1,480 residues)  
**Parser Engine**: Gemmi mmCIF Parser  
**Sequence Aligner**: Biopython PairwiseAligner (free end-gaps)  

---

## 1. Executive Summary

The Phase 3 structural mapping pipeline successfully ingested, parsed, and mapped all experimental mmCIF structures.
Crucially, the engine maintains strict separation between **sequence correspondence** (construct polymer alignment to UniProt P13569) and **experimental atomic coordinate observation** (presence of resolved atoms and Cα coordinates in the 3D model).

### Key Highlights:
1. **Full-Length vs. Domain vs. Non-CFTR Discrimination**:
   - **5UAK**, **6MSM**, and **8EIQ** represent authentic full-length human CFTR cryo-EM structures (>99.8% sequence identity to P13569).
   - **1AV5** (Human PKCI-1 / HINT1) and **7JQ4** (SARS-CoV-2 3C-like main protease) were identified as non-CFTR controls (<20% identity). Their metadata is preserved in the structure inventory, but their chains are excluded from the CFTR residue mapping table.
   - Low-identity fragments (e.g. Chain R / B in 5UAK, 6MSM, 8EIQ consisting of unassigned 'X' poly-peptides) are safely excluded.
2. **ΔF508 Deletion Confirmed in 8EIQ**:
   - In **8EIQ**, canonical position 508 (Phenylalanine) is authentically deleted in the experimental construct sequence.
   - Position 508 is explicitly mapped with `mapping_status = 'NOT_COVERED'` and warning `Construct deletion: Delta-F508`.
   - Canonical position 509 (Glycine) maps to construct `label_seq_id = 508` and author `auth_seq_id = 509`, demonstrating why author numbering must never be assumed equal to canonical or label numbering.
3. **Engineered Catalytic Mutations Detected (E1371Q)**:
   - In **6MSM** and **8EIQ**, canonical position 1371 (Glutamate, catalytic base in Walker B motif) is substituted by Glutamine (`E1371Q`) to trap ATP-bound states.
   - Position 1371 is mapped with `mapping_status = 'REFERENCE_MISMATCH'`, `reference_identity_match = False`, and `has_observed_residue = True`.
4. **Expression Tags Identified**:
   - In **5UAK** and **6MSM**, residues 1481..1489 (`SNSLEVLFQ`) represent a C-terminal PreScission protease expression tag outside the canonical 1480-residue sequence.

---

## 2. Structure Inventory & Provenance

| PDB ID | Source Filename | SHA-256 Digest | Method | Resolution (Å) | Deposition Date | Structure Type | CFTR Present |
| :--- | :--- | :--- | :--- | :---: | :---: | :--- | :---: |
| **1AV5** | `1AV5.cif` | `2e952dad865c...` | X-RAY DIFFRACTION | 2.00 | 1997-09-25 | `NON_CFTR` | No |
| **5UAK** | `5UAK(1).cif` | `9df8db7f7179...` | ELECTRON MICROSCOPY | 3.87 | 2016-12-19 | `FULL_LENGTH` | Yes |
| **6MSM** | `6MSM(1).cif` | `d387fdca7bf3...` | ELECTRON MICROSCOPY | 3.20 | 2018-10-16 | `FULL_LENGTH` | Yes |
| **7JQ4** | `7jq4.cif` | `c667dc1a232a...` | X-RAY DIFFRACTION | 1.65 | 2020-08-10 | `NON_CFTR` | No |
| **8EIQ** | `8EIQ(1).cif` | `541eadc473b4...` | ELECTRON MICROSCOPY | 3.00 | 2022-09-15 | `FULL_LENGTH` | Yes |

---

## 3. Structural Chain Mapping Metrics

> [!IMPORTANT]
> **Coverage Denominator Definitions**:
> - **Canonical Coverage**: $\frac{\text{Mapped Canonical Positions}}{\text{Canonical Sequence Length (1480)}} \times 100\%$
> - **Observed Coverage**: $\frac{\text{Positions with Observed Atomic Coordinates}}{\text{Canonical Sequence Length (1480)}} \times 100\%$
> - **Cα Coverage**: $\frac{\text{Positions with Resolved Cα Atom}}{\text{Canonical Sequence Length (1480)}} \times 100\%$

| Structure | Chain (Auth) | Entity Description | Polymer Len | Is CFTR | Identity | Canonical Cov | Observed Cov | Cα Cov | Observed Res | Cα Atoms | Unresolved | Not Covered | Mismatches |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1AV5** | `A` (`A`) | PROTEIN KINASE C INTERACTING P... | 126 | No | 19.8% | 0.0% | 0.0% | 0.0% | 0 | 0 | 0 | 0 | 0 |
| **1AV5** | `B` (`B`) | PROTEIN KINASE C INTERACTING P... | 126 | No | 19.8% | 0.0% | 0.0% | 0.0% | 0 | 0 | 0 | 0 | 0 |
| **5UAK** | `A` (`A`) | Cystic fibrosis transmembrane ... | 1489 | Yes | 100.0% | 100.0% | 77.0% | 77.0% | 1139 | 1139 | 341 | 0 | 0 |
| **5UAK** | `B` (`R`) | Cystic fibrosis transmembrane ... | 19 | No | 0.0% | 0.0% | 0.0% | 0.0% | 0 | 0 | 0 | 0 | 0 |
| **6MSM** | `A` (`A`) | Cystic fibrosis transmembrane ... | 1489 | Yes | 99.9% | 100.0% | 79.8% | 79.8% | 1181 | 1181 | 299 | 0 | 1 |
| **6MSM** | `B` (`B`) | Piece of Molecule-1 | 17 | No | 0.0% | 0.0% | 0.0% | 0.0% | 0 | 0 | 0 | 0 | 0 |
| **7JQ4** | `A` (`A`) | 3C-like proteinase | 306 | No | 15.1% | 0.0% | 0.0% | 0.0% | 0 | 0 | 0 | 0 | 0 |
| **8EIQ** | `A` (`A`) | Cystic fibrosis transmembrane ... | 1479 | Yes | 99.9% | 99.9% | 78.5% | 78.5% | 1162 | 1162 | 317 | 1 | 1 |
| **8EIQ** | `B` (`B`) | Cystic fibrosis transmembrane ... | 16 | No | 0.0% | 0.0% | 0.0% | 0.0% | 0 | 0 | 0 | 0 | 0 |

---

## 4. Residue Mapping Status Precedence & Distribution

Mapping statuses are assigned according to the following strict, unambiguous precedence hierarchy:
1. **`NOT_APPLICABLE`**: Chain or entity is not CFTR (e.g. 1AV5, 7JQ4, or unassigned peptide fragments).
2. **`ALIGNMENT_AMBIGUOUS`**: Canonical residue has multiple equally scoring alignment placements.
3. **`NOT_COVERED`**: Canonical residue is missing from the construct sequence (e.g. construct deletion like ΔF508, or terminal truncation).
4. **`REFERENCE_MISMATCH`**: Canonical residue aligns to a construct residue, but amino acid identity differs (e.g. E1371Q). Observability is tracked independently via `has_observed_residue`.
5. **`MAPPED_OBSERVED`**: Canonical residue aligns, identity matches, and coordinates are experimentally observed in the 3D model.
6. **`MAPPED_UNRESOLVED`**: Canonical residue aligns, identity matches, but coordinates are absent (disordered loop or flexible terminus).

### Mapping Breakdown for CFTR Chains:

| Structure ID | Chain | Total Canonical Residues | Mapped Observed | Mapped Unresolved | Reference Mismatch | Not Covered | Total Observed Coordinates | Total Cα Atoms |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **5UAK** | `A` | 1,480 | 1139 | 341 | 0 | 0 | 1139 | 1139 |
| **6MSM** | `A` | 1,480 | 1180 | 299 | 1 | 0 | 1181 | 1181 |
| **8EIQ** | `A` | 1,480 | 1161 | 317 | 1 | 1 | 1162 | 1162 |

---

## 5. Detailed Edge-Case Investigations

### 5.1. Authentic ΔF508 Patient Deletion in 8EIQ
- **Context**: 8EIQ is a 3.0 Å cryo-EM structure of human CFTR harboring the ΔF508 deletion in complex with Trikafta (elexacaftor, tezacaftor, ivacaftor) and ATP/Mg²⁺.
- **Sequence Verification**: In the 8EIQ construct sequence, residue 507 (`I`) is immediately followed by residue 508 (`G`), omitting Phenylalanine 508.
- **Canonical Mapping**: Position 508 is marked `NOT_COVERED`, with `has_observed_residue = False`, `has_ca_atom = False`, and warning: `Construct deletion: Delta-F508 (authentic CF-causing patient deletion)`.
- **Numbering Shift**: From position 509 onwards, construct `label_seq_id` is offset by $-1$ relative to canonical numbering, whereas author `auth_seq_id` retains standard UniProt numbering (e.g. Glycine is `label_seq_id = 508`, `auth_seq_id = '509'`). This validates the design requirement that author, label, and canonical numbering remain strictly decoupled.

### 5.2. Catalytic Inactivation Mutants (E1371Q in 6MSM and 8EIQ)
- **Context**: Cryo-EM studies of ATP-bound ABC transporters frequently introduce an engineered mutation at the catalytic base in the Walker B motif (`ILLLDEPSAHL`) to trap ATP-bound conformations by preventing hydrolysis.
- **Observation**: In both 6MSM and 8EIQ, canonical Glutamate 1371 (`E1371`) is replaced by Glutamine (`Q1371`).
- **Mapping Result**: Position 1371 is recorded with `mapping_status = 'REFERENCE_MISMATCH'`, `reference_identity_match = False`, `has_observed_residue = True`, and `has_ca_atom = True`.

### 5.3. Expression Tags and Non-Canonical Extensions
- **5UAK and 6MSM**: Both construct sequences contain 1,489 amino acids. Residues 1..1480 correspond to canonical CFTR, while residues 1481..1489 (`SNSLEVLFQ`) represent a C-terminal PreScission protease cleavage tag.
- **Audit Handling**: These 9 construct residues do not align to canonical CFTR and are logged under `structural_mapping_warnings.csv` with category `EXPRESSION_TAG`.

### 5.4. Excluded Negative-Control Structures (1AV5 and 7JQ4)
- **1AV5**: Human PKCI-1 / HINT1 (126 aa), an unrelated nucleotide-binding protein family. Sequence identity to CFTR is 19.8% (random background). Classified as `NON_CFTR`.
- **7JQ4**: SARS-CoV-2 3C-like main protease (306 aa). Sequence identity to CFTR is 15.1%. Classified as `NON_CFTR`.
- Both structures are preserved in `structure_inventory.csv` for audit traceability, with zero residues injected into the CFTR mapping table.

---

## 6. Generated Audit Artifacts

- [`structure_inventory.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/structure_inventory.csv): Metadata, SHA-256 hashes, resolution, and classifications for all structures.
- [`structural_chain_mapping_summary.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/structural_chain_mapping_summary.csv): Chain-by-chain coverage, observation, and mismatch counts.
- [`residue_structure_mapping.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/residue_structure_mapping.csv): Full residue-level table (4,440 rows across 3 CFTR chains).
- [`structural_mapping_warnings.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/structural_mapping_warnings.csv): Full log of construct deletions, substitutions, and excluded chains.

