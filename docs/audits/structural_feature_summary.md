# CFTR-StructMap: Structural Feature Extraction Summary

**Phase**: Phase 4 Structural Feature Extraction  
**Canonical Sequence**: UniProt [P13569](https://www.uniprot.org/uniprotkb/P13569/entry) (1,480 amino acids)  
**Target Structures**: 5UAK (Inward-facing apo), 6MSM (Outward-facing ATP-bound), 8EIQ (Trikafta-bound, Delta-F508)  
**Methods**: Shrake-Rupley SASA, Tien et al. (2013) Max SASA scale, KDTree Heavy-Atom Contacts (8.0 Å), Deposited mmCIF Secondary Structure  

---

## 1. Executive Summary

Phase 4 extracted quantitative, residue-level biophysical descriptors for all validated CFTR structures without predicting pathogenicity or generating machine-learning clinical scores. All features are anchored to canonical UniProt P13569 coordinates via Phase 3 mapping correspondence.

### Key Principles & Safeguards:
1. **Strict Coordinate Eligibility**: No structural features are calculated for unresolved residues or construct deletions. Missing values are recorded as `null` with explicit ineligibility reasons (`UNRESOLVED_COORDINATES`, `CONSTRUCT_DELETION`), never substituted by zero.
2. **Relative Solvent Accessibility (RSA)**: Raw SASA (in Å²) is normalized against the empirical theoretical maximum accessibility scale of **Tien et al. 2013** (Gly-X-Gly tripeptide model), distinguishing buried (< 20%), intermediate (20–50%), and surface (≥ 50%) residues.
3. **Local Residue Contacts**: Unique neighboring residues within 8.0 Å heavy-atom distance (excluding same-residue atoms) are enumerated, partitioned into intra-chain and inter-chain contacts.
4. **Deposited Secondary Structure**: Helices and sheets are extracted from official mmCIF deposited records (`_struct_conf`, `_struct_sheet_range`), while unannotated observed residues are assigned to coil with provenance method `Default Assignment (Unannotated in mmCIF)`. Unannotated residues are never described as explicitly deposited coil.
5. **Delta-F508 and E1371Q Preservation**: The authentic $\Delta\text{F508}$ deletion in 8EIQ is preserved as `INELIGIBLE_NOT_COVERED` with null values. 
The engineered catalytic mutant E1371Q in 6MSM and 8EIQ retains calculated construct features while being flagged as `CONSTRUCT_SUBSTITUTION` in warning logs.

---

## 2. Canonical CFTR Domain Architecture

The canonical sequence (1,480 residues) is partitioned into seven structurally and functionally defined domains 
(Liu et al. *Cell* 2017; Zhang et al. *Science* 2018; Fiedorczuk et al. *Cell* 2022):

| Domain Name | Boundary (UniProt) | Residue Count | Structural & Functional Description |
| :--- | :---: | :---: | :--- |
| **Lasso** | 1..79 | 79 | N-terminal Lasso motif (amphipathic helix regulating channel gating) |
| **TMD1** | 80..388 | 309 | Transmembrane Domain 1 (TM1-TM6, ICL1, ECL1, ICL2, ECL2) |
| **NBD1** | 389..678 | 290 | Nucleotide-Binding Domain 1 (ATP site 1 catalytic core, RI 405-436) |
| **R-domain** | 679..844 | 166 | Regulatory Domain (intrinsically disordered, multiple PKA sites) |
| **TMD2** | 845..1172 | 328 | Transmembrane Domain 2 (TM7-TM12, ICL3, ECL3, ICL4, ECL4) |
| **NBD2** | 1173..1436 | 264 | Nucleotide-Binding Domain 2 (ATP site 2 catalytic core) |
| **C-term** | 1437..1480 | 44 | C-terminal extension (flexible tail with PDZ-binding motif DTRL 1477-1480) |

---

## 3. Structural Feature Coverage & Extraction Summary

> [!IMPORTANT]
> **Coverage Denominators**:
> - **Total Positions**: 1,480 canonical residues.
> - **Eligible Positions**: Residues with observed atomic coordinates in the 3D model.
> - **Ineligible (Unresolved)**: Residues present in construct sequence but disordered/unmodeled in Cryo-EM map.
> - **Ineligible (Not Covered)**: Residues deleted from construct sequence (e.g. $\Delta\text{F508}$ in 8EIQ).

| Structure ID | Filename | Total Canonical | Eligible (Observed) | Unresolved | Not Covered | Mean SASA (Å²) | Mean RSA | Buried (<20%) | Intermediate (20-50%) | Surface (≥50%) | Mean Contacts (8 Å) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **5UAK** | `5UAK(1).cif` | 1,480 | **1139** | 341 | 0 | 54.3 | 0.270 | 535 (47.0%) | 426 (37.4%) | 178 (15.6%) | 22.6 |
| **6MSM** | `6MSM(1).cif` | 1,480 | **1181** | 299 | 0 | 42.5 | 0.211 | 689 (58.3%) | 365 (30.9%) | 127 (10.8%) | 25.3 |
| **8EIQ** | `8EIQ(1).cif` | 1,480 | **1162** | 317 | 1 | 43.2 | 0.216 | 662 (57.0%) | 386 (33.2%) | 114 (9.8%) | 24.4 |

### Secondary Structure Composition Breakdown:

| Structure ID | Total Eligible | Deposited Helix (`_struct_conf`) | Deposited Sheet (`_struct_sheet_range`) | Unannotated Coil (Default) |
| :--- | :---: | :---: | :---: | :---: |
| **5UAK** | **1139** | 807 (70.9%) | 91 (8.0%) | 241 (21.2%) |
| **6MSM** | **1181** | 805 (68.2%) | 100 (8.5%) | 276 (23.4%) |
| **8EIQ** | **1162** | 790 (68.0%) | 107 (9.2%) | 265 (22.8%) |

### Relative Solvent Accessibility (RSA) & Chain-Break Audit:
- **Theoretical Reference Scale**: Tien et al. (2013) maximum solvent accessibility scale derived from extended Gly-X-Gly tripeptide models.
- **Residues with RSA > 1.0**: Exactly 7 observed residues exhibit RSA > 1.0 across all three structures:
  - **6MSM**: Asn900 (RSA: 1.0277, SASA: 200.4 Å²), Asp1202 (RSA: 1.0972, SASA: 211.8 Å²), Pro1451 (RSA: 1.0308, SASA: 163.9 Å²).
  - **8EIQ**: Gly437 (RSA: 1.2108, SASA: 125.9 Å²), Leu636 (RSA: 1.0380, SASA: 208.6 Å²), Gly885 (RSA: 1.0031, SASA: 104.3 Å²), Ser902 (RSA: 1.0284, SASA: 159.4 Å²).
- **Physical Mechanism**: Every single one of these 7 residues directly borders an unobserved/disordered loop break or construct terminus. Lacking coordinates for the flanking residue in the atomic model, backbone atoms (carbonyl oxygen or amide nitrogen) normally shielded by covalent peptide bonds in an intact chain (and in the Gly-X-Gly tripeptide model) become accessible to the 1.4 Å solvent probe. This is an authentic physical consequence of missing flanking coordinates in Cryo-EM models, not an implementation defect.

---

## 4. Multi-Structure Per-Residue Comparative Analysis

Comparison of features across 5UAK (apo), 6MSM (ATP-bound), and 8EIQ (Trikafta/$\Delta$F508) reveals distinct observation patterns:

| Observation Pattern | Position Count | Fraction | Description |
| :--- | :---: | :---: | :--- |
| **ALL_THREE_OBSERVED** | **1127** | 76.1% | Observed coordinates present in 5UAK, 6MSM, and 8EIQ (eligible for 3-way comparative dynamics). |
| **ALL_UNRESOLVED** | **290** | 19.6% | Unmodeled in all three Cryo-EM reconstructions (disordered N-term, regulatory loop, R-domain, C-term). |
| **CONSTRUCT_DELETION** | **1** | 0.1% | Residue 508 (Phe) authentically deleted from 8EIQ construct sequence. |
| **PARTIALLY_OBSERVED** | **62** | 4.2% | Observed coordinates present in 1 or 2 structures (partially dynamic or partially resolved). |

### Multi-Structure Delta Metrics Definition & Non-Causal Usage:
- **Formulas**: $\Delta\text{RSA} = \max(\text{RSA}_{\text{obs}}) - \min(\text{RSA}_{\text{obs}})$ and $\Delta\text{Contacts} = \max(\text{Contacts}_{\text{obs}}) - \min(\text{Contacts}_{\text{obs}})$ computed across states with observed coordinates.
- **Interpretation**: Delta metrics represent unsigned maximal observational ranges of conformational plasticity across static experimental cryo-EM models.
- **Non-Causal Disclaimer**: These metrics must NOT be interpreted as directional responses, thermodynamic free energy shifts, or causal effects of drug binding, ATP ligation, or mutation, as structures represent distinct composite biochemical conditions.

### Key Comparative Observations (Reconciled Values):
1. **Fully Comparable Residues**: Exactly **1,126 canonical positions** have valid observed coordinates and identical wild-type amino acids across all three structures (1,127 all-observed minus position 1371).
2. **Residue 508 (Phe508)**:
   - **5UAK** (apo): SASA = 21.92 Å², RSA = 0.0913 (Buried, 9.1%), Total Contacts = 26 (intra: 26, inter: 0), Sec Struct = HELIX.
   - **6MSM** (ATP-bound): SASA = 5.90 Å², RSA = 0.0246 (Buried, 2.5%), Total Contacts = 28 (intra: 28, inter: 0), Sec Struct = HELIX.
   - **8EIQ** ($\Delta\text{F508}$): `INELIGIBLE_NOT_COVERED` with `failure_reason='CONSTRUCT_DELETION'` (null SASA, RSA, contacts).
   - **Observational Range**: $\Delta\text{RSA} = 0.0667$, $\Delta\text{Contacts} = 2$.
3. **Residue 1371 (Glu1371 / Catalytic Base)**: Located in Walker B catalytic motif (`ILLLDEPSAHL`):
   - **5UAK** (Canonical Glu): SASA = 44.47 Å², RSA = 0.1994 (Buried, 19.9%), Total Contacts = 21 (intra: 21, inter: 0), Sec Struct = HELIX.
   - **6MSM** (Catalytic mutant E1371Q): SASA = 4.29 Å², RSA = 0.0191 (Buried, 1.9%), Total Contacts = 35 (intra: 35 including coordination to ATP and Mg²⁺, inter: 0), Sec Struct = COIL (`Default Assignment (Unannotated in mmCIF)`), `failure_reason='CONSTRUCT_SUBSTITUTION'`.
   - **8EIQ** (Catalytic mutant E1371Q): SASA = 10.78 Å², RSA = 0.0479 (Buried, 4.8%), Total Contacts = 34 (intra: 34 including coordination to ATP and Mg²⁺, inter: 0), Sec Struct = COIL (`Default Assignment (Unannotated in mmCIF)`), `failure_reason='CONSTRUCT_SUBSTITUTION'`.
   - **Observational Range**: $\Delta\text{RSA} = 0.1803$, $\Delta\text{Contacts} = 14$, Secondary Structure Consistency = `DIVERGENT`.

---

## 5. Generated Phase 4 Audit Outputs

- [`structural_feature_summary.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/structural_feature_summary.csv): Aggregate feature metrics and eligibility counts per structure.
- [`residue_structural_features.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/residue_structural_features.csv): Full residue-level table (4,440 rows across 5UAK, 6MSM, 8EIQ).
- [`multi_structure_feature_comparison.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/multi_structure_feature_comparison.csv): 1,480-row cross-structure comparative analysis.
- [`structural_feature_warnings.csv`](file:///Users/amar/CFTR-StructMap/docs/audits/structural_feature_warnings.csv): Log of all feature ineligibility events and construct variations.

