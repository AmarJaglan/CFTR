# Phase 4 Scientific Validation and Audit Reconciliation Report

**Date**: October 11, 2026  
**Project**: CFTR-StructMap  
**Target Canonical Sequence**: Human CFTR (UniProt [P13569](https://www.uniprot.org/uniprotkb/P13569/entry), 1,480 amino acids)  
**Evaluated Structures**: 
- `5UAK` (Inward-facing, dephosphorylated apo state; PDB 5UAK)
- `6MSM` (Outward-facing, ATP-bound state with catalytic Walker B substitution E1371Q; PDB 6MSM)
- `8EIQ` (Inward-facing, Trikafta-bound $\Delta\text{F508}$ state with catalytic substitution E1371Q; PDB 8EIQ)

---

## 1. Executive Summary

This independent scientific audit evaluated the Phase 4 residue-level structural feature extraction implementation against its persisted CSV outputs and summary markdown documentation. 

### Key Findings & Corrections:
1. **Reconciliation of Residues 508 and 1371**: A numerical discrepancy was identified between the initial written summary text and the underlying CSV outputs. Investigation revealed that the markdown generator in `features.py` originally contained static draft placeholder strings for residues 508 and 1371 rather than values dynamically linked to the calculation pipeline. The underlying CSV values were verified as computationally correct and the markdown report has been updated with the audited values.
2. **Residues with RSA > 1.0**: Exactly 7 observed residues exhibit Relative Solvent Accessibility (RSA) greater than 1.0 (ranging from 1.003 to 1.211) across 6MSM and 8EIQ. All 7 residues directly border unmodeled loop breaks or construct termini. In the absence of coordinates for the flanking residue in the Cryo-EM model, backbone atoms that are covalently shielded in an unbroken chain (and in the intact Gly-X-Gly tripeptide reference model of Tien et al. 2013) become accessible to the 1.4 Å solvent probe. This is an authentic physical consequence of missing flanking coordinates in experimental structures, not an implementation defect.
3. **Secondary Structure Provenance**: mmCIF deposited records (`_struct_conf`, `_struct_sheet_range`) explicitly assign helices and sheets, but depositors do not deposit coil annotations. The implementation has been updated to assign `mmCIF Depositor (_struct_conf, _struct_sheet_range)` strictly to deposited helices and sheets, while observed unannotated residues are explicitly tagged with `Default Assignment (Unannotated in mmCIF)`.
4. **Construct Modification Warnings**: Construct substitutions (`E1371Q` in 6MSM and 8EIQ) were originally stored with an internal flag but omitted from `structural_feature_warnings.csv`. Both `E1371Q` (`CONSTRUCT_SUBSTITUTION`) and `\Delta\text{F508}` (`CONSTRUCT_DELETION`) are now logged in `structural_feature_warnings.csv`.
5. **Delta Metrics Framing**: Multi-structure delta metrics represent unsigned maximal observational ranges ($\max - \min$) across static Cryo-EM models. They are explicitly documented with non-causal warnings.
6. **Automated Verification**: Six new regression tests were added to `backend/tests/test_structural_features.py` (expanding the test suite to 14 Phase 4 tests, 54 tests total across the project), all passing with a 100% pass rate.

---

## 2. Investigation 1: Reconciliation of Canonical Residues 508 and 1371

### Discrepancy Analysis
In the initial Phase 4 summary report, placeholder draft text stated:
- Phe508: "5UAK (RSA: 0.175, Buried; Contacts: 28) and 6MSM (RSA: 0.088, Buried; Contacts: 30)"
- Glu1371: "6MSM (Contacts: 28) and 8EIQ (Contacts: 26)"

When audited against the computational pipeline and `docs/audits/residue_structural_features.csv`, the actual values were:
- Phe508:
  - 5UAK: RSA = **0.0913**, Contacts = **26**
  - 6MSM: RSA = **0.0246**, Contacts = **28**
- Glu1371:
  - 5UAK: RSA = **0.1994**, Contacts = **21**
  - 6MSM: RSA = **0.0191**, Contacts = **35**
  - 8EIQ: RSA = **0.0479**, Contacts = **34**

### Root Cause
In `backend/app/scientific/features.py`, the markdown report generator function `generate_feature_audit_markdown()` contained static text authored during initial testing. The calculation engine and CSV export functions had written the true values, but the markdown generator had not been updated.

### Reconciled Residue-Level Table

| Canonical Position | Structure ID | Amino Acid | Mapping Status | Eligibility Status | SASA (Å²) | RSA Ratio | Exposure Category | Total Contacts | Intra-Chain | Inter-Chain | Secondary Structure | Sec Struct Method |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **508** | **5UAK** | F | `MAPPED_OBSERVED` | `ELIGIBLE` | 21.92 | 0.0913 | BURIED | 26 | 26 | 0 | HELIX | mmCIF Depositor |
| **508** | **6MSM** | F | `MAPPED_OBSERVED` | `ELIGIBLE` | 5.90 | 0.0246 | BURIED | 28 | 28 | 0 | HELIX | mmCIF Depositor |
| **508** | **8EIQ** | F | `NOT_COVERED` | `INELIGIBLE_NOT_COVERED` | *null* | *null* | *null* | *null* | *null* | *null* | *null* | *null* |
| **1371** | **5UAK** | E | `MAPPED_OBSERVED` | `ELIGIBLE` | 44.47 | 0.1994 | BURIED | 21 | 21 | 0 | HELIX | mmCIF Depositor |
| **1371** | **6MSM** | Q | `REFERENCE_MISMATCH` | `ELIGIBLE` | 4.29 | 0.0191 | BURIED | 35 | 35 | 0 | COIL | Default Assignment |
| **1371** | **8EIQ** | Q | `REFERENCE_MISMATCH` | `ELIGIBLE` | 10.78 | 0.0479 | BURIED | 34 | 34 | 0 | COIL | Default Assignment |

### Multi-Structure Observational Deltas
- **Residue 508**:
  - Observation Pattern: `CONSTRUCT_DELETION`
  - $\Delta\text{RSA} = 0.0913 - 0.0246 = 0.0667$
  - $\Delta\text{Contacts} = 28 - 26 = 2$
  - Comparability: `is_fully_comparable = False` (absent in 8EIQ).
- **Residue 1371**:
  - Observation Pattern: `ALL_THREE_OBSERVED`
  - Construct Variation: `CONSTRUCT_SUBSTITUTION`
  - $\Delta\text{RSA} = \max(0.1994, 0.0191, 0.0479) - \min(0.1994, 0.0191, 0.0479) = 0.1994 - 0.0191 = 0.1803$
  - $\Delta\text{Contacts} = \max(21, 35, 34) - \min(21, 35, 34) = 35 - 21 = 14$
  - Secondary Structure Consistency: `DIVERGENT` (HELIX in 5UAK vs COIL in 6MSM and 8EIQ).
  - Comparability: `is_fully_comparable = False` (catalytic mutation E1371Q in 6MSM and 8EIQ).

---

## 3. Investigation 2: Validation of RSA, Tien et al. (2013) Scale, and RSA > 1.0 Audit

### Reference Scale Validation
The maximum SASA reference values in `TIEN_MAX_SASA_2013` were cross-checked against Table 1 ("Theoretical" column) of **Tien MZ et al. (2013) PLOS ONE 8(11): e80635**:
- All 20 canonical amino acids match the published theoretical maximum accessibility scale derived from extended Gly-X-Gly tripeptide models ($r_{\text{probe}} = 1.4\text{ \AA}$).
- Lookups correctly resolve both 1-letter and 3-letter codes. For construct variations (e.g. Q1371), the scale value for the modeled amino acid (Gln: 225.0 Å²) is used.

### Shrake-Rupley Implementation & Alternate Conformations
- **Algorithm**: `Bio.PDB.ShrakeRupley(probe_radius=1.4, n_points=100)`.
- **Atom Radii**: Chothia standard van der Waals radii.
- **Alternate Conformations**: Inspected `5UAK(1).cif`, `6MSM(1).cif`, and `8EIQ(1).cif` using Biopython. All three Cryo-EM models contain **0 disordered/alternate conformation atoms**.

### Investigation of Residues with RSA > 1.0
Querying `docs/audits/residue_structural_features.csv` identified exactly 7 residues with RSA > 1.0:

| Structure ID | Canonical Pos | Ref AA | Struct AA | Domain | Raw SASA (Å²) | Tien Max (Å²) | RSA Ratio | Preceding Position Modeled? | Succeeding Position Modeled? | Direct Chain Break Border? |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **6MSM** | 900 | N | N | TMD2 | 200.40 | 195.0 | 1.0277 | No (899 unmodeled) | Yes (901 modeled) | **Yes (N-term break)** |
| **6MSM** | 1202 | D | D | NBD2 | 211.75 | 193.0 | 1.0972 | No (1201 unmodeled) | Yes (1203 modeled) | **Yes (N-term break)** |
| **6MSM** | 1451 | P | P | C-term | 163.89 | 159.0 | 1.0308 | Yes (1450 modeled) | No (1452 unmodeled) | **Yes (C-term break)** |
| **8EIQ** | 437 | G | G | NBD1 | 125.92 | 104.0 | 1.2108 | No (436 unmodeled) | Yes (438 modeled) | **Yes (RI loop break)** |
| **8EIQ** | 636 | L | L | NBD1 | 208.64 | 201.0 | 1.0380 | Yes (635 modeled) | No (637 unmodeled) | **Yes (C-term break)** |
| **8EIQ** | 885 | G | G | TMD2 | 104.32 | 104.0 | 1.0031 | Yes (884 modeled) | No (886 unmodeled) | **Yes (C-term break)** |
| **8EIQ** | 902 | S | S | TMD2 | 159.40 | 155.0 | 1.0284 | No (901 unmodeled) | Yes (903 modeled) | **Yes (N-term break)** |

#### Physical Explanation
In 100% of cases (7 of 7), the residue directly borders a disordered region or construct terminus. In an intact protein chain (and in the Gly-X-Gly reference tripeptide), backbone carbonyl and amide groups are engaged in covalent peptide bonds that sterically occlude solvent. At coordinate chain breaks, these open backbone atoms become accessible to the 1.4 Å solvent probe, naturally yielding a raw SASA slightly larger than the Gly-X-Gly tripeptide maximum (by 0.3% to 21%). This represents an authentic physical property of experimental Cryo-EM models with missing loop coordinates, not an implementation defect.

---

## 4. Investigation 3: Local Residue Contacts & Chain/Entity Inclusion

### Methodological Audit
1. **Contact Sphere**: Heavy-atom Cartesian `KDTree` with configurable distance threshold (default $8.0\text{ \AA}$).
2. **Self-Exclusion**: Atoms belonging to the same residue are strictly excluded (`close_res != res`).
3. **Residue Counting**: Contacts count unique neighboring residues, not atom pairs.
4. **Partitioning**:
   - `intra_chain_contacts`: Neighboring residues belonging to author chain A.
   - `inter_chain_contacts`: Neighboring residues belonging to partner chains (chain B).
   - `total_contacts`: Strictly equals `intra_chain_contacts + inter_chain_contacts` across all 3,482 eligible residue records (0 discrepancies).

### Chain and Entity Inclusion
- **Excluded Atoms**: Hydrogen atoms (`atom.element == 'H'`) are strictly excluded.
- **Included Chains**: Author chain A (CFTR) and partner chain B (Fab, nanobody, or R-domain fragment).
  - Residues with `inter_chain_contacts > 0`: 48 in 5UAK, 45 in 6MSM, 43 in 8EIQ (representing structural contact with antibody/chaperone partners).
- **Hetero Entities / Ligands**: Biopython's `level="R"` search returns all entities modeled as residues. In 6MSM and 8EIQ, this includes co-modeled ligands (ATP, $\text{Mg}^{2+}$, cholesterol, lipids, elexacaftor, tezacaftor, ivacaftor). For example, residue 1371 in 6MSM has 35 total contacts, which includes 33 CFTR protein residues and 2 direct coordination contacts to ATP and $\text{Mg}^{2+}$.

---

## 5. Investigation 4: Secondary Structure Assignment Audit

### Findings on mmCIF Deposited Categories
- Depositors deposit helices in `_struct_conf` (`conf_type_id = 'HELX_P'`).
- Depositors deposit $\beta$-sheets in `_struct_sheet_range`.
- Depositors do **not** deposit coil annotations (`_struct_coil` does not exist in mmCIF).
- In the original implementation, unannotated observed residues were assigned `"COIL"`, but the provenance column `sec_struct_method` described all residues as `"mmCIF Depositor (_struct_conf, _struct_sheet_range)"`.

### Implemented Correction
The implementation now explicitly separates deposited records from default assignments:
- For residues in `_struct_conf` or `_struct_sheet_range`:
  - `secondary_structure`: `"HELIX"` or `"SHEET"`
  - `sec_struct_method`: `"mmCIF Depositor (_struct_conf, _struct_sheet_range)"`
- For observed residues absent from deposited categories:
  - `secondary_structure`: `"COIL"`
  - `sec_struct_method`: `"Default Assignment (Unannotated in mmCIF)"`

### Secondary Structure Breakdown

| Structure ID | Eligible Positions | Deposited Helices (`_struct_conf`) | Deposited Sheets (`_struct_sheet_range`) | Default Unannotated Coils |
| :---: | :---: | :---: | :---: | :---: |
| **5UAK** | **1,139** | 807 (70.9%) | 91 (8.0%) | 241 (21.2%) |
| **6MSM** | **1,181** | 805 (68.2%) | 100 (8.5%) | 276 (23.4%) |
| **8EIQ** | **1,162** | 790 (68.0%) | 107 (9.2%) | 265 (22.8%) |

`docs/audits/structural_feature_summary.csv` now records `deposited_helix_count`, `deposited_sheet_count`, and `unannotated_coil_count` columns.

---

## 6. Investigation 5: Construct Substitutions and Deletions in Feature Warnings

### Audit Finding
- The authentic construct deletion $\Delta\text{F508}$ in 8EIQ was logged in `structural_feature_warnings.csv`.
- However, the engineered catalytic substitution E1371Q in 6MSM and 8EIQ was stored in `ResidueFeature.failure_reason` but omitted from `structural_feature_warnings.csv`.

### Implemented Correction
`extract_structure_features()` was updated to append a structured warning whenever `mapping_status == STATUS_REFERENCE_MISMATCH` (`CONSTRUCT_SUBSTITUTION`).

### Reconciled `structural_feature_warnings.csv`
The regenerated `docs/audits/structural_feature_warnings.csv` contains exactly 3 entries:
1. `6MSM`, Pos 1371: `CONSTRUCT_SEQUENCE`, `CONSTRUCT_SUBSTITUTION`, "Construct substitution: canonical E1371 modeled as Q (CONSTRUCT_SUBSTITUTION)."
2. `8EIQ`, Pos 508: `ALL`, `CONSTRUCT_DELETION`, "Residue F508 absent from construct sequence (CONSTRUCT_DELETION)."
3. `8EIQ`, Pos 1371: `CONSTRUCT_SEQUENCE`, `CONSTRUCT_SUBSTITUTION`, "Construct substitution: canonical E1371 modeled as Q (CONSTRUCT_SUBSTITUTION)."

---

## 7. Investigation 6: Summary Metrics, Exposure Categories, and Denominators

### Accounting Integrity
Direct recomputation from `docs/audits/residue_structural_features.csv` (4,440 rows) verified:
- **Total Canonical Positions**: Exactly 1,480 for every structure.
- **Eligible + Unresolved + Not Covered = 1,480**:
  - 5UAK: $1,139 + 341 + 0 = 1,480$
  - 6MSM: $1,181 + 299 + 0 = 1,480$
  - 8EIQ: $1,162 + 317 + 1 = 1,480$
- **Exposure Categories Sum = Eligible**:
  - 5UAK: 535 Buried + 426 Intermediate + 178 Surface = 1,139.
  - 6MSM: 689 Buried + 365 Intermediate + 127 Surface = 1,181.
  - 8EIQ: 662 Buried + 386 Intermediate + 114 Surface = 1,162.
- **Null Handling**: All 957 ineligible residue instances across the three structures have strictly `null` values for SASA, RSA, exposure category, contacts, and secondary structure. No zero substitution occurred.

---

## 8. Investigation 7: Multi-Structure Delta Metrics Definition & Non-Causal Framing

### Mathematical Formulation
$$\Delta\text{RSA} = \max(\text{RSA}_{\text{observed}}) - \min(\text{RSA}_{\text{observed}})$$
$$\Delta\text{Contacts} = \max(\text{Contacts}_{\text{observed}}) - \min(\text{Contacts}_{\text{observed}})$$
Calculated whenever coordinates are observed in $\ge 2$ conformational states. If only 1 or 0 states are observed, delta values are `null`.

In `features.py`, `delta_rsa` calculation was updated to explicitly round to 4 decimal places (`round(..., 4)`) to eliminate floating-point precision artifacts.

### Non-Causal Usage Notice
These metrics quantify the unsigned maximal range of conformational variation across static Cryo-EM models. They must **not** be interpreted as causal effects of small-molecule modulator binding, phosphorylation, or mutation, because the compared structures differ across multiple experimental variables simultaneously.

---

## 9. Automated Regression Testing Results

Six dedicated regression tests were added to `backend/tests/test_structural_features.py`:
1. `test_reconciled_values_for_residues_508_and_1371`: Validates exact reconciled RSA, contacts, and deltas for 508 and 1371 across all three structures.
2. `test_rsa_greater_than_one_chain_break_audit`: Asserts that all residues with RSA > 1.0 border a disordered chain break.
3. `test_secondary_structure_deposited_vs_unannotated_distinction`: Asserts that helices/sheets are tagged with `mmCIF Depositor` while coils are tagged with `Default Assignment`.
4. `test_construct_variation_in_warnings`: Asserts that `structural_feature_warnings.csv` captures both E1371Q and $\Delta$F508.
5. `test_contact_counts_sum_and_self_exclusion`: Asserts that `total_contacts == intra_chain + inter_chain` and self-contacts are excluded.
6. `test_delta_metrics_are_strictly_non_causal_max_minus_min`: Asserts that delta metrics are unsigned maximal ranges.

### Test Execution Summary
```bash
.venv/bin/pytest backend/tests/ -v
```
- **Total Tests**: **54 passed**, 0 failed in 56.76s (100% pass rate).
  - Phase 1 (FastAPI backend, CORS, health): 6 passed
  - Phase 2 (Ingestion, normalization, reference validation): 18 passed
  - Phase 3 (mmCIF parsing, canonical mapping, coverage): 16 passed
  - Phase 4 & 4.1 (Feature extraction, Tien scale, contacts, audit reconciliation): 14 passed

---

## 10. Unresolved Limitations & Research Recommendations

1. **Local Cryo-EM Resolution Variation**: SASA and contact calculations use deposited Cartesian coordinates without weighting by local map resolution or B-factors. In regions of lower local Cryo-EM resolution (e.g. flexible loops), coordinate uncertainty is higher than in the transmembrane core.
2. **Absence of Hydrogen Atoms**: Experimental Cryo-EM models at 3.2–3.9 Å resolution lack modeled hydrogen coordinates. All contact searches and SASA calculations operate on heavy atoms only.
3. **Tripeptide Reference Model at Chain Termini**: The Tien et al. (2013) scale uses an intact Gly-X-Gly tripeptide. Terminal residues bordering unmodeled loops naturally display RSA > 1.0 because of open backbone faces. Downstream users should interpret RSA values for chain-break residues with this physical context.
