# CFTR-StructMap — Phase 5: Variant-Level Structural Evidence Integration Summary

**Generated:** 2026-10-11 16:48:16 UTC  
**Canonical Reference:** Human CFTR UniProt **P13569** (1,480 amino acids, SV=3)  
**Experimental Structures:** 5UAK (apo full-length), 6MSM (ATP-bound full-length), 8EIQ (Delta-F508 full-length)

---

## 1. Executive Summary & Purpose

Phase 5 establishes a reproducible, provenance-preserving structural evidence integration layer linking clinical and aggregate registry variant records from **CFTR2** and **ClinVar** to canonical residue coordinates and experimentally determined cryo-EM models.

### Key Scientific Distinctions
1. **Source Record Preservation**: Each raw source record is retained with its authentic identifiers (CFTR2 legacy name, ClinVar VariationID / VCV accession). Records sharing identical protein changes are **never merged** into single records in the primary evidence table.
2. **ClinVar VCV Aggregate Context**: Imported ClinVar records represent aggregate VariationID/VCV entities, not submitter-level SCV records.
3. **No Clinical Pathogenicity Prediction**: This module strictly extracts and contextualizes structural observations. It **does not** compute pathogenicity scores, damage probabilities, or machine-learning verdicts.
4. **Non-Causal Observational Divergence**: $\Delta\text{RSA}$ and $\Delta\text{Contacts}$ across 5UAK, 6MSM, and 8EIQ represent static structural divergence between distinct experimental models and conformations, **not** the causal biophysical effect of a patient mutation.

---

## 2. Source Record Accounting & Eligibility Partitioning

All source records from Phase 2 ingestion were processed without discarding any entry.

| Dataset | Total Source Records | Supported Missense | Synonymous | In-Frame Deletions | Reference Mismatches | Unsupported Classes | Unresolved Changes |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **CFTR2** | 2,092 | 779 | 40 | 16 | 1 | 696 | 560 |
| **ClinVar** | 2,664 | 2,643 | 0 | 0 | 0 | 8 | 13 |
| **Total** | **4,756** | **3,422** | **40** | **16** | **1** | **704** | **573** |

### Denominator Definitions
- **Supported Missense (`ELIGIBLE_MISSENSE`)**: Single-amino-acid substitutions where reported reference residue matches canonical UniProt P13569 (`MATCH`), alternate residue is a single amino acid, and `ref_aa != alt_aa`.
- **Synonymous (`SYNONYMOUS`)**: Single nucleotide alterations where amino acid identity is unaltered (`ref_aa == alt_aa`). Mapped to canonical position for contextual reference, but flagged as having no primary structural substitution.
- **In-Frame Deletion (`IN_FRAME_DELETION`)**: Deletion of one or more amino acids (e.g. F508del). Mapped to canonical coordinates with explicit construct deletion handling.
- **Reference Mismatch (`REFERENCE_MISMATCH`)**: Variant where the reported reference amino acid contradicts canonical UniProt P13569 (e.g. historical CFTR2 M470V). Preserved verbatim, flagged, and excluded from wild-type structural feature attachment.
- **Unsupported Classes (`UNSUPPORTED_VARIANT_CLASS`)**: Non-single-amino-acid molecular consequences, including nonsense (premature stops), frameshifts, complex haplotypes with non-coding modifiers, and delins.
- **Unresolved Changes (`UNRESOLVED_PROTEIN_CHANGE`)**: Variants where cDNA or protein change could not be resolved or normalized during ingestion (`p.?`, unparsed notation).

---

## 3. Structural Evidence Availability

For primary supported missense variants (3,422 total source records across CFTR2 and ClinVar):

| Structure | Model Identity | Conformational State | Variants with Features | Coverage of Eligible |
| :--- | :--- | :--- | :---: | :---: |
| **5UAK** | Human CFTR | Dephosphorylated, ATP-free (apo) | 2,830 | 82.7% |
| **6MSM** | Human CFTR (E1371Q) | Phosphorylated, ATP-bound dimer | 2,908 | 85.0% |
| **8EIQ** | Human CFTR (Delta-F508, E1371Q) | Phosphorylated, ATP-bound dimer | 2,877 | 84.1% |
| **All Three** | Multi-structure | Conformationally complete | **2,801** | **81.9%** |

### Observation Patterns (Supported Missense)
- **ALL_THREE_OBSERVED**: 2,801 variants (81.9%)
- **ALL_UNRESOLVED**: 492 variants (14.4%)
- **PARTIALLY_OBSERVED**: 124 variants (3.6%)
- **CONSTRUCT_DELETION**: 5 variants (0.1%)

---

## 4. Special Scientific Cases & Caveat Handling

### 1. F508del (Canonical Position 508)
- In **5UAK** (apo) and **6MSM** (ATP-bound), residue 508 is wild-type Phenylalanine with observed coordinates (5UAK RSA: 0.0913, contacts: 26; 6MSM RSA: 0.0246, contacts: 28).
- In **8EIQ**, residue 508 represents the authentic CF-causing **Delta-F508 patient deletion**. The residue is marked `NOT_COVERED` and `INELIGIBLE_NOT_COVERED` with failure reason `CONSTRUCT_DELETION`. Structural features for 8EIQ are strictly **null**, never fabricated or imputed.

### 2. E1371Q Catalytic Construct Substitution (Canonical Position 1371)
- Position 1371 in **5UAK** is wild-type Glutamate (E1371, RSA: 0.1994, contacts: 21).
- In **6MSM** and **8EIQ**, residue 1371 carries the engineered catalytic mutation **E1371Q** designed to arrest ATP hydrolysis. The observation is explicitly flagged as `CONSTRUCT_SUBSTITUTION` with `mapping_status = REFERENCE_MISMATCH`. Features are retained for structural reference but distinguished from wild-type observations.

### 3. M470V Canonical Reference Mismatch
- Canonical UniProt P13569 standardizes on **Valine at position 470** (rs1800089).
- CFTR2 historically recorded the legacy name `M470V` (presuming Met470 as reference).
- CFTR-StructMap preserves the historical CFTR2 record verbatim without altering its notation or clinical determination, while flagging it as `REFERENCE_MISMATCH`. Wild-type structural features of Valine 470 are **not** falsely attached to a putative wild-type Methionine.
- In ClinVar, VariationID 7130 correctly records `p.Val470Met` (`V470M`), which validates as `MATCH`.

### 4. Intrinsically Disordered R-Domain (Positions 679–844)
- Variants located within the intrinsically disordered regulatory domain (R-domain) have null structural coordinates across 5UAK, 6MSM, and 8EIQ (`MAPPED_UNRESOLVED`, `INELIGIBLE_UNRESOLVED`).
- The absence of structural features is documented strictly as an **evidence availability limitation**, never as an indicator of neutrality or benign impact.

---

## 5. CFTR2 and ClinVar Side-by-Side Classification Comparison

A total of **3,353** distinct normalized protein changes were identified across both datasets.
There are **655** distinct protein changes present in both CFTR2 and ClinVar, yielding **693** pairwise source comparisons.

| Concordance Category | Description | Pairwise Comparisons | Percentage |
| :--- | :--- | :---: | :---: |
| **CFTR2_NO_INTERPRETATION** | Classification concordance descriptor | 253 | 36.5% |
| **CONCORDANT_PATHOGENIC** | Classification concordance descriptor | 171 | 24.7% |
| **CLINVAR_CONFLICTING_CLASSIFICATIONS** | Classification concordance descriptor | 148 | 21.4% |
| **CFTR2_VARYING_CLINICAL_CONSEQUENCE** | Classification concordance descriptor | 80 | 11.5% |
| **CLINVAR_UNCERTAIN_SIGNIFICANCE** | Classification concordance descriptor | 26 | 3.8% |
| **OTHER_OR_DRUG_RESPONSE** | Classification concordance descriptor | 13 | 1.9% |
| **CONCORDANT_BENIGN** | Classification concordance descriptor | 1 | 0.1% |
| **DISCORDANT** | Classification concordance descriptor | 1 | 0.1% |

> [!IMPORTANT]
> CFTR-StructMap preserves clinical determinations verbatim and **never manufactures a consensus classification**. Inconsistencies between CFTR2 clinical verdicts and ClinVar submitter assertions are surfaced transparently to support expert review.

---

## 6. Audit Warnings & Anomaly Accounting

A total of **1,959** contextual warnings were logged during evidence integration:

| Warning Type | Total Occurrences | Scope & Interpretation |
| :--- | :---: | :--- |
| **UNSUPPORTED_VARIANT_CLASS** | 704 | Contextual caveat for downstream interpretation |
| **UNRESOLVED_COORDINATES** | 634 | Contextual caveat for downstream interpretation |
| **UNRESOLVED_PROTEIN_CHANGE** | 573 | Contextual caveat for downstream interpretation |
| **SYNONYMOUS_NO_STRUCTURAL_CHANGE** | 40 | Contextual caveat for downstream interpretation |
| **CONSTRUCT_DELETION** | 6 | Contextual caveat for downstream interpretation |
| **REFERENCE_MISMATCH** | 1 | Contextual caveat for downstream interpretation |
| **CONSTRUCT_SUBSTITUTION** | 1 | Contextual caveat for downstream interpretation |

---

## 7. Deliverables & Artifacts Generated

1. `docs/audits/variant_structural_evidence.csv`: Per-source-record integrated evidence table (4,756 rows).
2. `docs/audits/variant_structural_summary.csv`: Grouped view by normalized canonical protein change (3,353 rows).
3. `docs/audits/variant_classification_comparison.csv`: Pairwise CFTR2 vs ClinVar classifications (693 rows).
4. `docs/audits/variant_evidence_warnings.csv`: Complete audit log of construct and reference warnings.
5. `docs/audits/phase5_integration_summary.md`: This comprehensive audit report.
