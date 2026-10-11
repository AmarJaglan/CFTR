"""SQLAlchemy models for CFTR-StructMap persistence.

Maintains complete provenance, source immutability records, and distinct
ClinVar aggregate (VCV) and CFTR2 clinical classifications.
"""

from datetime import datetime, timezone
from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import declarative_base, relationship

Base = declarative_base()


class DataSource(Base):
    """Provenance tracking for ingested raw data files."""
    __tablename__ = "data_sources"

    id = Column(Integer, primary_key=True, autoincrement=True)
    source_name = Column(String(64), nullable=False, index=True)  # 'CFTR2', 'ClinVar'
    filename = Column(String(255), nullable=False)
    file_sha256 = Column(String(64), nullable=False, index=True)
    source_version = Column(String(64), nullable=True)
    imported_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)
    record_count = Column(Integer, default=0, nullable=False)
    summary_json = Column(Text, nullable=True)


class CFTR2Variant(Base):
    """Ingested CFTR2 variant records with exact interpretations and validation flags."""
    __tablename__ = "cftr2_variants"

    id = Column(Integer, primary_key=True, autoincrement=True)
    legacy_name = Column(String(128), nullable=False, unique=True, index=True)
    protein_name = Column(String(128), nullable=True)
    cdna_name = Column(String(128), nullable=True, index=True)
    
    # Normalized fields
    canonical_protein = Column(String(32), nullable=True, index=True)
    uniprot_pos = Column(Integer, nullable=True, index=True)
    ref_aa = Column(String(8), nullable=True)
    alt_aa = Column(String(8), nullable=True)
    variant_type = Column(String(32), nullable=False, default="unknown")
    
    # UniProt WT Sequence Validation
    is_valid_wt = Column(Boolean, nullable=True)
    validation_flag = Column(String(32), nullable=False, default="UNVALIDATED")
    validation_message = Column(Text, nullable=True)

    # CFTR2 Clinical and Phenotypic Annotations (Preserved verbatim)
    alleles_reported = Column(Integer, nullable=True)
    allele_frequency = Column(Float, nullable=True)
    determination_current = Column(String(128), nullable=True, index=True)
    determination_previous = Column(String(128), nullable=True)
    change_from_previous = Column(String(32), nullable=True)
    alternative_names = Column(Text, nullable=True)

    # Provenance
    source_id = Column(Integer, ForeignKey("data_sources.id"), nullable=True)
    raw_data_json = Column(Text, nullable=True)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)


class ClinVarVariant(Base):
    """Ingested ClinVar variants with VCV identifiers and verbatim assertions."""
    __tablename__ = "clinvar_variants"

    variation_id = Column(Integer, primary_key=True)  # Stable ClinVar VariationID
    vcv_accession = Column(String(32), nullable=True, index=True)
    name = Column(String(255), nullable=False)
    protein_change_raw = Column(String(64), nullable=True)
    
    # Normalized fields
    canonical_protein = Column(String(32), nullable=True, index=True)
    uniprot_pos = Column(Integer, nullable=True, index=True)
    ref_aa = Column(String(8), nullable=True)
    alt_aa = Column(String(8), nullable=True)
    variant_type = Column(String(64), nullable=False, default="unknown")
    molecular_consequence = Column(String(128), nullable=True)

    # UniProt WT Sequence Validation
    is_valid_wt = Column(Boolean, nullable=True)
    validation_flag = Column(String(32), nullable=False, default="UNVALIDATED")
    validation_message = Column(Text, nullable=True)

    # Clinical Significance & Review Status (Preserved verbatim)
    germline_classification = Column(String(255), nullable=True, index=True)
    review_status = Column(String(255), nullable=True)
    date_last_evaluated = Column(String(64), nullable=True)
    conditions = Column(Text, nullable=True)
    
    # Genomic & Cross-references
    grch38_location = Column(String(64), nullable=True)
    grch37_location = Column(String(64), nullable=True)
    dbsnp_id = Column(String(32), nullable=True)
    allele_id = Column(String(32), nullable=True)

    # Provenance
    source_id = Column(Integer, ForeignKey("data_sources.id"), nullable=True)
    raw_data_json = Column(Text, nullable=True)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)


class Structure(Base):
    """Experimentally determined macromolecular structure record."""
    __tablename__ = "structures"

    id = Column(String(32), primary_key=True)  # PDB identifier (e.g. 5UAK, 6MSM, 8EIQ)
    source_filename = Column(String(255), nullable=False)
    file_sha256 = Column(String(64), nullable=False, index=True)
    title = Column(Text, nullable=True)
    experimental_method = Column(String(64), nullable=True)
    resolution = Column(Float, nullable=True)
    structure_type = Column(String(32), nullable=False, default="FULL_LENGTH")  # FULL_LENGTH, DOMAIN, NON_CFTR
    has_cftr = Column(Boolean, nullable=False, default=True)
    deposition_date = Column(String(32), nullable=True)
    summary_json = Column(Text, nullable=True)
    imported_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)

    chains = relationship("StructuralChain", back_populates="structure", cascade="all, delete-orphan")
    residue_mappings = relationship("ResidueMapping", back_populates="structure", cascade="all, delete-orphan")
    features = relationship("ResidueFeature", back_populates="structure", cascade="all, delete-orphan")


class StructuralChain(Base):
    """Metadata and mapping summary for an individual polymer chain."""
    __tablename__ = "structural_chains"

    id = Column(Integer, primary_key=True, autoincrement=True)
    structure_id = Column(String(32), ForeignKey("structures.id"), nullable=False, index=True)
    entity_id = Column(String(32), nullable=False)
    label_asym_id = Column(String(32), nullable=False)
    auth_asym_id = Column(String(32), nullable=False)
    entity_description = Column(Text, nullable=True)
    organism = Column(String(128), nullable=True)
    polymer_seq_len = Column(Integer, nullable=False, default=0)
    is_cftr = Column(Boolean, nullable=False, default=False, index=True)
    alignment_identity = Column(Float, nullable=True)
    canonical_coverage = Column(Float, nullable=True)
    observed_coverage = Column(Float, nullable=True)
    ca_coverage = Column(Float, nullable=True)
    mapped_positions_count = Column(Integer, nullable=False, default=0)
    observed_residues_count = Column(Integer, nullable=False, default=0)
    ca_atoms_count = Column(Integer, nullable=False, default=0)
    unresolved_count = Column(Integer, nullable=False, default=0)
    not_covered_count = Column(Integer, nullable=False, default=0)
    mismatch_count = Column(Integer, nullable=False, default=0)
    warning_count = Column(Integer, nullable=False, default=0)

    structure = relationship("Structure", back_populates="chains")

    __table_args__ = (
        UniqueConstraint("structure_id", "label_asym_id", name="uq_structure_chain"),
    )


class ResidueMapping(Base):
    """Residue-level correspondence between canonical UniProt sequence and experimental structure."""
    __tablename__ = "residue_mappings"

    id = Column(Integer, primary_key=True, autoincrement=True)
    structure_id = Column(String(32), ForeignKey("structures.id"), nullable=False, index=True)
    source_filename = Column(String(255), nullable=False)
    entity_id = Column(String(32), nullable=False)
    label_asym_id = Column(String(32), nullable=False)
    auth_asym_id = Column(String(32), nullable=False)
    uniprot_accession = Column(String(32), nullable=False, default="P13569")
    canonical_position = Column(Integer, nullable=False, index=True)
    canonical_amino_acid = Column(String(8), nullable=False)
    entity_sequence_position = Column(Integer, nullable=True)
    label_seq_id = Column(Integer, nullable=True)
    auth_seq_id = Column(String(16), nullable=True)
    insertion_code = Column(String(8), nullable=True)
    structure_amino_acid = Column(String(8), nullable=True)
    has_observed_residue = Column(Boolean, nullable=False, default=False)
    has_ca_atom = Column(Boolean, nullable=False, default=False)
    reference_identity_match = Column(Boolean, nullable=True)
    mapping_status = Column(String(32), nullable=False, index=True)
    mapping_warning = Column(Text, nullable=True)

    structure = relationship("Structure", back_populates="residue_mappings")

    __table_args__ = (
        UniqueConstraint("structure_id", "label_asym_id", "canonical_position", name="uq_structure_chain_canonical_pos"),
    )


class ResidueFeature(Base):
    """Calculated residue-level structural features for an individual experimental structure."""
    __tablename__ = "residue_features"

    id = Column(Integer, primary_key=True, autoincrement=True)
    structure_id = Column(String(32), ForeignKey("structures.id"), nullable=False, index=True)
    source_filename = Column(String(255), nullable=False)
    file_sha256 = Column(String(64), nullable=False)
    label_asym_id = Column(String(32), nullable=False)
    auth_asym_id = Column(String(32), nullable=False)
    canonical_position = Column(Integer, nullable=False, index=True)
    canonical_amino_acid = Column(String(8), nullable=False)
    domain_name = Column(String(64), nullable=False, index=True)
    mapping_status = Column(String(32), nullable=False)
    structure_amino_acid = Column(String(8), nullable=True)
    has_observed_residue = Column(Boolean, nullable=False, default=False)
    has_ca_atom = Column(Boolean, nullable=False, default=False)
    eligibility_status = Column(String(32), nullable=False, default="ELIGIBLE")
    failure_reason = Column(String(128), nullable=True)

    # Solvent Accessibility (SASA & RSA)
    sasa_raw_angstrom2 = Column(Float, nullable=True)
    rsa_ratio = Column(Float, nullable=True)
    exposure_category = Column(String(32), nullable=True)
    rsa_method = Column(String(64), nullable=True)
    rsa_scale = Column(String(64), nullable=True)

    # Local Residue Contacts
    contact_distance_threshold = Column(Float, nullable=False, default=8.0)
    total_contacts = Column(Integer, nullable=True)
    intra_chain_contacts = Column(Integer, nullable=True)
    inter_chain_contacts = Column(Integer, nullable=True)

    # Secondary Structure
    secondary_structure = Column(String(32), nullable=True)
    sec_struct_method = Column(String(64), nullable=True)

    # Metadata & Provenance
    software_version = Column(String(64), nullable=True)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)

    structure = relationship("Structure", back_populates="features")

    __table_args__ = (
        UniqueConstraint("structure_id", "label_asym_id", "canonical_position", name="uq_structure_chain_canonical_feature"),
    )


class MultiStructureComparison(Base):
    """Per-residue structural feature comparison across 5UAK, 6MSM, and 8EIQ."""
    __tablename__ = "multi_structure_comparisons"

    canonical_position = Column(Integer, primary_key=True)
    canonical_amino_acid = Column(String(8), nullable=False)
    domain_name = Column(String(64), nullable=False, index=True)
    observation_pattern = Column(String(64), nullable=False)
    construct_variation = Column(String(64), nullable=False)
    is_fully_comparable = Column(Boolean, nullable=False, default=False)

    # 5UAK Features
    sasa_5uak = Column(Float, nullable=True)
    rsa_5uak = Column(Float, nullable=True)
    contacts_5uak = Column(Integer, nullable=True)
    sec_struct_5uak = Column(String(32), nullable=True)

    # 6MSM Features
    sasa_6msm = Column(Float, nullable=True)
    rsa_6msm = Column(Float, nullable=True)
    contacts_6msm = Column(Integer, nullable=True)
    sec_struct_6msm = Column(String(32), nullable=True)

    # 8EIQ Features
    sasa_8eiq = Column(Float, nullable=True)
    rsa_8eiq = Column(Float, nullable=True)
    contacts_8eiq = Column(Integer, nullable=True)
    sec_struct_8eiq = Column(String(32), nullable=True)

    # Cross-structure Deltas & Dynamics
    delta_rsa = Column(Float, nullable=True)
    delta_contacts = Column(Integer, nullable=True)
    sec_struct_consistency = Column(String(32), nullable=True)
    comparison_notes = Column(Text, nullable=True)


class VariantStructuralEvidence(Base):
    """Integrated variant-level structural evidence linking source records to experimental structures.

    Preserves source-record identity (CFTR2 legacy name / ClinVar VariationID) and verbatim
    determinations, attaching per-structure observations and cross-structure metrics.
    """
    __tablename__ = "variant_structural_evidence"

    id = Column(Integer, primary_key=True, autoincrement=True)
    source_db = Column(String(32), nullable=False, index=True)  # 'CFTR2' or 'ClinVar'
    source_record_id = Column(String(128), nullable=False, index=True)  # legacy_name or variation_id
    source_variant_name = Column(String(255), nullable=True)
    protein_change_raw = Column(String(128), nullable=True)
    canonical_protein = Column(String(64), nullable=True, index=True)
    canonical_position = Column(Integer, nullable=True, index=True)
    ref_aa = Column(String(8), nullable=True)
    alt_aa = Column(String(8), nullable=True)
    variant_type = Column(String(64), nullable=False)

    # Verbatim Source Classifications & Review (No consensus synthesis)
    source_determination = Column(String(255), nullable=True, index=True)
    clinvar_review_status = Column(String(255), nullable=True)
    clinvar_vcv_accession = Column(String(32), nullable=True, index=True)
    cftr2_alleles_reported = Column(Integer, nullable=True)

    # Reference Validation & Structural Eligibility
    validation_flag = Column(String(32), nullable=False)
    is_valid_wt = Column(Boolean, nullable=True)
    eligibility_status = Column(String(64), nullable=False, index=True)
    ineligibility_reason = Column(Text, nullable=True)
    domain_name = Column(String(64), nullable=True, index=True)

    # 5UAK (apo full-length) Evidence
    mapping_status_5uak = Column(String(32), nullable=True)
    structure_aa_5uak = Column(String(8), nullable=True)
    feature_status_5uak = Column(String(32), nullable=True)
    failure_reason_5uak = Column(String(128), nullable=True)
    sasa_5uak = Column(Float, nullable=True)
    rsa_5uak = Column(Float, nullable=True)
    exposure_cat_5uak = Column(String(32), nullable=True)
    contacts_5uak = Column(Integer, nullable=True)
    intra_contacts_5uak = Column(Integer, nullable=True)
    inter_contacts_5uak = Column(Integer, nullable=True)
    sec_struct_5uak = Column(String(32), nullable=True)
    construct_note_5uak = Column(Text, nullable=True)

    # 6MSM (ATP-bound full-length) Evidence
    mapping_status_6msm = Column(String(32), nullable=True)
    structure_aa_6msm = Column(String(8), nullable=True)
    feature_status_6msm = Column(String(32), nullable=True)
    failure_reason_6msm = Column(String(128), nullable=True)
    sasa_6msm = Column(Float, nullable=True)
    rsa_6msm = Column(Float, nullable=True)
    exposure_cat_6msm = Column(String(32), nullable=True)
    contacts_6msm = Column(Integer, nullable=True)
    intra_contacts_6msm = Column(Integer, nullable=True)
    inter_contacts_6msm = Column(Integer, nullable=True)
    sec_struct_6msm = Column(String(32), nullable=True)
    construct_note_6msm = Column(Text, nullable=True)

    # 8EIQ (Delta-F508 full-length) Evidence
    mapping_status_8eiq = Column(String(32), nullable=True)
    structure_aa_8eiq = Column(String(8), nullable=True)
    feature_status_8eiq = Column(String(32), nullable=True)
    failure_reason_8eiq = Column(String(128), nullable=True)
    sasa_8eiq = Column(Float, nullable=True)
    rsa_8eiq = Column(Float, nullable=True)
    exposure_cat_8eiq = Column(String(32), nullable=True)
    contacts_8eiq = Column(Integer, nullable=True)
    intra_contacts_8eiq = Column(Integer, nullable=True)
    inter_contacts_8eiq = Column(Integer, nullable=True)
    sec_struct_8eiq = Column(String(32), nullable=True)
    construct_note_8eiq = Column(Text, nullable=True)

    # Multi-structure Comparative Metrics (Observational divergence, not causal impact)
    observation_pattern = Column(String(64), nullable=True)
    construct_variation = Column(String(64), nullable=True)
    is_fully_comparable = Column(Boolean, nullable=False, default=False)
    delta_rsa = Column(Float, nullable=True)
    delta_contacts = Column(Integer, nullable=True)
    sec_struct_consistency = Column(String(32), nullable=True)
    comparison_notes = Column(Text, nullable=True)

    # Provenance
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)

    __table_args__ = (
        UniqueConstraint("source_db", "source_record_id", name="uq_variant_source_record"),
    )


class VariantStructuralSummary(Base):
    """Grouped view of structural evidence by normalized canonical protein change.

    Links back to all contributing source records in CFTR2 and ClinVar without
    synthesizing artificial consensus classifications.
    """
    __tablename__ = "variant_structural_summaries"

    canonical_protein = Column(String(64), primary_key=True)
    canonical_position = Column(Integer, nullable=True, index=True)
    ref_aa = Column(String(8), nullable=True)
    alt_aa = Column(String(8), nullable=True)
    variant_type = Column(String(64), nullable=True)
    domain_name = Column(String(64), nullable=True, index=True)

    # Source Associations
    cftr2_record_count = Column(Integer, default=0, nullable=False)
    cftr2_record_ids = Column(Text, nullable=True)
    cftr2_determinations = Column(Text, nullable=True)

    clinvar_record_count = Column(Integer, default=0, nullable=False)
    clinvar_variation_ids = Column(Text, nullable=True)
    clinvar_vcv_accessions = Column(Text, nullable=True)
    clinvar_classifications = Column(Text, nullable=True)
    clinvar_review_statuses = Column(Text, nullable=True)

    # Structural Metrics
    observation_pattern = Column(String(64), nullable=True)
    construct_variation = Column(String(64), nullable=True)
    is_fully_comparable = Column(Boolean, default=False, nullable=False)

    rsa_5uak = Column(Float, nullable=True)
    contacts_5uak = Column(Integer, nullable=True)
    sec_struct_5uak = Column(String(32), nullable=True)

    rsa_6msm = Column(Float, nullable=True)
    contacts_6msm = Column(Integer, nullable=True)
    sec_struct_6msm = Column(String(32), nullable=True)

    rsa_8eiq = Column(Float, nullable=True)
    contacts_8eiq = Column(Integer, nullable=True)
    sec_struct_8eiq = Column(String(32), nullable=True)

    delta_rsa = Column(Float, nullable=True)
    delta_contacts = Column(Integer, nullable=True)
    sec_struct_consistency = Column(String(32), nullable=True)
    caveats_or_warnings = Column(Text, nullable=True)

    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)

