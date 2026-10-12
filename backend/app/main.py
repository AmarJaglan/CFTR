from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import settings

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.VERSION,
    description=settings.DESCRIPTION,
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/openapi.json",
)

# CORS configuration restricted to configured development origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
)


@app.get("/", tags=["General"])
async def root():
    """Root endpoint returning basic metadata and service information."""
    return {
        "name": settings.PROJECT_NAME,
        "version": settings.VERSION,
        "description": settings.DESCRIPTION,
        "docs_url": "/docs",
        "status": "online",
    }


import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, List, Optional
from pydantic import BaseModel, Field
from fastapi import Depends, HTTPException, Query
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from backend.app.db.models import (
    CFTR2Variant,
    ClinVarVariant,
    MultiStructureComparison,
    ResidueFeature,
    ResidueMapping,
    StructuralChain,
    Structure,
    VariantStructuralEvidence,
    VariantStructuralSummary,
)
from backend.app.db.session import get_db
from backend.app.scientific.normalization import normalize_protein_change



@app.get("/health", tags=["Health"])
async def health_check():
    """Health check endpoint for container orchestration and uptime monitoring."""
    return {
        "status": "healthy",
        "service": settings.PROJECT_NAME,
        "version": settings.VERSION,
    }


@app.get("/stats", tags=["Statistics"])
async def get_database_statistics(db: Session = Depends(get_db)):
    """Retrieve verified database aggregate statistics and data-loading status."""
    try:
        total_evidence = db.query(VariantStructuralEvidence).count()
        total_summaries = db.query(VariantStructuralSummary).count()
        cftr2_count = db.query(CFTR2Variant).count()
        clinvar_count = db.query(ClinVarVariant).count()
        structures_count = db.query(Structure).filter_by(has_cftr=True).count()
        control_count = db.query(Structure).filter_by(has_cftr=False).count()
        comparable_count = db.query(MultiStructureComparison).filter_by(is_fully_comparable=True).count()
        return {
            "database_connected": True,
            "evidence_loaded": total_evidence > 0,
            "total_evidence_records": total_evidence,
            "grouped_protein_changes": total_summaries,
            "cftr2_records": cftr2_count,
            "clinvar_records": clinvar_count,
            "cftr_structures_count": structures_count,
            "control_structures_count": control_count,
            "comparable_positions": comparable_count,
        }
    except Exception as e:
        return {
            "database_connected": False,
            "evidence_loaded": False,
            "error": str(e),
            "total_evidence_records": 0,
            "grouped_protein_changes": 0,
            "cftr2_records": 0,
            "clinvar_records": 0,
            "cftr_structures_count": 0,
            "control_structures_count": 0,
            "comparable_positions": 0,
        }


@app.get("/variants", tags=["Variants"])
async def list_variants(
    position: int | None = Query(None, description="Filter by 1-based canonical CFTR position (1..1480)"),
    domain: str | None = Query(None, description="Filter by canonical domain name (e.g. NBD1, TMD1)"),
    eligibility: str | None = Query(None, description="Filter by eligibility status"),
    source_db: str | None = Query(None, description="Filter by source database ('CFTR2' or 'ClinVar')"),
    search: str | None = Query(None, description="Search term for variant name, protein change, or source ID"),
    limit: int = Query(50, ge=1, le=500, description="Page limit"),
    offset: int = Query(0, ge=0, description="Page offset"),
    db: Session = Depends(get_db),
):
    """Retrieve paginated variant structural evidence records with optional filtering."""
    query = db.query(VariantStructuralEvidence)
    if search:
        s = search.strip()
        if s.isdigit():
            pos = int(s)
            query = query.filter(
                (VariantStructuralEvidence.canonical_position == pos)
                | (VariantStructuralEvidence.source_record_id.ilike(f"%{s}%"))
            )
        else:
            query = query.filter(
                (VariantStructuralEvidence.canonical_protein.ilike(f"%{s}%"))
                | (VariantStructuralEvidence.source_record_id.ilike(f"%{s}%"))
                | (VariantStructuralEvidence.source_variant_name.ilike(f"%{s}%"))
            )
    if position is not None:
        query = query.filter(VariantStructuralEvidence.canonical_position == position)
    if domain:
        query = query.filter(VariantStructuralEvidence.domain_name.ilike(f"%{domain}%"))
    if eligibility:
        query = query.filter(VariantStructuralEvidence.eligibility_status == eligibility)
    if source_db:
        query = query.filter(VariantStructuralEvidence.source_db.ilike(source_db))

    total = query.count()
    records = query.order_by(VariantStructuralEvidence.canonical_position, VariantStructuralEvidence.id).offset(offset).limit(limit).all()

    return {
        "total": total,
        "limit": limit,
        "offset": offset,
        "results": [
            {
                "source_db": r.source_db,
                "source_record_id": r.source_record_id,
                "source_variant_name": r.source_variant_name,
                "canonical_protein": r.canonical_protein,
                "canonical_position": r.canonical_position,
                "variant_type": r.variant_type,
                "source_determination": r.source_determination,
                "clinvar_review_status": r.clinvar_review_status,
                "eligibility_status": r.eligibility_status,
                "domain_name": r.domain_name,
                "observation_pattern": r.observation_pattern,
                "is_fully_comparable": r.is_fully_comparable,
            }
            for r in records
        ],
    }


class BatchVariantRequest(BaseModel):
    identifiers: List[str] = Field(
        ...,
        description="List of variant identifiers to analyze (up to 500 rows)",
        examples=[["G551D", "F508del", "p.Gly551Asp", "E1371Q", "M470V"]],
    )


class BatchRowResult(BaseModel):
    row_number: int
    original_input: str
    normalized_identifier: Optional[str] = None
    status: str  # 'RESOLVED', 'NOT_FOUND', 'UNSUPPORTED_FORMAT', 'EMPTY'
    canonical_protein: Optional[str] = None
    canonical_position: Optional[int] = None
    ref_aa: Optional[str] = None
    alt_aa: Optional[str] = None
    domain_name: Optional[str] = None
    source_record_count: int = 0
    source_determinations: List[str] = []
    has_structural_evidence: bool = False
    structures: Optional[dict[str, Any]] = None
    comparative_metrics: Optional[dict[str, Any]] = None
    warnings: List[dict[str, str]] = []
    failure_reason: Optional[str] = None


class BatchVariantResponse(BaseModel):
    batch_id: str
    timestamp: str
    total_input_rows: int
    unique_identifiers: int
    resolved_count: int
    unresolved_count: int
    evidence_available_count: int
    results: List[BatchRowResult]


ENGINEERED_CONSTRUCT_SUBSTITUTIONS = {
    "E1371Q": {
        "canonical_protein": "E1371Q",
        "canonical_position": 1371,
        "ref_aa": "E",
        "alt_aa": "Q",
        "domain_name": "NBD2",
        "warning_type": "CONSTRUCT_SUBSTITUTION",
        "warning_message": "Engineered catalytic mutation present at target position in 6MSM and 8EIQ (E1371Q).",
        "construct_note": "Construct substitution: E1371Q (engineered catalytic mutation to trap ATP-bound state)",
        "structures_with_substitution": ["6MSM", "8EIQ"],
    }
}


def resolve_construct_substitution_identifier(clean_id: str) -> Optional[str]:
    """Check if identifier explicitly refers to a known engineered construct substitution."""
    upper = clean_id.strip().upper()
    if upper in ("E1371Q", "P.GLU1371GLN", "GLU1371GLN", "P.E1371Q"):
        return "E1371Q"
    norm = normalize_protein_change(clean_id)
    if norm.success and norm.canonical_protein in ENGINEERED_CONSTRUCT_SUBSTITUTIONS:
        return norm.canonical_protein
    return None


def build_construct_substitution_payload(canonical_key: str, db: Session) -> Optional[dict]:
    """Build structural evidence payload for engineered construct substitutions directly from ResidueFeature."""
    config = ENGINEERED_CONSTRUCT_SUBSTITUTIONS.get(canonical_key)
    if not config:
        return None

    pos = config["canonical_position"]
    msc = db.query(MultiStructureComparison).filter(MultiStructureComparison.canonical_position == pos).first()
    rfs = db.query(ResidueFeature).filter(ResidueFeature.canonical_position == pos).all()
    rf_by_struct = {r.structure_id: r for r in rfs}

    rf_5uak = rf_by_struct.get("5UAK")
    rf_6msm = rf_by_struct.get("6MSM")
    rf_8eiq = rf_by_struct.get("8EIQ")

    contact_def = "Heavy-atom distance <= 8.0 Å (includes neighboring protein residues and bound ligands if present)"

    structures_payload = {
        "5UAK": {
            "state": "apo, dephosphorylated",
            "mapping_status": rf_5uak.mapping_status if rf_5uak else None,
            "structure_aa": rf_5uak.structure_amino_acid if rf_5uak else None,
            "feature_status": rf_5uak.eligibility_status if rf_5uak else None,
            "failure_reason": rf_5uak.failure_reason if rf_5uak else None,
            "sasa_angstrom2": rf_5uak.sasa_raw_angstrom2 if rf_5uak else None,
            "rsa_ratio": rf_5uak.rsa_ratio if rf_5uak else None,
            "exposure_category": rf_5uak.exposure_category if rf_5uak else None,
            "total_contacts": rf_5uak.total_contacts if rf_5uak else None,
            "intra_contacts": rf_5uak.intra_chain_contacts if rf_5uak else None,
            "inter_contacts": rf_5uak.inter_chain_contacts if rf_5uak else None,
            "secondary_structure": rf_5uak.secondary_structure if rf_5uak else None,
            "construct_note": None,
            "contact_definition": contact_def,
        },
        "6MSM": {
            "state": "ATP-bound dimer (E1371Q)",
            "mapping_status": rf_6msm.mapping_status if rf_6msm else None,
            "structure_aa": rf_6msm.structure_amino_acid if rf_6msm else None,
            "feature_status": rf_6msm.eligibility_status if rf_6msm else None,
            "failure_reason": rf_6msm.failure_reason if rf_6msm else None,
            "sasa_angstrom2": rf_6msm.sasa_raw_angstrom2 if rf_6msm else None,
            "rsa_ratio": rf_6msm.rsa_ratio if rf_6msm else None,
            "exposure_category": rf_6msm.exposure_category if rf_6msm else None,
            "total_contacts": rf_6msm.total_contacts if rf_6msm else None,
            "intra_contacts": rf_6msm.intra_chain_contacts if rf_6msm else None,
            "inter_contacts": rf_6msm.inter_chain_contacts if rf_6msm else None,
            "secondary_structure": rf_6msm.secondary_structure if rf_6msm else None,
            "construct_note": config["construct_note"] if "6MSM" in config["structures_with_substitution"] else None,
            "contact_definition": contact_def,
        },
        "8EIQ": {
            "state": "ATP-bound dimer (Delta-F508, E1371Q)",
            "mapping_status": rf_8eiq.mapping_status if rf_8eiq else None,
            "structure_aa": rf_8eiq.structure_amino_acid if rf_8eiq else None,
            "feature_status": rf_8eiq.eligibility_status if rf_8eiq else None,
            "failure_reason": rf_8eiq.failure_reason if rf_8eiq else None,
            "sasa_angstrom2": rf_8eiq.sasa_raw_angstrom2 if rf_8eiq else None,
            "rsa_ratio": rf_8eiq.rsa_ratio if rf_8eiq else None,
            "exposure_category": rf_8eiq.exposure_category if rf_8eiq else None,
            "total_contacts": rf_8eiq.total_contacts if rf_8eiq else None,
            "intra_contacts": rf_8eiq.intra_chain_contacts if rf_8eiq else None,
            "inter_contacts": rf_8eiq.inter_chain_contacts if rf_8eiq else None,
            "secondary_structure": rf_8eiq.secondary_structure if rf_8eiq else None,
            "construct_note": config["construct_note"] if "8EIQ" in config["structures_with_substitution"] else None,
            "contact_definition": contact_def,
        },
    }

    all_structures_obs = (msc.observation_pattern == "ALL_THREE_OBSERVED") if msc else True

    comparative_metrics = {
        "observation_pattern": msc.observation_pattern if msc else "ALL_THREE_OBSERVED",
        "construct_variation": msc.construct_variation if msc else "CONSTRUCT_SUBSTITUTION",
        "is_fully_comparable": msc.is_fully_comparable if msc else False,
        "all_structures_observed": all_structures_obs,
        "canonical_identity_match": False,
        "construct_context_caveat": True,
        "environment_comparability": "INCOMPARABLE_CONSTRUCT_SUBSTITUTION",
        "comparison_caveats": [
            "6MSM and 8EIQ carry engineered catalytic substitution E1371Q."
        ],
        "delta_rsa": msc.delta_rsa if msc else None,
        "delta_contacts": msc.delta_contacts if msc else None,
        "sec_struct_consistency": msc.sec_struct_consistency if msc else None,
        "comparison_notes": msc.comparison_notes if msc else None,
    }

    warnings = [
        {
            "type": config["warning_type"],
            "message": config["warning_message"],
        }
    ]

    return {
        "canonical_protein": config["canonical_protein"],
        "canonical_position": pos,
        "ref_aa": config["ref_aa"],
        "alt_aa": config["alt_aa"],
        "domain_name": msc.domain_name if msc else config["domain_name"],
        "source_associations": {
            "cftr2_count": 0,
            "cftr2_records": [],
            "clinvar_count": 0,
            "clinvar_records": [],
        },
        "structures": structures_payload,
        "comparative_metrics": comparative_metrics,
        "warnings": warnings,
        "guardrail_notice": (
            "This endpoint provides structural evidence integration observations only. "
            "It does NOT predict clinical pathogenicity, compute weighted damage scores, "
            "or synthesize manufactured consensus classifications."
        ),
    }


def find_variant_records(clean_id: str, db: Session) -> list[VariantStructuralEvidence]:
    """Retrieve VariantStructuralEvidence records using canonical, legacy, or HGVS lookups."""
    clean_id = clean_id.strip()
    if not clean_id:
        return []

    # 1. Exact match on canonical_protein, source_record_id, or source_variant_name
    query = db.query(VariantStructuralEvidence).filter(
        (VariantStructuralEvidence.canonical_protein == clean_id)
        | (VariantStructuralEvidence.source_record_id == clean_id)
        | (VariantStructuralEvidence.source_variant_name == clean_id)
    )

    if clean_id.isdigit():
        pos = int(clean_id)
        query = db.query(VariantStructuralEvidence).filter(
            (VariantStructuralEvidence.canonical_position == pos)
            | (VariantStructuralEvidence.source_record_id == clean_id)
        )

    records = query.all()

    # Fallback 1: Normalized protein notation (e.g. 'p.Gly551Asp' -> 'G551D', 'p.Phe508del' -> 'F508del')
    if not records:
        norm = normalize_protein_change(clean_id)
        if norm.success and norm.canonical_protein:
            records = db.query(VariantStructuralEvidence).filter(
                (VariantStructuralEvidence.canonical_protein == norm.canonical_protein)
                | (VariantStructuralEvidence.source_record_id == norm.canonical_protein)
            ).all()

    # Fallback 2: Raw protein change substring
    if not records:
        records = db.query(VariantStructuralEvidence).filter(
            VariantStructuralEvidence.protein_change_raw.ilike(f"%{clean_id}%")
        ).all()

    return records


def format_variant_evidence_payload(records: list[VariantStructuralEvidence]) -> dict:
    """Format matching evidence records into structured API payload."""
    if not records:
        return {}

    first = records[0]

    cftr2_records = [
        {
            "legacy_name": r.source_record_id,
            "determination_verbatim": r.source_determination,
            "alleles_reported": r.cftr2_alleles_reported,
            "validation_flag": r.validation_flag,
            "eligibility_status": r.eligibility_status,
            "ineligibility_reason": r.ineligibility_reason,
        }
        for r in records
        if r.source_db == "CFTR2"
    ]

    clinvar_records = [
        {
            "variation_id": r.source_record_id,
            "vcv_accession": r.clinvar_vcv_accession,
            "name": r.source_variant_name,
            "classification_verbatim": r.source_determination,
            "review_status": r.clinvar_review_status,
            "validation_flag": r.validation_flag,
            "eligibility_status": r.eligibility_status,
            "ineligibility_reason": r.ineligibility_reason,
        }
        for r in records
        if r.source_db == "ClinVar"
    ]

    warnings = []
    if any(r.eligibility_status == "REFERENCE_MISMATCH" for r in records):
        warnings.append({
            "type": "REFERENCE_MISMATCH",
            "message": "Canonical reference mismatch: reported reference residue does not match canonical CFTR P13569.",
        })
    if first.construct_variation == "CONSTRUCT_SUBSTITUTION" or first.canonical_position == 1371:
        warnings.append({
            "type": "CONSTRUCT_SUBSTITUTION",
            "message": "Engineered catalytic mutation present at target position in 6MSM and 8EIQ (E1371Q).",
        })
    if first.construct_variation == "CONSTRUCT_DELETION" or (
        first.canonical_position == 508 and first.failure_reason_8eiq == "CONSTRUCT_DELETION"
    ):
        warnings.append({
            "type": "CONSTRUCT_DELETION",
            "message": "Residue authentically deleted in 8EIQ (Delta-F508); coordinates absent.",
        })
    if first.observation_pattern == "ALL_UNRESOLVED" or first.domain_name == "R-domain":
        warnings.append({
            "type": "UNRESOLVED_COORDINATES",
            "message": f"Coordinates unmodeled in experimental models at position {first.canonical_position} ({first.domain_name}) due to disorder.",
        })
    if any(r.eligibility_status == "SYNONYMOUS" for r in records):
        warnings.append({
            "type": "SYNONYMOUS",
            "message": "Synonymous alteration at nucleotide level; no structural amino-acid substitution.",
        })
    if first.canonical_position == 551:
        warnings.append({
            "type": "DIMER_INTERFACE_LIGANDS",
            "message": "Residue G551 is in NBD1 LSGGQ signature motif and forms direct contacts with bound ATP1504, Mg1502, and engineered mutant Gln1371 in 6MSM and 8EIQ.",
        })

    all_structures_obs = (first.observation_pattern == "ALL_THREE_OBSERVED")
    can_id_match = (first.construct_variation == "CANONICAL_WT")
    has_construct_caveat = False
    comp_caveats = []

    if first.canonical_position == 508 or first.construct_variation == "CONSTRUCT_DELETION":
        has_construct_caveat = True
        comp_caveats.append("8EIQ construct carries authentic patient deletion Delta-F508; coordinates absent.")
        env_comp = "INCOMPARABLE_CONSTRUCT_DELETION"
    elif first.canonical_position == 1371 or first.construct_variation == "CONSTRUCT_SUBSTITUTION":
        has_construct_caveat = True
        comp_caveats.append("6MSM and 8EIQ carry engineered catalytic substitution E1371Q.")
        env_comp = "INCOMPARABLE_CONSTRUCT_SUBSTITUTION"
    else:
        if all_structures_obs:
            comp_caveats.append("5UAK represents unphosphorylated apo state, whereas 6MSM and 8EIQ represent phosphorylated ATP/Mg-bound dimers.")
            comp_caveats.append("6MSM and 8EIQ carry engineered catalytic substitution E1371Q to trap the ATP-bound dimer.")
            if first.canonical_position == 551:
                has_construct_caveat = True
                comp_caveats.append("Residue G551 directly contacts bound ATP1504, Mg1502, and engineered Gln1371 in 6MSM and 8EIQ.")
            env_comp = "QUALIFIED_DIFFERENT_CONFORMATIONAL_STATES"
        elif not all_structures_obs:
            env_comp = "PARTIALLY_OBSERVED"
            comp_caveats.append("Residue coordinates unmodeled in one or more conformational states.")
        else:
            env_comp = "CONSTRUCT_DIVERGENT"

    contact_def = "Heavy-atom distance <= 8.0 Å (includes neighboring protein residues and bound ligands if present)"

    return {
        "canonical_protein": first.canonical_protein,
        "canonical_position": first.canonical_position,
        "ref_aa": first.ref_aa,
        "alt_aa": first.alt_aa,
        "domain_name": first.domain_name,
        "source_associations": {
            "cftr2_count": len(cftr2_records),
            "cftr2_records": cftr2_records,
            "clinvar_count": len(clinvar_records),
            "clinvar_records": clinvar_records,
        },
        "structures": {
            "5UAK": {
                "state": "apo, dephosphorylated",
                "mapping_status": first.mapping_status_5uak,
                "structure_aa": first.structure_aa_5uak,
                "feature_status": first.feature_status_5uak,
                "failure_reason": first.failure_reason_5uak,
                "sasa_angstrom2": first.sasa_5uak,
                "rsa_ratio": first.rsa_5uak,
                "exposure_category": first.exposure_cat_5uak,
                "total_contacts": first.contacts_5uak,
                "intra_contacts": first.intra_contacts_5uak,
                "inter_contacts": first.inter_contacts_5uak,
                "secondary_structure": first.sec_struct_5uak,
                "construct_note": first.construct_note_5uak,
                "contact_definition": contact_def,
            },
            "6MSM": {
                "state": "ATP-bound dimer (E1371Q)",
                "mapping_status": first.mapping_status_6msm,
                "structure_aa": first.structure_aa_6msm,
                "feature_status": first.feature_status_6msm,
                "failure_reason": first.failure_reason_6msm,
                "sasa_angstrom2": first.sasa_6msm,
                "rsa_ratio": first.rsa_6msm,
                "exposure_category": first.exposure_cat_6msm,
                "total_contacts": first.contacts_6msm,
                "intra_contacts": first.intra_contacts_6msm,
                "inter_contacts": first.inter_contacts_6msm,
                "secondary_structure": first.sec_struct_6msm,
                "construct_note": first.construct_note_6msm,
                "contact_definition": contact_def,
            },
            "8EIQ": {
                "state": "ATP-bound dimer (Delta-F508, E1371Q)",
                "mapping_status": first.mapping_status_8eiq,
                "structure_aa": first.structure_aa_8eiq,
                "feature_status": first.feature_status_8eiq,
                "failure_reason": first.failure_reason_8eiq,
                "sasa_angstrom2": first.sasa_8eiq,
                "rsa_ratio": first.rsa_8eiq,
                "exposure_category": first.exposure_cat_8eiq,
                "total_contacts": first.contacts_8eiq,
                "intra_contacts": first.intra_contacts_8eiq,
                "inter_contacts": first.inter_contacts_8eiq,
                "secondary_structure": first.sec_struct_8eiq,
                "construct_note": first.construct_note_8eiq,
                "contact_definition": contact_def,
            },
        },
        "comparative_metrics": {
            "observation_pattern": first.observation_pattern,
            "construct_variation": first.construct_variation,
            "is_fully_comparable": first.is_fully_comparable,
            "all_structures_observed": all_structures_obs,
            "canonical_identity_match": can_id_match,
            "construct_context_caveat": has_construct_caveat,
            "environment_comparability": env_comp,
            "comparison_caveats": comp_caveats,
            "delta_rsa": first.delta_rsa,
            "delta_contacts": first.delta_contacts,
            "sec_struct_consistency": first.sec_struct_consistency,
            "comparison_notes": first.comparison_notes,
        },
        "warnings": warnings,
        "guardrail_notice": (
            "This endpoint provides structural evidence integration observations only. "
            "It does NOT predict clinical pathogenicity, compute weighted damage scores, "
            "or synthesize manufactured consensus classifications."
        ),
    }


@app.get("/variants/{identifier}", tags=["Variants"])
async def get_variant_evidence(
    identifier: str,
    db: Session = Depends(get_db),
):
    """Retrieve comprehensive structural evidence and source classifications for a variant."""
    clean_id = identifier.strip()

    # Check for engineered construct substitutions first (e.g. E1371Q)
    construct_key = resolve_construct_substitution_identifier(clean_id)
    if construct_key:
        payload = build_construct_substitution_payload(construct_key, db)
        if payload:
            return payload

    records = find_variant_records(clean_id, db)
    if not records:
        raise HTTPException(
            status_code=404,
            detail=f"Variant '{identifier}' not found in CFTR-StructMap evidence database.",
        )
    return format_variant_evidence_payload(records)


@app.post("/variants/batch", response_model=BatchVariantResponse, tags=["Variants"])
async def batch_variant_evidence(
    request: BatchVariantRequest,
    db: Session = Depends(get_db),
):
    """Analyze a batch of variant identifiers (up to 500 rows).

    Preserves input order, row numbers, and original input text.
    Deduplicates database lookups and returns traceable structural evidence.
    """
    if not request.identifiers:
        raise HTTPException(status_code=400, detail="Identifiers list cannot be empty.")

    if len(request.identifiers) > 500:
        raise HTTPException(
            status_code=400,
            detail=f"Batch size of {len(request.identifiers)} exceeds maximum limit of 500 identifiers.",
        )

    # Cache lookups for unique non-empty identifiers to avoid redundant DB queries
    lookup_cache: dict[str, tuple[list[VariantStructuralEvidence], Optional[dict], Optional[str], str]] = {}
    unique_identifiers_set = set()

    for raw in request.identifiers:
        cleaned = raw.strip()
        if cleaned:
            unique_identifiers_set.add(cleaned)

    # Perform lookups once per unique identifier
    for clean_id in unique_identifiers_set:
        construct_key = resolve_construct_substitution_identifier(clean_id)
        if construct_key:
            payload = build_construct_substitution_payload(construct_key, db)
            if payload:
                lookup_cache[clean_id] = ([], payload, payload.get("canonical_protein"), "RESOLVED")
                continue

        records = find_variant_records(clean_id, db)
        if records:
            payload = format_variant_evidence_payload(records)
            lookup_cache[clean_id] = (records, payload, payload.get("canonical_protein"), "RESOLVED")
        else:
            # Check syntax normalization
            norm = normalize_protein_change(clean_id)
            if norm.success and norm.canonical_protein:
                lookup_cache[clean_id] = (
                    [],
                    None,
                    norm.canonical_protein,
                    "NOT_FOUND",
                )
            else:
                lookup_cache[clean_id] = (
                    [],
                    None,
                    None,
                    "UNSUPPORTED_FORMAT",
                )

    results: list[BatchRowResult] = []
    resolved_count = 0
    unresolved_count = 0
    evidence_available_count = 0

    for idx, raw in enumerate(request.identifiers, start=1):
        cleaned = raw.strip()
        if not cleaned:
            results.append(
                BatchRowResult(
                    row_number=idx,
                    original_input=raw,
                    status="EMPTY",
                    failure_reason="Blank input row.",
                )
            )
            unresolved_count += 1
            continue

        records, payload, norm_key, status = lookup_cache.get(cleaned, ([], None, None, "UNSUPPORTED_FORMAT"))

        if status == "RESOLVED" and payload:
            determinations = (
                list(dict.fromkeys([r.source_determination for r in records if r.source_determination]))
                if records else []
            )

            results.append(
                BatchRowResult(
                    row_number=idx,
                    original_input=raw,
                    normalized_identifier=payload.get("canonical_protein") or norm_key,
                    status="RESOLVED",
                    canonical_protein=payload.get("canonical_protein"),
                    canonical_position=payload.get("canonical_position"),
                    ref_aa=payload.get("ref_aa"),
                    alt_aa=payload.get("alt_aa"),
                    domain_name=payload.get("domain_name"),
                    source_record_count=len(records),
                    source_determinations=determinations,
                    has_structural_evidence=True,
                    structures=payload.get("structures"),
                    comparative_metrics=payload.get("comparative_metrics"),
                    warnings=payload.get("warnings", []),
                )
            )
            resolved_count += 1
            evidence_available_count += 1
        elif status == "NOT_FOUND":
            results.append(
                BatchRowResult(
                    row_number=idx,
                    original_input=raw,
                    normalized_identifier=norm_key,
                    status="NOT_FOUND",
                    failure_reason=f"Recognized notation '{norm_key}' but no records exist in CFTR2 or ClinVar database.",
                )
            )
            unresolved_count += 1
        else:
            results.append(
                BatchRowResult(
                    row_number=idx,
                    original_input=raw,
                    status="UNSUPPORTED_FORMAT",
                    failure_reason=f"Unable to resolve '{cleaned}' as a supported CFTR variant identifier.",
                )
            )
            unresolved_count += 1

    return BatchVariantResponse(
        batch_id=f"batch-{uuid.uuid4().hex[:12]}",
        timestamp=datetime.now(timezone.utc).isoformat(),
        total_input_rows=len(request.identifiers),
        unique_identifiers=len(unique_identifiers_set),
        resolved_count=resolved_count,
        unresolved_count=unresolved_count,
        evidence_available_count=evidence_available_count,
        results=results,
    )


# ==============================================================================
# Experimental Macromolecular Structure Delivery & Coordinate Mapping
# ==============================================================================

SUPPORTED_STRUCTURES = {
    "5UAK": "5UAK.cif",
    "6MSM": "6MSM.cif",
    "8EIQ": "8EIQ.cif",
    "1AV5": "1AV5.cif",
    "7JQ4": "7jq4.cif",
}


@app.get("/structures/{structure_id}/coordinates", tags=["Structures"])
async def get_structure_coordinates(structure_id: str):
    """Securely stream mmCIF coordinate files for supported experimental structures."""
    sid = structure_id.strip().upper()
    if sid not in SUPPORTED_STRUCTURES:
        raise HTTPException(
            status_code=404,
            detail=f"Structure '{structure_id}' is not in the supported structure allowlist. Supported: {list(SUPPORTED_STRUCTURES.keys())}.",
        )

    filename = SUPPORTED_STRUCTURES[sid]
    raw_dir = Path("data/raw").resolve()
    target_path = (raw_dir / filename).resolve()

    # Prevent path traversal and verify local file existence
    if not str(target_path).startswith(str(raw_dir)) or not target_path.is_file():
        raise HTTPException(
            status_code=404,
            detail=f"Coordinate file for structure '{sid}' not found on server.",
        )

    return FileResponse(
        path=target_path,
        media_type="chemical/x-mmcif",
        filename=filename,
        headers={"Cache-Control": "public, max-age=86400"},
    )


@app.get("/structures/{structure_id}/residues/{canonical_position}", tags=["Structures"])
async def get_structure_residue_mapping(
    structure_id: str,
    canonical_position: int,
    db: Session = Depends(get_db),
):
    """Retrieve structure-specific chain and residue numbering for a canonical CFTR position."""
    sid = structure_id.strip().upper()
    if sid not in SUPPORTED_STRUCTURES:
        raise HTTPException(
            status_code=404,
            detail=f"Structure '{structure_id}' is not in the supported structure allowlist.",
        )

    mapping = (
        db.query(ResidueMapping)
        .filter(
            ResidueMapping.structure_id == sid,
            ResidueMapping.canonical_position == canonical_position,
        )
        .first()
    )

    if not mapping:
        raise HTTPException(
            status_code=404,
            detail=f"Residue position {canonical_position} not mapped in structure '{sid}'.",
        )

    return {
        "structure_id": sid,
        "canonical_position": canonical_position,
        "canonical_amino_acid": mapping.canonical_amino_acid,
        "entity_id": mapping.entity_id,
        "label_asym_id": mapping.label_asym_id,
        "auth_asym_id": mapping.auth_asym_id,
        "label_seq_id": mapping.label_seq_id,
        "auth_seq_id": mapping.auth_seq_id,
        "insertion_code": mapping.insertion_code,
        "structure_amino_acid": mapping.structure_amino_acid,
        "has_observed_residue": mapping.has_observed_residue,
        "has_ca_atom": mapping.has_ca_atom,
        "reference_identity_match": mapping.reference_identity_match,
        "mapping_status": mapping.mapping_status,
        "mapping_warning": mapping.mapping_warning,
    }


@app.get("/structures/{structure_id}/metadata", tags=["Structures"])
async def get_structure_metadata(
    structure_id: str,
    db: Session = Depends(get_db),
):
    """Retrieve metadata and polymer chains for an experimental structure."""
    sid = structure_id.strip().upper()
    if sid not in SUPPORTED_STRUCTURES:
        raise HTTPException(
            status_code=404,
            detail=f"Structure '{structure_id}' is not supported.",
        )

    struct = db.query(Structure).filter(Structure.id == sid).first()
    if not struct:
        raise HTTPException(
            status_code=404,
            detail=f"Structure '{sid}' metadata not found in database.",
        )

    chains = (
        db.query(StructuralChain)
        .filter(StructuralChain.structure_id == sid)
        .all()
    )

    return {
        "id": struct.id,
        "title": struct.title,
        "experimental_method": struct.experimental_method,
        "resolution": struct.resolution,
        "structure_type": struct.structure_type,
        "has_cftr": struct.has_cftr,
        "deposition_date": struct.deposition_date,
        "chains": [
            {
                "label_asym_id": c.label_asym_id,
                "auth_asym_id": c.auth_asym_id,
                "is_cftr": c.is_cftr,
                "entity_description": c.entity_description,
                "polymer_seq_len": c.polymer_seq_len,
                "canonical_coverage": c.canonical_coverage,
                "observed_coverage": c.observed_coverage,
            }
            for c in chains
        ],
    }


