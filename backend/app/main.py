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


@app.get("/health", tags=["Health"])
async def health_check():
    """Health check endpoint for container orchestration and uptime monitoring."""
    return {
        "status": "healthy",
        "service": settings.PROJECT_NAME,
        "version": settings.VERSION,
    }


from fastapi import Depends, HTTPException, Query
from sqlalchemy.orm import Session
from backend.app.db.models import VariantStructuralEvidence, VariantStructuralSummary
from backend.app.db.session import get_db


@app.get("/variants", tags=["Variants"])
async def list_variants(
    position: int | None = Query(None, description="Filter by 1-based canonical CFTR position (1..1480)"),
    domain: str | None = Query(None, description="Filter by canonical domain name (e.g. NBD1, TMD1)"),
    eligibility: str | None = Query(None, description="Filter by eligibility status"),
    source_db: str | None = Query(None, description="Filter by source database ('CFTR2' or 'ClinVar')"),
    limit: int = Query(50, ge=1, le=500, description="Page limit"),
    offset: int = Query(0, ge=0, description="Page offset"),
    db: Session = Depends(get_db),
):
    """Retrieve paginated variant structural evidence records with optional filtering."""
    query = db.query(VariantStructuralEvidence)
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


@app.get("/variants/{identifier}", tags=["Variants"])
async def get_variant_evidence(
    identifier: str,
    db: Session = Depends(get_db),
):
    """Retrieve comprehensive structural evidence and source classifications for a variant.

    Supports lookup by:
    - Normalized protein notation (e.g. 'G551D', 'F508del', 'M470V')
    - Canonical CFTR residue position (e.g. '508', '551', '1371')
    - ClinVar VariationID (e.g. '7120')
    - CFTR2 legacy name (e.g. 'F508del')
    """
    clean_id = identifier.strip()

    # Query evidence matching identifier as canonical_protein, source_record_id, or position
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
    if not records:
        raise HTTPException(
            status_code=404,
            detail=f"Variant '{identifier}' not found in CFTR-StructMap evidence database.",
        )

    first = records[0]

    # Collect source records
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

    # Collect active warnings
    warnings = []
    if any(r.eligibility_status == "REFERENCE_MISMATCH" for r in records):
        warnings.append({
            "type": "REFERENCE_MISMATCH",
            "message": "Canonical reference mismatch: reported reference residue does not match canonical CFTR P13569.",
        })
    if first.construct_variation == "CONSTRUCT_SUBSTITUTION":
        warnings.append({
            "type": "CONSTRUCT_SUBSTITUTION",
            "message": "Engineered catalytic mutation present at target position in 6MSM and 8EIQ (E1371Q).",
        })
    if first.construct_variation == "CONSTRUCT_DELETION" or (first.canonical_position == 508 and first.failure_reason_8eiq == "CONSTRUCT_DELETION"):
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
            },
        },
        "comparative_metrics": {
            "observation_pattern": first.observation_pattern,
            "construct_variation": first.construct_variation,
            "is_fully_comparable": first.is_fully_comparable,
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
