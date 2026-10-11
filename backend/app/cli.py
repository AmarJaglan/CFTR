"""Command line interface for reproducible CFTR-StructMap data ingestion."""

import argparse
from pathlib import Path
import sys

from backend.app.db.models import (
    CFTR2Variant,
    ClinVarVariant,
    DataSource,
    MultiStructureComparison,
    ResidueFeature,
    ResidueMapping,
    StructuralChain,
    Structure,
    VariantStructuralEvidence,
    VariantStructuralSummary,
)
from backend.app.db.session import get_engine, get_session_factory, init_db
from backend.app.ingestion.cftr2 import ingest_cftr2
from backend.app.ingestion.clinvar import ingest_clinvar
from backend.app.scientific.features import (
    export_feature_audit_csvs,
    generate_feature_audit_markdown,
    persist_feature_results,
    run_all_feature_extractions,
)
from backend.app.scientific.mapping import (
    export_audit_csvs,
    generate_audit_markdown,
    map_all_structures,
    persist_mapping_results,
)
from backend.app.scientific.variant_evidence import (
    build_variant_evidence,
    export_variant_audit_artifacts,
    persist_variant_evidence,
)

DEFAULT_CFTR2_FILE = Path("data/raw/CFTR2_30January2026.xlsx")
DEFAULT_CLINVAR_FILE = Path("data/raw/clinvar_download_20261009_143804.tsv")


def cmd_import_cftr2(args):
    file_path = Path(args.file)
    if not file_path.exists():
        print(f"Error: CFTR2 source file not found at: {file_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Starting CFTR2 ingestion from: {file_path}")
    engine = get_engine(args.db)
    init_db(engine)
    session_factory = get_session_factory(engine)
    
    with session_factory() as session:
        summary = ingest_cftr2(file_path, db_session=session, dry_run=args.dry_run)
    
    print("\n" + summary.to_text_report())


def cmd_import_clinvar(args):
    file_path = Path(args.file)
    if not file_path.exists():
        print(f"Error: ClinVar source file not found at: {file_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Starting ClinVar ingestion from: {file_path}")
    engine = get_engine(args.db)
    init_db(engine)
    session_factory = get_session_factory(engine)
    
    with session_factory() as session:
        summary = ingest_clinvar(file_path, db_session=session, dry_run=args.dry_run)
    
    print("\n" + summary.to_text_report())


def cmd_import_all(args):
    print("=== Running Complete CFTR-StructMap Data Ingestion ===")
    args.file = args.cftr2_file
    cmd_import_cftr2(args)
    print("\n" + "=" * 50 + "\n")
    args.file = args.clinvar_file
    cmd_import_clinvar(args)


def cmd_map_structures(args):
    structures_dir = Path(args.structures_dir)
    output_dir = Path(args.output_dir)

    target_names = [
        "1AV5.cif",
        "5UAK(1).cif",
        "6MSM(1).cif",
        "7jq4.cif",
        "8EIQ(1).cif",
    ]

    cif_paths = []
    missing_files = []
    for name in target_names:
        p = structures_dir / name
        if not p.exists():
            # Try fallback without (1)
            alt_name = name.replace("(1)", "")
            p_alt = structures_dir / alt_name
            if p_alt.exists():
                p = p_alt
            else:
                missing_files.append(name)
                continue
        cif_paths.append(p)

    if missing_files:
        print(f"Warning: The following expected structure files were not found in {structures_dir}: {missing_files}", file=sys.stderr)

    if not cif_paths:
        print(f"Error: No valid structure files found in {structures_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"=== Mapping {len(cif_paths)} Experimental Structures ===")
    for p in cif_paths:
        print(f"  - {p.name}")

    results = map_all_structures(cif_paths)

    if not args.dry_run:
        engine = get_engine(args.db)
        init_db(engine)
        session_factory = get_session_factory(engine)
        with session_factory() as session:
            persist_mapping_results(results, session)
        print("Successfully persisted structural mappings to database.")

    # Export audit artifacts
    csv_paths = export_audit_csvs(results, output_dir)
    md_path = generate_audit_markdown(results, output_dir / "structural_mapping_summary.md")
    print(f"Exported audit CSVs to: {output_dir}")
    print(f"Generated audit summary report: {md_path}")

    # Print summary table
    print("\n" + "=" * 70)
    print("STRUCTURE MAPPING SUMMARY")
    print("=" * 70)
    for r in results:
        print(f"Structure: {r.structure_id} ({r.source_filename}) | Type: {r.structure_type} | CFTR: {r.has_cftr} | Method: {r.experimental_method or 'N/A'} ({r.resolution or 'N/A'} Å)")
        for cs in r.chain_summaries:
            cftr_flag = "CFTR" if cs.is_cftr else "Non-CFTR"
            cov = f"{cs.canonical_coverage * 100:.1f}%"
            print(f"  Chain {cs.label_asym_id} (auth {cs.auth_asym_id}) [{cftr_flag}]: Len={cs.polymer_seq_len}, Cov={cov}, Obs={cs.observed_residues_count}, CA={cs.ca_atoms_count}, Unobs={cs.unresolved_count}, NotCovered={cs.not_covered_count}, Mismatches={cs.reference_mismatch_count}")
    print("=" * 70 + "\n")


def cmd_extract_features(args):
    structures_dir = Path(args.structures_dir)
    output_dir = Path(args.output_dir)
    contact_threshold = args.contact_threshold

    if not structures_dir.exists():
        print(f"Error: Structures directory not found: {structures_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"=== Extracting Structural Features (Contact Threshold: {contact_threshold} Å) ===")
    features_by_structure, comparisons, all_warnings = run_all_feature_extractions(
        structures_dir=structures_dir,
        contact_threshold=contact_threshold,
    )

    if not args.dry_run:
        engine = get_engine(args.db)
        init_db(engine)
        session_factory = get_session_factory(engine)
        with session_factory() as session:
            persist_feature_results(features_by_structure, comparisons, session)
        print("Successfully persisted structural features to database.")

    # Export audit artifacts
    csv_paths = export_feature_audit_csvs(
        features_by_structure=features_by_structure,
        comparisons=comparisons,
        all_warnings=all_warnings,
        output_dir=output_dir,
    )
    md_path = generate_feature_audit_markdown(
        features_by_structure=features_by_structure,
        comparisons=comparisons,
        output_path=output_dir / "structural_feature_summary.md",
    )
    print(f"Exported feature audit CSVs to: {output_dir}")
    print(f"Generated feature audit summary report: {md_path}")

    # Print summary table
    print("\n" + "=" * 75)
    print("STRUCTURAL FEATURE EXTRACTION SUMMARY")
    print("=" * 75)
    for sid, feats in features_by_structure.items():
        fname = feats[0].source_filename if feats else "N/A"
        el = [r for r in feats if r.eligibility_status == "ELIGIBLE"]
        un = [r for r in feats if r.eligibility_status == "INELIGIBLE_UNRESOLVED"]
        nc = [r for r in feats if r.eligibility_status == "INELIGIBLE_NOT_COVERED"]
        rsas = [r.rsa_ratio for r in el if r.rsa_ratio is not None]
        cnts = [r.total_contacts for r in el if r.total_contacts is not None]
        mean_rsa = f"{sum(rsas)/len(rsas):.3f}" if rsas else "N/A"
        mean_cnt = f"{sum(cnts)/len(cnts):.1f}" if cnts else "N/A"
        print(f"Structure {sid} ({fname}): Total={len(feats)}, Eligible={len(el)}, Unresolved={len(un)}, NotCovered={len(nc)}, Mean RSA={mean_rsa}, Mean Contacts={mean_cnt}")

    comp_pat_counts = {}
    for c in comparisons:
        comp_pat_counts[c.observation_pattern] = comp_pat_counts.get(c.observation_pattern, 0) + 1
    print("\nMulti-Structure Observation Patterns (1,480 positions):")
    for pat, count in sorted(comp_pat_counts.items()):
        print(f"  - {pat}: {count} ({count/1480*100:.1f}%)")
    print("=" * 75 + "\n")


def cmd_integrate_variants(args):
    output_dir = Path(args.output_dir)

    print("=== Integrating Variant-Level Structural Evidence (Phase 5) ===")
    engine = get_engine(args.db)
    init_db(engine)
    session_factory = get_session_factory(engine)

    with session_factory() as session:
        evidence_list, summary_list, comparisons_list, warnings_list = build_variant_evidence(session)

        if not args.dry_run:
            persisted_ev, persisted_sum = persist_variant_evidence(evidence_list, summary_list, session)
            print(f"Successfully persisted {persisted_ev} variant structural evidence records.")
            print(f"Successfully persisted {persisted_sum} grouped variant structural summary records.")
        else:
            print(f"[Dry-run] Built {len(evidence_list)} evidence records and {len(summary_list)} summary records (no DB write).")

    # Export audit CSVs and markdown summary
    csv_paths = export_variant_audit_artifacts(
        evidence_list=evidence_list,
        summary_list=summary_list,
        comparisons_list=comparisons_list,
        warnings_list=warnings_list,
        output_dir=output_dir,
    )
    print(f"Exported audit CSVs to: {output_dir}")
    print(f"Generated audit summary report: {csv_paths['markdown']}")

    # Print summary breakdown
    cftr2_cnt = sum(1 for e in evidence_list if e.source_db == "CFTR2")
    cv_cnt = sum(1 for e in evidence_list if e.source_db == "ClinVar")
    el_cnt = sum(1 for e in evidence_list if e.eligibility_status == "ELIGIBLE_MISSENSE")
    syn_cnt = sum(1 for e in evidence_list if e.eligibility_status == "SYNONYMOUS")
    indel_cnt = sum(1 for e in evidence_list if e.eligibility_status == "IN_FRAME_DELETION")
    ref_mismatch_cnt = sum(1 for e in evidence_list if e.eligibility_status == "REFERENCE_MISMATCH")

    print("\n" + "=" * 75)
    print("VARIANT STRUCTURAL EVIDENCE INTEGRATION SUMMARY")
    print("=" * 75)
    print(f"Total Source Records: {len(evidence_list)} (CFTR2: {cftr2_cnt}, ClinVar VCV: {cv_cnt})")
    print(f"Unique Normalized Protein Changes: {len(summary_list)}")
    print(f"Eligible Missense Substitutions: {el_cnt}")
    print(f"Synonymous Variants (Separately Handled): {syn_cnt}")
    print(f"In-Frame Deletions (e.g. F508del): {indel_cnt}")
    print(f"Canonical Reference Mismatches (e.g. M470V): {ref_mismatch_cnt}")
    print(f"Pairwise CFTR2/ClinVar Overlaps: {len(comparisons_list)}")
    print(f"Total Audit Warnings Logged: {len(warnings_list)}")
    print("=" * 75 + "\n")


def cmd_status(args):
    engine = get_engine(args.db)
    init_db(engine)
    session_factory = get_session_factory(engine)
    with session_factory() as session:
        sources = session.query(DataSource).all()
        cftr2_count = session.query(CFTR2Variant).count()
        clinvar_count = session.query(ClinVarVariant).count()
        struct_count = session.query(Structure).count()
        chain_count = session.query(StructuralChain).count()
        cftr_chain_count = session.query(StructuralChain).filter_by(is_cftr=True).count()
        res_count = session.query(ResidueMapping).count()
        feat_count = session.query(ResidueFeature).count()
        comp_count = session.query(MultiStructureComparison).count()
        ev_count = session.query(VariantStructuralEvidence).count()
        sum_count = session.query(VariantStructuralSummary).count()

        print("=== CFTR-StructMap Database Status ===")
        print(f"Database target: {args.db or 'default'}")
        print(f"CFTR2 Variants: {cftr2_count}")
        print(f"ClinVar Variants: {clinvar_count}")
        print(f"Ingested Sources: {len(sources)}")
        for src in sources:
            print(f"  - [{src.source_name}] {src.filename} (SHA256: {src.file_sha256[:12]}..., Records: {src.record_count})")
        print(f"Structures: {struct_count}")
        print(f"Structural Chains: {chain_count} (CFTR chains: {cftr_chain_count})")
        print(f"Residue Mappings: {res_count}")
        print(f"Residue Features: {feat_count}")
        print(f"Multi-Structure Comparisons: {comp_count}")
        print(f"Variant Structural Evidence: {ev_count}")
        print(f"Variant Structural Summaries: {sum_count}")


def main():
    parser = argparse.ArgumentParser(
        prog="cftr-structmap-cli",
        description="CFTR-StructMap scientific data ingestion, structural mapping, and provenance CLI",
    )
    parser.add_argument("--db", type=str, default=None, help="SQLite database URL or path")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # import-cftr2
    p_cftr2 = subparsers.add_parser("import-cftr2", help="Ingest CFTR2 dataset")
    p_cftr2.add_argument("--file", type=str, default=str(DEFAULT_CFTR2_FILE), help="Path to CFTR2 Excel file")
    p_cftr2.add_argument("--dry-run", action="store_true", help="Parse and validate without saving to database")
    p_cftr2.set_defaults(func=cmd_import_cftr2)

    # import-clinvar
    p_clinvar = subparsers.add_parser("import-clinvar", help="Ingest ClinVar dataset")
    p_clinvar.add_argument("--file", type=str, default=str(DEFAULT_CLINVAR_FILE), help="Path to ClinVar TSV file")
    p_clinvar.add_argument("--dry-run", action="store_true", help="Parse and validate without saving to database")
    p_clinvar.set_defaults(func=cmd_import_clinvar)

    # import-all
    p_all = subparsers.add_parser("import-all", help="Ingest both CFTR2 and ClinVar datasets")
    p_all.add_argument("--cftr2-file", type=str, default=str(DEFAULT_CFTR2_FILE), help="Path to CFTR2 Excel file")
    p_all.add_argument("--clinvar-file", type=str, default=str(DEFAULT_CLINVAR_FILE), help="Path to ClinVar TSV file")
    p_all.add_argument("--dry-run", action="store_true", help="Parse and validate without saving to database")
    p_all.set_defaults(func=cmd_import_all)

    # map-structures
    p_map = subparsers.add_parser("map-structures", help="Parse mmCIF structures and map to canonical CFTR")
    p_map.add_argument("--structures-dir", type=str, default="data/raw", help="Directory containing .cif structure files")
    p_map.add_argument("--output-dir", type=str, default="docs/audits", help="Directory for audit CSV and markdown outputs")
    p_map.add_argument("--dry-run", action="store_true", help="Map structures without saving to database")
    p_map.set_defaults(func=cmd_map_structures)

    # extract-features
    p_feat = subparsers.add_parser("extract-features", help="Extract residue-level biophysical features from validated structures")
    p_feat.add_argument("--structures-dir", type=str, default="data/raw", help="Directory containing .cif structure files")
    p_feat.add_argument("--output-dir", type=str, default="docs/audits", help="Directory for audit CSV and markdown outputs")
    p_feat.add_argument("--contact-threshold", type=float, default=8.0, help="Heavy-atom distance threshold in Angstroms (default 8.0)")
    p_feat.add_argument("--dry-run", action="store_true", help="Extract features without saving to database")
    p_feat.set_defaults(func=cmd_extract_features)

    # integrate-variants
    p_int = subparsers.add_parser("integrate-variants", help="Integrate variant records with structural features (Phase 5)")
    p_int.add_argument("--output-dir", type=str, default="docs/audits", help="Directory for audit CSV and markdown outputs")
    p_int.add_argument("--dry-run", action="store_true", help="Integrate variants without saving to database")
    p_int.set_defaults(func=cmd_integrate_variants)

    # status
    p_status = subparsers.add_parser("status", help="Show database record counts and provenance status")
    p_status.set_defaults(func=cmd_status)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
