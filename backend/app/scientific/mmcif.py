"""Gemmi-based mmCIF parsing and structural metadata extraction.

Extracts PDB metadata, polymer entity descriptions, author vs label numbering,
observed coordinate status, C-alpha atom presence, and construct sequence differences.
"""

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import gemmi


@dataclass
class PolymerResidueInfo:
    """Individual residue in a polymer chain's full construct sequence."""
    seq_id: int  # label_seq_id (1-based index in entity polymer)
    mon_id: str  # 3-letter chemical component code (e.g., 'PHE')
    one_letter_code: str  # Single-letter amino acid code (e.g., 'F')
    auth_seq_id: Optional[str]  # Author residue identifier (e.g., '508')
    auth_comp_id: Optional[str]  # Author component code
    auth_asym_id: str  # Author chain identifier (e.g., 'A')
    label_asym_id: str  # mmCIF label chain identifier (e.g., 'A')
    insertion_code: Optional[str] = None  # Insertion code if any
    is_observed: bool = False  # Has atomic coordinates in model
    has_ca: bool = False  # Has C-alpha atom in model
    is_hetero: bool = False  # Flagged as hetero residue


@dataclass
class ChainMetadata:
    """Polymer chain sequence, numbering systems, and coordinate observability."""
    label_asym_id: str
    auth_asym_id: str
    entity_id: str
    entity_description: Optional[str]
    organism: Optional[str]
    polymer_sequence: str  # 1-letter code sequence of the construct
    polymer_length: int
    residues: List[PolymerResidueInfo] = field(default_factory=list)
    observed_count: int = 0
    ca_count: int = 0
    unobserved_count: int = 0
    unobserved_ranges: List[Tuple[int, int]] = field(default_factory=list)


@dataclass
class StructureMetadata:
    """Parsed mmCIF macromolecular structure metadata."""
    structure_id: str  # Uppercase PDB ID
    source_filename: str
    file_sha256: str
    title: Optional[str] = None
    experimental_method: Optional[str] = None
    resolution: Optional[float] = None
    deposition_date: Optional[str] = None
    chains: Dict[str, ChainMetadata] = field(default_factory=dict)  # keyed by label_asym_id
    polymer_entities: Dict[str, dict] = field(default_factory=dict)  # entity_id -> info
    non_polymer_entities: List[dict] = field(default_factory=list)
    ref_seq_differences: List[dict] = field(default_factory=list)
    uniprot_accessions: Dict[str, str] = field(default_factory=dict)  # entity_id -> UniProt accession
    parsing_warnings: List[str] = field(default_factory=list)


def calculate_file_sha256(file_path: Path) -> str:
    """Compute SHA-256 digest of raw binary file."""
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def clean_cif_value(val: Optional[object]) -> Optional[str]:
    """Clean string extracted from CIF block, treating CIF nulls as None."""
    if val is None or isinstance(val, bool):
        return None
    if not isinstance(val, str):
        val = str(val)
    val = val.strip().strip("'\"").strip(";")
    val = val.strip()
    if val in ("?", ".", "", "None"):
        return None
    return val


def parse_mmcif_file(file_path: Path) -> StructureMetadata:
    """Parse an mmCIF file using Gemmi and extract structured metadata.
    
    Args:
        file_path: Path to the .cif file.
        
    Returns:
        StructureMetadata containing parsed entities, chains, residue schemes,
        and coordinate observability.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"mmCIF file not found: {path}")

    file_sha256 = calculate_file_sha256(path)

    try:
        cif_doc = gemmi.cif.read(str(path))
        block = cif_doc.sole_block()
    except Exception as e:
        raise ValueError(f"Failed to parse CIF syntax in {path.name}: {e}") from e

    try:
        st = gemmi.read_structure(str(path))
    except Exception as e:
        raise ValueError(f"Gemmi failed to build structure model from {path.name}: {e}") from e

    warnings: List[str] = []

    # 1. PDB ID
    entry_id = clean_cif_value(block.find_value("_entry.id"))
    if not entry_id:
        entry_id = clean_cif_value(st.name)
    if not entry_id:
        # Fallback to file stem
        stem = path.stem
        # Clean any (1) suffix
        if "(" in stem:
            stem = stem.split("(")[0]
        entry_id = stem
    structure_id = entry_id.upper()

    # 2. Title
    raw_title = block.find_value("_struct.title")
    title = clean_cif_value(raw_title)

    # 3. Experimental Method
    exptl = block.get_mmcif_category("_exptl")
    experimental_method = None
    if exptl and "method" in exptl:
        methods = [clean_cif_value(m) for m in exptl.get("method", []) if clean_cif_value(m)]
        if methods:
            experimental_method = "; ".join(methods)

    # 4. Reported Resolution
    resolution: Optional[float] = None
    if st.resolution and st.resolution > 0:
        resolution = float(round(st.resolution, 2))
    else:
        # Check EM nominal resolution or refinement resolution
        for tag in (
            "_em_3d_reconstruction.nominal_resolution",
            "_refine.ls_d_res_high",
            "_reflns.d_resolution_high",
        ):
            res_val = clean_cif_value(block.find_value(tag))
            if res_val:
                try:
                    resolution = float(round(float(res_val), 2))
                    break
                except ValueError:
                    pass

    # 5. Deposition Date
    deposition_date = clean_cif_value(block.find_value("_pdbx_database_status.recvd_initial_deposition_date"))
    if not deposition_date:
        deposition_date = clean_cif_value(block.find_value("_database_PDB_rev.date_original"))

    # 6. Entity Descriptions & Organism
    entities_cat = block.get_mmcif_category("_entity")
    polymer_entities: Dict[str, dict] = {}
    non_polymer_entities: List[dict] = []

    if entities_cat and "id" in entities_cat:
        for i in range(len(entities_cat["id"])):
            eid = clean_cif_value(entities_cat["id"][i])
            if not eid:
                continue
            etype = clean_cif_value(entities_cat.get("type", [None])[i])
            desc = clean_cif_value(entities_cat.get("pdbx_description", [None])[i])
            if not desc:
                desc = clean_cif_value(entities_cat.get("description", [None])[i])

            if etype == "polymer":
                polymer_entities[eid] = {
                    "entity_id": eid,
                    "type": etype,
                    "description": desc,
                    "organism": None,
                }
            else:
                non_polymer_entities.append({
                    "entity_id": eid,
                    "type": etype,
                    "description": desc,
                })

    # Add organism metadata
    src_gen = block.get_mmcif_category("_entity_src_gen")
    if src_gen and "entity_id" in src_gen:
        for i in range(len(src_gen["entity_id"])):
            eid = clean_cif_value(src_gen["entity_id"][i])
            org = clean_cif_value(src_gen.get("pdbx_gene_src_scientific_name", [None])[i])
            if eid in polymer_entities and org:
                polymer_entities[eid]["organism"] = org

    src_nat = block.get_mmcif_category("_entity_src_nat")
    if src_nat and "entity_id" in src_nat:
        for i in range(len(src_nat["entity_id"])):
            eid = clean_cif_value(src_nat["entity_id"][i])
            org = clean_cif_value(src_nat.get("pdbx_organism_scientific_name", [None])[i])
            if eid in polymer_entities and org and not polymer_entities[eid].get("organism"):
                polymer_entities[eid]["organism"] = org

    # 7. UniProt Cross-references (_struct_ref)
    uniprot_accessions: Dict[str, str] = {}
    struct_ref = block.get_mmcif_category("_struct_ref")
    if struct_ref and "entity_id" in struct_ref:
        for i in range(len(struct_ref["entity_id"])):
            eid = clean_cif_value(struct_ref["entity_id"][i])
            db_name = clean_cif_value(struct_ref.get("db_name", [None])[i])
            db_acc = clean_cif_value(struct_ref.get("pdbx_db_accession", [None])[i])
            if eid and db_name in ("UNP", "UniProt") and db_acc:
                uniprot_accessions[eid] = db_acc

    # 8. Sequence Differences (_struct_ref_seq_dif)
    ref_seq_differences: List[dict] = []
    dif_cat = block.get_mmcif_category("_struct_ref_seq_dif")
    if dif_cat and "align_id" in dif_cat:
        for i in range(len(dif_cat["align_id"])):
            dif_row = {
                "align_id": clean_cif_value(dif_cat.get("align_id", [None])[i]),
                "mon_id": clean_cif_value(dif_cat.get("mon_id", [None])[i]),
                "db_mon_id": clean_cif_value(dif_cat.get("db_mon_id", [None])[i]),
                "seq_num": clean_cif_value(dif_cat.get("seq_num", [None])[i]),
                "details": clean_cif_value(dif_cat.get("details", [None])[i]),
                "strand_id": clean_cif_value(dif_cat.get("pdbx_pdb_strand_id", [None])[i]),
            }
            ref_seq_differences.append(dif_row)

    # 9. Coordinate Model Residues (Model 0)
    model0 = st[0] if len(st) > 0 else None
    observed_residues_by_chain: Dict[str, Dict[int, gemmi.Residue]] = {}
    # Keyed by (auth_asym_id, label_seq) or (subchain_id, label_seq)
    if model0:
        for ch in model0:
            for res in ch:
                if res.label_seq is not None:
                    # Map both by author chain and label asym (subchain)
                    subchain = clean_cif_value(res.subchain) or ch.name
                    if subchain not in observed_residues_by_chain:
                        observed_residues_by_chain[subchain] = {}
                    observed_residues_by_chain[subchain][res.label_seq] = res

                    if ch.name not in observed_residues_by_chain:
                        observed_residues_by_chain[ch.name] = {}
                    observed_residues_by_chain[ch.name][res.label_seq] = res

    # 10. Extract Polymer Chains from _pdbx_poly_seq_scheme
    scheme = block.get_mmcif_category("_pdbx_poly_seq_scheme")
    chains: Dict[str, ChainMetadata] = {}

    if scheme and "asym_id" in scheme and "seq_id" in scheme:
        # Group rows by asym_id (label_asym_id)
        rows_by_asym: Dict[str, List[int]] = {}
        for idx, aid in enumerate(scheme["asym_id"]):
            clean_aid = clean_cif_value(aid)
            if clean_aid:
                rows_by_asym.setdefault(clean_aid, []).append(idx)

        for label_asym_id, indices in rows_by_asym.items():
            first_idx = indices[0]
            entity_id = clean_cif_value(scheme["entity_id"][first_idx]) or "1"
            auth_asym_id = clean_cif_value(scheme.get("pdb_strand_id", [label_asym_id])[first_idx]) or label_asym_id

            ent_info = polymer_entities.get(entity_id, {})
            ent_desc = ent_info.get("description")
            organism = ent_info.get("organism")

            residues_list: List[PolymerResidueInfo] = []
            seq_letters: List[str] = []
            observed_count = 0
            ca_count = 0
            unobserved_count = 0

            # Find coordinates map for this chain
            obs_map = observed_residues_by_chain.get(label_asym_id)
            if not obs_map:
                obs_map = observed_residues_by_chain.get(auth_asym_id, {})

            for idx in indices:
                seq_id_str = clean_cif_value(scheme["seq_id"][idx])
                if not seq_id_str:
                    continue
                try:
                    seq_id = int(seq_id_str)
                except ValueError:
                    continue

                mon_id = clean_cif_value(scheme["mon_id"][idx]) or "UNK"
                tab_res = gemmi.find_tabulated_residue(mon_id)
                one_letter = tab_res.one_letter_code if tab_res.is_amino_acid() else (
                    "X" if mon_id in ("UNK", "MSE") else tab_res.one_letter_code or "X"
                )
                seq_letters.append(one_letter)

                # Author identifiers
                auth_seq_id = clean_cif_value(scheme.get("pdb_seq_num", [None])[idx])
                if not auth_seq_id:
                    auth_seq_id = clean_cif_value(scheme.get("auth_seq_num", [None])[idx])

                auth_comp_id = clean_cif_value(scheme.get("auth_mon_id", [None])[idx])
                if not auth_comp_id:
                    auth_comp_id = clean_cif_value(scheme.get("pdb_mon_id", [None])[idx])

                ins_code = clean_cif_value(scheme.get("pdb_ins_code", [None])[idx])
                if ins_code in ("?", ".", "False", "None", ""):
                    ins_code = None

                hetero_flag = clean_cif_value(scheme.get("hetero", [None])[idx])
                is_hetero = (hetero_flag == "y")

                # Check observation in coordinates
                obs_res = obs_map.get(seq_id)
                is_obs = (obs_res is not None)
                has_ca = (obs_res.get_ca() is not None) if obs_res else False

                if is_obs:
                    observed_count += 1
                else:
                    unobserved_count += 1

                if has_ca:
                    ca_count += 1

                residues_list.append(
                    PolymerResidueInfo(
                        seq_id=seq_id,
                        mon_id=mon_id,
                        one_letter_code=one_letter,
                        auth_seq_id=auth_seq_id,
                        auth_comp_id=auth_comp_id,
                        auth_asym_id=auth_asym_id,
                        label_asym_id=label_asym_id,
                        insertion_code=ins_code,
                        is_observed=is_obs,
                        has_ca=has_ca,
                        is_hetero=is_hetero,
                    )
                )

            # Determine contiguous unobserved ranges
            unobserved_ranges: List[Tuple[int, int]] = []
            range_start: Optional[int] = None
            for r in residues_list:
                if not r.is_observed:
                    if range_start is None:
                        range_start = r.seq_id
                else:
                    if range_start is not None:
                        unobserved_ranges.append((range_start, r.seq_id - 1))
                        range_start = None
            if range_start is not None and residues_list:
                unobserved_ranges.append((range_start, residues_list[-1].seq_id))

            chain_meta = ChainMetadata(
                label_asym_id=label_asym_id,
                auth_asym_id=auth_asym_id,
                entity_id=entity_id,
                entity_description=ent_desc,
                organism=organism,
                polymer_sequence="".join(seq_letters),
                polymer_length=len(residues_list),
                residues=residues_list,
                observed_count=observed_count,
                ca_count=ca_count,
                unobserved_count=unobserved_count,
                unobserved_ranges=unobserved_ranges,
            )
            chains[label_asym_id] = chain_meta

    else:
        # Fallback if _pdbx_poly_seq_scheme is missing (e.g., minimal test CIF)
        warnings.append("Missing _pdbx_poly_seq_scheme; falling back to model chains.")
        if model0:
            for ch in model0:
                amino_acids = [r for r in ch if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
                if not amino_acids:
                    continue
                label_asym_id = clean_cif_value(amino_acids[0].subchain) or ch.name
                auth_asym_id = ch.name
                entity_id = clean_cif_value(amino_acids[0].entity_id) or "1"
                ent_desc = polymer_entities.get(entity_id, {}).get("description")
                organism = polymer_entities.get(entity_id, {}).get("organism")

                res_list: List[PolymerResidueInfo] = []
                seq_letters = []
                ca_cnt = 0
                for idx, r in enumerate(amino_acids, start=1):
                    seq_id = r.label_seq if r.label_seq is not None else idx
                    tab = gemmi.find_tabulated_residue(r.name)
                    olc = tab.one_letter_code or "X"
                    seq_letters.append(olc)
                    has_ca = (r.get_ca() is not None)
                    if has_ca:
                        ca_cnt += 1
                    res_list.append(
                        PolymerResidueInfo(
                            seq_id=seq_id,
                            mon_id=r.name,
                            one_letter_code=olc,
                            auth_seq_id=str(r.seqid.num),
                            auth_comp_id=r.name,
                            auth_asym_id=auth_asym_id,
                            label_asym_id=label_asym_id,
                            insertion_code=clean_cif_value(str(r.seqid.icode)) if r.seqid.icode != ' ' else None,
                            is_observed=True,
                            has_ca=has_ca,
                            is_hetero=(r.het_flag == 'H'),
                        )
                    )
                chain_meta = ChainMetadata(
                    label_asym_id=label_asym_id,
                    auth_asym_id=auth_asym_id,
                    entity_id=entity_id,
                    entity_description=ent_desc,
                    organism=organism,
                    polymer_sequence="".join(seq_letters),
                    polymer_length=len(res_list),
                    residues=res_list,
                    observed_count=len(res_list),
                    ca_count=ca_cnt,
                    unobserved_count=0,
                    unobserved_ranges=[],
                )
                chains[label_asym_id] = chain_meta

    return StructureMetadata(
        structure_id=structure_id,
        source_filename=path.name,
        file_sha256=file_sha256,
        title=title,
        experimental_method=experimental_method,
        resolution=resolution,
        deposition_date=deposition_date,
        chains=chains,
        polymer_entities=polymer_entities,
        non_polymer_entities=non_polymer_entities,
        ref_seq_differences=ref_seq_differences,
        uniprot_accessions=uniprot_accessions,
        parsing_warnings=warnings,
    )
