"""Variant protein change normalization.

Converts diverse protein change notations (HGVS 3-letter, 1-letter, legacy)
into canonical short forms (e.g., 'N1303K', 'G551D', 'F508del') while
strictly preserving original source notation.
"""

from dataclasses import dataclass
import re
from typing import Optional

from .reference import AA3_TO_AA1


@dataclass(frozen=True)
class NormalizedVariant:
    """Standardized representation of a parsed CFTR protein variant."""
    canonical_protein: Optional[str]
    source_protein: str
    ref_aa: Optional[str]
    pos: Optional[int]
    alt_aa: Optional[str]
    variant_type: str
    is_missense: bool
    success: bool
    failure_reason: Optional[str] = None


# Regex patterns
RE_HGVS_3LETTER_MISSENSE = re.compile(
    r"^(?:p\.)?(?:\()?(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2})(?:\))?$"
)
RE_1LETTER_MISSENSE = re.compile(
    r"^(?:p\.)?(?P<ref>[A-Z])(?P<pos>\d+)(?P<alt>[A-Z])$"
)
RE_3LETTER_DEL = re.compile(
    r"^(?:p\.)?(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)del(?:[A-Za-z0-9_]*)?$"
)
RE_1LETTER_DEL = re.compile(
    r"^(?:p\.)?(?P<ref>[A-Z])(?P<pos>\d+)del(?:[A-Za-z0-9_]*)?$"
)
RE_NONSENSE = re.compile(
    r"^(?:p\.)?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?P<alt>X|\*|Ter)$"
)
RE_FRAMESHIFT = re.compile(
    r"^(?:p\.)?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)(?:[A-Za-z0-9_]*)?fs(?:[A-Za-z0-9_*]*)?$"
)
RE_SYNONYMOUS = re.compile(
    r"^(?:p\.)?(?:\()?(?P<ref>[A-Z][a-z]{2}|[A-Z])(?P<pos>\d+)=(?:\))?$"
)
RE_EXTRACT_P_FROM_NAME = re.compile(r"\(p\.([^\)]+)\)")


def normalize_protein_change(protein_raw: Optional[str], fallback_name: Optional[str] = None) -> NormalizedVariant:
    """Normalize protein notation to canonical representation.
    
    Args:
        protein_raw: Raw protein notation (e.g. 'M1V', 'p.Gly551Asp', 'F508del', 'p.?').
        fallback_name: Optional full variant name to extract protein change from if raw is empty
                       (e.g., 'NM_000492.4(CFTR):c.1A>G (p.Met1Val)').
                       
    Returns:
        NormalizedVariant dataclass preserving source notation and parsed fields.
    """
    raw_str = (protein_raw or "").strip()
    source_notation = raw_str

    if not raw_str and fallback_name:
        match_p = RE_EXTRACT_P_FROM_NAME.search(fallback_name)
        if match_p:
            raw_str = f"p.{match_p.group(1)}"
            source_notation = raw_str

    if not raw_str or raw_str in {"p.?", "?", "None", ".", "-"}:
        return NormalizedVariant(
            canonical_protein=None,
            source_protein=source_notation,
            ref_aa=None,
            pos=None,
            alt_aa=None,
            variant_type="unknown",
            is_missense=False,
            success=False,
            failure_reason="MISSING_OR_UNKNOWN_PROTEIN_CHANGE",
        )

    # 1. Check nonsense / stop-gain (e.g., G542X, p.Gly542X, p.Trp1282*, p.Trp1282Ter)
    match_non = RE_NONSENSE.match(raw_str)
    if match_non:
        ref = match_non.group("ref")
        pos = int(match_non.group("pos"))
        if len(ref) == 3:
            ref = AA3_TO_AA1.get(ref, ref)
        canonical = f"{ref}{pos}X"
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref,
            pos=pos,
            alt_aa="X",
            variant_type="nonsense",
            is_missense=False,
            success=True,
        )

    # 2. Check 3-letter HGVS missense (e.g., p.Gly551Asp, Gly551Asp)
    match_3m = RE_HGVS_3LETTER_MISSENSE.match(raw_str)
    if match_3m:
        ref3 = match_3m.group("ref")
        pos = int(match_3m.group("pos"))
        alt3 = match_3m.group("alt")
        ref1 = AA3_TO_AA1.get(ref3)
        alt1 = AA3_TO_AA1.get(alt3)
        if not ref1 or not alt1:
            return NormalizedVariant(
                canonical_protein=None,
                source_protein=source_notation,
                ref_aa=None,
                pos=pos,
                alt_aa=None,
                variant_type="missense",
                is_missense=True,
                success=False,
                failure_reason=f"UNKNOWN_AMINO_ACID_CODE_{ref3}_{alt3}",
            )
        variant_type = "synonymous" if ref1 == alt1 else "missense"
        canonical = f"{ref1}{pos}{alt1}"
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref1,
            pos=pos,
            alt_aa=alt1,
            variant_type=variant_type,
            is_missense=(variant_type == "missense"),
            success=True,
        )

    # 2. Check 1-letter missense (e.g., M1V, G551D, N1303K)
    match_1m = RE_1LETTER_MISSENSE.match(raw_str)
    if match_1m:
        ref1 = match_1m.group("ref")
        pos = int(match_1m.group("pos"))
        alt1 = match_1m.group("alt")
        variant_type = "synonymous" if ref1 == alt1 else "missense"
        canonical = f"{ref1}{pos}{alt1}"
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref1,
            pos=pos,
            alt_aa=alt1,
            variant_type=variant_type,
            is_missense=(variant_type == "missense"),
            success=True,
        )

    # 3. Check 3-letter deletion (e.g., p.Phe508del)
    match_3d = RE_3LETTER_DEL.match(raw_str)
    if match_3d:
        ref3 = match_3d.group("ref")
        pos = int(match_3d.group("pos"))
        ref1 = AA3_TO_AA1.get(ref3, ref3)
        canonical = f"{ref1}{pos}del"
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref1,
            pos=pos,
            alt_aa=None,
            variant_type="inframe_deletion",
            is_missense=False,
            success=True,
        )

    # 4. Check 1-letter deletion (e.g., F508del)
    match_1d = RE_1LETTER_DEL.match(raw_str)
    if match_1d:
        ref1 = match_1d.group("ref")
        pos = int(match_1d.group("pos"))
        canonical = f"{ref1}{pos}del"
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref1,
            pos=pos,
            alt_aa=None,
            variant_type="inframe_deletion",
            is_missense=False,
            success=True,
        )

    # 5. Check frameshift (e.g., p.Leu167fs)
    match_fs = RE_FRAMESHIFT.match(raw_str)
    if match_fs:
        ref = match_fs.group("ref")
        pos = int(match_fs.group("pos"))
        if len(ref) == 3:
            ref = AA3_TO_AA1.get(ref, ref)
        canonical = f"{ref}{pos}fs"
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref,
            pos=pos,
            alt_aa=None,
            variant_type="frameshift",
            is_missense=False,
            success=True,
        )

    # 7. Check synonymous (e.g., p.His139=, p.Thr854=, p.H139=)
    match_syn = RE_SYNONYMOUS.match(raw_str)
    if match_syn:
        ref = match_syn.group("ref")
        pos = int(match_syn.group("pos"))
        if len(ref) == 3:
            ref1 = AA3_TO_AA1.get(ref)
            if not ref1:
                return NormalizedVariant(
                    canonical_protein=None,
                    source_protein=source_notation,
                    ref_aa=None,
                    pos=pos,
                    alt_aa=None,
                    variant_type="synonymous",
                    is_missense=False,
                    success=False,
                    failure_reason=f"UNKNOWN_AMINO_ACID_CODE_{ref}",
                )
        else:
            ref1 = ref
        canonical = f"{ref1}{pos}="
        return NormalizedVariant(
            canonical_protein=canonical,
            source_protein=source_notation,
            ref_aa=ref1,
            pos=pos,
            alt_aa=ref1,
            variant_type="synonymous",
            is_missense=False,
            success=True,
        )

    # 8. Unparsed complex variant
    return NormalizedVariant(
        canonical_protein=None,
        source_protein=source_notation,
        ref_aa=None,
        pos=None,
        alt_aa=None,
        variant_type="complex_unparsed",
        is_missense=False,
        success=False,
        failure_reason="UNRECOGNIZED_PROTEIN_FORMAT",
    )
