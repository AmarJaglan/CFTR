import pytest

from backend.app.scientific.normalization import normalize_protein_change


def test_normalize_1letter_missense():
    v = normalize_protein_change("N1303K")
    assert v.success is True
    assert v.canonical_protein == "N1303K"
    assert v.source_protein == "N1303K"
    assert v.ref_aa == "N"
    assert v.pos == 1303
    assert v.alt_aa == "K"
    assert v.is_missense is True
    assert v.variant_type == "missense"


def test_normalize_3letter_missense():
    v = normalize_protein_change("p.Gly551Asp")
    assert v.success is True
    assert v.canonical_protein == "G551D"
    assert v.source_protein == "p.Gly551Asp"
    assert v.ref_aa == "G"
    assert v.pos == 551
    assert v.alt_aa == "D"
    assert v.is_missense is True
    assert v.variant_type == "missense"


def test_normalize_with_fallback_name():
    v = normalize_protein_change("", fallback_name="NM_000492.4(CFTR):c.1A>G (p.Met1Val)")
    assert v.success is True
    assert v.canonical_protein == "M1V"
    assert v.source_protein == "p.Met1Val"
    assert v.ref_aa == "M"
    assert v.pos == 1
    assert v.alt_aa == "V"
    assert v.is_missense is True


def test_normalize_deletion():
    # 3-letter deletion
    v3 = normalize_protein_change("p.Phe508del")
    assert v3.success is True
    assert v3.canonical_protein == "F508del"
    assert v3.source_protein == "p.Phe508del"
    assert v3.ref_aa == "F"
    assert v3.pos == 508
    assert v3.alt_aa is None
    assert v3.is_missense is False
    assert v3.variant_type == "inframe_deletion"

    # 1-letter deletion
    v1 = normalize_protein_change("F508del")
    assert v1.success is True
    assert v1.canonical_protein == "F508del"
    assert v1.source_protein == "F508del"
    assert v1.ref_aa == "F"
    assert v1.pos == 508


def test_normalize_nonsense():
    v = normalize_protein_change("p.Gly542X")
    assert v.success is True
    assert v.canonical_protein == "G542X"
    assert v.ref_aa == "G"
    assert v.pos == 542
    assert v.alt_aa == "X"
    assert v.variant_type == "nonsense"


def test_normalize_unknown_or_empty():
    v_empty = normalize_protein_change("")
    assert v_empty.success is False
    assert v_empty.canonical_protein is None
    assert v_empty.failure_reason == "MISSING_OR_UNKNOWN_PROTEIN_CHANGE"

    v_question = normalize_protein_change("p.?")
    assert v_question.success is False
    assert v_question.canonical_protein is None
    assert v_question.failure_reason == "MISSING_OR_UNKNOWN_PROTEIN_CHANGE"


def test_normalize_synonymous():
    # 3-letter synonymous HGVS notation
    v3 = normalize_protein_change("p.His139=")
    assert v3.success is True
    assert v3.canonical_protein == "H139="
    assert v3.source_protein == "p.His139="
    assert v3.ref_aa == "H"
    assert v3.pos == 139
    assert v3.alt_aa == "H"
    assert v3.variant_type == "synonymous"
    assert v3.is_missense is False

    # Another 3-letter synonymous
    vt = normalize_protein_change("p.Thr854=")
    assert vt.success is True
    assert vt.canonical_protein == "T854="
    assert vt.pos == 854

    # 1-letter synonymous
    v1 = normalize_protein_change("p.H139=")
    assert v1.success is True
    assert v1.canonical_protein == "H139="


def test_normalize_nonsense_star_and_ter():
    v_star = normalize_protein_change("p.Trp1282*")
    assert v_star.success is True
    assert v_star.canonical_protein == "W1282X"
    assert v_star.pos == 1282
    assert v_star.alt_aa == "X"

    v_ter = normalize_protein_change("p.Trp1282Ter")
    assert v_ter.success is True
    assert v_ter.canonical_protein == "W1282X"


def test_complex_and_haplotypes_not_reduced_to_single_missense():
    """Ensure complex variants and haplotypes fail normalization rather than being silently reduced."""
    # Multi-site haplotype
    v_hap = normalize_protein_change("S912L, G1244V")
    assert v_hap.success is False
    assert v_hap.canonical_protein is None
    assert v_hap.variant_type == "complex_unparsed"

    # Bracketed complex allele
    v_bracket = normalize_protein_change("p.[Arg117His;7T]")
    assert v_bracket.success is False
    assert v_bracket.canonical_protein is None

    # Multi-residue deletion
    v_multi_del = normalize_protein_change("p.Ile1023_Val1024del")
    assert v_multi_del.success is False
    assert v_multi_del.canonical_protein is None

