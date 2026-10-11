import pytest

from backend.app.scientific.reference import (
    CFTR_CANONICAL_SEQUENCE,
    SEQUENCE_LENGTH,
    get_canonical_residue,
    validate_wt_residue,
)


def test_sequence_length_and_endpoints():
    assert len(CFTR_CANONICAL_SEQUENCE) == 1480
    assert SEQUENCE_LENGTH == 1480
    assert get_canonical_residue(1) == "M"
    assert get_canonical_residue(1480) == "L"
    assert get_canonical_residue(0) is None
    assert get_canonical_residue(1481) is None


def test_validate_wt_canonical_matches():
    # Position 508: Phe
    res_508 = validate_wt_residue(508, "F")
    assert res_508.is_valid is True
    assert res_508.canonical_ref == "F"
    assert res_508.flag == "MATCH"

    # With 3-letter code
    res_508_3 = validate_wt_residue(508, "Phe")
    assert res_508_3.is_valid is True
    assert res_508_3.canonical_ref == "F"
    assert res_508_3.flag == "MATCH"

    # Position 551: Gly
    res_551 = validate_wt_residue(551, "G")
    assert res_551.is_valid is True
    assert res_551.canonical_ref == "G"
    assert res_551.flag == "MATCH"

    # Position 1303: Asn
    res_1303 = validate_wt_residue(1303, "N")
    assert res_1303.is_valid is True
    assert res_1303.canonical_ref == "N"
    assert res_1303.flag == "MATCH"


def test_validate_wt_mismatches_flagged_not_corrected():
    # If someone reports G508 instead of F508
    res = validate_wt_residue(508, "G")
    assert res.is_valid is False
    assert res.flag == "WT_MISMATCH"
    assert res.reported_ref == "G"
    assert res.canonical_ref == "F"
    assert "Reference mismatch at position 508" in res.message


def test_validate_wt_out_of_bounds():
    res_neg = validate_wt_residue(-5, "M")
    assert res_neg.is_valid is False
    assert res_neg.flag == "OUT_OF_BOUNDS"

    res_too_high = validate_wt_residue(2000, "M")
    assert res_too_high.is_valid is False
    assert res_too_high.flag == "OUT_OF_BOUNDS"


def test_validate_wt_m470v_polymorphism():
    """Verify that position 470 canonical Val is enforced, and Met is flagged as WT_MISMATCH with historical context."""
    # Reported reference V (as in V470M / modern HGVS) -> Matches
    res_v = validate_wt_residue(470, "V")
    assert res_v.is_valid is True
    assert res_v.canonical_ref == "V"
    assert res_v.flag == "MATCH"

    res_val = validate_wt_residue(470, "Val")
    assert res_val.is_valid is True
    assert res_val.canonical_ref == "V"
    assert res_val.flag == "MATCH"

    # Reported reference M (as in M470V / legacy 1989 cDNA) -> Mismatch flagged with context
    res_m = validate_wt_residue(470, "M")
    assert res_m.is_valid is False
    assert res_m.canonical_ref == "V"
    assert res_m.reported_ref == "M"
    assert res_m.flag == "WT_MISMATCH"
    assert "Historical CFTR cDNAs" in res_m.message
    assert "rs1800089" in res_m.message

    res_met = validate_wt_residue(470, "Met")
    assert res_met.is_valid is False
    assert res_met.canonical_ref == "V"
    assert res_met.flag == "WT_MISMATCH"

