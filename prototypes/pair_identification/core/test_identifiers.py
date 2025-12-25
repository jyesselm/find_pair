"""Test identifier conversion utilities."""

from prototypes.pair_identification.core.identifiers import (
    dssr_to_res_id,
    res_id_to_dssr,
    parse_res_id,
    parse_dssr_id,
    make_res_id,
    normalize_chain,
    extract_sequence,
    is_standard_wc_sequence,
    is_canonical_pair,
)


def test_dssr_to_res_id():
    """Test DSSR to res_id conversion."""
    assert dssr_to_res_id("A.G1") == "A-G-1"
    assert dssr_to_res_id("0.C530") == "0-C-530"
    assert dssr_to_res_id("D.DG10") == "D-DG-10"
    assert dssr_to_res_id("A.G1^A") == "A-G-1"  # Insertion code stripped
    assert dssr_to_res_id("") is None
    assert dssr_to_res_id("invalid") is None


def test_res_id_to_dssr():
    """Test res_id to DSSR conversion."""
    assert res_id_to_dssr("A-G-1") == "A.G1"
    assert res_id_to_dssr("0-C-530") == "0.C530"
    assert res_id_to_dssr("D-DG-10") == "D.DG10"
    assert res_id_to_dssr("A-G--5") == "A.G-5"  # Negative residue number
    assert res_id_to_dssr("") is None
    assert res_id_to_dssr("A-G") is None  # Too few parts


def test_parse_res_id():
    """Test res_id parsing."""
    assert parse_res_id("A-G-1") == ("A", "G", "1")
    assert parse_res_id("0-C-530") == ("0", "C", "530")
    assert parse_res_id("D-DG-10") == ("D", "DG", "10")
    assert parse_res_id("A-G--5") == ("A", "G", "-5")  # Negative
    assert parse_res_id("A-G-1A") == ("A", "G", "1A")  # Insertion
    assert parse_res_id("") is None
    assert parse_res_id("A-G") is None


def test_parse_dssr_id():
    """Test DSSR ID parsing."""
    assert parse_dssr_id("A.G1") == ("A", "G", "1")
    assert parse_dssr_id("0.C530") == ("0", "C", "530")
    assert parse_dssr_id("D.DG10") == ("D", "DG", "10")
    assert parse_dssr_id("A.G1^A") == ("A", "G", "1")  # Alt location stripped
    assert parse_dssr_id("A.G-5") == ("A", "G", "-5")  # Negative
    assert parse_dssr_id("") is None
    assert parse_dssr_id("invalid") is None


def test_make_res_id():
    """Test res_id construction."""
    assert make_res_id("A", "G", "1") == "A-G-1"
    assert make_res_id("0", "C", "530") == "0-C-530"
    assert make_res_id("D", "DG", "10") == "D-DG-10"


def test_normalize_chain():
    """Test chain normalization."""
    assert normalize_chain("A") == "A"
    assert normalize_chain(" A ") == "A"
    assert normalize_chain("") == "A"  # Empty -> "A"
    assert normalize_chain("  ") == "A"


def test_extract_sequence():
    """Test sequence extraction from pair."""
    assert extract_sequence("A-G-1", "A-C-10") == "GC"
    assert extract_sequence("A-A-1", "A-U-10") == "AU"
    assert extract_sequence("invalid", "A-G-1") is None
    assert extract_sequence("A-G-1", "invalid") is None


def test_is_standard_wc_sequence():
    """Test standard WC sequence check."""
    assert is_standard_wc_sequence("GC")
    assert is_standard_wc_sequence("CG")
    assert is_standard_wc_sequence("AU")
    assert is_standard_wc_sequence("UA")
    assert is_standard_wc_sequence("AT")
    assert is_standard_wc_sequence("TA")
    assert is_standard_wc_sequence("gc")  # Case insensitive
    assert not is_standard_wc_sequence("GU")
    assert not is_standard_wc_sequence("AA")


def test_is_canonical_pair():
    """Test canonical pair check."""
    # Standard WC
    assert is_canonical_pair("GC")
    assert is_canonical_pair("CG")
    assert is_canonical_pair("AU")
    assert is_canonical_pair("UA")
    assert is_canonical_pair("AT")
    assert is_canonical_pair("TA")

    # Wobble
    assert is_canonical_pair("GU")
    assert is_canonical_pair("UG")
    assert is_canonical_pair("GT")
    assert is_canonical_pair("TG")

    # Non-canonical
    assert not is_canonical_pair("AA")
    assert not is_canonical_pair("GG")
    assert not is_canonical_pair("AC")


def test_round_trip_conversion():
    """Test round-trip DSSR <-> res_id conversion."""
    test_cases = [
        "A.G1",
        "0.C530",
        "D.DG10",
        "B.U25",
    ]

    for dssr_id in test_cases:
        res_id = dssr_to_res_id(dssr_id)
        assert res_id is not None

        back_to_dssr = res_id_to_dssr(res_id)
        assert back_to_dssr == dssr_id
