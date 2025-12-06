"""
Stages 6, 7, 9, 10: Pair-related comparisons.

Stage 6: pair_validation - is_valid, bp_type_id
Stage 7: distance_checks - dorg, dNN, plane_angle, d_v
Stage 9: base_pair - bp_type, orien_i, org_i
Stage 10: find_bestpair_selection - num_bp, pairs
"""

from typing import Any, Dict, List, Tuple

from .base import (
    CompareResult, check_required_field,
    compare_scalar, compare_vector, compare_matrix
)
from ..config import Tolerance


# Type alias for pair keys: ((chain1, seq1, ins1), (chain2, seq2, ins2))
PairKey = Tuple[Tuple[str, int, str], Tuple[str, int, str]]


def _make_pair_key(rec: Dict[str, Any]) -> PairKey:
    """Create a pair key from a record."""
    return (
        (rec.get("chain_id_i", ""), rec.get("residue_seq_i", 0), rec.get("insertion_i", " ")),
        (rec.get("chain_id_j", ""), rec.get("residue_seq_j", 0), rec.get("insertion_j", " "))
    )


def _build_pair_lookup(records: List[Dict[str, Any]]) -> Dict[PairKey, Dict[str, Any]]:
    """Build lookup by pair identity."""
    return {_make_pair_key(rec): rec for rec in records}


def compare_pair_validation(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 0.0
) -> CompareResult:
    """
    Compare pair_validation records (Stage 6).
    
    Args:
        legacy_records: Legacy pair_validation records
        modern_records: Modern pair_validation records
        tolerance: Unused
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        
        leg_valid, mod_valid = check_required_field(leg, mod, "is_valid", key, errors)
        compare_scalar(leg_valid, mod_valid, "is_valid", key, errors)
        
        leg_type, mod_type = check_required_field(leg, mod, "bp_type_id", key, errors)
        compare_scalar(leg_type, mod_type, "bp_type_id", key, errors)
    
    _report_pair_diffs(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def compare_distance_checks(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.DISTANCE
) -> CompareResult:
    """
    Compare distance_checks records (Stage 7).
    
    Args:
        legacy_records: Legacy distance_checks records
        modern_records: Modern distance_checks records
        tolerance: Distance tolerance (default 1e-6)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_distances(key, leg, mod, errors, tolerance)
    
    return len(errors) == 0, errors


def _compare_distances(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare distance fields for a pair."""
    for field in ["dorg", "dNN", "plane_angle", "d_v"]:
        leg_val, mod_val = check_required_field(leg, mod, field, key, errors)
        compare_scalar(leg_val, mod_val, field, key, errors, tolerance)


def compare_base_pair(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.MATRIX
) -> CompareResult:
    """
    Compare base_pair records (Stage 9).
    
    Args:
        legacy_records: Legacy base_pair records
        modern_records: Modern base_pair records
        tolerance: Matrix tolerance (default 1e-4)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_base_pair_fields(key, leg, mod, errors, tolerance)
    
    _report_pair_diffs(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _compare_base_pair_fields(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare fields for a base_pair record."""
    # bp_type
    leg_type, mod_type = check_required_field(leg, mod, "bp_type", key, errors)
    compare_scalar(leg_type, mod_type, "bp_type", key, errors)
    
    # orien_i (rotation matrix)
    leg_orien, mod_orien = check_required_field(leg, mod, "orien_i", key, errors)
    compare_matrix(leg_orien, mod_orien, "orien_i", key, errors, tolerance)
    
    # org_i (origin)
    leg_org, mod_org = check_required_field(leg, mod, "org_i", key, errors)
    compare_vector(leg_org, mod_org, "org_i", key, errors, Tolerance.COORDINATE)


def compare_best_pair_selection(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 0.0
) -> CompareResult:
    """
    Compare find_bestpair_selection records (Stage 10).
    
    Args:
        legacy_records: Legacy selection records
        modern_records: Modern selection records
        tolerance: Unused
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    # Expect exactly one record
    if len(legacy_records) != 1 or len(modern_records) != 1:
        errors.append(f"Expected 1 record each, got {len(legacy_records)} vs {len(modern_records)}")
        return False, errors
    
    leg, mod = legacy_records[0], modern_records[0]
    
    # num_bp
    leg_num, mod_num = check_required_field(leg, mod, "num_bp", "selection", errors)
    compare_scalar(leg_num, mod_num, "num_bp", "selection", errors)
    
    # pairs (list of pair indices)
    leg_pairs = leg.get("pairs")
    mod_pairs = mod.get("pairs")
    if leg_pairs is None:
        errors.append("Legacy missing pairs field")
    elif mod_pairs is None:
        errors.append("Modern missing pairs field")
    elif set(map(tuple, leg_pairs)) != set(map(tuple, mod_pairs)):
        errors.append(f"pairs differ: {len(leg_pairs)} vs {len(mod_pairs)}")
    
    return len(errors) == 0, errors


def _report_pair_diffs(
    leg_lookup: Dict[PairKey, Any],
    mod_lookup: Dict[PairKey, Any],
    errors: List[str]
) -> None:
    """Report pairs missing from either dataset."""
    missing = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for key in list(missing)[:3]:
        errors.append(f"Missing in modern: {key}")
    
    extra = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for key in list(extra)[:3]:
        errors.append(f"Extra in modern: {key}")

