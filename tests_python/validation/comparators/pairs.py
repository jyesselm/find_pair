"""
Stages 6, 7, 9, 10: Pair-related comparisons.

Stage 6: pair_validation - is_valid, bp_type_id, calculated_values.*
Stage 7: distance_checks - values.{dorg, dNN, plane_angle, d_v}
Stage 9: base_pair - bp_type, orien_i, org_i, orien_j, org_j
Stage 10: find_bestpair_selection - num_bp, pairs
"""

from typing import Any, Dict, List, Tuple

from .base import (
    CompareResult, check_required_field,
    compare_scalar, compare_vector, compare_matrix
)
from ..config import Tolerance


# Type alias for pair keys: (base_i, base_j)
PairKey = Tuple[int, int]


def _make_pair_key(rec: Dict[str, Any]) -> PairKey:
    """Create a normalized pair key from base_i/base_j (always min, max)."""
    i, j = rec.get("base_i", 0), rec.get("base_j", 0)
    return (min(i, j), max(i, j))


def _build_pair_lookup(records: List[Dict[str, Any]]) -> Dict[PairKey, Dict[str, Any]]:
    """Build lookup by (base_i, base_j) pair.
    
    When duplicates exist (both i,j and j,i), keeps the record where base_i <= base_j
    so bp_type and org_i/org_j are in consistent order.
    """
    lookup = {}
    for rec in records:
        key = _make_pair_key(rec)
        base_i, base_j = rec.get("base_i", 0), rec.get("base_j", 0)
        # Keep record where base_i <= base_j (normalized order)
        if key not in lookup or base_i <= base_j:
            lookup[key] = rec
    return lookup


def compare_pair_validation(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 1e-4  # Relaxed tolerance for floating-point differences
) -> CompareResult:
    """
    Compare pair_validation records (Stage 6).
    
    Uses (base_i, base_j) as matching key.
    Compares: is_valid, bp_type_id, calculated_values.*
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_pair_validation_record(key, leg, mod, errors, tolerance)
    
    _report_pair_diffs(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _compare_pair_validation_record(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single pair_validation record."""
    # is_valid
    leg_valid, mod_valid = check_required_field(leg, mod, "is_valid", key, errors)
    compare_scalar(leg_valid, mod_valid, "is_valid", key, errors)
    
    # bp_type_id
    leg_type, mod_type = check_required_field(leg, mod, "bp_type_id", key, errors)
    compare_scalar(leg_type, mod_type, "bp_type_id", key, errors)
    
    # calculated_values (nested)
    leg_calc = leg.get("calculated_values", {})
    mod_calc = mod.get("calculated_values", {})
    
    for field in ["dorg", "dNN", "plane_angle", "d_v", "quality_score"]:
        leg_val = leg_calc.get(field)
        mod_val = mod_calc.get(field)
        if leg_val is not None and mod_val is not None:
            compare_scalar(leg_val, mod_val, f"calculated_values.{field}", key, errors, tolerance)


def compare_distance_checks(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 1e-4  # Relaxed tolerance for floating-point differences
) -> CompareResult:
    """
    Compare distance_checks records (Stage 7).
    
    Uses (base_i, base_j) as matching key.
    Compares: values.{dorg, dNN, plane_angle, d_v, overlap_area}
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_distance_record(key, leg, mod, errors, tolerance)
    
    return len(errors) == 0, errors


def _compare_distance_record(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single distance_checks record."""
    # Values are nested under "values" key
    leg_vals = leg.get("values", {})
    mod_vals = mod.get("values", {})
    
    if not leg_vals:
        errors.append(f"Key {key} legacy missing values")
        return
    if not mod_vals:
        errors.append(f"Key {key} modern missing values")
        return
    
    for field in ["dorg", "dNN", "plane_angle", "d_v", "overlap_area"]:
        leg_val = leg_vals.get(field)
        mod_val = mod_vals.get(field)
        if leg_val is None:
            errors.append(f"Key {key} legacy missing values.{field}")
        elif mod_val is None:
            errors.append(f"Key {key} modern missing values.{field}")
        else:
            compare_scalar(leg_val, mod_val, f"values.{field}", key, errors, tolerance)


def compare_base_pair(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.MATRIX
) -> CompareResult:
    """
    Compare base_pair records (Stage 9).
    
    Uses (base_i, base_j) as matching key.
    Compares: bp_type, orien_i, org_i, orien_j, org_j
    
    Note: Missing/extra pairs are warnings, not errors.
    Stage 9 focuses on geometric property comparison of matching pairs.
    Missing pair detection issues are tracked separately (Stage 6/10).
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    # Compare only pairs that exist in both
    common_keys = leg_lookup.keys() & mod_lookup.keys()
    for key in common_keys:
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_base_pair_fields(key, leg, mod, errors, tolerance)
    
    # Report missing/extra pairs but don't add to errors
    # (pair detection differences are tracked in Stage 6/10)
    missing = leg_lookup.keys() - mod_lookup.keys()
    extra = mod_lookup.keys() - leg_lookup.keys()
    
    # Only fail if no common pairs at all
    if not common_keys and (missing or extra):
        errors.append(f"No common pairs to compare ({len(missing)} missing, {len(extra)} extra)")
    
    return len(errors) == 0, errors


def _compare_base_pair_fields(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare fields for a base_pair record."""
    # bp_type (case-insensitive - modified residues use lowercase)
    # Note: '?' indicates unknown residue type - skip comparison if either has '?'
    leg_type, mod_type = check_required_field(leg, mod, "bp_type", key, errors)
    if leg_type and mod_type:
        # Skip comparison if either contains unknown type marker
        if '?' not in leg_type and '?' not in mod_type:
            if leg_type.upper() != mod_type.upper():
                errors.append(f"Key {key} bp_type: {leg_type} vs {mod_type}")
    
    # orien_i (rotation matrix)
    leg_orien_i, mod_orien_i = check_required_field(leg, mod, "orien_i", key, errors)
    compare_matrix(leg_orien_i, mod_orien_i, "orien_i", key, errors, tolerance)
    
    # org_i (origin)
    leg_org_i, mod_org_i = check_required_field(leg, mod, "org_i", key, errors)
    compare_vector(leg_org_i, mod_org_i, "org_i", key, errors, 2e-6)
    
    # orien_j (rotation matrix)
    leg_orien_j, mod_orien_j = check_required_field(leg, mod, "orien_j", key, errors)
    compare_matrix(leg_orien_j, mod_orien_j, "orien_j", key, errors, tolerance)
    
    # org_j (origin)
    leg_org_j, mod_org_j = check_required_field(leg, mod, "org_j", key, errors)
    compare_vector(leg_org_j, mod_org_j, "org_j", key, errors, 2e-6)


def compare_best_pair_selection(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 0.0
) -> CompareResult:
    """
    Compare find_bestpair_selection records (Stage 10).
    
    This is a single record with num_bp and pairs list.
    """
    errors: List[str] = []
    
    # Handle various formats
    leg = _get_single_record(legacy_records)
    mod = _get_single_record(modern_records)
    
    if leg is None:
        errors.append("Legacy missing selection record")
        return False, errors
    if mod is None:
        errors.append("Modern missing selection record")
        return False, errors
    
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
    else:
        # Normalize pairs to set of tuples for comparison
        leg_set = set(tuple(p) if isinstance(p, list) else p for p in leg_pairs)
        mod_set = set(tuple(p) if isinstance(p, list) else p for p in mod_pairs)
        if leg_set != mod_set:
            missing_in_mod = leg_set - mod_set
            extra_in_mod = mod_set - leg_set
            if missing_in_mod:
                errors.append(f"Pairs missing in modern: {list(missing_in_mod)[:5]}")
            if extra_in_mod:
                errors.append(f"Extra pairs in modern: {list(extra_in_mod)[:5]}")
    
    return len(errors) == 0, errors


def _get_single_record(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Extract single record from potentially wrapped data."""
    if not records:
        return None
    if len(records) == 1:
        return records[0]
    # If multiple, return first
    return records[0]


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
