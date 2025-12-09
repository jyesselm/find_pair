"""
Stages 11, 12: Parameter comparisons.

Stage 11: bpstep_params - shift, slide, rise, tilt, roll, twist
Stage 12: helical_params - x_displacement, y_displacement, rise, inclination, tip, twist

Legacy format uses "params" dict/array, modern uses individual fields.
"""

from typing import Any, Dict, List, Tuple

from .base import CompareResult, compare_scalar
from ..config import Tolerance


# Step parameter fields (Stage 11)
# Order matches legacy "params" dict keys: Shift, Slide, Rise, Tilt, Roll, Twist
STEP_PARAM_FIELDS = ["shift", "slide", "rise", "tilt", "roll", "twist"]
STEP_PARAM_LEGACY_KEYS = ["Shift", "Slide", "Rise", "Tilt", "Roll", "Twist"]

# Helical parameter fields (Stage 12)
# Order matches legacy "params" array: x_displacement, y_displacement, rise, inclination, tip, twist
HELICAL_PARAM_FIELDS = ["x_displacement", "y_displacement", "rise", "inclination", "tip", "twist"]


def compare_step_params(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.PARAMETER
) -> CompareResult:
    """
    Compare bpstep_params records (Stage 11).
    
    Args:
        legacy_records: Legacy step parameter records
        modern_records: Modern step parameter records
        tolerance: Parameter tolerance (default 1e-6)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Count mismatch: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_step_params_record(key, leg, mod, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def compare_helical_params(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.PARAMETER
) -> CompareResult:
    """
    Compare helical_params records (Stage 12).
    
    Args:
        legacy_records: Legacy helical parameter records
        modern_records: Modern helical parameter records
        tolerance: Parameter tolerance (default 1e-6)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = _build_pair_lookup(legacy_records)
    mod_lookup = _build_pair_lookup(modern_records)
    
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Count mismatch: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_helical_params_record(key, leg, mod, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _build_pair_lookup(records: List[Dict[str, Any]]) -> Dict[Tuple[int, int], Dict[str, Any]]:
    """Build lookup by (bp_idx1, bp_idx2) pair.
    
    Only keeps the first occurrence of each key (for duplex 1 in legacy data
    which contains both duplex 1 and duplex 2 with same bp_idx pairs).
    """
    lookup = {}
    for rec in records:
        bp_idx1 = rec.get("bp_idx1")
        bp_idx2 = rec.get("bp_idx2")
        if bp_idx1 is not None and bp_idx2 is not None:
            key = (bp_idx1, bp_idx2)
            # Only keep first occurrence (duplex 1)
            if key not in lookup:
                lookup[key] = rec
    return lookup


def _compare_step_params_record(
    key: Tuple[int, int],
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare step parameter fields between legacy and modern."""
    # Legacy has params as dict with capitalized keys: {"Shift": ..., "Slide": ...}
    # Modern has lowercase individual fields: {"shift": ..., "slide": ...}
    leg_params = leg.get("params", {})
    
    for field, leg_key in zip(STEP_PARAM_FIELDS, STEP_PARAM_LEGACY_KEYS):
        leg_val = leg_params.get(leg_key) if isinstance(leg_params, dict) else None
        mod_val = mod.get(field)
        
        if leg_val is None:
            errors.append(f"Key {key} missing {leg_key} in legacy")
            continue
        if mod_val is None:
            errors.append(f"Key {key} missing {field} in modern")
            continue
            
        compare_scalar(leg_val, mod_val, field, key, errors, tolerance)


def _compare_helical_params_record(
    key: Tuple[int, int],
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare helical parameter fields between legacy and modern."""
    # Legacy has params as array: [x_disp, y_disp, rise, incl, tip, twist]
    # Modern has individual fields: {"x_displacement": ..., "y_displacement": ...}
    leg_params = leg.get("params", [])
    
    for i, field in enumerate(HELICAL_PARAM_FIELDS):
        leg_val = leg_params[i] if isinstance(leg_params, list) and i < len(leg_params) else None
        mod_val = mod.get(field)
        
        if leg_val is None:
            errors.append(f"Key {key} missing param[{i}] in legacy")
            continue
        if mod_val is None:
            errors.append(f"Key {key} missing {field} in modern")
            continue
            
        compare_scalar(leg_val, mod_val, field, key, errors, tolerance)


def _report_missing(
    leg_lookup: Dict[Tuple[int, int], Any],
    mod_lookup: Dict[Tuple[int, int], Any],
    errors: List[str]
) -> None:
    """Report missing pair indices."""
    missing = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for key in list(missing)[:5]:
        errors.append(f"Missing pair {key} in modern")
    
    extra = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for key in list(extra)[:5]:
        errors.append(f"Extra pair {key} in modern")
