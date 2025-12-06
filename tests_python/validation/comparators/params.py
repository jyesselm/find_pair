"""
Stages 11, 12: Parameter comparisons.

Stage 11: bpstep_params - shift, slide, rise, tilt, roll, twist
Stage 12: helical_params - x_displacement, y_displacement, rise, inclination, tip, twist
"""

from typing import Any, Dict, List

from .base import CompareResult, check_required_field, compare_scalar
from ..config import Tolerance


# Step parameter fields (Stage 11)
STEP_PARAM_FIELDS = ["shift", "slide", "rise", "tilt", "roll", "twist"]

# Helical parameter fields (Stage 12)
HELICAL_PARAM_FIELDS = ["x_displacement", "y_displacement", "rise", "inclination", "tip", "twist"]


def compare_step_params(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.ANGLE
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
    
    leg_lookup = _build_step_lookup(legacy_records)
    mod_lookup = _build_step_lookup(modern_records)
    
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Count mismatch: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_params(key, leg, mod, STEP_PARAM_FIELDS, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def compare_helical_params(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.ANGLE
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
    
    leg_lookup = _build_step_lookup(legacy_records)
    mod_lookup = _build_step_lookup(modern_records)
    
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Count mismatch: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_params(key, leg, mod, HELICAL_PARAM_FIELDS, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _build_step_lookup(records: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    """Build lookup by step index."""
    lookup = {}
    for rec in records:
        # Use explicit None check - 0 is a valid index!
        idx = rec.get("step_idx")
        if idx is None:
            idx = rec.get("idx")
        if idx is not None:
            lookup[idx] = rec
    return lookup


def _compare_params(
    key: Any,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    fields: List[str],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare parameter fields."""
    for field in fields:
        leg_val, mod_val = check_required_field(leg, mod, field, key, errors)
        compare_scalar(leg_val, mod_val, field, key, errors, tolerance)


def _report_missing(
    leg_lookup: Dict[int, Any],
    mod_lookup: Dict[int, Any],
    errors: List[str]
) -> None:
    """Report missing step indices."""
    missing = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for key in list(missing)[:5]:
        errors.append(f"Missing step {key} in modern")
    
    extra = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for key in list(extra)[:5]:
        errors.append(f"Extra step {key} in modern")

