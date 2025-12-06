"""
Stage 2: Residue indices comparison.

Compares: legacy_residue_idx, start_atom_idx, end_atom_idx
"""

from typing import Any, Dict, List

from .base import CompareResult, build_residue_lookup, check_required_field, compare_scalar


def compare_residues(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 0.0  # Not used, but kept for consistent signature
) -> CompareResult:
    """
    Compare residue_indices records (Stage 2).
    
    Args:
        legacy_records: Legacy residue records
        modern_records: Modern residue records
        tolerance: Unused (all comparisons are exact)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = build_residue_lookup(modern_records, errors, "modern")
    
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Count mismatch: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    _compare_common_residues(leg_lookup, mod_lookup, errors)
    _report_missing_residues(leg_lookup, mod_lookup, errors)
    
    return len(errors) == 0, errors


def _compare_common_residues(
    leg_lookup: Dict,
    mod_lookup: Dict,
    errors: List[str]
) -> None:
    """Compare residues present in both datasets."""
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg = leg_lookup[key]
        mod = mod_lookup[key]
        
        _compare_residue_idx(key, leg, mod, errors)
        _compare_atom_range(key, leg, mod, errors)


def _compare_residue_idx(key: Any, leg: Dict, mod: Dict, errors: List[str]) -> None:
    """Compare legacy_residue_idx field."""
    leg_idx = leg.get("residue_idx") or leg.get("legacy_residue_idx")
    mod_idx = mod.get("legacy_residue_idx")
    
    if leg_idx is None:
        errors.append(f"Key {key} legacy missing residue_idx")
    elif mod_idx is None:
        errors.append(f"Key {key} modern missing legacy_residue_idx")
    elif leg_idx != mod_idx:
        errors.append(f"Key {key} residue_idx: {leg_idx} vs {mod_idx}")


def _compare_atom_range(key: Any, leg: Dict, mod: Dict, errors: List[str]) -> None:
    """Compare start_atom_idx and end_atom_idx fields."""
    # start_atom_idx
    leg_start = leg.get("start_atom_idx") or leg.get("legacy_start_atom_idx")
    mod_start = mod.get("legacy_start_atom_idx")
    if leg_start is None:
        errors.append(f"Key {key} legacy missing start_atom_idx")
    elif mod_start is None:
        errors.append(f"Key {key} modern missing start_atom_idx")
    elif leg_start != mod_start:
        errors.append(f"Key {key} start_atom_idx: {leg_start} vs {mod_start}")
    
    # end_atom_idx
    leg_end = leg.get("end_atom_idx") or leg.get("legacy_end_atom_idx")
    mod_end = mod.get("legacy_end_atom_idx")
    if leg_end is None:
        errors.append(f"Key {key} legacy missing end_atom_idx")
    elif mod_end is None:
        errors.append(f"Key {key} modern missing end_atom_idx")
    elif leg_end != mod_end:
        errors.append(f"Key {key} end_atom_idx: {leg_end} vs {mod_end}")


def _report_missing_residues(
    leg_lookup: Dict,
    mod_lookup: Dict,
    errors: List[str]
) -> None:
    """Report residues missing from either dataset."""
    missing_in_mod = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for key in list(missing_in_mod)[:5]:
        errors.append(f"Missing in modern: {key}")
    
    extra_in_mod = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for key in list(extra_in_mod)[:5]:
        errors.append(f"Extra in modern: {key}")

