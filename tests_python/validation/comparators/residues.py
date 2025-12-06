"""
Stage 2: Residue indices comparison.

Compares: residue_idx, start_atom, end_atom

Note: The actual JSON format uses simple index-based format without
chain_id/residue_seq fields. Records are matched by residue_idx.
"""

from typing import Any, Dict, List

from .base import CompareResult


def compare_residues(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 0.0  # Not used, all comparisons are exact
) -> CompareResult:
    """
    Compare residue_indices records (Stage 2).
    
    Args:
        legacy_records: Legacy residue records (from seidx array)
        modern_records: Modern residue records (from seidx array)
        tolerance: Unused (all comparisons are exact)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Count mismatch: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    _compare_common_residues(leg_lookup, mod_lookup, errors)
    _report_missing_residues(leg_lookup, mod_lookup, errors)
    
    return len(errors) == 0, errors


def _build_residue_lookup(
    records: List[Dict[str, Any]],
    errors: List[str],
    label: str
) -> Dict[int, Dict[str, Any]]:
    """Build lookup dict keyed by residue_idx."""
    lookup = {}
    for rec in records:
        idx = rec.get("residue_idx")
        if idx is None:
            errors.append(f"{label} record missing residue_idx")
            continue
        lookup[idx] = rec
    return lookup


def _compare_common_residues(
    leg_lookup: Dict[int, Dict[str, Any]],
    mod_lookup: Dict[int, Dict[str, Any]],
    errors: List[str]
) -> None:
    """Compare residues present in both datasets."""
    for idx in leg_lookup.keys() & mod_lookup.keys():
        leg = leg_lookup[idx]
        mod = mod_lookup[idx]
        
        _compare_atom_range(idx, leg, mod, errors)


def _compare_atom_range(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str]
) -> None:
    """Compare start_atom and end_atom fields."""
    # start_atom
    leg_start = leg.get("start_atom")
    mod_start = mod.get("start_atom")
    if leg_start is None:
        errors.append(f"Residue {idx} legacy missing start_atom")
    elif mod_start is None:
        errors.append(f"Residue {idx} modern missing start_atom")
    elif leg_start != mod_start:
        errors.append(f"Residue {idx} start_atom: {leg_start} vs {mod_start}")
    
    # end_atom
    leg_end = leg.get("end_atom")
    mod_end = mod.get("end_atom")
    if leg_end is None:
        errors.append(f"Residue {idx} legacy missing end_atom")
    elif mod_end is None:
        errors.append(f"Residue {idx} modern missing end_atom")
    elif leg_end != mod_end:
        errors.append(f"Residue {idx} end_atom: {leg_end} vs {mod_end}")


def _report_missing_residues(
    leg_lookup: Dict[int, Any],
    mod_lookup: Dict[int, Any],
    errors: List[str]
) -> None:
    """Report residues missing from either dataset."""
    missing_in_mod = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for idx in list(missing_in_mod)[:5]:
        errors.append(f"Missing in modern: residue {idx}")
    
    extra_in_mod = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for idx in list(extra_in_mod)[:5]:
        errors.append(f"Extra in modern: residue {idx}")
