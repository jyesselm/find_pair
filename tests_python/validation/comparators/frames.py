"""
Stages 3, 4, 5: Frame calculation comparisons.

Stage 3: base_frame_calc - base_type, rms_fit, num_matched_atoms, matched_atoms
Stage 4: ls_fitting - rms_fit, num_points, rotation_matrix, translation
Stage 5: frame_calc - base_type, rms_fit, num_matched_atoms, matched_coordinates
"""

from typing import Any, Dict, List, Set

from .base import (
    CompareResult, check_required_field,
    compare_scalar, compare_vector, compare_matrix
)
from ..config import Tolerance


def compare_base_frame_calc(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.RMS
) -> CompareResult:
    """
    Compare base_frame_calc records (Stage 3).
    
    Uses residue_idx as matching key.
    Compares: base_type, rms_fit, num_matched_atoms, matched_atoms
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    for idx in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[idx], mod_lookup[idx]
        _compare_base_frame(idx, leg, mod, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _build_residue_lookup(
    records: List[Dict[str, Any]],
    errors: List[str],
    label: str
) -> Dict[int, Dict[str, Any]]:
    """Build lookup by residue_idx."""
    lookup = {}
    for rec in records:
        # Legacy uses residue_idx, modern uses legacy_residue_idx
        idx = rec.get("residue_idx") or rec.get("legacy_residue_idx")
        if idx is None:
            errors.append(f"{label} record missing residue_idx")
            continue
        lookup[idx] = rec
    return lookup


def _compare_base_frame(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single base_frame_calc record."""
    key = f"residue {idx}"
    
    # base_type
    leg_type, mod_type = check_required_field(leg, mod, "base_type", key, errors)
    compare_scalar(leg_type, mod_type, "base_type", key, errors)
    
    # rms_fit
    leg_rms, mod_rms = check_required_field(leg, mod, "rms_fit", key, errors)
    compare_scalar(leg_rms, mod_rms, "rms_fit", key, errors, tolerance)
    
    # num_matched_atoms
    leg_num, mod_num = check_required_field(leg, mod, "num_matched_atoms", key, errors)
    compare_scalar(leg_num, mod_num, "num_matched_atoms", key, errors)
    
    # matched_atoms (set comparison)
    leg_atoms, mod_atoms = check_required_field(leg, mod, "matched_atoms", key, errors)
    if leg_atoms and mod_atoms:
        _compare_atom_sets(key, set(leg_atoms), set(mod_atoms), errors)


def _compare_atom_sets(key: Any, leg: Set[str], mod: Set[str], errors: List[str]) -> None:
    """Compare matched_atoms sets."""
    if leg != mod:
        only_leg = leg - mod
        only_mod = mod - leg
        errors.append(f"Key {key} matched_atoms differ: only_legacy={only_leg}, only_modern={only_mod}")


def compare_ls_fitting(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.RMS
) -> CompareResult:
    """
    Compare ls_fitting records (Stage 4).
    
    Uses residue_idx as matching key.
    Compares: rms_fit, num_points, rotation_matrix, translation
    Note: Legacy does NOT have base_type field.
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    for idx in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[idx], mod_lookup[idx]
        _compare_ls_fit(idx, leg, mod, errors, tolerance)
    
    return len(errors) == 0, errors


def _compare_ls_fit(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single ls_fitting record."""
    key = f"residue {idx}"
    
    # rms_fit
    leg_rms, mod_rms = check_required_field(leg, mod, "rms_fit", key, errors)
    compare_scalar(leg_rms, mod_rms, "rms_fit", key, errors, tolerance)
    
    # num_points
    leg_num, mod_num = check_required_field(leg, mod, "num_points", key, errors)
    compare_scalar(leg_num, mod_num, "num_points", key, errors)
    
    # rotation_matrix
    leg_rot, mod_rot = check_required_field(leg, mod, "rotation_matrix", key, errors)
    compare_matrix(leg_rot, mod_rot, "rotation_matrix", key, errors, Tolerance.MATRIX)
    
    # translation
    leg_trans, mod_trans = check_required_field(leg, mod, "translation", key, errors)
    compare_vector(leg_trans, mod_trans, "translation", key, errors, Tolerance.COORDINATE)


def compare_frame_calc(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.MATRIX
) -> CompareResult:
    """
    Compare frame_calc records (Stage 5).
    
    Uses residue_idx as matching key.
    Compares: base_type, rms_fit, num_matched_atoms
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    for idx in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[idx], mod_lookup[idx]
        _compare_frame(idx, leg, mod, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _compare_frame(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single frame_calc record."""
    key = f"residue {idx}"
    
    # base_type
    leg_type, mod_type = check_required_field(leg, mod, "base_type", key, errors)
    compare_scalar(leg_type, mod_type, "base_type", key, errors)
    
    # rms_fit
    leg_rms, mod_rms = check_required_field(leg, mod, "rms_fit", key, errors)
    compare_scalar(leg_rms, mod_rms, "rms_fit", key, errors, Tolerance.RMS)
    
    # num_matched_atoms
    leg_num, mod_num = check_required_field(leg, mod, "num_matched_atoms", key, errors)
    compare_scalar(leg_num, mod_num, "num_matched_atoms", key, errors)


def _report_missing(leg_lookup: Dict, mod_lookup: Dict, errors: List[str]) -> None:
    """Report missing keys."""
    missing = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for idx in list(missing)[:5]:
        errors.append(f"Missing in modern: residue {idx}")
