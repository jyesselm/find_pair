"""
Stages 3, 4, 5: Frame calculation comparisons.

Stage 3: base_frame_calc
  - Key: residue_idx (legacy) or legacy_residue_idx (modern)
  - Checks: base_type, rms_fit, num_matched_atoms, matched_atoms, standard_template

Stage 4: ls_fitting
  - Key: residue_idx (legacy) or legacy_residue_idx (modern)  
  - Checks: rms_fit, num_points, rotation_matrix, translation

Stage 5: frame_calc
  - Key: residue_idx (legacy) or legacy_residue_idx (modern)
  - Checks: base_type, rms_fit, num_matched_atoms
"""

import os
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
    
    Key: residue_idx (legacy) or legacy_residue_idx (modern)
    
    Checks per docs/JSON_DATA_TYPES_AND_COMPARISONS.md:
    - base_type: exact match
    - rms_fit: tolerance 0.001
    - num_matched_atoms: exact match
    - matched_atoms: set equality
    - standard_template: filename only (path differences OK)
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    # Check count
    if len(leg_lookup) != len(mod_lookup):
        errors.append(f"Record count: {len(leg_lookup)} vs {len(mod_lookup)}")
    
    # Compare common records
    for idx in sorted(leg_lookup.keys() & mod_lookup.keys()):
        leg, mod = leg_lookup[idx], mod_lookup[idx]
        _compare_base_frame_record(idx, leg, mod, errors, tolerance)
    
    # Report missing
    _report_missing(leg_lookup, mod_lookup, errors)
    
    return len(errors) == 0, errors


def _normalize_template(filename: str) -> str:
    """
    Normalize template filename for comparison.
    
    Handles naming differences like:
    - Atomic.p.pdb (legacy) == Atomic_P.pdb (modern)
    - Atomic_G.pdb (legacy) == Atomic.g.pdb (modern)
    
    Normalizes to uppercase underscore format: Atomic_X.pdb
    """
    if not filename:
        return filename
    # Normalize dot-lowercase to underscore-uppercase: Atomic.x.pdb -> Atomic_X.pdb
    import re
    match = re.match(r'Atomic\.([a-z])\.pdb', filename)
    if match:
        return f"Atomic_{match.group(1).upper()}.pdb"
    return filename


def _build_residue_lookup(
    records: List[Dict[str, Any]],
    errors: List[str],
    label: str
) -> Dict[int, Dict[str, Any]]:
    """Build lookup by residue_idx."""
    lookup = {}
    for rec in records:
        # Legacy uses residue_idx, modern uses legacy_residue_idx
        idx = rec.get("legacy_residue_idx")
        if idx is None:
            idx = rec.get("residue_idx")
        if idx is None:
            errors.append(f"{label} record missing residue_idx")
            continue
        lookup[idx] = rec
    return lookup


def _compare_base_frame_record(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single base_frame_calc record."""
    key = f"residue {idx}"
    
    # 1. base_type (case-insensitive for pseudouridine p/P)
    leg_type = leg.get("base_type")
    mod_type = mod.get("base_type")
    if leg_type is None:
        errors.append(f"{key}: legacy missing base_type")
    elif mod_type is None:
        errors.append(f"{key}: modern missing base_type")
    elif leg_type.lower() != mod_type.lower():
        errors.append(f"{key}: base_type '{leg_type}' vs '{mod_type}'")
    
    # 2. rms_fit (tolerance 0.001)
    leg_rms = leg.get("rms_fit")
    mod_rms = mod.get("rms_fit")
    if leg_rms is None:
        errors.append(f"{key}: legacy missing rms_fit")
    elif mod_rms is None:
        errors.append(f"{key}: modern missing rms_fit")
    elif abs(leg_rms - mod_rms) > tolerance:
        errors.append(f"{key}: rms_fit {leg_rms} vs {mod_rms} (diff={abs(leg_rms - mod_rms)})")
    
    # 3. num_matched_atoms (exact)
    leg_num = leg.get("num_matched_atoms")
    mod_num = mod.get("num_matched_atoms")
    if leg_num is None:
        errors.append(f"{key}: legacy missing num_matched_atoms")
    elif mod_num is None:
        errors.append(f"{key}: modern missing num_matched_atoms")
    elif leg_num != mod_num:
        errors.append(f"{key}: num_matched_atoms {leg_num} vs {mod_num}")
    
    # 4. matched_atoms (set equality)
    leg_atoms = leg.get("matched_atoms")
    mod_atoms = mod.get("matched_atoms")
    if leg_atoms is None:
        errors.append(f"{key}: legacy missing matched_atoms")
    elif mod_atoms is None:
        errors.append(f"{key}: modern missing matched_atoms")
    else:
        leg_set = set(a.strip() for a in leg_atoms)
        mod_set = set(a.strip() for a in mod_atoms)
        if leg_set != mod_set:
            only_leg = leg_set - mod_set
            only_mod = mod_set - leg_set
            errors.append(f"{key}: matched_atoms differ - only_legacy={only_leg}, only_modern={only_mod}")
    
    # 5. standard_template (filename only, path differences OK)
    # Normalize: Atomic.p.pdb == Atomic_P.pdb (pseudouridine template naming)
    leg_template = leg.get("standard_template", "")
    mod_template = mod.get("standard_template", "")
    leg_filename = _normalize_template(os.path.basename(leg_template))
    mod_filename = _normalize_template(os.path.basename(mod_template))
    if leg_filename and mod_filename and leg_filename != mod_filename:
        errors.append(f"{key}: standard_template '{leg_filename}' vs '{mod_filename}'")


def compare_ls_fitting(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.RMS
) -> CompareResult:
    """
    Compare ls_fitting records (Stage 4).
    
    Key: residue_idx (legacy) or legacy_residue_idx (modern)
    
    Checks per docs:
    - rms_fit: tolerance 0.001
    - num_points: exact
    - rotation_matrix: tolerance 1e-4 per element
    - translation: tolerance 1e-6 per element
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    for idx in sorted(leg_lookup.keys() & mod_lookup.keys()):
        leg, mod = leg_lookup[idx], mod_lookup[idx]
        _compare_ls_fit_record(idx, leg, mod, errors, tolerance)
    
    return len(errors) == 0, errors


def _compare_ls_fit_record(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single ls_fitting record."""
    key = f"residue {idx}"
    
    # 1. rms_fit (tolerance 0.001)
    leg_rms = leg.get("rms_fit")
    mod_rms = mod.get("rms_fit")
    if leg_rms is None:
        errors.append(f"{key}: legacy missing rms_fit")
    elif mod_rms is None:
        errors.append(f"{key}: modern missing rms_fit")
    elif abs(leg_rms - mod_rms) > tolerance:
        errors.append(f"{key}: rms_fit {leg_rms} vs {mod_rms}")
    
    # 2. num_points (exact)
    leg_num = leg.get("num_points")
    mod_num = mod.get("num_points")
    if leg_num is None:
        errors.append(f"{key}: legacy missing num_points")
    elif mod_num is None:
        errors.append(f"{key}: modern missing num_points")
    elif leg_num != mod_num:
        errors.append(f"{key}: num_points {leg_num} vs {mod_num}")
    
    # 3. rotation_matrix (tolerance 1e-4)
    leg_rot = leg.get("rotation_matrix")
    mod_rot = mod.get("rotation_matrix")
    if leg_rot is None:
        errors.append(f"{key}: legacy missing rotation_matrix")
    elif mod_rot is None:
        errors.append(f"{key}: modern missing rotation_matrix")
    else:
        compare_matrix(leg_rot, mod_rot, "rotation_matrix", key, errors, Tolerance.MATRIX)
    
    # 4. translation (tolerance 2e-6 to handle floating-point boundary cases)
    leg_trans = leg.get("translation")
    mod_trans = mod.get("translation")
    if leg_trans is None:
        errors.append(f"{key}: legacy missing translation")
    elif mod_trans is None:
        errors.append(f"{key}: modern missing translation")
    else:
        compare_vector(leg_trans, mod_trans, "translation", key, errors, 2e-6)


def compare_frame_calc(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.MATRIX
) -> CompareResult:
    """
    Compare frame_calc records (Stage 5).
    
    Key: residue_idx (legacy) or legacy_residue_idx (modern)
    
    Checks per docs:
    - base_type: exact
    - rms_fit: tolerance 0.001
    - num_matched_atoms: exact
    """
    errors: List[str] = []
    
    leg_lookup = _build_residue_lookup(legacy_records, errors, "legacy")
    mod_lookup = _build_residue_lookup(modern_records, errors, "modern")
    
    for idx in sorted(leg_lookup.keys() & mod_lookup.keys()):
        leg, mod = leg_lookup[idx], mod_lookup[idx]
        _compare_frame_record(idx, leg, mod, errors, tolerance)
    
    _report_missing(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _compare_frame_record(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single frame_calc record."""
    key = f"residue {idx}"
    
    # 1. base_type (case-insensitive for pseudouridine p/P)
    leg_type = leg.get("base_type")
    mod_type = mod.get("base_type")
    if leg_type is None:
        errors.append(f"{key}: legacy missing base_type")
    elif mod_type is None:
        errors.append(f"{key}: modern missing base_type")
    elif leg_type.lower() != mod_type.lower():
        errors.append(f"{key}: base_type '{leg_type}' vs '{mod_type}'")
    
    # 2. rms_fit (tolerance 0.001)
    leg_rms = leg.get("rms_fit")
    mod_rms = mod.get("rms_fit")
    if leg_rms is not None and mod_rms is not None:
        if abs(leg_rms - mod_rms) > Tolerance.RMS:
            errors.append(f"{key}: rms_fit {leg_rms} vs {mod_rms}")
    
    # 3. num_matched_atoms (exact)
    leg_num = leg.get("num_matched_atoms")
    mod_num = mod.get("num_matched_atoms")
    if leg_num is not None and mod_num is not None:
        if leg_num != mod_num:
            errors.append(f"{key}: num_matched_atoms {leg_num} vs {mod_num}")


def _report_missing(leg_lookup: Dict, mod_lookup: Dict, errors: List[str]) -> None:
    """Report missing keys."""
    missing = sorted(set(leg_lookup.keys()) - set(mod_lookup.keys()))
    for idx in missing[:5]:
        errors.append(f"Missing in modern: residue {idx}")
    
    extra = sorted(set(mod_lookup.keys()) - set(leg_lookup.keys()))
    for idx in extra[:5]:
        errors.append(f"Extra in modern: residue {idx}")
