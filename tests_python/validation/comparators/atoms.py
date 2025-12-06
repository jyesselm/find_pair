"""
Stage 1: Atom comparison (pdb_atoms).

Compares: xyz, atom_name, residue_name, chain_id
"""

from typing import Any, Dict, List

from .base import CompareResult, check_required_field, compare_vector
from ..config import Tolerance


def compare_atoms(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.COORDINATE
) -> CompareResult:
    """
    Compare pdb_atoms records (Stage 1).
    
    Args:
        legacy_records: Legacy atom records
        modern_records: Modern atom records
        tolerance: Coordinate tolerance (default 1e-6)
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    legacy_by_idx = _build_atom_index(legacy_records, "legacy", errors)
    modern_by_idx = _build_atom_index(modern_records, "modern", errors)
    
    if len(legacy_by_idx) != len(modern_by_idx):
        errors.append(f"Count mismatch: {len(legacy_by_idx)} vs {len(modern_by_idx)}")
    
    for idx, legacy_atom in legacy_by_idx.items():
        if idx not in modern_by_idx:
            errors.append(f"Atom {idx} missing in modern")
            continue
        _compare_single_atom(idx, legacy_atom, modern_by_idx[idx], errors, tolerance)
    
    return len(errors) == 0, errors


def _build_atom_index(
    records: List[Dict[str, Any]],
    label: str,
    errors: List[str]
) -> Dict[int, Dict[str, Any]]:
    """Build lookup by legacy atom index."""
    lookup = {}
    for rec in records:
        idx = rec.get("legacy_atom_idx") or rec.get("atom_idx")
        if idx is None:
            errors.append(f"{label} atom missing index")
            continue
        lookup[idx] = rec
    return lookup


def _compare_single_atom(
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single atom's fields."""
    # xyz coordinates
    leg_xyz, mod_xyz = check_required_field(leg, mod, "xyz", idx, errors)
    compare_vector(leg_xyz, mod_xyz, "xyz", idx, errors, tolerance)
    
    # atom_name
    leg_name, mod_name = check_required_field(leg, mod, "atom_name", idx, errors)
    if leg_name and mod_name and leg_name.strip() != mod_name.strip():
        errors.append(f"Atom {idx} atom_name: '{leg_name}' vs '{mod_name}'")
    
    # residue_name
    leg_res, mod_res = check_required_field(leg, mod, "residue_name", idx, errors)
    if leg_res and mod_res and leg_res.strip() != mod_res.strip():
        errors.append(f"Atom {idx} residue_name: '{leg_res}' vs '{mod_res}'")
    
    # chain_id
    leg_chain, mod_chain = check_required_field(leg, mod, "chain_id", idx, errors)
    if leg_chain and mod_chain and leg_chain != mod_chain:
        errors.append(f"Atom {idx} chain_id: '{leg_chain}' vs '{mod_chain}'")

