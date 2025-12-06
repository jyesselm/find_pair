#!/usr/bin/env python3
"""
Unified stage-by-stage validation for legacy vs modern JSON comparison.

This script validates each stage of the base pair finding algorithm by comparing
legacy and modern JSON outputs. It supports stop-on-first-failure, parallel
processing, and automatic cleanup of temp files.

Stage Reference:
    1. pdb_atoms           - Atom parsing
    2. residue_indices     - Residue index mapping
    3. base_frame_calc     - Base frame calculation
    4. ls_fitting          - Least squares fitting
    5. frame_calc          - Reference frame calculation
    6. pair_validation     - Pair validation
    7. distance_checks     - Distance measurements
    8. hbond_list          - Hydrogen bond list
    9. base_pair           - Base pair records
    10. find_bestpair_selection - Final pair selection (PRIMARY OUTPUT)
    11. bpstep_params      - Step parameters
    12. helical_params     - Helical parameters

Usage:
    # CLI mode
    python test_stage_validation.py 1 --max-pdbs 100        # Stage 1
    python test_stage_validation.py frames --max-pdbs 100   # Stages 3,4,5
    python test_stage_validation.py all --pdb 1EHZ          # All stages, single PDB

    # pytest mode
    pytest test_stage_validation.py -v -k "stage1"
    pytest test_stage_validation.py -v --max-pdbs=100
"""

import json
import pytest
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scripts.test_utils import find_executables, load_valid_pdbs_fast


# ============================================================================
# Stage Configuration
# ============================================================================

@dataclass
class StageConfig:
    """Configuration for a validation stage."""
    stage_num: int
    stage_id: str
    name: str
    description: str
    json_type: str
    key_fields: List[str]
    dependencies: List[int] = field(default_factory=list)
    tolerance: float = 1e-6


# Stage definitions aligned with docs/JSON_DATA_TYPES_AND_COMPARISONS.md
STAGES = {
    1: StageConfig(
        stage_num=1,
        stage_id="pdb_atoms",
        name="Atom Parsing",
        description="Parse PDB file atoms and coordinates",
        json_type="pdb_atoms",
        key_fields=["atom_idx", "xyz", "atom_name", "residue_name", "chain_id"],
        dependencies=[]
    ),
    2: StageConfig(
        stage_num=2,
        stage_id="residue_indices",
        name="Residue Indices",
        description="Map residues to atom index ranges",
        json_type="residue_indices",
        key_fields=["chain_id", "residue_seq", "insertion", "legacy_residue_idx"],
        dependencies=[1]
    ),
    3: StageConfig(
        stage_num=3,
        stage_id="base_frame_calc",
        name="Base Frame Calc",
        description="Calculate base frame via template matching",
        json_type="base_frame_calc",
        key_fields=["legacy_residue_idx", "rms_fit", "num_matched_atoms", "matched_atoms"],
        dependencies=[1, 2],
        tolerance=0.001
    ),
    4: StageConfig(
        stage_num=4,
        stage_id="ls_fitting",
        name="LS Fitting",
        description="Least squares fitting results",
        json_type="ls_fitting",
        key_fields=["legacy_residue_idx", "rms_fit", "num_points", "rotation_matrix", "translation"],
        dependencies=[1, 2],
        tolerance=0.001
    ),
    5: StageConfig(
        stage_num=5,
        stage_id="frame_calc",
        name="Frame Calc",
        description="Reference frame (rotation matrix + origin)",
        json_type="frame_calc",
        key_fields=["legacy_residue_idx", "orien", "org"],
        dependencies=[3, 4],
        tolerance=1e-4
    ),
    6: StageConfig(
        stage_num=6,
        stage_id="pair_validation",
        name="Pair Validation",
        description="Validate potential base pairs",
        json_type="pair_validation",
        key_fields=["base_i", "base_j", "is_valid", "bp_type_id", "quality_score"],
        dependencies=[5]
    ),
    7: StageConfig(
        stage_num=7,
        stage_id="distance_checks",
        name="Distance Checks",
        description="Geometric measurements between pairs",
        json_type="distance_checks",
        key_fields=["base_i", "base_j", "dorg", "dNN", "plane_angle", "d_v"],
        dependencies=[5]
    ),
    8: StageConfig(
        stage_num=8,
        stage_id="hbond_list",
        name="H-bond List",
        description="Hydrogen bond detection",
        json_type="hbond_list",
        key_fields=["base_i", "base_j", "num_hbonds", "hbonds"],
        dependencies=[6, 7]
    ),
    9: StageConfig(
        stage_num=9,
        stage_id="base_pair",
        name="Base Pair",
        description="Identified base pair records",
        json_type="base_pair",
        key_fields=["base_i", "base_j", "bp_type", "orien_i", "org_i"],
        dependencies=[6, 7, 8]
    ),
    10: StageConfig(
        stage_num=10,
        stage_id="find_bestpair_selection",
        name="Best Pair Selection",
        description="Final selected base pairs (PRIMARY OUTPUT)",
        json_type="find_bestpair_selection",
        key_fields=["num_bp", "pairs"],
        dependencies=[9]
    ),
    11: StageConfig(
        stage_num=11,
        stage_id="bpstep_params",
        name="Step Parameters",
        description="Base pair step parameters",
        json_type="bpstep_params",
        key_fields=["bp_idx1", "bp_idx2", "shift", "slide", "rise", "tilt", "roll", "twist"],
        dependencies=[10]
    ),
    12: StageConfig(
        stage_num=12,
        stage_id="helical_params",
        name="Helical Parameters",
        description="Helical axis parameters",
        json_type="helical_params",
        key_fields=["bp_idx1", "bp_idx2", "x_displacement", "y_displacement", "rise", "inclination", "tip", "twist"],
        dependencies=[10, 11]
    ),
}

# Stage groups for convenience
STAGE_GROUPS = {
    "atoms": [1],
    "residue": [2],
    "frames": [3, 4, 5],
    "pairs": [6, 7, 9, 10],
    "hbonds": [8],
    "steps": [11, 12],
    "all": list(range(1, 13)),
}

# Map stage IDs to numbers for lookup
STAGE_ID_TO_NUM = {cfg.stage_id: cfg.stage_num for cfg in STAGES.values()}


def resolve_stages(stage_args: List[str]) -> List[int]:
    """Resolve stage arguments to list of stage numbers."""
    if not stage_args:
        return list(range(1, 13))  # All stages
    
    stages = []
    for arg in stage_args:
        # Try as group name
        if arg.lower() in STAGE_GROUPS:
            stages.extend(STAGE_GROUPS[arg.lower()])
        # Try as stage number
        elif arg.isdigit():
            num = int(arg)
            if 1 <= num <= 12:
                stages.append(num)
        # Try as stage ID
        elif arg.lower() in STAGE_ID_TO_NUM:
            stages.append(STAGE_ID_TO_NUM[arg.lower()])
    
    return sorted(set(stages))


# ============================================================================
# Comparison Functions
# ============================================================================

def normalize_pair(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize a pair to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


def compare_values(legacy_val: Any, modern_val: Any, tolerance: float = 1e-6) -> Tuple[bool, str]:
    """Compare two values with tolerance for floats."""
    if legacy_val is None and modern_val is None:
        return True, ""
    if legacy_val is None or modern_val is None:
        return False, f"One value is None: legacy={legacy_val}, modern={modern_val}"
    
    if isinstance(legacy_val, (int, float)) and isinstance(modern_val, (int, float)):
        if abs(legacy_val - modern_val) <= tolerance:
            return True, ""
        return False, f"Value mismatch: legacy={legacy_val}, modern={modern_val}, diff={abs(legacy_val - modern_val)}"
    
    if isinstance(legacy_val, str) and isinstance(modern_val, str):
        if legacy_val.strip() == modern_val.strip():
            return True, ""
        return False, f"String mismatch: legacy='{legacy_val}', modern='{modern_val}'"
    
    if isinstance(legacy_val, list) and isinstance(modern_val, list):
        if len(legacy_val) != len(modern_val):
            return False, f"List length mismatch: legacy={len(legacy_val)}, modern={len(modern_val)}"
        for i, (l, m) in enumerate(zip(legacy_val, modern_val)):
            match, msg = compare_values(l, m, tolerance)
            if not match:
                return False, f"List element {i}: {msg}"
        return True, ""
    
    # Direct comparison for other types
    if legacy_val == modern_val:
        return True, ""
    return False, f"Type/value mismatch: legacy={legacy_val} ({type(legacy_val).__name__}), modern={modern_val} ({type(modern_val).__name__})"


def compare_matrix(legacy_mat: List[List[float]], modern_mat: List[List[float]], 
                   tolerance: float = 1e-4) -> Tuple[bool, str]:
    """Compare 3x3 rotation matrices."""
    if len(legacy_mat) != len(modern_mat):
        return False, f"Matrix dimension mismatch: {len(legacy_mat)} vs {len(modern_mat)}"
    
    max_diff = 0.0
    for i, (l_row, m_row) in enumerate(zip(legacy_mat, modern_mat)):
        for j, (l_val, m_val) in enumerate(zip(l_row, m_row)):
            diff = abs(l_val - m_val)
            max_diff = max(max_diff, diff)
            if diff > tolerance:
                return False, f"Matrix element [{i}][{j}] differs: {l_val} vs {m_val} (diff={diff})"
    
    return True, ""


# ============================================================================
# Stage-Specific Comparisons
# ============================================================================

def compare_pdb_atoms(legacy_records: List[Dict], modern_records: List[Dict],
                      tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare pdb_atoms records (Stage 1)."""
    errors = []
    
    # Handle nested structure (atoms array inside record)
    if legacy_records and isinstance(legacy_records[0], dict):
        if 'atoms' in legacy_records[0]:
            legacy_atoms = legacy_records[0].get('atoms', [])
        elif legacy_records[0].get('type') == 'pdb_atoms':
            legacy_atoms = legacy_records[0].get('atoms', [])
        else:
            legacy_atoms = legacy_records
    else:
        legacy_atoms = legacy_records
    
    if modern_records and isinstance(modern_records[0], dict):
        if 'atoms' in modern_records[0]:
            modern_atoms = modern_records[0].get('atoms', [])
        elif modern_records[0].get('type') == 'pdb_atoms':
            modern_atoms = modern_records[0].get('atoms', [])
        else:
            modern_atoms = modern_records
    else:
        modern_atoms = modern_records
    
    # Check count
    if len(legacy_atoms) != len(modern_atoms):
        errors.append(f"Atom count mismatch: legacy={len(legacy_atoms)}, modern={len(modern_atoms)}")
        return False, errors
    
    # Build lookup by legacy_atom_idx
    legacy_by_idx = {}
    for atom in legacy_atoms:
        idx = atom.get('atom_idx') or atom.get('legacy_atom_idx')
        if idx is not None:
            legacy_by_idx[idx] = atom
    
    modern_by_legacy_idx = {}
    for atom in modern_atoms:
        idx = atom.get('legacy_atom_idx')
        if idx is not None:
            modern_by_legacy_idx[idx] = atom
    
    # Compare atoms by index
    for idx, legacy_atom in legacy_by_idx.items():
        if idx not in modern_by_legacy_idx:
            errors.append(f"Atom {idx} missing in modern")
            continue
        
        modern_atom = modern_by_legacy_idx[idx]
        
        # Compare xyz (required, tolerance 1e-6)
        leg_xyz = legacy_atom.get('xyz')
        mod_xyz = modern_atom.get('xyz')
        if leg_xyz is None:
            errors.append(f"Atom {idx} legacy missing xyz")
        elif mod_xyz is None:
            errors.append(f"Atom {idx} modern missing xyz")
        elif len(leg_xyz) != 3 or len(mod_xyz) != 3:
            errors.append(f"Atom {idx} xyz wrong length")
        else:
            for i, (l, m) in enumerate(zip(leg_xyz, mod_xyz)):
                if abs(l - m) > tolerance:
                    errors.append(f"Atom {idx} xyz[{i}]: {l} vs {m}")
        
        # Compare atom_name (required, exact)
        leg_name = legacy_atom.get('atom_name')
        mod_name = modern_atom.get('atom_name')
        if leg_name is None:
            errors.append(f"Atom {idx} legacy missing atom_name")
        elif mod_name is None:
            errors.append(f"Atom {idx} modern missing atom_name")
        elif leg_name.strip() != mod_name.strip():
            errors.append(f"Atom {idx} atom_name: '{leg_name}' vs '{mod_name}'")
        
        # Compare residue_name (required, exact)
        leg_resname = legacy_atom.get('residue_name')
        mod_resname = modern_atom.get('residue_name')
        if leg_resname is None:
            errors.append(f"Atom {idx} legacy missing residue_name")
        elif mod_resname is None:
            errors.append(f"Atom {idx} modern missing residue_name")
        elif leg_resname.strip() != mod_resname.strip():
            errors.append(f"Atom {idx} residue_name: '{leg_resname}' vs '{mod_resname}'")
        
        # Compare chain_id (required, exact)
        leg_chain = legacy_atom.get('chain_id')
        mod_chain = modern_atom.get('chain_id')
        if leg_chain is None:
            errors.append(f"Atom {idx} legacy missing chain_id")
        elif mod_chain is None:
            errors.append(f"Atom {idx} modern missing chain_id")
        elif leg_chain != mod_chain:
            errors.append(f"Atom {idx} chain_id: '{leg_chain}' vs '{mod_chain}'")
    
    return len(errors) == 0, errors


def compare_residue_indices(legacy_records: List[Dict], modern_records: List[Dict], 
                            tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare residue indices records (Stage 2).
    
    Required fields: legacy_residue_idx, start_atom_idx, end_atom_idx
    """
    errors = []
    
    if len(legacy_records) != len(modern_records):
        errors.append(f"Count mismatch: legacy={len(legacy_records)}, modern={len(modern_records)}")
        # Don't return early - continue to find what's different
    
    # Build lookup by (chain_id, residue_seq, insertion)
    legacy_by_key = {}
    for rec in legacy_records:
        chain = rec.get('chain_id')
        seq = rec.get('residue_seq')
        if chain is None or seq is None:
            errors.append(f"Legacy record missing chain_id or residue_seq")
            continue
        key = (chain, seq, rec.get('insertion', ' '))
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        chain = rec.get('chain_id')
        seq = rec.get('residue_seq')
        if chain is None or seq is None:
            errors.append(f"Modern record missing chain_id or residue_seq")
            continue
        key = (chain, seq, rec.get('insertion', ' '))
        modern_by_key[key] = rec
    
    # Compare
    for key, legacy_rec in legacy_by_key.items():
        if key not in modern_by_key:
            errors.append(f"Missing in modern: {key}")
            continue
        
        modern_rec = modern_by_key[key]
        
        # Compare legacy_residue_idx (required, exact)
        leg_idx = legacy_rec.get('residue_idx') or legacy_rec.get('legacy_residue_idx')
        mod_idx = modern_rec.get('legacy_residue_idx')
        if leg_idx is None:
            errors.append(f"Key {key} legacy missing residue_idx")
        elif mod_idx is None:
            errors.append(f"Key {key} modern missing legacy_residue_idx")
        elif leg_idx != mod_idx:
            errors.append(f"Key {key} legacy_residue_idx: {leg_idx} vs {mod_idx}")
        
        # Compare start_atom_idx (required, exact)
        leg_start = legacy_rec.get('start_atom_idx') or legacy_rec.get('legacy_start_atom_idx')
        mod_start = modern_rec.get('legacy_start_atom_idx')
        if leg_start is None:
            errors.append(f"Key {key} legacy missing start_atom_idx")
        elif mod_start is None:
            errors.append(f"Key {key} modern missing legacy_start_atom_idx")
        elif leg_start != mod_start:
            errors.append(f"Key {key} start_atom_idx: {leg_start} vs {mod_start}")
        
        # Compare end_atom_idx (required, exact)
        leg_end = legacy_rec.get('end_atom_idx') or legacy_rec.get('legacy_end_atom_idx')
        mod_end = modern_rec.get('legacy_end_atom_idx')
        if leg_end is None:
            errors.append(f"Key {key} legacy missing end_atom_idx")
        elif mod_end is None:
            errors.append(f"Key {key} modern missing legacy_end_atom_idx")
        elif leg_end != mod_end:
            errors.append(f"Key {key} end_atom_idx: {leg_end} vs {mod_end}")
    
    # Check for extra keys in modern
    extra_in_modern = set(modern_by_key.keys()) - set(legacy_by_key.keys())
    if extra_in_modern:
        for key in list(extra_in_modern)[:5]:
            errors.append(f"Extra in modern: {key}")
    
    return len(errors) == 0, errors


def compare_base_frame_calc(legacy_records: List[Dict], modern_records: List[Dict],
                            tolerance: float = 0.001) -> Tuple[bool, List[str]]:
    """Compare base_frame_calc records (Stage 3).
    
    Required fields: base_type, rms_fit, num_matched_atoms, matched_atoms
    """
    errors = []
    
    # Build lookup by (chain_id, residue_seq, insertion)
    legacy_by_key = {}
    for rec in legacy_records:
        key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion', ' '))
        if None in key[:2]:
            errors.append(f"Legacy record missing chain_id or residue_seq: {rec}")
            continue
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion', ' '))
        if None in key[:2]:
            errors.append(f"Modern record missing chain_id or residue_seq: {rec}")
            continue
        modern_by_key[key] = rec
    
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())
    
    for key in common_keys:
        leg_rec = legacy_by_key[key]
        mod_rec = modern_by_key[key]
        
        # Compare base_type (required, critical for modified nucleotides)
        leg_type = leg_rec.get('base_type')
        mod_type = mod_rec.get('base_type')
        if leg_type is None:
            errors.append(f"Key {key} legacy missing base_type")
        elif mod_type is None:
            errors.append(f"Key {key} modern missing base_type")
        elif leg_type != mod_type:
            errors.append(f"Key {key} base_type: '{leg_type}' vs '{mod_type}'")
        
        # Compare rms_fit (required, tolerance 0.001)
        leg_rms = leg_rec.get('rms_fit')
        mod_rms = mod_rec.get('rms_fit')
        if leg_rms is None:
            errors.append(f"Key {key} legacy missing rms_fit")
        elif mod_rms is None:
            errors.append(f"Key {key} modern missing rms_fit")
        elif abs(leg_rms - mod_rms) > tolerance:
            errors.append(f"Key {key} rms_fit: {leg_rms} vs {mod_rms} (diff={abs(leg_rms - mod_rms)})")
        
        # Compare num_matched_atoms (required, exact)
        leg_num = leg_rec.get('num_matched_atoms')
        mod_num = mod_rec.get('num_matched_atoms')
        if leg_num is None:
            errors.append(f"Key {key} legacy missing num_matched_atoms")
        elif mod_num is None:
            errors.append(f"Key {key} modern missing num_matched_atoms")
        elif leg_num != mod_num:
            errors.append(f"Key {key} num_matched_atoms: {leg_num} vs {mod_num}")
        
        # Compare matched_atoms (required, set equality)
        leg_atoms = leg_rec.get('matched_atoms')
        mod_atoms = mod_rec.get('matched_atoms')
        if leg_atoms is None:
            errors.append(f"Key {key} legacy missing matched_atoms")
        elif mod_atoms is None:
            errors.append(f"Key {key} modern missing matched_atoms")
        elif set(leg_atoms) != set(mod_atoms):
            only_leg = set(leg_atoms) - set(mod_atoms)
            only_mod = set(mod_atoms) - set(leg_atoms)
            errors.append(f"Key {key} matched_atoms differ: only_legacy={only_leg}, only_modern={only_mod}")
    
    # Check for missing keys
    missing_in_modern = set(legacy_by_key.keys()) - common_keys
    if missing_in_modern:
        for key in list(missing_in_modern)[:5]:
            errors.append(f"Missing in modern: {key}")
    
    return len(errors) == 0, errors


def compare_ls_fitting(legacy_records: List[Dict], modern_records: List[Dict],
                       tolerance: float = 0.001) -> Tuple[bool, List[str]]:
    """Compare ls_fitting records (Stage 4).
    
    Required fields: base_type, rms_fit, num_points, rotation_matrix, translation
    """
    errors = []
    
    # Build lookup by (chain_id, residue_seq, insertion)
    legacy_by_key = {}
    for rec in legacy_records:
        key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion', ' '))
        if None in key[:2]:
            errors.append(f"Legacy record missing chain_id or residue_seq")
            continue
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion', ' '))
        if None in key[:2]:
            errors.append(f"Modern record missing chain_id or residue_seq")
            continue
        modern_by_key[key] = rec
    
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())
    
    for key in common_keys:
        leg_rec = legacy_by_key[key]
        mod_rec = modern_by_key[key]
        
        # Compare base_type (required)
        leg_type = leg_rec.get('base_type')
        mod_type = mod_rec.get('base_type')
        if leg_type is None:
            errors.append(f"Key {key} legacy missing base_type")
        elif mod_type is None:
            errors.append(f"Key {key} modern missing base_type")
        elif leg_type != mod_type:
            errors.append(f"Key {key} base_type: '{leg_type}' vs '{mod_type}'")
        
        # Compare rms_fit (required, tolerance 0.001)
        leg_rms = leg_rec.get('rms_fit')
        mod_rms = mod_rec.get('rms_fit')
        if leg_rms is None:
            errors.append(f"Key {key} legacy missing rms_fit")
        elif mod_rms is None:
            errors.append(f"Key {key} modern missing rms_fit")
        elif abs(leg_rms - mod_rms) > tolerance:
            errors.append(f"Key {key} rms_fit: {leg_rms} vs {mod_rms}")
        
        # Compare num_points (required, exact)
        leg_num = leg_rec.get('num_points')
        mod_num = mod_rec.get('num_points')
        if leg_num is None:
            errors.append(f"Key {key} legacy missing num_points")
        elif mod_num is None:
            errors.append(f"Key {key} modern missing num_points")
        elif leg_num != mod_num:
            errors.append(f"Key {key} num_points: {leg_num} vs {mod_num}")
        
        # Compare rotation_matrix (required, tolerance 1e-4 per element)
        leg_rot = leg_rec.get('rotation_matrix')
        mod_rot = mod_rec.get('rotation_matrix')
        if leg_rot is None:
            errors.append(f"Key {key} legacy missing rotation_matrix")
        elif mod_rot is None:
            errors.append(f"Key {key} modern missing rotation_matrix")
        else:
            match, msg = compare_matrix(leg_rot, mod_rot, tolerance=1e-4)
            if not match:
                errors.append(f"Key {key} rotation_matrix: {msg}")
        
        # Compare translation (required, tolerance 1e-6 per element)
        leg_trans = leg_rec.get('translation')
        mod_trans = mod_rec.get('translation')
        if leg_trans is None:
            errors.append(f"Key {key} legacy missing translation")
        elif mod_trans is None:
            errors.append(f"Key {key} modern missing translation")
        elif len(leg_trans) != 3 or len(mod_trans) != 3:
            errors.append(f"Key {key} translation wrong length: {len(leg_trans)} vs {len(mod_trans)}")
        else:
            max_diff = max(abs(leg_trans[i] - mod_trans[i]) for i in range(3))
            if max_diff > 1e-6:
                errors.append(f"Key {key} translation max_diff: {max_diff}")
    
    return len(errors) == 0, errors


def compare_frame_calc(legacy_records: List[Dict], modern_records: List[Dict],
                       tolerance: float = 1e-4) -> Tuple[bool, List[str]]:
    """Compare frame_calc records (Stage 5).
    
    Required fields: base_type, orien, org
    """
    errors = []
    
    # Build lookup by (chain_id, residue_seq, insertion)
    legacy_by_key = {}
    for rec in legacy_records:
        key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion', ' '))
        if None in key[:2]:
            errors.append(f"Legacy record missing chain_id or residue_seq")
            continue
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion', ' '))
        if None in key[:2]:
            errors.append(f"Modern record missing chain_id or residue_seq")
            continue
        modern_by_key[key] = rec
    
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())
    
    for key in common_keys:
        leg_rec = legacy_by_key[key]
        mod_rec = modern_by_key[key]
        
        # Compare base_type (required, critical for modified nucleotides)
        leg_type = leg_rec.get('base_type')
        mod_type = mod_rec.get('base_type')
        if leg_type is None:
            errors.append(f"Key {key} legacy missing base_type")
        elif mod_type is None:
            errors.append(f"Key {key} modern missing base_type")
        elif leg_type != mod_type:
            errors.append(f"Key {key} base_type: '{leg_type}' vs '{mod_type}'")
        
        # Compare orien (required, rotation matrix, tolerance 1e-4)
        leg_orien = leg_rec.get('orien')
        mod_orien = mod_rec.get('orien')
        if leg_orien is None:
            errors.append(f"Key {key} legacy missing orien")
        elif mod_orien is None:
            errors.append(f"Key {key} modern missing orien")
        else:
            match, msg = compare_matrix(leg_orien, mod_orien, tolerance)
            if not match:
                errors.append(f"Key {key} orien: {msg}")
        
        # Compare org (required, origin coordinates, tolerance 1e-6)
        leg_org = leg_rec.get('org')
        mod_org = mod_rec.get('org')
        if leg_org is None:
            errors.append(f"Key {key} legacy missing org")
        elif mod_org is None:
            errors.append(f"Key {key} modern missing org")
        elif len(leg_org) != 3 or len(mod_org) != 3:
            errors.append(f"Key {key} org wrong length: {len(leg_org)} vs {len(mod_org)}")
        else:
            max_diff = max(abs(leg_org[i] - mod_org[i]) for i in range(3))
            if max_diff > 1e-6:
                errors.append(f"Key {key} org max_diff: {max_diff}")
    
    # Check for missing keys
    missing_in_modern = set(legacy_by_key.keys()) - common_keys
    if missing_in_modern:
        for key in list(missing_in_modern)[:5]:
            errors.append(f"Missing in modern: {key}")
    
    return len(errors) == 0, errors


def compare_distance_checks(legacy_records: List[Dict], modern_records: List[Dict],
                           tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare distance checks records (Stage 7)."""
    errors = []
    
    # Build lookup by normalized pair
    legacy_by_pair = {}
    for rec in legacy_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        legacy_by_pair[pair] = rec
    
    modern_by_pair = {}
    for rec in modern_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        modern_by_pair[pair] = rec
    
    common_pairs = set(legacy_by_pair.keys()) & set(modern_by_pair.keys())
    
    if not common_pairs and (legacy_by_pair or modern_by_pair):
        errors.append(f"No common pairs found (legacy={len(legacy_by_pair)}, modern={len(modern_by_pair)})")
        return False, errors
    
    fields_to_compare = ['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area']
    
    for pair in common_pairs:
        legacy_rec = legacy_by_pair[pair]
        modern_rec = modern_by_pair[pair]
        
        for field in fields_to_compare:
            leg_val = legacy_rec.get(field)
            mod_val = modern_rec.get(field)
            
            if leg_val is None or mod_val is None:
                continue
            
            match, msg = compare_values(leg_val, mod_val, tolerance)
            if not match:
                errors.append(f"Pair {pair} field {field}: {msg}")
    
    return len(errors) == 0, errors


def compare_hbond_list(legacy_records: List[Dict], modern_records: List[Dict],
                       tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare H-bond list records (Stage 8)."""
    errors = []
    
    # Build lookup by normalized pair
    legacy_by_pair = {}
    for rec in legacy_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        legacy_by_pair[pair] = rec
    
    modern_by_pair = {}
    for rec in modern_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        modern_by_pair[pair] = rec
    
    common_pairs = set(legacy_by_pair.keys()) & set(modern_by_pair.keys())
    
    for pair in common_pairs:
        legacy_rec = legacy_by_pair[pair]
        modern_rec = modern_by_pair[pair]
        
        # Compare num_hbonds
        leg_count = legacy_rec.get('num_hbonds', 0)
        mod_count = modern_rec.get('num_hbonds', 0)
        
        if leg_count != mod_count:
            errors.append(f"Pair {pair} num_hbonds: legacy={leg_count}, modern={mod_count}")
            continue
        
        # Compare individual H-bonds
        leg_hbonds = legacy_rec.get('hbonds', [])
        mod_hbonds = modern_rec.get('hbonds', [])
        
        if len(leg_hbonds) != len(mod_hbonds):
            errors.append(f"Pair {pair} hbonds count: legacy={len(leg_hbonds)}, modern={len(mod_hbonds)}")
            continue
        
        for i, (leg_hb, mod_hb) in enumerate(zip(leg_hbonds, mod_hbonds)):
            # Compare donor, acceptor, distance
            if leg_hb.get('donor_atom', '').strip() != mod_hb.get('donor_atom', '').strip():
                errors.append(f"Pair {pair} hbond {i} donor: {leg_hb.get('donor_atom')} vs {mod_hb.get('donor_atom')}")
            if leg_hb.get('acceptor_atom', '').strip() != mod_hb.get('acceptor_atom', '').strip():
                errors.append(f"Pair {pair} hbond {i} acceptor: {leg_hb.get('acceptor_atom')} vs {mod_hb.get('acceptor_atom')}")
            
            leg_dist = leg_hb.get('distance', 0.0)
            mod_dist = mod_hb.get('distance', 0.0)
            if abs(leg_dist - mod_dist) > tolerance:
                errors.append(f"Pair {pair} hbond {i} distance: {leg_dist} vs {mod_dist}")
    
    return len(errors) == 0, errors


def compare_pair_validation(legacy_records: List[Dict], modern_records: List[Dict],
                            tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare pair validation records (Stage 6)."""
    errors = []
    
    # Build lookup by normalized pair
    legacy_by_pair = {}
    for rec in legacy_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        legacy_by_pair[pair] = rec
    
    modern_by_pair = {}
    for rec in modern_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        modern_by_pair[pair] = rec
    
    common_pairs = set(legacy_by_pair.keys()) & set(modern_by_pair.keys())
    
    for pair in common_pairs:
        legacy_rec = legacy_by_pair[pair]
        modern_rec = modern_by_pair[pair]
        
        # Compare is_valid
        if legacy_rec.get('is_valid') != modern_rec.get('is_valid'):
            errors.append(f"Pair {pair} is_valid: legacy={legacy_rec.get('is_valid')}, modern={modern_rec.get('is_valid')}")
        
        # Compare bp_type_id
        if legacy_rec.get('bp_type_id') != modern_rec.get('bp_type_id'):
            errors.append(f"Pair {pair} bp_type_id: legacy={legacy_rec.get('bp_type_id')}, modern={modern_rec.get('bp_type_id')}")
        
        # Compare calculated values
        leg_vals = legacy_rec.get('calculated_values', {})
        mod_vals = modern_rec.get('calculated_values', {})
        
        for field in ['dorg', 'd_v', 'plane_angle', 'dNN', 'quality_score']:
            leg_val = leg_vals.get(field)
            mod_val = mod_vals.get(field)
            if leg_val is not None and mod_val is not None:
                match, msg = compare_values(leg_val, mod_val, tolerance)
                if not match:
                    errors.append(f"Pair {pair} {field}: {msg}")
    
    return len(errors) == 0, errors


def compare_find_bestpair_selection(legacy_records: List[Dict], modern_records: List[Dict],
                                    tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare find_bestpair_selection records - THE PRIMARY OUTPUT (Stage 10)."""
    errors = []
    
    if not legacy_records or not modern_records:
        if not legacy_records and not modern_records:
            return True, []
        errors.append(f"One set is empty: legacy={len(legacy_records)}, modern={len(modern_records)}")
        return False, errors
    
    # Get the selection record (usually single record per PDB)
    legacy_rec = legacy_records[0] if legacy_records else {}
    modern_rec = modern_records[0] if modern_records else {}
    
    # Compare num_bp
    leg_count = legacy_rec.get('num_bp', 0)
    mod_count = modern_rec.get('num_bp', 0)
    
    if leg_count != mod_count:
        errors.append(f"num_bp mismatch: legacy={leg_count}, modern={mod_count}")
    
    # Compare pairs
    leg_pairs = set(tuple(normalize_pair(p[0], p[1])) for p in legacy_rec.get('pairs', []))
    mod_pairs = set(tuple(normalize_pair(p[0], p[1])) for p in modern_rec.get('pairs', []))
    
    missing_in_modern = leg_pairs - mod_pairs
    extra_in_modern = mod_pairs - leg_pairs
    
    if missing_in_modern:
        errors.append(f"Missing in modern: {sorted(missing_in_modern)}")
    if extra_in_modern:
        errors.append(f"Extra in modern: {sorted(extra_in_modern)}")
    
    return len(errors) == 0, errors


def compare_base_pair(legacy_records: List[Dict], modern_records: List[Dict],
                      tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare base pair records (Stage 9)."""
    errors = []
    
    # Build lookup by normalized pair
    legacy_by_pair = {}
    for rec in legacy_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        legacy_by_pair[pair] = rec
    
    modern_by_pair = {}
    for rec in modern_records:
        pair = normalize_pair(rec.get('base_i', 0), rec.get('base_j', 0))
        modern_by_pair[pair] = rec
    
    common_pairs = set(legacy_by_pair.keys()) & set(modern_by_pair.keys())
    
    for pair in common_pairs:
        legacy_rec = legacy_by_pair[pair]
        modern_rec = modern_by_pair[pair]
        
        # Compare bp_type
        if legacy_rec.get('bp_type') != modern_rec.get('bp_type'):
            errors.append(f"Pair {pair} bp_type: legacy={legacy_rec.get('bp_type')}, modern={modern_rec.get('bp_type')}")
        
        # Compare rotation matrices
        for matrix_name in ['orien_i', 'orien_j']:
            leg_mat = legacy_rec.get(matrix_name, [])
            mod_mat = modern_rec.get(matrix_name, [])
            if leg_mat and mod_mat:
                match, msg = compare_matrix(leg_mat, mod_mat, tolerance=1e-4)
                if not match:
                    errors.append(f"Pair {pair} {matrix_name}: {msg}")
        
        # Compare origins
        for org_name in ['org_i', 'org_j']:
            leg_org = legacy_rec.get(org_name, [])
            mod_org = modern_rec.get(org_name, [])
            if leg_org and mod_org:
                match, msg = compare_values(leg_org, mod_org, tolerance)
                if not match:
                    errors.append(f"Pair {pair} {org_name}: {msg}")
    
    return len(errors) == 0, errors


def compare_step_params(legacy_records: List[Dict], modern_records: List[Dict],
                        tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare step parameter records (Stage 11)."""
    errors = []
    
    # Build lookup by (bp_idx1, bp_idx2)
    legacy_by_key = {}
    for rec in legacy_records:
        key = (rec.get('bp_idx1', 0), rec.get('bp_idx2', 0))
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        key = (rec.get('bp_idx1', 0), rec.get('bp_idx2', 0))
        modern_by_key[key] = rec
    
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())
    
    if not common_keys and (legacy_by_key or modern_by_key):
        errors.append(f"No common steps (legacy={len(legacy_by_key)}, modern={len(modern_by_key)})")
        return False, errors
    
    param_fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
    
    for key in common_keys:
        legacy_rec = legacy_by_key[key]
        modern_rec = modern_by_key[key]
        
        for field in param_fields:
            leg_val = legacy_rec.get(field)
            mod_val = modern_rec.get(field)
            
            if leg_val is not None and mod_val is not None:
                match, msg = compare_values(leg_val, mod_val, tolerance)
                if not match:
                    errors.append(f"Step {key} {field}: {msg}")
    
    return len(errors) == 0, errors


def compare_helical_params(legacy_records: List[Dict], modern_records: List[Dict],
                           tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare helical parameter records (Stage 12)."""
    errors = []
    
    # Build lookup by (bp_idx1, bp_idx2)
    legacy_by_key = {}
    for rec in legacy_records:
        key = (rec.get('bp_idx1', 0), rec.get('bp_idx2', 0))
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        key = (rec.get('bp_idx1', 0), rec.get('bp_idx2', 0))
        modern_by_key[key] = rec
    
    common_keys = set(legacy_by_key.keys()) & set(modern_by_key.keys())
    
    param_fields = ['x_displacement', 'y_displacement', 'rise', 'inclination', 'tip', 'twist']
    
    for key in common_keys:
        legacy_rec = legacy_by_key[key]
        modern_rec = modern_by_key[key]
        
        for field in param_fields:
            leg_val = legacy_rec.get(field)
            mod_val = modern_rec.get(field)
            
            if leg_val is not None and mod_val is not None:
                match, msg = compare_values(leg_val, mod_val, tolerance)
                if not match:
                    errors.append(f"Helical {key} {field}: {msg}")
    
    return len(errors) == 0, errors


# Map stage number to comparison function
COMPARISON_FUNCTIONS = {
    1: compare_pdb_atoms,
    2: compare_residue_indices,
    3: compare_base_frame_calc,
    4: compare_ls_fitting,
    5: compare_frame_calc,
    6: compare_pair_validation,
    7: compare_distance_checks,
    8: compare_hbond_list,
    9: compare_base_pair,
    10: compare_find_bestpair_selection,
    11: compare_step_params,
    12: compare_helical_params,
}


# ============================================================================
# Test Infrastructure
# ============================================================================

def load_json_records(json_path: Path) -> List[Dict]:
    """Load JSON records from file."""
    if not json_path.exists():
        return []
    
    with open(json_path) as f:
        data = json.load(f)
    
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return [data]
    return []


def test_single_pdb_stage(pdb_id: str, stage_num: int, project_root: Path) -> Dict:
    """Test a single PDB for a specific stage."""
    stage_config = STAGES.get(stage_num)
    if not stage_config:
        return {"pdb_id": pdb_id, "stage": stage_num, "passed": False, "errors": ["Unknown stage"]}
    
    result = {
        "pdb_id": pdb_id,
        "stage": stage_num,
        "stage_name": stage_config.name,
        "passed": False,
        "errors": [],
        "details": {}
    }
    
    legacy_base = project_root / "data" / "json_legacy"
    modern_base = project_root / "data" / "json"
    
    json_type = stage_config.json_type
    legacy_file = legacy_base / json_type / f"{pdb_id}.json"
    modern_file = modern_base / json_type / f"{pdb_id}.json"
    
    # Check if legacy exists
    if not legacy_file.exists():
        result["errors"].append(f"No legacy JSON for {json_type}")
        result["details"]["legacy_exists"] = False
        return result
    
    # Check if modern exists
    if not modern_file.exists():
        result["errors"].append(f"No modern JSON for {json_type}")
        result["details"]["legacy_exists"] = True
        result["details"]["modern_exists"] = False
        return result
    
    # Load records
    legacy_records = load_json_records(legacy_file)
    modern_records = load_json_records(modern_file)
    
    result["details"] = {
        "legacy_exists": True,
        "modern_exists": True,
        "legacy_count": len(legacy_records),
        "modern_count": len(modern_records)
    }
    
    # Get comparison function
    compare_fn = COMPARISON_FUNCTIONS.get(stage_num)
    if compare_fn:
        passed, errors = compare_fn(legacy_records, modern_records, stage_config.tolerance)
        result["details"]["passed"] = passed
        result["errors"].extend(errors)
    else:
        # Generic comparison: just check counts match
        if len(legacy_records) != len(modern_records):
            result["errors"].append(f"{json_type} count mismatch: {len(legacy_records)} vs {len(modern_records)}")
        else:
            result["details"]["passed"] = True
    
    result["passed"] = len(result["errors"]) == 0
    return result


def run_stage_validation(stage_num: int, pdb_ids: List[str], project_root: Path,
                         stop_on_failure: bool = True, verbose: bool = False) -> Dict:
    """Run validation for a stage across all PDBs."""
    stage_config = STAGES.get(stage_num)
    if not stage_config:
        raise ValueError(f"Unknown stage: {stage_num}")
    
    results = {
        "stage_num": stage_num,
        "stage_id": stage_config.stage_id,
        "stage_name": stage_config.name,
        "test_date": datetime.now().isoformat(),
        "total_pdbs": len(pdb_ids),
        "passed": 0,
        "failed": 0,
        "skipped": 0,
        "failed_pdbs": [],
        "elapsed_seconds": 0.0
    }
    
    start_time = datetime.now()
    
    for i, pdb_id in enumerate(pdb_ids):
        result = test_single_pdb_stage(pdb_id, stage_num, project_root)
        
        if result["passed"]:
            results["passed"] += 1
            if verbose:
                print(f"  ✅ [{i+1}/{len(pdb_ids)}] {pdb_id}")
        elif not result["details"]:
            results["skipped"] += 1
            if verbose:
                print(f"  ⏭️ [{i+1}/{len(pdb_ids)}] {pdb_id} (skipped)")
        else:
            results["failed"] += 1
            results["failed_pdbs"].append({
                "pdb_id": pdb_id,
                "errors": result["errors"][:5],  # Limit to first 5 errors
                "details": result["details"]
            })
            
            if verbose:
                print(f"  ❌ [{i+1}/{len(pdb_ids)}] {pdb_id}")
                for error in result["errors"][:3]:
                    print(f"       {error}")
            
            if stop_on_failure:
                print(f"\n❌ FAILED at PDB {pdb_id}")
                for error in result["errors"][:5]:
                    print(f"  - {error}")
                break
        
        # Progress update (if not verbose)
        if not verbose and (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{len(pdb_ids)} ({results['passed']} passed, {results['failed']} failed)")
    
    results["elapsed_seconds"] = (datetime.now() - start_time).total_seconds()
    results["pass_rate"] = (results["passed"] / results["total_pdbs"] * 100) if results["total_pdbs"] > 0 else 0
    
    return results


# ============================================================================
# Pytest Tests
# ============================================================================

@pytest.fixture(scope="session")
def valid_pdbs():
    """Load valid fast PDBs."""
    return load_valid_pdbs_fast(project_root)


def _run_stage_test(stage_num: int, valid_pdbs: List[str], request) -> Dict:
    """Helper to run a stage test with pytest options."""
    max_pdbs = request.config.getoption("--max-pdbs", default=None)
    verbose = request.config.getoption("--verbose", default=False)
    pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
    
    results = run_stage_validation(stage_num, pdb_list, project_root, 
                                   stop_on_failure=True, verbose=verbose)
    
    # Save results
    stage_config = STAGES[stage_num]
    results_file = project_root / "data" / "validation_results" / f"stage{stage_num}_{stage_config.stage_id}_pytest.json"
    results_file.parent.mkdir(parents=True, exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    return results


class TestStageValidation:
    """Test class for stage validation."""
    
    def test_stage1_pdb_atoms(self, valid_pdbs, request):
        """Stage 1: Validate atom parsing."""
        results = _run_stage_test(1, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 1 (pdb_atoms) failed: {results['failed']} PDBs"
    
    def test_stage2_residue_indices(self, valid_pdbs, request):
        """Stage 2: Validate residue indices."""
        results = _run_stage_test(2, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 2 (residue_indices) failed: {results['failed']} PDBs"
    
    def test_stage3_base_frame_calc(self, valid_pdbs, request):
        """Stage 3: Validate base frame calculation."""
        results = _run_stage_test(3, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 3 (base_frame_calc) failed: {results['failed']} PDBs"
    
    def test_stage4_ls_fitting(self, valid_pdbs, request):
        """Stage 4: Validate least squares fitting."""
        results = _run_stage_test(4, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 4 (ls_fitting) failed: {results['failed']} PDBs"
    
    def test_stage5_frame_calc(self, valid_pdbs, request):
        """Stage 5: Validate frame calculation."""
        results = _run_stage_test(5, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 5 (frame_calc) failed: {results['failed']} PDBs"
    
    def test_stage6_pair_validation(self, valid_pdbs, request):
        """Stage 6: Validate pair validation."""
        results = _run_stage_test(6, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 6 (pair_validation) failed: {results['failed']} PDBs"
    
    def test_stage7_distance_checks(self, valid_pdbs, request):
        """Stage 7: Validate distance checks."""
        results = _run_stage_test(7, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 7 (distance_checks) failed: {results['failed']} PDBs"
    
    def test_stage8_hbond_list(self, valid_pdbs, request):
        """Stage 8: Validate H-bond list."""
        results = _run_stage_test(8, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 8 (hbond_list) failed: {results['failed']} PDBs"
    
    def test_stage9_base_pair(self, valid_pdbs, request):
        """Stage 9: Validate base pair records."""
        results = _run_stage_test(9, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 9 (base_pair) failed: {results['failed']} PDBs"
    
    def test_stage10_find_bestpair_selection(self, valid_pdbs, request):
        """Stage 10: Validate find_bestpair_selection - THE PRIMARY OUTPUT."""
        results = _run_stage_test(10, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 10 (PRIMARY OUTPUT) failed: {results['failed']} PDBs"
    
    def test_stage11_bpstep_params(self, valid_pdbs, request):
        """Stage 11: Validate step parameters."""
        results = _run_stage_test(11, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 11 (bpstep_params) failed: {results['failed']} PDBs"
    
    def test_stage12_helical_params(self, valid_pdbs, request):
        """Stage 12: Validate helical parameters."""
        results = _run_stage_test(12, valid_pdbs, request)
        assert results["failed"] == 0, f"Stage 12 (helical_params) failed: {results['failed']} PDBs"


def pytest_addoption(parser):
    """Add custom pytest options."""
    parser.addoption(
        "--max-pdbs",
        action="store",
        default=None,
        help="Maximum number of PDBs to test"
    )
    parser.addoption(
        "--pdb",
        action="store",
        default=None,
        help="Test a specific PDB ID"
    )


# ============================================================================
# CLI Interface
# ============================================================================

def main():
    """CLI interface for running validation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Stage validation CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python test_stage_validation.py 1 --max-pdbs 100        # Stage 1
    python test_stage_validation.py frames --max-pdbs 100   # Stages 3,4,5
    python test_stage_validation.py all --pdb 1EHZ          # All stages, single PDB
    python test_stage_validation.py 3 -v -s                 # Stage 3, verbose, stop on first

Stage Groups:
    atoms   = 1         (atom parsing)
    residue = 2         (residue indices)
    frames  = 3,4,5     (frame calculations)
    pairs   = 6,7,9,10  (pair validation/selection)
    hbonds  = 8         (hydrogen bonds)
    steps   = 11,12     (step/helical parameters)
    all     = 1-12      (all stages)
"""
    )
    parser.add_argument("stages", nargs="*", default=["all"],
                        help="Stage(s) to validate: number (1-12), name (pdb_atoms), or group (frames)")
    parser.add_argument("--max-pdbs", type=int, help="Maximum number of PDBs to test")
    parser.add_argument("--pdb", help="Test a specific PDB ID")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("-s", "--stop-on-failure", action="store_true", default=True,
                        help="Stop on first failure (default)")
    parser.add_argument("--no-stop-on-failure", dest="stop_on_failure", action="store_false",
                        help="Don't stop on first failure")
    
    args = parser.parse_args()
    
    # Resolve stages
    stage_nums = resolve_stages(args.stages)
    if not stage_nums:
        print(f"Error: No valid stages found in: {args.stages}")
        sys.exit(1)
    
    # Load PDB list
    if args.pdb:
        pdb_ids = [args.pdb.upper()]
    else:
        try:
            pdb_ids = load_valid_pdbs_fast(project_root)
        except (FileNotFoundError, ValueError) as e:
            print(f"Error: {e}")
            sys.exit(1)
        if args.max_pdbs:
            pdb_ids = pdb_ids[:args.max_pdbs]
    
    print(f"Running validation for stages: {stage_nums}")
    print(f"PDBs: {len(pdb_ids)}")
    print()
    
    all_results = []
    any_failures = False
    
    for stage_num in stage_nums:
        stage_config = STAGES[stage_num]
        print(f"{'='*60}")
        print(f"STAGE {stage_num}: {stage_config.name}")
        print(f"{'='*60}")
        
        results = run_stage_validation(
            stage_num, pdb_ids, project_root,
            stop_on_failure=args.stop_on_failure,
            verbose=args.verbose
        )
        all_results.append(results)
        
        # Print summary
        print(f"\nResults:")
        print(f"  Total: {results['total_pdbs']}")
        print(f"  Passed: {results['passed']} ({results['pass_rate']:.1f}%)")
        print(f"  Failed: {results['failed']}")
        print(f"  Skipped: {results['skipped']}")
        print(f"  Time: {results['elapsed_seconds']:.2f}s")
        
        if results['failed_pdbs']:
            print(f"\nFailed PDBs (first 5):")
            for fail in results['failed_pdbs'][:5]:
                print(f"  {fail['pdb_id']}: {fail['errors'][0] if fail['errors'] else 'Unknown'}")
        
        # Save results
        results_file = project_root / "data" / "validation_results" / f"stage{stage_num}_{stage_config.stage_id}_cli.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        if results['failed'] > 0:
            any_failures = True
            if args.stop_on_failure:
                print(f"\nStopping due to failures in stage {stage_num}")
                break
        
        print()
    
    # Overall summary
    print(f"\n{'='*60}")
    print("OVERALL SUMMARY")
    print(f"{'='*60}")
    for r in all_results:
        status = "✅" if r['failed'] == 0 else "❌"
        print(f"  {status} Stage {r['stage_num']} ({r['stage_id']}): {r['passed']}/{r['total_pdbs']} passed")
    
    return 1 if any_failures else 0


if __name__ == "__main__":
    sys.exit(main())
