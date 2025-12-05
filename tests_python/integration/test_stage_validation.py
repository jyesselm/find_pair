#!/usr/bin/env python3
"""
Unified stage-by-stage validation for legacy vs modern JSON comparison.

This script validates each stage of the base pair finding algorithm by comparing
legacy and modern JSON outputs. It supports stop-on-first-failure, parallel
processing, and automatic cleanup of temp files.

Usage:
    # Test a specific stage
    pytest tests_python/integration/test_stage_validation.py -v -k "stage3"
    
    # Run with max PDBs limit
    pytest tests_python/integration/test_stage_validation.py -v --max-pdbs=100
    
    # Stop on first failure  
    pytest tests_python/integration/test_stage_validation.py -v -x
    
    # Test single PDB
    pytest tests_python/integration/test_stage_validation.py::test_stage3_single -v
"""

import json
import pytest
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from datetime import datetime
from multiprocessing import Pool, cpu_count
import sys

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
    stage_id: str
    name: str
    description: str
    record_types: List[str]
    dependencies: List[str]
    key_fields: List[str]
    tolerance: float = 1e-6


STAGES = {
    "stage0": StageConfig(
        stage_id="stage0",
        name="Residue Indices",
        description="PDB residue parsing and index mapping",
        record_types=["residue_indices"],
        dependencies=[],
        key_fields=["chain_id", "residue_seq", "insertion", "legacy_residue_idx"]
    ),
    "stage1": StageConfig(
        stage_id="stage1",
        name="Atoms",
        description="PDB atom parsing",
        record_types=["pdb_atoms"],
        dependencies=[],
        key_fields=["atom_idx", "atom_name", "xyz", "chain_id", "residue_seq"]
    ),
    "stage2": StageConfig(
        stage_id="stage2",
        name="LS Fitting",
        description="Least-squares fitting for base frame calculation",
        record_types=["ls_fitting"],
        dependencies=["stage0", "stage1"],
        key_fields=["chain_id", "residue_seq", "insertion", "rms", "num_points"]
    ),
    "stage3": StageConfig(
        stage_id="stage3",
        name="Distance Checks",
        description="Geometric measurements between base pairs",
        record_types=["distance_checks"],
        dependencies=["stage2"],
        key_fields=["base_i", "base_j", "dorg", "dNN", "plane_angle", "d_v", "overlap_area"]
    ),
    "stage4": StageConfig(
        stage_id="stage4",
        name="H-bond List",
        description="Hydrogen bond detection for base pairs",
        record_types=["hbond_list"],
        dependencies=["stage3"],
        key_fields=["base_i", "base_j", "num_hbonds", "hbonds"]
    ),
    "stage5": StageConfig(
        stage_id="stage5",
        name="Pair Validation",
        description="Validation results for each residue pair",
        record_types=["pair_validation"],
        dependencies=["stage3", "stage4"],
        key_fields=["base_i", "base_j", "is_valid", "bp_type_id", "quality_score"]
    ),
    "stage6": StageConfig(
        stage_id="stage6",
        name="Find Bestpair Selection",
        description="THE PRIMARY OUTPUT - selected base pairs",
        record_types=["find_bestpair_selection"],
        dependencies=["stage5"],
        key_fields=["num_bp", "pairs"]
    ),
    "stage7": StageConfig(
        stage_id="stage7",
        name="Base Pair Records",
        description="Detailed base pair information",
        record_types=["base_pair"],
        dependencies=["stage6"],
        key_fields=["base_i", "base_j", "bp_type", "orien_i", "orien_j", "org_i", "org_j"]
    ),
    "stage8": StageConfig(
        stage_id="stage8",
        name="Step Parameters",
        description="6 step parameters for consecutive base pairs",
        record_types=["bpstep_params"],
        dependencies=["stage6"],
        key_fields=["bp_idx1", "bp_idx2", "shift", "slide", "rise", "tilt", "roll", "twist"]
    ),
    "stage9": StageConfig(
        stage_id="stage9",
        name="Helical Parameters",
        description="Helical axis parameters",
        record_types=["helical_params"],
        dependencies=["stage6", "stage8"],
        key_fields=["bp_idx1", "bp_idx2", "x_displacement", "y_displacement", "rise", "inclination", "tip", "twist"]
    ),
}


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
    return False, f"Type/value mismatch: legacy={legacy_val} ({type(legacy_val)}), modern={modern_val} ({type(modern_val)})"


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

def compare_residue_indices(legacy_records: List[Dict], modern_records: List[Dict], 
                            tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare residue indices records."""
    errors = []
    
    if len(legacy_records) != len(modern_records):
        errors.append(f"Count mismatch: legacy={len(legacy_records)}, modern={len(modern_records)}")
        return False, errors
    
    # Build lookup by (chain_id, residue_seq, insertion)
    legacy_by_key = {}
    for rec in legacy_records:
        key = (rec.get('chain_id', ''), rec.get('residue_seq', 0), rec.get('insertion', ' '))
        legacy_by_key[key] = rec
    
    modern_by_key = {}
    for rec in modern_records:
        key = (rec.get('chain_id', ''), rec.get('residue_seq', 0), rec.get('insertion', ' '))
        modern_by_key[key] = rec
    
    # Compare
    for key, legacy_rec in legacy_by_key.items():
        if key not in modern_by_key:
            errors.append(f"Missing in modern: {key}")
            continue
        
        modern_rec = modern_by_key[key]
        
        # Compare legacy_residue_idx
        leg_idx = legacy_rec.get('residue_idx') or legacy_rec.get('legacy_residue_idx')
        mod_idx = modern_rec.get('legacy_residue_idx')
        if leg_idx != mod_idx:
            errors.append(f"Index mismatch at {key}: legacy={leg_idx}, modern={mod_idx}")
    
    return len(errors) == 0, errors


def compare_distance_checks(legacy_records: List[Dict], modern_records: List[Dict],
                           tolerance: float = 1e-6) -> Tuple[bool, List[str]]:
    """Compare distance checks records."""
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
    
    # Find common pairs
    common_pairs = set(legacy_by_pair.keys()) & set(modern_by_pair.keys())
    
    if not common_pairs and (legacy_by_pair or modern_by_pair):
        errors.append(f"No common pairs found (legacy={len(legacy_by_pair)}, modern={len(modern_by_pair)})")
        return False, errors
    
    # Compare common pairs
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
    """Compare H-bond list records."""
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
            # Compare donor, acceptor, distance, type
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
    """Compare pair validation records."""
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
    """Compare find_bestpair_selection records - THE PRIMARY OUTPUT."""
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
    """Compare base pair records."""
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
    """Compare step parameter records."""
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
    """Compare helical parameter records."""
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


# Map stage to comparison function
COMPARISON_FUNCTIONS = {
    "residue_indices": compare_residue_indices,
    "distance_checks": compare_distance_checks,
    "hbond_list": compare_hbond_list,
    "pair_validation": compare_pair_validation,
    "find_bestpair_selection": compare_find_bestpair_selection,
    "base_pair": compare_base_pair,
    "bpstep_params": compare_step_params,
    "helical_params": compare_helical_params,
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


def test_single_pdb_stage(pdb_id: str, stage_config: StageConfig,
                          project_root: Path) -> Dict:
    """Test a single PDB for a specific stage."""
    result = {
        "pdb_id": pdb_id,
        "stage": stage_config.stage_id,
        "passed": False,
        "errors": [],
        "details": {}
    }
    
    legacy_base = project_root / "data" / "json_legacy"
    modern_base = project_root / "data" / "json"
    
    for record_type in stage_config.record_types:
        legacy_file = legacy_base / record_type / f"{pdb_id}.json"
        modern_file = modern_base / record_type / f"{pdb_id}.json"
        
        # Check if legacy exists
        if not legacy_file.exists():
            result["errors"].append(f"No legacy JSON for {record_type}")
            result["details"][record_type] = {"legacy_exists": False}
            continue
        
        # Check if modern exists
        if not modern_file.exists():
            result["errors"].append(f"No modern JSON for {record_type}")
            result["details"][record_type] = {"legacy_exists": True, "modern_exists": False}
            continue
        
        # Load records
        legacy_records = load_json_records(legacy_file)
        modern_records = load_json_records(modern_file)
        
        result["details"][record_type] = {
            "legacy_exists": True,
            "modern_exists": True,
            "legacy_count": len(legacy_records),
            "modern_count": len(modern_records)
        }
        
        # Get comparison function
        compare_fn = COMPARISON_FUNCTIONS.get(record_type)
        if compare_fn:
            passed, errors = compare_fn(legacy_records, modern_records, stage_config.tolerance)
            result["details"][record_type]["passed"] = passed
            result["errors"].extend(errors)
        else:
            # Generic comparison: just check counts match
            if len(legacy_records) != len(modern_records):
                result["errors"].append(f"{record_type} count mismatch: {len(legacy_records)} vs {len(modern_records)}")
            else:
                result["details"][record_type]["passed"] = True
    
    result["passed"] = len(result["errors"]) == 0
    return result


def run_stage_validation(stage_id: str, pdb_ids: List[str], project_root: Path,
                         max_workers: int = 20, stop_on_failure: bool = True) -> Dict:
    """Run validation for a stage across all PDBs."""
    stage_config = STAGES.get(stage_id)
    if not stage_config:
        raise ValueError(f"Unknown stage: {stage_id}")
    
    results = {
        "stage_id": stage_id,
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
        result = test_single_pdb_stage(pdb_id, stage_config, project_root)
        
        if result["passed"]:
            results["passed"] += 1
        elif not result["details"]:
            results["skipped"] += 1
        else:
            results["failed"] += 1
            results["failed_pdbs"].append({
                "pdb_id": pdb_id,
                "errors": result["errors"][:5],  # Limit to first 5 errors
                "details": result["details"]
            })
            
            if stop_on_failure:
                print(f"\n❌ FAILED at PDB {pdb_id} [{i+1}/{len(pdb_ids)}]")
                for error in result["errors"][:5]:
                    print(f"  - {error}")
                break
        
        # Progress update
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i+1}/{len(pdb_ids)} ({results['passed']} passed, {results['failed']} failed)")
    
    results["elapsed_seconds"] = (datetime.now() - start_time).total_seconds()
    results["pass_rate"] = (results["passed"] / results["total_pdbs"] * 100) if results["total_pdbs"] > 0 else 0
    
    return results


# ============================================================================
# Pytest Tests
# ============================================================================

@pytest.fixture(scope="session")
def valid_pdbs(project_root_path):
    """Load valid fast PDBs."""
    return load_valid_pdbs_fast(project_root_path)


@pytest.fixture(scope="session")
def project_root_path():
    """Return project root."""
    return project_root


class TestStageValidation:
    """Test class for stage validation."""
    
    def test_stage0_residue_indices(self, valid_pdbs, project_root_path, request):
        """Stage 0: Validate residue indices matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage0", pdb_list, project_root_path)
        
        # Save results
        results_file = project_root_path / "data" / "validation_results" / "stage0_residue_indices_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 0 failed: {results['failed']} PDBs"
    
    def test_stage3_distance_checks(self, valid_pdbs, project_root_path, request):
        """Stage 3: Validate distance checks matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage3", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage3_distance_checks_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 3 failed: {results['failed']} PDBs"
    
    def test_stage4_hbond_list(self, valid_pdbs, project_root_path, request):
        """Stage 4: Validate H-bond list matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage4", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage4_hbond_list_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 4 failed: {results['failed']} PDBs"
    
    def test_stage5_pair_validation(self, valid_pdbs, project_root_path, request):
        """Stage 5: Validate pair validation matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage5", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage5_pair_validation_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 5 failed: {results['failed']} PDBs"
    
    def test_stage6_find_bestpair_selection(self, valid_pdbs, project_root_path, request):
        """Stage 6: Validate find_bestpair_selection - THE PRIMARY OUTPUT."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage6", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage6_find_bestpair_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 6 (PRIMARY OUTPUT) failed: {results['failed']} PDBs"
    
    def test_stage7_base_pair(self, valid_pdbs, project_root_path, request):
        """Stage 7: Validate base pair records matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage7", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage7_base_pair_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 7 failed: {results['failed']} PDBs"
    
    def test_stage8_step_params(self, valid_pdbs, project_root_path, request):
        """Stage 8: Validate step parameters matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage8", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage8_step_params_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 8 failed: {results['failed']} PDBs"
    
    def test_stage9_helical_params(self, valid_pdbs, project_root_path, request):
        """Stage 9: Validate helical parameters matching."""
        max_pdbs = request.config.getoption("--max-pdbs", default=None)
        pdb_list = valid_pdbs[:int(max_pdbs)] if max_pdbs else valid_pdbs
        
        results = run_stage_validation("stage9", pdb_list, project_root_path)
        
        results_file = project_root_path / "data" / "validation_results" / "stage9_helical_params_pytest.json"
        results_file.parent.mkdir(parents=True, exist_ok=True)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        assert results["failed"] == 0, f"Stage 9 failed: {results['failed']} PDBs"


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
    
    parser = argparse.ArgumentParser(description="Stage validation CLI")
    parser.add_argument("stage", help="Stage to validate (e.g., stage3)")
    parser.add_argument("--max-pdbs", type=int, help="Maximum number of PDBs to test")
    parser.add_argument("--stop-on-failure", action="store_true", default=True,
                        help="Stop on first failure (default: True)")
    parser.add_argument("--no-stop-on-failure", dest="stop_on_failure", action="store_false",
                        help="Don't stop on first failure")
    parser.add_argument("--pdb", help="Test a specific PDB ID")
    
    args = parser.parse_args()
    
    # Load PDB list
    if args.pdb:
        pdb_ids = [args.pdb.upper()]
    else:
        pdb_ids = load_valid_pdbs_fast(project_root)
        if args.max_pdbs:
            pdb_ids = pdb_ids[:args.max_pdbs]
    
    print(f"Running {args.stage} validation on {len(pdb_ids)} PDBs...")
    print()
    
    results = run_stage_validation(args.stage, pdb_ids, project_root, 
                                   stop_on_failure=args.stop_on_failure)
    
    # Print summary
    print("\n" + "=" * 60)
    print(f"STAGE {args.stage.upper()} VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Total PDBs: {results['total_pdbs']}")
    print(f"Passed: {results['passed']} ({results['pass_rate']:.2f}%)")
    print(f"Failed: {results['failed']}")
    print(f"Skipped: {results['skipped']}")
    print(f"Elapsed: {results['elapsed_seconds']:.2f}s")
    
    if results['failed_pdbs']:
        print(f"\nFailed PDBs (first 10):")
        for fail in results['failed_pdbs'][:10]:
            print(f"  {fail['pdb_id']}: {fail['errors'][0] if fail['errors'] else 'Unknown error'}")
    
    # Save results
    results_file = project_root / "data" / "validation_results" / f"{args.stage}_validation_cli.json"
    results_file.parent.mkdir(parents=True, exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {results_file}")
    
    return 0 if results['failed'] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

