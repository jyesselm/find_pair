"""
Validation runner for executing comparisons.

This module handles:
- Loading JSON files
- Running comparators
- Collecting and reporting results
- Parallel processing for large batches
"""

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .config import STAGES, StageConfig
from .comparators import COMPARATORS


@dataclass
class ValidationResult:
    """Result of validating a single PDB for a stage."""
    
    pdb_id: str
    stage_num: int
    passed: bool
    errors: List[str] = field(default_factory=list)
    skipped: bool = False
    skip_reason: str = ""


@dataclass
class StageResult:
    """Aggregated result for a stage across all PDBs."""
    
    stage_num: int
    stage_name: str
    total: int = 0
    passed: int = 0
    failed: int = 0
    skipped: int = 0
    failed_pdbs: List[str] = field(default_factory=list)
    first_errors: Dict[str, List[str]] = field(default_factory=dict)


def validate_pdb(
    pdb_id: str,
    stage_num: int,
    json_dir: Path,
    verbose: bool = False
) -> ValidationResult:
    """
    Validate a single PDB for a specific stage.
    
    Args:
        pdb_id: PDB identifier (e.g., "100D")
        stage_num: Stage number (1-12)
        json_dir: Path to JSON directory (contains type subdirs)
        verbose: Print verbose output
        
    Returns:
        ValidationResult with pass/fail status and any errors
    """
    config = STAGES.get(stage_num)
    if config is None:
        return ValidationResult(pdb_id, stage_num, False, [f"Unknown stage {stage_num}"])
    
    legacy_data = _load_json(json_dir, pdb_id, config.json_type, legacy=True)
    modern_data = _load_json(json_dir, pdb_id, config.json_type, legacy=False)
    
    if legacy_data is None:
        return ValidationResult(
            pdb_id, stage_num, False, [],
            skipped=True, skip_reason="Missing legacy JSON"
        )
    if modern_data is None:
        return ValidationResult(
            pdb_id, stage_num, False, [],
            skipped=True, skip_reason="Missing modern JSON"
        )
    
    comparator = COMPARATORS.get(stage_num)
    if comparator is None:
        return ValidationResult(pdb_id, stage_num, False, [f"No comparator for stage {stage_num}"])
    
    passed, errors = comparator(legacy_data, modern_data, config.tolerance)
    
    if verbose and errors:
        print(f"  {pdb_id}: {len(errors)} error(s)")
        for err in errors[:3]:
            print(f"    - {err}")
    
    return ValidationResult(pdb_id, stage_num, passed, errors)


def _validate_pdb_worker(args: Tuple[str, int, str, bool]) -> ValidationResult:
    """Worker function for parallel validation (must be at module level for pickling)."""
    pdb_id, stage_num, json_dir_str, verbose = args
    return validate_pdb(pdb_id, stage_num, Path(json_dir_str), verbose)


def validate_stage(
    stage_num: int,
    pdb_ids: List[str],
    json_dir: Path,
    verbose: bool = False,
    stop_on_failure: bool = False,
    num_workers: int = 1
) -> StageResult:
    """
    Validate all PDBs for a specific stage.
    
    Args:
        stage_num: Stage number (1-12)
        pdb_ids: List of PDB identifiers
        json_dir: Path to JSON directory
        verbose: Print verbose output
        stop_on_failure: Stop on first failure
        num_workers: Number of parallel workers (1 = sequential)
        
    Returns:
        StageResult with aggregated statistics
    """
    config = STAGES.get(stage_num)
    if config is None:
        return StageResult(stage_num, f"Unknown ({stage_num})")
    
    result = StageResult(stage_num, config.name)
    
    if num_workers <= 1:
        # Sequential processing
        for pdb_id in pdb_ids:
            vr = validate_pdb(pdb_id, stage_num, json_dir, verbose)
            result.total += 1
            
            if vr.skipped:
                result.skipped += 1
            elif vr.passed:
                result.passed += 1
            else:
                result.failed += 1
                result.failed_pdbs.append(pdb_id)
                if vr.errors:
                    result.first_errors[pdb_id] = vr.errors[:3]
                
                if stop_on_failure:
                    break
    else:
        # Parallel processing
        json_dir_str = str(json_dir)
        args_list = [(pdb_id, stage_num, json_dir_str, verbose) for pdb_id in pdb_ids]
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(_validate_pdb_worker, args): args[0] 
                       for args in args_list}
            
            for future in as_completed(futures):
                pdb_id = futures[future]
                try:
                    vr = future.result()
                except Exception as e:
                    vr = ValidationResult(pdb_id, stage_num, False, [str(e)])
                
                result.total += 1
                
                if vr.skipped:
                    result.skipped += 1
                elif vr.passed:
                    result.passed += 1
                else:
                    result.failed += 1
                    result.failed_pdbs.append(pdb_id)
                    if vr.errors:
                        result.first_errors[pdb_id] = vr.errors[:3]
    
    return result


def _load_json(
    json_dir: Path,
    pdb_id: str,
    json_type: str,
    legacy: bool = True
) -> Optional[List[Dict[str, Any]]]:
    """
    Load JSON data for a PDB and record type.
    
    Directory structure: json_dir/<json_type>/<pdb_id>.json
    For legacy comparison, also checks json_dir/../json_legacy/<json_type>/<pdb_id>.json
    
    Args:
        json_dir: Base JSON directory
        pdb_id: PDB identifier
        json_type: JSON record type (e.g., "pdb_atoms")
        legacy: If True, load legacy data; if False, load modern data
        
    Returns:
        List of records or None if file not found/invalid
    """
    if legacy:
        # Try legacy directory first, then fall back to main
        legacy_dir = json_dir.parent / "json_legacy" / json_type
        if legacy_dir.exists():
            json_path = legacy_dir / f"{pdb_id}.json"
            if json_path.exists():
                return _read_json_file(json_path)
    
    # Modern data or legacy fallback
    json_path = json_dir / json_type / f"{pdb_id}.json"
    return _read_json_file(json_path)


def _read_json_file(path: Path) -> Optional[List[Dict[str, Any]]]:
    """Read and parse a JSON file."""
    if not path.exists():
        return None
    
    try:
        with open(path) as f:
            data = json.load(f)
        return _extract_records(data)
    except (json.JSONDecodeError, IOError):
        return None


def _extract_records(data: Any) -> Optional[List[Dict[str, Any]]]:
    """
    Extract list of records from various JSON formats.
    
    Handles:
    - List of records: [{...}, {...}]
    - Wrapped in dict: {atoms: [{...}]} or {records: [{...}]} or {seidx: [{...}]}
    - List with wrapped dict: [{atoms: [{...}]}] (legacy format)
    - Single record: {...}
    """
    # Direct list of records
    if isinstance(data, list):
        if not data:
            return []
        # Check if it's a list of wrapped dicts (legacy format)
        first = data[0]
        if isinstance(first, dict):
            if "atoms" in first:
                return first["atoms"]
            if "records" in first:
                return first["records"]
            if "seidx" in first:
                return first["seidx"]
        # Regular list of records
        return data
    
    # Dict with atoms/records/seidx key
    if isinstance(data, dict):
        if "atoms" in data:
            return data["atoms"]
        if "records" in data:
            return data["records"]
        if "seidx" in data:
            return data["seidx"]
        # Single record
        return [data]
    
    return None


def print_stage_result(result: StageResult, show_errors: bool = True) -> None:
    """
    Print a formatted stage result.
    
    Args:
        result: StageResult to print
        show_errors: Show first few errors for failed PDBs
    """
    print(f"\n{'='*60}")
    print(f"STAGE {result.stage_num}: {result.stage_name}")
    print(f"{'='*60}")
    
    pct = (result.passed / result.total * 100) if result.total > 0 else 0
    print(f"\nResults:")
    print(f"  Total: {result.total}")
    print(f"  Passed: {result.passed} ({pct:.1f}%)")
    print(f"  Failed: {result.failed}")
    print(f"  Skipped: {result.skipped}")
    
    if result.failed_pdbs and show_errors:
        print(f"\nFailed PDBs (first 5):")
        for pdb_id in result.failed_pdbs[:5]:
            errs = result.first_errors.get(pdb_id, [])
            first_err = errs[0] if errs else "Unknown error"
            print(f"  {pdb_id}: {first_err}")
