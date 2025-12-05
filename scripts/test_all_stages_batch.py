#!/usr/bin/env python3
"""
Unified batch validation for all stages of the base pair finding algorithm.

This script validates each stage by comparing legacy and modern JSON outputs.
Features:
- Stop-on-first-failure
- Parallel processing
- Cleanup temp files after each batch (keep only mismatches)
- Uses x3dna_json_compare module for comparisons

Usage:
    # Run Stage 3 (distance_checks) 
    python scripts/test_all_stages_batch.py stage3
    
    # Run all stages starting from stage 3
    python scripts/test_all_stages_batch.py stage3 --continue-stages
    
    # Run with max PDBs limit
    python scripts/test_all_stages_batch.py stage3 --max-pdbs=100
    
    # Single PDB test
    python scripts/test_all_stages_batch.py stage3 --pdb=1H4S
    
    # Continue from where we left off (if stopped early)
    python scripts/test_all_stages_batch.py stage3 --resume
"""

import argparse
import json
import subprocess
import shutil
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field

# Import comparison modules
from x3dna_json_compare.distance_comparison import compare_distance_checks
from x3dna_json_compare.hbond_comparison import compare_hbond_lists
from x3dna_json_compare.pair_validation_comparison import compare_pair_validation
from x3dna_json_compare.selection_comparison import compare_pair_selection
from x3dna_json_compare.base_pair_comparison import compare_base_pairs
from x3dna_json_compare.step_comparison import compare_step_parameters


# Stage definitions
STAGES = {
    "stage3": {
        "name": "Distance Checks",
        "legacy_dir": "distance_checks",
        "modern_dir": "distance_checks",
        "stage_flag": "distances",
        "compare_func": compare_distance_checks,
        "dependencies": ["stage2"],  # Requires frames
    },
    "stage4": {
        "name": "H-Bonds",
        "legacy_dir": "hbond_list",
        "modern_dir": "hbond_list",
        "stage_flag": "hbonds",
        "compare_func": compare_hbond_lists,
        "dependencies": ["stage3"],
    },
    "stage5": {
        "name": "Pair Validation",
        "legacy_dir": "pair_validation",
        "modern_dir": "pair_validation",
        "stage_flag": "validation",
        "compare_func": compare_pair_validation,
        "dependencies": ["stage4"],
    },
    "stage6": {
        "name": "Find Bestpair Selection (PRIMARY)",
        "legacy_dir": "find_bestpair_selection",
        "modern_dir": "find_bestpair_selection",
        "stage_flag": "selection",
        "compare_func": compare_pair_selection,
        "dependencies": ["stage5"],
    },
    "stage7": {
        "name": "Base Pair Records",
        "legacy_dir": "base_pair",
        "modern_dir": "base_pair",
        "stage_flag": "all",  # Generated with all
        "compare_func": compare_base_pairs,
        "dependencies": ["stage6"],
    },
    "stage8": {
        "name": "Step Parameters",
        "legacy_dir": "bpstep_params",
        "modern_dir": "step_params",
        "stage_flag": "steps",
        "compare_func": compare_step_parameters,
        "dependencies": ["stage7"],
    },
}


@dataclass
class StageResult:
    """Result of validating one stage for one PDB."""
    pdb_id: str
    stage_id: str
    passed: bool
    error: str = ""
    legacy_count: int = 0
    modern_count: int = 0
    matched: int = 0
    mismatched: int = 0
    missing_in_modern: int = 0
    extra_in_modern: int = 0
    skipped: bool = False  # True if skipped due to corrupted legacy JSON


def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    modern_exe = project_root / "build" / "generate_modern_json"
    
    if not legacy_exe.exists():
        legacy_exe = None
    if not modern_exe.exists():
        modern_exe = None
    
    return legacy_exe, modern_exe


def generate_modern_json(
    pdb_id: str,
    pdb_file: Path,
    output_dir: Path,
    modern_exe: Path,
    stage_flag: str
) -> Tuple[bool, str]:
    """Generate modern JSON for a specific stage."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [str(modern_exe), str(pdb_file), str(output_dir), f"--stage={stage_flag}"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            return False, f"Generation failed: {result.stderr[:200]}"
        
        return True, "OK"
        
    except subprocess.TimeoutExpired:
        return False, "Generation timeout"
    except Exception as e:
        return False, f"Error: {str(e)}"


def load_legacy_json(legacy_dir: Path, pdb_id: str) -> Tuple[Optional[List], str]:
    """Load legacy JSON records for a PDB."""
    json_file = legacy_dir / f"{pdb_id}.json"
    
    if not json_file.exists():
        return None, f"Legacy JSON not found: {json_file}"
    
    try:
        with open(json_file) as f:
            data = json.load(f)
        
        # Handle different JSON formats
        if isinstance(data, list):
            return data, "OK"
        elif isinstance(data, dict):
            # Look for records array
            if "records" in data:
                return data["records"], "OK"
            elif "calculations" in data:
                return data["calculations"], "OK"
            else:
                # Maybe the dict is a single record, wrap it
                return [data], "OK"
        else:
            return None, "Unexpected JSON format"
            
    except json.JSONDecodeError as e:
        # Corrupted JSON - mark as skip
        return None, f"CORRUPT_LEGACY: {str(e)[:50]}"
    except Exception as e:
        return None, f"Error loading: {str(e)}"


def load_modern_json(modern_dir: Path, pdb_id: str, subdir: str) -> Tuple[Optional[List], str]:
    """Load modern JSON records for a PDB."""
    json_file = modern_dir / subdir / f"{pdb_id}.json"
    
    if not json_file.exists():
        return None, f"Modern JSON not found: {json_file}"
    
    try:
        with open(json_file) as f:
            data = json.load(f)
        
        if isinstance(data, list):
            return data, "OK"
        elif isinstance(data, dict):
            if "records" in data:
                return data["records"], "OK"
            else:
                return [data], "OK"
        else:
            return None, "Unexpected JSON format"
            
    except Exception as e:
        return None, f"Error loading: {str(e)}"


def validate_stage(
    pdb_id: str,
    stage_id: str,
    project_root: Path,
    modern_exe: Path,
    temp_dir: Path,
) -> StageResult:
    """
    Validate a single stage for a single PDB.
    
    Args:
        pdb_id: PDB ID to validate
        stage_id: Stage identifier (e.g., "stage3")
        project_root: Project root path
        modern_exe: Path to generate_modern_json executable
        temp_dir: Temporary directory for modern JSON output
    
    Returns:
        StageResult with validation results
    """
    stage_config = STAGES[stage_id]
    result = StageResult(pdb_id=pdb_id, stage_id=stage_id, passed=False)
    
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        result.error = f"PDB file not found: {pdb_file}"
        return result
    
    # Load legacy JSON
    legacy_dir = project_root / "data" / "json_legacy" / stage_config["legacy_dir"]
    legacy_records, legacy_msg = load_legacy_json(legacy_dir, pdb_id)
    
    if legacy_records is None:
        result.error = f"Legacy: {legacy_msg}"
        # Mark as skipped if legacy is corrupt (not a real failure)
        if "CORRUPT_LEGACY" in legacy_msg:
            result.skipped = True
        return result
    
    result.legacy_count = len(legacy_records)
    
    # Generate modern JSON to temp directory
    output_dir = temp_dir / pdb_id
    gen_ok, gen_msg = generate_modern_json(
        pdb_id, pdb_file, output_dir, modern_exe, stage_config["stage_flag"]
    )
    
    if not gen_ok:
        result.error = f"Generation: {gen_msg}"
        return result
    
    # Load modern JSON
    modern_records, modern_msg = load_modern_json(
        output_dir, pdb_id, stage_config["modern_dir"]
    )
    
    if modern_records is None:
        result.error = f"Modern: {modern_msg}"
        return result
    
    result.modern_count = len(modern_records)
    
    # Compare using appropriate comparison function
    compare_func = stage_config["compare_func"]
    try:
        # Use tolerance of 1e-5 for floating point comparisons
        # (1e-6 can trigger false positives from floating point precision)
        if stage_id == "stage8":  # step_params doesn't have tolerance param
            comparison = compare_func(legacy_records, modern_records)
        elif stage_id == "stage6":  # selection has different signature
            comparison = compare_func(legacy_records, modern_records, tolerance=1e-5)
        else:
            comparison = compare_func(legacy_records, modern_records, tolerance=1e-5)
        
        # Extract results from comparison object
        # Different comparators may have different attribute names
        if hasattr(comparison, 'matched'):
            result.matched = comparison.matched
        if hasattr(comparison, 'mismatched_values'):
            result.mismatched = len(comparison.mismatched_values)
        elif hasattr(comparison, 'mismatches'):
            result.mismatched = len(comparison.mismatches)
        if hasattr(comparison, 'missing_in_modern'):
            result.missing_in_modern = len(comparison.missing_in_modern)
        if hasattr(comparison, 'extra_in_modern'):
            result.extra_in_modern = len(comparison.extra_in_modern)
        
        # Determine if passed
        # A pass means no mismatches and no missing/extra records
        has_mismatches = result.mismatched > 0
        has_missing = result.missing_in_modern > 0
        has_extra = result.extra_in_modern > 0
        
        if has_mismatches or has_missing or has_extra:
            result.passed = False
            errors = []
            if has_mismatches:
                errors.append(f"{result.mismatched} mismatched values")
            if has_missing:
                errors.append(f"{result.missing_in_modern} missing in modern")
            if has_extra:
                errors.append(f"{result.extra_in_modern} extra in modern")
            result.error = "; ".join(errors)
        else:
            result.passed = True
            
    except Exception as e:
        result.error = f"Comparison error: {str(e)}"
    
    return result


def validate_stage_wrapper(args: Tuple) -> StageResult:
    """Wrapper for parallel execution - unpacks arguments."""
    return validate_stage(*args)


def run_batch_validation(
    stage_id: str,
    pdb_ids: List[str],
    project_root: Path,
    modern_exe: Path,
    temp_base: Path,
    num_workers: int = 10,
    batch_size: int = 10,
    stop_on_failure: bool = True,
    results_file: Optional[Path] = None
) -> Tuple[int, int, List[StageResult]]:
    """
    Run batch validation for a stage across multiple PDBs.
    
    Returns:
        Tuple of (passed_count, failed_count, all_results)
    """
    all_results = []
    passed = 0
    failed = 0
    
    num_batches = (len(pdb_ids) + batch_size - 1) // batch_size
    start_time = time.time()
    
    print(f"\nValidating Stage: {STAGES[stage_id]['name']}")
    print(f"  PDBs: {len(pdb_ids)}")
    print(f"  Batches: {num_batches} (batch size: {batch_size})")
    print(f"  Workers: {num_workers}")
    print()
    
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(pdb_ids))
        batch_pdb_ids = pdb_ids[start_idx:end_idx]
        
        print(f"Batch {batch_num + 1}/{num_batches} (PDBs {start_idx + 1}-{end_idx})...")
        
        # Create batch temp directory
        batch_temp = temp_base / f"batch_{batch_num + 1}"
        batch_temp.mkdir(parents=True, exist_ok=True)
        
        batch_results = []
        has_failure = False
        
        skipped = 0
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {
                executor.submit(
                    validate_stage,
                    pdb_id,
                    stage_id,
                    project_root,
                    modern_exe,
                    batch_temp
                ): pdb_id
                for pdb_id in batch_pdb_ids
            }
            
            for future in as_completed(futures):
                pdb_id = futures[future]
                try:
                    result = future.result()
                    batch_results.append(result)
                    
                    if result.skipped:
                        status = "⏭️ "
                        skipped += 1
                        print(f"  {status} {pdb_id} - SKIPPED (corrupt legacy)", end="")
                    elif result.passed:
                        status = "✅"
                        passed += 1
                        print(f"  {status} {pdb_id}", end="")
                    else:
                        status = "❌"
                        failed += 1
                        has_failure = True
                        print(f"  {status} {pdb_id} - {result.error[:80]}", end="")
                    print()
                        
                except Exception as e:
                    print(f"  ❌ {pdb_id} - Exception: {str(e)[:80]}")
                    batch_results.append(StageResult(
                        pdb_id=pdb_id,
                        stage_id=stage_id,
                        passed=False,
                        error=f"Exception: {str(e)}"
                    ))
                    failed += 1
                    has_failure = True
        
        all_results.extend(batch_results)
        
        # Cleanup batch temp if all passed
        if not has_failure:
            try:
                if batch_temp.exists():
                    shutil.rmtree(batch_temp)
            except Exception as e:
                print(f"  ⚠️  Warning: Could not clean up: {e}")
        else:
            print(f"  ⚠️  Output preserved: {batch_temp}")
            if stop_on_failure:
                print(f"\n❌ FAILURE DETECTED - Stopping")
                break
        
        # Save progress after each batch
        if results_file:
            _save_results(results_file, stage_id, all_results, start_time, pdb_ids)
        
        batch_passed = sum(1 for r in batch_results if r.passed)
        batch_skipped = sum(1 for r in batch_results if r.skipped)
        batch_failed = len(batch_results) - batch_passed - batch_skipped
        print(f"  Batch summary: {batch_passed}/{len(batch_results)} passed, {batch_skipped} skipped, {batch_failed} failed")
        print()
    
    # Final save
    if results_file:
        _save_results(results_file, stage_id, all_results, start_time, pdb_ids)
    
    return passed, failed, all_results


def _save_results(
    results_file: Path,
    stage_id: str,
    results: List[StageResult],
    start_time: float,
    all_pdb_ids: List[str]
):
    """Save results to JSON file."""
    elapsed = time.time() - start_time
    
    passed_results = [r for r in results if r.passed]
    skipped_results = [r for r in results if r.skipped]
    failed_results = [r for r in results if not r.passed and not r.skipped]
    
    tested_count = len(passed_results) + len(failed_results)
    
    data = {
        "stage_id": stage_id,
        "stage_name": STAGES[stage_id]["name"],
        "test_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_pdbs": len(all_pdb_ids),
        "processed": len(results),
        "passed": len(passed_results),
        "failed": len(failed_results),
        "skipped": len(skipped_results),
        "elapsed_seconds": round(elapsed, 2),
        "stopped_early": len(results) < len(all_pdb_ids),
        "summary": {
            "pass_rate": round(100.0 * len(passed_results) / tested_count, 2) if tested_count else 0,
            "note": f"Pass rate excludes {len(skipped_results)} skipped PDBs (corrupt legacy JSON)"
        },
        "results": {
            "passed": [{"pdb_id": r.pdb_id, "matched": r.matched} for r in passed_results],
            "skipped": [{"pdb_id": r.pdb_id, "reason": r.error} for r in skipped_results],
            "failed": [
                {
                    "pdb_id": r.pdb_id,
                    "error": r.error,
                    "legacy_count": r.legacy_count,
                    "modern_count": r.modern_count,
                    "mismatched": r.mismatched,
                    "missing_in_modern": r.missing_in_modern,
                    "extra_in_modern": r.extra_in_modern,
                }
                for r in failed_results
            ]
        }
    }
    
    results_file.parent.mkdir(parents=True, exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(data, f, indent=2)


def load_valid_pdb_ids(project_root: Path) -> List[str]:
    """Load PDB IDs from valid_pdbs_fast.json."""
    fast_pdbs_file = project_root / "data" / "valid_pdbs_fast.json"
    
    if not fast_pdbs_file.exists():
        print(f"❌ {fast_pdbs_file} not found!")
        return []
    
    with open(fast_pdbs_file) as f:
        data = json.load(f)
    
    if isinstance(data, list):
        return data
    elif isinstance(data, dict):
        # Try common keys
        for key in ["valid_pdbs_with_atoms_and_frames", "valid_pdbs"]:
            if key in data:
                return data[key]
        # Assume keys are PDB IDs
        return list(data.keys())
    
    return []


def main():
    parser = argparse.ArgumentParser(description="Batch stage validation")
    parser.add_argument("stage", choices=list(STAGES.keys()), help="Stage to validate")
    parser.add_argument("--max-pdbs", type=int, default=None, help="Max PDBs to process")
    parser.add_argument("--pdb", type=str, default=None, help="Single PDB to test")
    parser.add_argument("--workers", type=int, default=10, help="Number of parallel workers")
    parser.add_argument("--batch-size", type=int, default=10, help="Batch size")
    parser.add_argument("--no-stop", action="store_true", help="Don't stop on failure")
    parser.add_argument("--continue-stages", action="store_true", help="Continue to next stages")
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    
    # Find executables
    _, modern_exe = find_executables(project_root)
    
    if not modern_exe:
        print("❌ Modern executable (generate_modern_json) not found!")
        return 1
    
    print(f"✅ Modern executable: {modern_exe}")
    
    # Get PDB list
    if args.pdb:
        pdb_ids = [args.pdb]
        print(f"Testing single PDB: {args.pdb}")
    else:
        pdb_ids = load_valid_pdb_ids(project_root)
        if not pdb_ids:
            return 1
        print(f"Loaded {len(pdb_ids)} PDBs from valid_pdbs_fast.json")
        
        if args.max_pdbs:
            pdb_ids = pdb_ids[:args.max_pdbs]
            print(f"Limited to {len(pdb_ids)} PDBs")
    
    # Determine stages to run
    stages_to_run = [args.stage]
    if args.continue_stages:
        stage_list = list(STAGES.keys())
        start_idx = stage_list.index(args.stage)
        stages_to_run = stage_list[start_idx:]
        print(f"Will run stages: {stages_to_run}")
    
    temp_base = Path("/tmp/stage_validation")
    
    total_passed = 0
    total_failed = 0
    
    for stage_id in stages_to_run:
        results_file = project_root / "data" / "validation_results" / f"{stage_id}_{STAGES[stage_id]['legacy_dir']}_results.json"
        
        passed, failed, results = run_batch_validation(
            stage_id=stage_id,
            pdb_ids=pdb_ids,
            project_root=project_root,
            modern_exe=modern_exe,
            temp_base=temp_base / stage_id,
            num_workers=args.workers,
            batch_size=args.batch_size,
            stop_on_failure=not args.no_stop,
            results_file=results_file
        )
        
        total_passed += passed
        total_failed += failed
        
        print(f"\n{'='*60}")
        print(f"Stage {stage_id} Results:")
        print(f"  Passed: {passed}")
        print(f"  Failed: {failed}")
        print(f"  Results saved to: {results_file}")
        
        if failed > 0 and not args.no_stop:
            print(f"\n❌ Stage {stage_id} has failures - stopping")
            break
    
    # Final cleanup
    try:
        if temp_base.exists() and total_failed == 0:
            shutil.rmtree(temp_base)
            print(f"\n✅ Cleaned up temp directory")
    except Exception as e:
        print(f"\n⚠️  Warning: Could not clean up temp directory: {e}")
    
    print(f"\n{'='*60}")
    print(f"Final Summary:")
    print(f"  Total Passed: {total_passed}")
    print(f"  Total Failed: {total_failed}")
    
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    exit(main())

