#!/usr/bin/env python3
"""
Batch test ls_fitting JSON generation and comparison across all fast PDBs.

Runs ls_fitting comparison for all PDBs in valid_pdbs_fast.json.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List
from datetime import datetime
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from tests_python.integration.test_ls_fitting import test_single_pdb
from scripts.test_utils import find_executables, load_valid_pdbs_fast


def test_single_pdb_wrapper(args):
    """Wrapper for test_single_pdb to work with multiprocessing.Pool.map"""
    pdb_id, output_dir, legacy_exe, modern_exe, project_root, index, total = args
    
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    
    if not pdb_file.exists():
        print(f"[{index}/{total}] {pdb_id}: SKIP (PDB file not found)", flush=True)
        return pdb_id, {"status": "skip", "error": "PDB file not found"}
    
    try:
        result = test_single_pdb(
            pdb_id, pdb_file, output_dir, 
            legacy_exe, modern_exe, project_root
        )
        
        if result["compare_ok"]:
            print(f"[{index}/{total}] {pdb_id}: PASS ({result.get('matched', 0)}/{result.get('total', 0)} matched)", flush=True)
            return pdb_id, {
                "status": "pass",
                "matched": result.get("matched", 0),
                "total": result.get("total", 0)
            }
        else:
            error_msg = "; ".join(result.get("errors", ["Unknown error"]))
            print(f"[{index}/{total}] {pdb_id}: FAIL - {error_msg}", flush=True)
            return pdb_id, {
                "status": "fail",
                "errors": result.get("errors", []),
                "matched": result.get("matched", 0),
                "total": result.get("total", 0)
            }
            
    except Exception as e:
        print(f"[{index}/{total}] {pdb_id}: ERROR - {str(e)}", flush=True)
        return pdb_id, {"status": "error", "error": str(e)}


def run_batch_test(pdb_ids: List[str], output_dir: Path, 
                   legacy_exe: Path, modern_exe: Path, 
                   project_root: Path, max_pdbs: int = None,
                   num_workers: int = None) -> Dict:
    """Run ls_fitting tests on a batch of PDBs in parallel."""
    
    if max_pdbs:
        pdb_ids = pdb_ids[:max_pdbs]
    
    if num_workers is None:
        num_workers = min(cpu_count(), 20)  # Default to CPU count or 20, whichever is smaller
    
    print(f"Testing ls_fitting for {len(pdb_ids)} PDBs...")
    print(f"Using {num_workers} parallel workers")
    print(f"Output directory: {output_dir}")
    print()
    
    results = {
        "total": len(pdb_ids),
        "passed": 0,
        "failed": 0,
        "errors": 0,
        "pdb_results": {},
        "start_time": datetime.now().isoformat(),
        "num_workers": num_workers,
    }
    
    # Prepare arguments for parallel processing
    args_list = [
        (pdb_id, output_dir, legacy_exe, modern_exe, project_root, i+1, len(pdb_ids))
        for i, pdb_id in enumerate(pdb_ids)
    ]
    
    # Run tests in parallel
    with Pool(processes=num_workers) as pool:
        pdb_results = pool.map(test_single_pdb_wrapper, args_list)
    
    # Aggregate results
    for pdb_id, result in pdb_results:
        results["pdb_results"][pdb_id] = result
        
        if result["status"] == "pass":
            results["passed"] += 1
        elif result["status"] == "fail":
            results["failed"] += 1
        else:  # skip or error
            results["errors"] += 1
    
    results["end_time"] = datetime.now().isoformat()
    return results


def print_summary(results: Dict):
    """Print summary statistics."""
    print("\n" + "="*80)
    print("LS_FITTING BATCH TEST SUMMARY")
    print("="*80)
    print(f"Total PDBs tested: {results['total']}")
    print(f"Passed: {results['passed']} ({results['passed']/results['total']*100:.1f}%)")
    print(f"Failed: {results['failed']} ({results['failed']/results['total']*100:.1f}%)")
    print(f"Errors: {results['errors']} ({results['errors']/results['total']*100:.1f}%)")
    print(f"Start time: {results['start_time']}")
    print(f"End time: {results['end_time']}")
    print("="*80)
    
    # Show some failed examples
    if results['failed'] > 0:
        print("\nFailed PDBs (first 10):")
        failed_pdbs = [pdb for pdb, res in results['pdb_results'].items() 
                       if res['status'] == 'fail']
        for pdb in failed_pdbs[:10]:
            res = results['pdb_results'][pdb]
            matched = res.get('matched', 0)
            total = res.get('total', 0)
            errors = "; ".join(res.get('errors', []))[:100]
            print(f"  {pdb}: {matched}/{total} matched - {errors}")


def main():
    parser = argparse.ArgumentParser(description="Batch test ls_fitting comparison")
    parser.add_argument("--max-pdbs", type=int, help="Maximum number of PDBs to test")
    parser.add_argument("--output", type=str, 
                        default="data/validation_results/ls_fitting_batch_results.json",
                        help="Output file for results")
    parser.add_argument("--workers", type=int, default=20,
                        help="Number of parallel workers (default: 20)")
    args = parser.parse_args()
    
    # Setup
    legacy_exe, modern_exe = find_executables(project_root)
    
    if not modern_exe:
        print("Error: Modern executable not found")
        print(f"Expected: {project_root / 'build' / 'generate_modern_json'}")
        sys.exit(1)
    
    if not legacy_exe:
        print("Warning: Legacy executable not found - tests will fail")
    
    # Load fast PDBs
    try:
        pdb_ids = load_valid_pdbs_fast(project_root)
        print(f"Loaded {len(pdb_ids)} fast PDBs from valid_pdbs_fast.json")
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Create output directory
    output_dir = project_root / "tmp" / "ls_fitting_batch"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run batch test
    results = run_batch_test(
        pdb_ids, output_dir, legacy_exe, modern_exe, 
        project_root, max_pdbs=args.max_pdbs,
        num_workers=args.workers
    )
    
    # Print summary
    print_summary(results)
    
    # Save results
    output_file = project_root / args.output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_file}")
    
    # Return exit code
    if results['failed'] > 0 or results['errors'] > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

