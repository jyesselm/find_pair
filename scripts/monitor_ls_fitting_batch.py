#!/usr/bin/env python3
"""
Monitor and analyze ls_fitting batch test results.
"""

import json
import sys
from pathlib import Path
from datetime import datetime


def analyze_results(results_file: Path):
    """Analyze batch test results."""
    
    if not results_file.exists():
        print(f"Results file not found: {results_file}")
        return
    
    with open(results_file) as f:
        results = json.load(f)
    
    total = results['total']
    passed = results['passed']
    failed = results['failed']
    errors = results['errors']
    
    print("="*80)
    print("LS_FITTING BATCH TEST ANALYSIS")
    print("="*80)
    print(f"Total PDBs: {total}")
    print(f"Passed: {passed} ({passed/total*100:.1f}%)")
    print(f"Failed: {failed} ({failed/total*100:.1f}%)")
    print(f"Errors: {errors} ({errors/total*100:.1f}%)")
    print()
    
    # Analyze mismatches
    mismatch_counts = {}
    total_matched = 0
    total_records = 0
    
    for pdb_id, res in results['pdb_results'].items():
        if res['status'] == 'fail':
            matched = res.get('matched', 0)
            total = res.get('total', 0)
            if total > 0:
                mismatches = total - matched
                if mismatches not in mismatch_counts:
                    mismatch_counts[mismatches] = []
                mismatch_counts[mismatches].append(pdb_id)
                total_matched += matched
                total_records += total
        elif res['status'] == 'pass':
            matched = res.get('matched', 0)
            total = res.get('total', 0)
            total_matched += matched
            total_records += total
    
    if total_records > 0:
        print(f"Overall match rate: {total_matched}/{total_records} ({total_matched/total_records*100:.2f}%)")
        print()
    
    # Show mismatch distribution
    if mismatch_counts:
        print("Mismatch distribution:")
        for count in sorted(mismatch_counts.keys()):
            pdbs = mismatch_counts[count]
            print(f"  {count} mismatch(es): {len(pdbs)} PDBs")
            if len(pdbs) <= 5:
                print(f"    {', '.join(pdbs)}")
        print()
    
    # Show worst cases
    failed_pdbs = [(pdb, res) for pdb, res in results['pdb_results'].items() 
                   if res['status'] == 'fail']
    failed_pdbs.sort(key=lambda x: x[1].get('total', 0) - x[1].get('matched', 0), reverse=True)
    
    if failed_pdbs:
        print("Worst failures (most mismatches):")
        for pdb, res in failed_pdbs[:20]:
            matched = res.get('matched', 0)
            total = res.get('total', 0)
            mismatches = total - matched
            pct = (matched/total*100) if total > 0 else 0
            print(f"  {pdb}: {matched}/{total} matched ({pct:.1f}%), {mismatches} mismatches")


def monitor_log(log_file: Path):
    """Monitor the log file for progress."""
    
    if not log_file.exists():
        print(f"Log file not found: {log_file}")
        return
    
    with open(log_file) as f:
        lines = f.readlines()
    
    # Find summary if available
    summary_start = -1
    for i, line in enumerate(lines):
        if "LS_FITTING BATCH TEST SUMMARY" in line:
            summary_start = i
            break
    
    if summary_start >= 0:
        # Test is complete, show summary
        print("Test complete! Showing summary:")
        print("".join(lines[summary_start:]))
    else:
        # Test is still running, show progress
        test_lines = [l for l in lines if "] " in l and ("PASS" in l or "FAIL" in l or "SKIP" in l)]
        if test_lines:
            print(f"Progress: {len(test_lines)} PDBs tested")
            print("Last 10 results:")
            for line in test_lines[-10:]:
                print(line.strip())
        else:
            print("No test results yet...")


def main():
    project_root = Path(".")
    
    results_file = project_root / "data" / "validation_results" / "ls_fitting_batch_results.json"
    log_file = project_root / "data" / "validation_results" / "ls_fitting_batch.log"
    
    print(f"Checking results at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    if results_file.exists():
        analyze_results(results_file)
    elif log_file.exists():
        monitor_log(log_file)
    else:
        print("No results or log file found yet")


if __name__ == "__main__":
    main()

