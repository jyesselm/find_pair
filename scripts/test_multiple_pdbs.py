#!/usr/bin/env python3
"""
Test multiple PDBs for mismatches after the is_nucleotide() fix

Usage:
    python3 scripts/test_multiple_pdbs.py [--limit N] [--pdbs PDB1 PDB2 ...]
"""

import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse

def test_pdb(pdb_id, project_root):
    """Test a single PDB for mismatches"""
    try:
        # Use longer timeout for large PDBs (check file size)
        pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
        timeout = 300  # 5 minutes default
        if pdb_file.exists():
            # Check file size - large files need more time
            file_size = pdb_file.stat().st_size
            if file_size > 10_000_000:  # > 10MB
                timeout = 600  # 10 minutes for very large files
            elif file_size > 5_000_000:  # > 5MB
                timeout = 450  # 7.5 minutes for large files
        
        result = subprocess.run(
            [sys.executable, str(project_root / "scripts" / "analyze_mismatched_pairs.py"), pdb_id],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=project_root
        )
        
        output = result.stdout + result.stderr
        
        # Check for perfect match
        if "Perfect match" in output:
            return (pdb_id, "PASS", None, None)
        
        # Check for mismatches
        missing_count = 0
        extra_count = 0
        
        for line in output.split('\n'):
            if 'Missing Pairs' in line or 'missing pairs' in line.lower():
                # Extract number
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        missing_count = int(parts[1].strip().split()[0])
                    except:
                        pass
            if 'Extra Pairs' in line or 'extra pairs' in line.lower():
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        extra_count = int(parts[1].strip().split()[0])
                    except:
                        pass
        
        if missing_count > 0 or extra_count > 0:
            return (pdb_id, "FAIL", missing_count, extra_count)
        else:
            return (pdb_id, "PASS", None, None)
            
    except subprocess.TimeoutExpired:
        return (pdb_id, "TIMEOUT", None, None)
    except Exception as e:
        return (pdb_id, "ERROR", str(e), None)

def main():
    parser = argparse.ArgumentParser(description='Test multiple PDBs for mismatches')
    parser.add_argument('--limit', type=int, help='Limit number of PDBs to test')
    parser.add_argument('--pdbs', nargs='+', help='Specific PDB IDs to test')
    parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers')
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    
    # Get list of PDBs to test
    if args.pdbs:
        pdb_ids = args.pdbs
    else:
        legacy_dir = project_root / "data" / "json_legacy" / "find_bestpair_selection"
        pdb_ids = sorted([f.stem for f in legacy_dir.glob("*.json")])
    
    if args.limit:
        pdb_ids = pdb_ids[:args.limit]
    
    print(f"Testing {len(pdb_ids)} PDBs...")
    print("=" * 70)
    
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(test_pdb, pdb_id, project_root): pdb_id 
                   for pdb_id in pdb_ids}
        
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            pdb_id, status, missing, extra = result
            if status == "PASS":
                print(f"✅ {pdb_id}: Perfect match")
            elif status == "FAIL":
                print(f"❌ {pdb_id}: {missing} missing, {extra} extra")
            elif status == "TIMEOUT":
                print(f"⏱️  {pdb_id}: Timeout")
            else:
                print(f"⚠️  {pdb_id}: {status}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    passed = [r for r in results if r[1] == "PASS"]
    failed = [r for r in results if r[1] == "FAIL"]
    timeouts = [r for r in results if r[1] == "TIMEOUT"]
    errors = [r for r in results if r[1] == "ERROR"]
    
    print(f"Total tested: {len(results)}")
    print(f"✅ Perfect matches: {len(passed)} ({100*len(passed)/len(results):.1f}%)")
    print(f"❌ Mismatches: {len(failed)}")
    print(f"⏱️  Timeouts: {len(timeouts)}")
    print(f"⚠️  Errors: {len(errors)}")
    
    if failed:
        print("\nFailed PDBs:")
        for pdb_id, status, missing, extra in failed:
            print(f"  {pdb_id}: {missing} missing, {extra} extra")
    
    if timeouts:
        print("\nTimeout PDBs:")
        for pdb_id, status, _, _ in timeouts:
            print(f"  {pdb_id}")
    
    if errors:
        print("\nError PDBs:")
        for pdb_id, status, error, _ in errors:
            print(f"  {pdb_id}: {error}")

if __name__ == "__main__":
    main()

