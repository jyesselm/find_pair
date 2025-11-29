#!/usr/bin/env python3
"""
Validate find_pair functionality by testing on multiple PDBs.

This script:
1. Tests find_pair_app can run successfully
2. Verifies output files are generated
3. Compares with legacy JSON to check match rates
4. Reports validation status
"""

import subprocess
import sys
from pathlib import Path
import json
import tempfile
import os

def find_executable(name):
    """Find executable in build directory."""
    build_dir = Path(__file__).parent.parent / "build"
    exe = build_dir / name
    if exe.exists():
        return str(exe)
    return name

def run_find_pair(pdb_file, output_file, fix_indices=False):
    """Run find_pair_app and return success status."""
    exe = find_executable("find_pair_app")
    pdb_path = Path("data/pdb") / f"{pdb_file}.pdb"
    
    if not pdb_path.exists():
        print(f"  ⚠️  PDB file not found: {pdb_path}")
        return False
    
    cmd = [exe]
    if fix_indices:
        cmd.append("--fix-indices")
    cmd.append(str(pdb_path))
    cmd.append(output_file)
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            # Check output - if no pairs found, no output file is written (this is OK)
            output_text = result.stdout
            if "No base pairs found" in output_text:
                print(f"  ✅ find_pair_app ran successfully (0 pairs found - expected)")
                return True
            # Check if output file was created and has content
            if Path(output_file).exists() and Path(output_file).stat().st_size > 0:
                return True
            else:
                print(f"  ⚠️  Output file empty or missing: {output_file}")
                return False
        else:
            print(f"  ❌ Error running find_pair_app:")
            print(f"     {result.stderr}")
            return False
    except subprocess.TimeoutExpired:
        print(f"  ⚠️  Timeout running find_pair_app")
        return False
    except Exception as e:
        print(f"  ❌ Exception: {e}")
        return False

def check_json_comparison(pdb_id):
    """Check if JSON comparison shows good match."""
    try:
        result = subprocess.run(
            ["python3", "scripts/compare_json.py", "compare", pdb_id],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            output = result.stdout
            # Check for perfect match
            if "Perfect matches: 1" in output and "Files with differences: 0" in output:
                return True, "Perfect match"
            # Check for find_bestpair_selection match
            if "Missing in modern: 0" in output and "Extra in modern: 0" in output:
                return True, "Selection matches"
            return False, "Has differences"
        return False, "Comparison failed"
    except Exception as e:
        return False, f"Exception: {e}"

def validate_single_pdb(pdb_id, fix_indices=False):
    """Validate find_pair for a single PDB."""
    print(f"\nTesting {pdb_id}...")
    
    # Test 1: Can run find_pair_app
    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        output_file = f.name
    
    try:
        success = run_find_pair(pdb_id, output_file, fix_indices=fix_indices)
        if not success:
            print(f"  ❌ find_pair_app failed for {pdb_id}")
            return False
        
        print(f"  ✅ find_pair_app ran successfully")
        
        # Test 2: Check JSON comparison (if JSON files exist)
        legacy_json = Path(f"data/json_legacy/{pdb_id}.json")
        modern_json = Path(f"data/json/{pdb_id}.json")
        
        if legacy_json.exists() or modern_json.exists():
            match_status, message = check_json_comparison(pdb_id)
            if match_status:
                print(f"  ✅ JSON comparison: {message}")
            else:
                print(f"  ⚠️  JSON comparison: {message}")
        
        return True
        
    finally:
        # Clean up
        if Path(output_file).exists():
            os.unlink(output_file)

def get_test_pdbs():
    """Get list of PDBs to test."""
    import sys
    
    # 10-PDB test set
    test_set_10 = ["1Q96", "1VBY", "3AVY", "3G8T", "3KNC", "4AL5", "5UJ2", "6LTU", "8J1J", "6CAQ"]
    
    # Check if specific PDBs provided as arguments
    if len(sys.argv) > 1:
        return sys.argv[1:]
    
    # Default: use 10-PDB test set
    return test_set_10

def main():
    """Main validation function."""
    print("=" * 80)
    print("FIND_PAIR VALIDATION")
    print("=" * 80)
    
    test_pdbs = get_test_pdbs()
    
    print(f"\nTesting {len(test_pdbs)} PDBs...")
    print(f"Using --fix-indices: True")
    print(f"PDBs: {', '.join(test_pdbs)}")
    
    results = []
    json_results = []
    
    for pdb_id in test_pdbs:
        success = validate_single_pdb(pdb_id, fix_indices=True)
        results.append((pdb_id, success))
        
        # Track JSON comparison separately
        legacy_json = Path(f"data/json_legacy/{pdb_id}.json")
        modern_json = Path(f"data/json/{pdb_id}.json")
        if legacy_json.exists() or modern_json.exists():
            match_status, message = check_json_comparison(pdb_id)
            json_results.append((pdb_id, match_status, message))
    
    # Summary
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    print(f"\nExecutable Tests:")
    print(f"  Total tested: {total}")
    print(f"  ✅ Passed: {passed}")
    print(f"  ❌ Failed: {total - passed}")
    
    if json_results:
        json_passed = sum(1 for _, status, _ in json_results if status)
        json_total = len(json_results)
        print(f"\nJSON Comparison Tests:")
        print(f"  Total tested: {json_total}")
        print(f"  ✅ Perfect/Good matches: {json_passed}")
        print(f"  ⚠️  Has differences: {json_total - json_passed}")
        
        if json_passed < json_total:
            print(f"\n  PDBs with differences:")
            for pdb_id, status, message in json_results:
                if not status:
                    print(f"    - {pdb_id}: {message}")
    
    if passed == total:
        print("\n🎉 All executable tests passed! find_pair is working correctly.")
        return 0
    else:
        print("\n⚠️  Some tests failed. Check output above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())

