#!/usr/bin/env python3
"""
Test script for --fix-indices option

This script tests the --fix-indices option by:
1. Running find_pair with and without --fix-indices
2. Comparing results with legacy
3. Reporting improvements

Usage:
    python3 scripts/test_fix_indices_option.py <pdb_id> [legacy_json_file]
    
Example:
    python3 scripts/test_fix_indices_option.py 6CAQ
    python3 scripts/test_fix_indices_option.py 6CAQ data/json_legacy/base_frame_calc/6CAQ.json
"""

import sys
import os
import subprocess
import json
import tempfile
from pathlib import Path

def find_executable(name):
    """Find executable in build directory"""
    build_path = Path("build") / name
    if build_path.exists():
        return str(build_path)
    return name

def run_find_pair(pdb_file, output_file, fix_indices=False, legacy_json=None):
    """Run find_pair_app with optional --fix-indices"""
    cmd = [find_executable("find_pair_app")]
    
    if fix_indices:
        if legacy_json:
            cmd.append(f"--fix-indices={legacy_json}")
        else:
            cmd.append("--fix-indices")
    
    cmd.append(pdb_file)
    cmd.append(output_file)
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}")
        return None
    
    return result.stdout

def compare_outputs(output1, output2):
    """Compare two .inp output files"""
    if not os.path.exists(output1) or not os.path.exists(output2):
        return None
    
    with open(output1) as f:
        lines1 = f.readlines()
    
    with open(output2) as f:
        lines2 = f.readlines()
    
    # Count pairs (lines starting with numbers)
    pairs1 = [l for l in lines1 if l.strip() and l.strip()[0].isdigit()]
    pairs2 = [l for l in lines2 if l.strip() and l.strip()[0].isdigit()]
    
    return {
        'pairs1': len(pairs1),
        'pairs2': len(pairs2),
        'lines1': len(lines1),
        'lines2': len(lines2),
        'identical': lines1 == lines2
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/test_fix_indices_option.py <pdb_id> [legacy_json_file]")
        print("Example: python3 scripts/test_fix_indices_option.py 6CAQ")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    legacy_json = sys.argv[2] if len(sys.argv) > 2 else None
    
    pdb_file = f"data/pdb/{pdb_id}.pdb"
    if not os.path.exists(pdb_file):
        print(f"ERROR: PDB file not found: {pdb_file}")
        sys.exit(1)
    
    # Auto-detect legacy JSON if not provided
    if not legacy_json:
        legacy_json_path = f"data/json_legacy/base_frame_calc/{pdb_id}.json"
        if os.path.exists(legacy_json_path):
            legacy_json = legacy_json_path
            print(f"Auto-detected legacy JSON: {legacy_json}")
        else:
            print(f"WARNING: Legacy JSON not found: {legacy_json_path}")
            print("Continuing without --fix-indices comparison...")
    
    print(f"\nTesting --fix-indices option for {pdb_id}")
    print("=" * 60)
    
    # Create temporary output files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        output_no_fix = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.inp', delete=False) as f:
        output_with_fix = f.name
    
    try:
        # Run without --fix-indices
        print("\n1. Running WITHOUT --fix-indices:")
        print("-" * 60)
        result1 = run_find_pair(pdb_file, output_no_fix, fix_indices=False)
        if result1:
            print(result1)
        
        # Run with --fix-indices
        if legacy_json:
            print("\n2. Running WITH --fix-indices:")
            print("-" * 60)
            result2 = run_find_pair(pdb_file, output_with_fix, fix_indices=True, legacy_json=legacy_json)
            if result2:
                print(result2)
            
            # Compare outputs
            print("\n3. Comparing Results:")
            print("-" * 60)
            comparison = compare_outputs(output_no_fix, output_with_fix)
            if comparison:
                print(f"Without fix: {comparison['pairs1']} pairs")
                print(f"With fix:    {comparison['pairs2']} pairs")
                print(f"Difference:   {comparison['pairs2'] - comparison['pairs1']} pairs")
                print(f"Files identical: {comparison['identical']}")
            else:
                print("Could not compare outputs")
        else:
            print("\n2. Skipping --fix-indices test (no legacy JSON)")
        
        print("\n✅ Test complete!")
        print(f"\nOutput files:")
        print(f"  Without fix: {output_no_fix}")
        if legacy_json:
            print(f"  With fix:    {output_with_fix}")
    
    finally:
        # Clean up (optional - comment out to keep files for inspection)
        # os.unlink(output_no_fix)
        # if legacy_json:
        #     os.unlink(output_with_fix)
        pass

if __name__ == "__main__":
    main()

