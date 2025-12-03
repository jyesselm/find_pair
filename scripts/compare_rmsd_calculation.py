#!/usr/bin/env python3
"""
Compare RMSD calculations between legacy and modern for specific residues.
"""

import subprocess
import json
import re

def test_residue(pdb_id, chain, seq):
    """Test a specific residue's RMSD calculation."""
    
    print(f"\n=== Testing {pdb_id} chain {chain} seq {seq} ===")
    
    # Run legacy
    print("\nLegacy output:")
    result = subprocess.run(
        [f"./org/build/bin/find_pair_analyze", f"data/pdb/{pdb_id}.pdb", "tmp/legacy_test"],
        capture_output=True,
        text=True
    )
    
    # Look for the residue in stderr
    for line in result.stderr.split('\n'):
        if f"chain {chain}" in line and f"#{seq}" in line:
            print(f"  {line}")
    
    # Run modern with structural variant debug (need to rebuild with debug first)
    print("\nModern output:")
    result = subprocess.run(
        [f"./build/generate_modern_json", f"data/pdb/{pdb_id}.pdb", "tmp/modern_test", "--stage=frames"],
        capture_output=True,
        text=True
    )
    
    # Look for RMSD debug output
    for line in result.stderr.split('\n'):
        if f"chain {chain}" in line and f"seq {seq}" in line:
            print(f"  {line}")

if __name__ == "__main__":
    # Test 2YR in 9CJI
    test_residue("9CJI", "C", "7")
    
    # Test other problematic cases
    test_residue("7S36", "T", "6")  # 2YR with lower RMSD

