#!/usr/bin/env python3
"""
Debug RMSD calculation for a specific residue to understand differences.
"""

import subprocess
import re

def extract_atoms_from_pdb(pdb_id, res_name, chain, seq):
    """Extract ring atom coordinates from PDB file."""
    
    ring_atoms = ["C4", "N3", "C2", "N1", "C6", "C5", "N7", "C8", "N9"]
    
    print(f"\n=== Ring atoms in PDB for {res_name} {chain} {seq} ===")
    
    with open(f"data/pdb/{pdb_id}.pdb") as f:
        for line in f:
            if "HETATM" in line or "ATOM" in line:
                # Parse PDB line
                atom_name = line[12:16].strip()
                res = line[17:20].strip()
                ch = line[21]
                res_seq = int(line[22:26])
                
                if res == res_name and ch == chain and res_seq == seq:
                    if atom_name in ring_atoms:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        print(f"  {atom_name:4s}: ({x:8.3f}, {y:8.3f}, {z:8.3f})")
    
    # Also check for C1' or C1R
    print(f"\n=== Sugar atoms ===")
    with open(f"data/pdb/{pdb_id}.pdb") as f:
        for line in f:
            if "HETATM" in line or "ATOM" in line:
                atom_name = line[12:16].strip()
                res = line[17:20].strip()
                ch = line[21]
                res_seq = int(line[22:26])
                
                if res == res_name and ch == chain and res_seq == seq:
                    if atom_name in ["C1'", "C1R", "O4'", "O4R"]:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        print(f"  {atom_name:4s}: ({x:8.3f}, {y:8.3f}, {z:8.3f})")

def get_legacy_rmsd(pdb_id, res_name, chain, seq):
    """Get RMSD from legacy code."""
    
    print(f"\n=== Legacy RMSD ===")
    
    # Run legacy and capture stderr
    result = subprocess.run(
        ["./org/build/bin/find_pair_analyze", f"data/pdb/{pdb_id}.pdb", "tmp/legacy_debug"],
        capture_output=True,
        text=True
    )
    
    # Look for the residue in output
    for line in result.stderr.split('\n'):
        if res_name in line and f"chain {chain}" in line:
            print(f"  {line}")

def get_modern_rmsd(pdb_id, res_name, chain, seq):
    """Get RMSD from modern code with debug output."""
    
    print(f"\n=== Modern RMSD ===")
    
    # Run modern
    result = subprocess.run(
        ["./build/generate_modern_json", f"data/pdb/{pdb_id}.pdb", "tmp/modern_debug", "--stage=frames"],
        capture_output=True,
        text=True
    )
    
    # Look for debug output
    in_target = False
    for line in result.stderr.split('\n'):
        if f"residue: {res_name}" in line and f"chain {chain}" in line and f"seq {seq}" in line:
            in_target = True
        
        if in_target:
            print(f"  {line}")
            if "FAILED" in line or "PASSED" in line or "Skipping" in line:
                break

if __name__ == "__main__":
    # Test 2YR in 9CJI
    pdb_id = "9CJI"
    res_name = "2YR"
    chain = "C"
    seq = 7
    
    print(f"Debugging {res_name} in {pdb_id} chain {chain} seq {seq}")
    print("=" * 70)
    
    extract_atoms_from_pdb(pdb_id, res_name, chain, seq)
    get_legacy_rmsd(pdb_id, res_name, chain, seq)
    get_modern_rmsd(pdb_id, res_name, chain, seq)

