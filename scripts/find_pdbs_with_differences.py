#!/usr/bin/env python3
"""
Quick script to find PDBs that still have differences between legacy and modern JSON.
"""

import json
import sys
from pathlib import Path
from collections import defaultdict

def get_atom_count(json_file):
    """Get atom count from JSON file."""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        for calc in data.get('calculations', []):
            if calc.get('type') == 'pdb_atoms':
                return len(calc.get('atoms', []))
    except:
        pass
    return None

def main():
    project_root = Path(__file__).parent.parent.absolute()
    
    # Load problematic PDBs
    problem_file = project_root / 'docs' / 'problematic_pdbs.txt'
    if not problem_file.exists():
        print("No problematic_pdbs.txt found")
        return
    
    pdbs = []
    with open(problem_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            pdb_id = line.split()[0] if line.split() else None
            if pdb_id:
                pdbs.append(pdb_id)
    
    print(f"Checking {len(pdbs)} problematic PDBs...\n")
    
    different = []
    for pdb_id in pdbs:
        legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
        modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
        
        if not legacy_file.exists() or not modern_file.exists():
            continue
        
        legacy_count = get_atom_count(str(legacy_file))
        modern_count = get_atom_count(str(modern_file))
        
        if legacy_count is None or modern_count is None:
            continue
        
        if legacy_count != modern_count:
            different.append((pdb_id, legacy_count, modern_count, legacy_count - modern_count))
    
    if different:
        print(f"Found {len(different)} PDBs with differences:\n")
        print(f"{'PDB':<10} {'Legacy':<10} {'Modern':<10} {'Diff':<10}")
        print("-" * 40)
        for pdb_id, leg, mod, diff in sorted(different, key=lambda x: abs(x[3]), reverse=True):
            print(f"{pdb_id:<10} {leg:<10} {mod:<10} {diff:<10}")
    else:
        print("✓ All checked PDBs match perfectly!")

if __name__ == '__main__':
    main()

