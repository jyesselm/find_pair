#!/usr/bin/env python3
"""
Validation script that uses removed atom information from legacy JSON
to verify that the modern parser correctly filters out atoms.

Usage:
    python3 scripts/validate_removed_atoms.py <pdb_id>
"""

import json
import sys
import os

def load_removed_atoms(json_file):
    """Load removed atoms from legacy JSON file."""
    if not os.path.exists(json_file):
        print(f"Error: Legacy JSON file not found: {json_file}")
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    removed = [x for x in data.get('calculations', []) if x.get('type') == 'removed_atom']
    
    # Create a map of removed atoms by (chain, residue_seq, atom_name, alt_loc)
    removed_map = {}
    for r in removed:
        line = r.get('pdb_line', '')
        if len(line) > 16:
            alt_loc = line[16]
        else:
            alt_loc = ' '
        
        key = (r.get('chain_id', ''), r.get('residue_seq', 0), r.get('atom_name', ''), alt_loc)
        removed_map[key] = r
    
    # Also get summary
    summary = [x for x in data.get('calculations', []) if x.get('type') == 'removed_atoms_summary']
    num_removed = summary[0].get('num_removed', 0) if summary else 0
    
    return {
        'removed_atoms': removed,
        'removed_map': removed_map,
        'num_removed': num_removed
    }

def load_modern_atoms(json_file):
    """Load atoms from modern JSON file."""
    if not os.path.exists(json_file):
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    for calc in data.get('calculations', []):
        if calc.get('type') == 'pdb_atoms':
            return calc.get('atoms', [])
    
    return []

def validate_filtering(pdb_id, legacy_json_dir='data/json_legacy', modern_json_dir='data/json'):
    """Validate that modern parser correctly filters removed atoms."""
    legacy_json = os.path.join(legacy_json_dir, f'{pdb_id}.json')
    modern_json = os.path.join(modern_json_dir, f'{pdb_id}.json')
    
    print(f"Validating filtering for {pdb_id}...")
    print("=" * 70)
    
    # Load removed atoms
    removed_info = load_removed_atoms(legacy_json)
    if removed_info is None:
        return False
    
    num_removed = removed_info['num_removed']
    removed_map = removed_info['removed_map']
    
    print(f"\nRemoved atoms in legacy: {num_removed}")
    print(f"Unique removed atom keys: {len(removed_map)}")
    
    if num_removed == 0:
        print("\n✓ No atoms were removed - nothing to validate")
        return True
    
    # Group by reason
    reasons = {}
    for r in removed_info['removed_atoms']:
        reason = r.get('reason', 'unknown')
        if reason not in reasons:
            reasons[reason] = []
        reasons[reason].append(r)
    
    print("\nRemoved by reason:")
    for reason, atoms in reasons.items():
        print(f"  {reason}: {len(atoms)}")
    
    # Load modern atoms
    modern_atoms = load_modern_atoms(modern_json)
    if modern_atoms is None:
        print(f"\n⚠ Warning: Modern JSON not found: {modern_json}")
        print("  Cannot validate - please generate modern JSON first")
        return False
    
    print(f"\nModern JSON has {len(modern_atoms)} atoms")
    
    # Check if any removed atoms are in modern JSON
    found_removed = []
    for atom in modern_atoms:
        chain = atom.get('chain_id', '')
        resseq = atom.get('residue_seq', 0)
        atom_name = atom.get('atom_name', '')
        alt_loc = atom.get('alt_loc', ' ')
        
        key = (chain, resseq, atom_name, alt_loc)
        if key in removed_map:
            found_removed.append((atom, removed_map[key]))
    
    print(f"\n{'=' * 70}")
    if len(found_removed) == 0:
        print("✓ SUCCESS: No removed atoms found in modern JSON")
        print("  The modern parser correctly filters out all removed atoms!")
        return True
    else:
        print(f"✗ ERROR: Found {len(found_removed)} removed atoms in modern JSON:")
        print("\nFirst 10 violations:")
        for i, (atom, removed_info) in enumerate(found_removed[:10], 1):
            print(f"\n  {i}. {atom.get('atom_name')} in chain {atom.get('chain_id')} "
                  f"residue {atom.get('residue_seq')} alt_loc='{atom.get('alt_loc')}'")
            print(f"     Should be removed because: {removed_info.get('reason')}")
            print(f"     PDB line: {removed_info.get('pdb_line', '')[:70]}...")
        
        if len(found_removed) > 10:
            print(f"\n  ... and {len(found_removed) - 10} more violations")
        
        return False

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/validate_removed_atoms.py <pdb_id>")
        print("Example: python3 scripts/validate_removed_atoms.py 4P9R")
        sys.exit(1)
    
    pdb_id = sys.argv[1].upper()
    success = validate_filtering(pdb_id)
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()

