#!/usr/bin/env python3
"""
Create a detailed summary of atom differences between legacy and modern JSON outputs.
Uses removed atom information to explain differences.

Usage:
    python3 scripts/atom_diff_summary.py <pdb_id>
"""

import json
import sys
import os
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional

def load_json_atoms(json_file: str) -> List[Dict]:
    """Load atoms from JSON file."""
    if not os.path.exists(json_file):
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    for calc in data.get('calculations', []):
        if calc.get('type') == 'pdb_atoms':
            return calc.get('atoms', [])
    
    return []

def load_removed_atoms(json_file: str) -> Dict:
    """Load removed atoms from legacy JSON."""
    if not os.path.exists(json_file):
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    removed = [x for x in data.get('calculations', []) if x.get('type') == 'removed_atom']
    
    # Create map by (chain, residue_seq, atom_name, alt_loc)
    removed_map = {}
    for r in removed:
        line = r.get('pdb_line', '')
        alt_loc = line[16] if len(line) > 16 else ' '
        key = (r.get('chain_id', ''), r.get('residue_seq', 0), r.get('atom_name', ''), alt_loc)
        removed_map[key] = r
    
    # Get summary
    summary = [x for x in data.get('calculations', []) if x.get('type') == 'removed_atoms_summary']
    num_removed = summary[0].get('num_removed', 0) if summary else 0
    
    # Group by reason
    reasons = defaultdict(list)
    for r in removed:
        reason = r.get('reason', 'unknown')
        reasons[reason].append(r)
    
    return {
        'removed_atoms': removed,
        'removed_map': removed_map,
        'num_removed': num_removed,
        'by_reason': dict(reasons)
    }

def create_atom_key(atom: Dict) -> Tuple:
    """Create a unique key for an atom."""
    return (
        atom.get('chain_id', ''),
        atom.get('residue_seq', 0),
        atom.get('atom_name', ''),
        atom.get('alt_loc', ' ')
    )

def create_atom_dict(atoms: List[Dict]) -> Dict[Tuple, Dict]:
    """Create a dictionary keyed by atom identifier."""
    atom_dict = {}
    for atom in atoms:
        key = create_atom_key(atom)
        atom_dict[key] = atom
    return atom_dict

def compare_atoms(pdb_id: str, legacy_json: str = 'data/json_legacy', modern_json: str = 'data/json') -> Dict:
    """Compare atoms between legacy and modern JSON files."""
    legacy_file = os.path.join(legacy_json, f'{pdb_id}.json')
    modern_file = os.path.join(modern_json, f'{pdb_id}.json')
    
    # Load atoms
    legacy_atoms = load_json_atoms(legacy_file)
    modern_atoms = load_json_atoms(modern_file)
    
    if legacy_atoms is None:
        return {'error': f'Legacy JSON not found: {legacy_file}'}
    if modern_atoms is None:
        return {'error': f'Modern JSON not found: {modern_file}'}
    
    # Create dictionaries
    legacy_dict = create_atom_dict(legacy_atoms)
    modern_dict = create_atom_dict(modern_atoms)
    
    legacy_keys = set(legacy_dict.keys())
    modern_keys = set(modern_dict.keys())
    
    # Find differences
    missing = legacy_keys - modern_keys  # In legacy but not modern
    extra = modern_keys - legacy_keys    # In modern but not legacy
    common = legacy_keys & modern_keys
    
    # Check for coordinate/field differences in common atoms
    field_diffs = []
    for key in common:
        leg_atom = legacy_dict[key]
        mod_atom = modern_dict[key]
        
        diffs = []
        if leg_atom.get('xyz') != mod_atom.get('xyz'):
            leg_xyz = leg_atom.get('xyz', [])
            mod_xyz = mod_atom.get('xyz', [])
            if len(leg_xyz) == 3 and len(mod_xyz) == 3:
                max_diff = max(abs(leg_xyz[i] - mod_xyz[i]) for i in range(3))
                if max_diff > 0.001:
                    diffs.append(f'xyz diff: {max_diff:.6f}')
        
        if leg_atom.get('residue_name') != mod_atom.get('residue_name'):
            diffs.append(f'residue_name: {leg_atom.get("residue_name")} != {mod_atom.get("residue_name")}')
        
        if leg_atom.get('chain_id') != mod_atom.get('chain_id'):
            diffs.append(f'chain_id: {leg_atom.get("chain_id")} != {mod_atom.get("chain_id")}')
        
        if leg_atom.get('alt_loc') != mod_atom.get('alt_loc'):
            diffs.append(f'alt_loc: {leg_atom.get("alt_loc")} != {mod_atom.get("alt_loc")}')
        
        if diffs:
            field_diffs.append((key, leg_atom, mod_atom, diffs))
    
    # Load removed atoms information
    removed_info = load_removed_atoms(legacy_file)
    
    return {
        'pdb_id': pdb_id,
        'legacy_count': len(legacy_atoms),
        'modern_count': len(modern_atoms),
        'missing': list(missing),
        'extra': list(extra),
        'common': len(common),
        'field_diffs': field_diffs,
        'missing_atoms': [legacy_dict[k] for k in missing],
        'extra_atoms': [modern_dict[k] for k in extra],
        'removed_info': removed_info
    }

def print_summary(comparison: Dict):
    """Print a detailed summary of atom differences."""
    pdb_id = comparison.get('pdb_id', 'UNKNOWN')
    
    print("=" * 80)
    print(f"ATOM DIFFERENCE SUMMARY: {pdb_id}")
    print("=" * 80)
    print()
    
    # Basic counts
    print(f"Legacy atoms:     {comparison.get('legacy_count', 0)}")
    print(f"Modern atoms:     {comparison.get('modern_count', 0)}")
    print(f"Difference:       {comparison.get('modern_count', 0) - comparison.get('legacy_count', 0)}")
    print()
    
    # Missing atoms
    missing = comparison.get('missing', [])
    print(f"Missing atoms (in legacy but not modern): {len(missing)}")
    
    # Extra atoms
    extra = comparison.get('extra', [])
    print(f"Extra atoms (in modern but not legacy):   {len(extra)}")
    
    # Common atoms
    print(f"Common atoms:                              {comparison.get('common', 0)}")
    print()
    
    # Removed atoms info
    removed_info = comparison.get('removed_info')
    if removed_info:
        print("=" * 80)
        print("REMOVED ATOMS (from legacy JSON)")
        print("=" * 80)
        print(f"Total removed: {removed_info.get('num_removed', 0)}")
        
        by_reason = removed_info.get('by_reason', {})
        if by_reason:
            print("\nRemoved by reason:")
            for reason, atoms in sorted(by_reason.items()):
                print(f"  {reason}: {len(atoms)}")
        print()
    
    # Missing atoms details
    if missing:
        print("=" * 80)
        print(f"MISSING ATOMS ({len(missing)} total)")
        print("=" * 80)
        
        missing_atoms = comparison.get('missing_atoms', [])
        
        # Group by chain and residue
        by_residue = defaultdict(list)
        for atom in missing_atoms:
            chain = atom.get('chain_id', '?')
            resseq = atom.get('residue_seq', 0)
            by_residue[(chain, resseq)].append(atom)
        
        print(f"\nGrouped by residue ({len(by_residue)} residues):")
        for (chain, resseq), atoms in sorted(by_residue.items())[:20]:
            print(f"\n  Chain {chain}, Residue {resseq}: {len(atoms)} atoms")
            for atom in atoms[:5]:
                alt = atom.get('alt_loc', ' ')
                alt_str = f" alt_loc='{alt}'" if alt != ' ' else ''
                print(f"    - {atom.get('atom_name', '?')}{alt_str}")
            if len(atoms) > 5:
                print(f"    ... and {len(atoms) - 5} more")
        
        if len(by_residue) > 20:
            print(f"\n  ... and {len(by_residue) - 20} more residues")
        
        # Check if missing atoms match removed atoms
        if removed_info:
            removed_map = removed_info.get('removed_map', {})
            matched_removed = []
            unmatched_missing = []
            
            for atom in missing_atoms:
                key = create_atom_key(atom)
                # Try with different alt_loc values
                found = False
                for alt_loc in [' ', 'A', '1', 'B', 'C', 'D']:
                    test_key = (key[0], key[1], key[2], alt_loc)
                    if test_key in removed_map:
                        matched_removed.append((atom, removed_map[test_key]))
                        found = True
                        break
                if not found:
                    unmatched_missing.append(atom)
            
            print(f"\n  Matched with removed atoms: {len(matched_removed)}")
            if matched_removed:
                print("\n  First 5 matched with removed atoms:")
                for atom, removed in matched_removed[:5]:
                    reason = removed.get('reason', 'unknown')
                    print(f"    - {atom.get('atom_name')} in chain {atom.get('chain_id')} "
                          f"res {atom.get('residue_seq')}: {reason}")
            
            if unmatched_missing:
                print(f"\n  ⚠ Unmatched missing atoms: {len(unmatched_missing)}")
                print("    (These are missing but not in removed atoms list)")
                for atom in unmatched_missing[:5]:
                    print(f"    - {atom.get('atom_name')} in chain {atom.get('chain_id')} "
                          f"res {atom.get('residue_seq')}")
        print()
    
    # Extra atoms details
    if extra:
        print("=" * 80)
        print(f"EXTRA ATOMS ({len(extra)} total)")
        print("=" * 80)
        
        extra_atoms = comparison.get('extra_atoms', [])
        
        # Group by chain and residue
        by_residue = defaultdict(list)
        for atom in extra_atoms:
            chain = atom.get('chain_id', '?')
            resseq = atom.get('residue_seq', 0)
            by_residue[(chain, resseq)].append(atom)
        
        print(f"\nGrouped by residue ({len(by_residue)} residues):")
        for (chain, resseq), atoms in sorted(by_residue.items())[:20]:
            print(f"\n  Chain {chain}, Residue {resseq}: {len(atoms)} atoms")
            for atom in atoms[:5]:
                alt = atom.get('alt_loc', ' ')
                alt_str = f" alt_loc='{alt}'" if alt != ' ' else ''
                print(f"    - {atom.get('atom_name', '?')}{alt_str}")
            if len(atoms) > 5:
                print(f"    ... and {len(atoms) - 5} more")
        
        if len(by_residue) > 20:
            print(f"\n  ... and {len(by_residue) - 20} more residues")
        print()
    
    # Field differences
    field_diffs = comparison.get('field_diffs', [])
    if field_diffs:
        print("=" * 80)
        print(f"FIELD DIFFERENCES ({len(field_diffs)} atoms with field differences)")
        print("=" * 80)
        print("\nFirst 10 atoms with field differences:")
        for key, leg_atom, mod_atom, diffs in field_diffs[:10]:
            chain, resseq, atom_name, alt_loc = key
            alt_str = f" alt_loc='{alt_loc}'" if alt_loc != ' ' else ''
            print(f"\n  {atom_name} in chain {chain} res {resseq}{alt_str}:")
            for diff in diffs:
                print(f"    - {diff}")
        
        if len(field_diffs) > 10:
            print(f"\n  ... and {len(field_diffs) - 10} more atoms with field differences")
        print()
    
    # Summary conclusion
    print("=" * 80)
    if len(missing) == 0 and len(extra) == 0 and len(field_diffs) == 0:
        print("✓ PERFECT MATCH - No differences found!")
    else:
        print("⚠ DIFFERENCES FOUND")
        if removed_info and removed_info.get('num_removed', 0) > 0:
            print(f"  Note: {removed_info.get('num_removed', 0)} atoms were removed during legacy parsing")
            print("  Check if missing atoms match removed atoms above.")
    print("=" * 80)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/atom_diff_summary.py <pdb_id>")
        print("Example: python3 scripts/atom_diff_summary.py 4P9R")
        sys.exit(1)
    
    pdb_id = sys.argv[1].upper()
    comparison = compare_atoms(pdb_id)
    
    if 'error' in comparison:
        print(f"Error: {comparison['error']}")
        sys.exit(1)
    
    print_summary(comparison)

if __name__ == '__main__':
    main()

