#!/usr/bin/env python3
"""
Deep analysis of a small sample of problematic PDBs using removed atom tracking.

This script analyzes a few PDBs in detail to understand why they differ between
legacy and modern JSON outputs.

Usage:
    python3 scripts/analyze_sample_pdbs.py [pdb_id1] [pdb_id2] ...
    If no PDBs provided, analyzes first 5 from problematic list.
"""

import json
import sys
import os
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set, Optional

def load_json_atoms(json_file: str) -> Optional[List[Dict]]:
    """Load atoms from JSON file."""
    if not os.path.exists(json_file):
        return None
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        for calc in data.get('calculations', []):
            if calc.get('type') == 'pdb_atoms':
                return calc.get('atoms', [])
    except json.JSONDecodeError as e:
        print(f"  ERROR: JSON decode error in {json_file}: {e}")
        return None
    
    return []

def load_removed_atoms_from_legacy(json_file: str) -> Dict:
    """Load removed atom information from legacy JSON."""
    if not os.path.exists(json_file):
        return {}
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        removed_atoms = []
        removed_summary = {}
        
        for calc in data.get('calculations', []):
            if calc.get('type') == 'removed_atom':
                removed_atoms.append(calc)
            elif calc.get('type') == 'removed_atoms_summary':
                removed_summary = calc
        
        # Create lookup by atom key
        removed_map = {}
        for atom in removed_atoms:
            # Create key: (chain_id, residue_seq, residue_name, atom_name, alt_loc)
            key = (
                atom.get('chain_id', ' '),
                atom.get('residue_seq', 0),
                atom.get('residue_name', ''),
                atom.get('atom_name', ''),
                atom.get('alt_loc', ' ') if 'alt_loc' in atom else ' '
            )
            removed_map[key] = atom
        
        return {
            'removed_atoms': removed_atoms,
            'removed_summary': removed_summary,
            'removed_map': removed_map,
            'count': len(removed_atoms),
            'by_reason': Counter(a.get('reason', 'unknown') for a in removed_atoms)
        }
    except json.JSONDecodeError as e:
        print(f"  ERROR: JSON decode error in {json_file}: {e}")
        return {}
    except Exception as e:
        print(f"  ERROR loading removed atoms from {json_file}: {e}")
        return {}

def create_atom_key(atom: Dict) -> Tuple:
    """Create unique key for an atom."""
    return (
        atom.get('chain_id', ' '),
        atom.get('residue_seq', 0),
        atom.get('residue_name', ''),
        atom.get('atom_name', ''),
        atom.get('alt_loc', ' ') if 'alt_loc' in atom else ' '
    )

def analyze_single_pdb(pdb_id: str, project_root: Path) -> Dict:
    """Analyze a single PDB in detail."""
    legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    print(f"\n{'='*80}")
    print(f"ANALYZING: {pdb_id}")
    print(f"{'='*80}")
    
    # Load atoms
    legacy_atoms = load_json_atoms(str(legacy_file))
    modern_atoms = load_json_atoms(str(modern_file))
    
    if legacy_atoms is None:
        return {'error': f'Legacy JSON not found or invalid: {legacy_file}'}
    if modern_atoms is None:
        return {'error': f'Modern JSON not found or invalid: {modern_file}'}
    
    # Load removed atoms info
    removed_info = load_removed_atoms_from_legacy(str(legacy_file))
    
    # Create atom dictionaries
    legacy_dict = {create_atom_key(a): a for a in legacy_atoms}
    modern_dict = {create_atom_key(a): a for a in modern_atoms}
    
    legacy_keys = set(legacy_dict.keys())
    modern_keys = set(modern_dict.keys())
    
    # Find differences
    missing_in_modern = legacy_keys - modern_keys  # In legacy but not modern
    extra_in_modern = modern_keys - legacy_keys    # In modern but not legacy
    common = legacy_keys & modern_keys
    
    # Analyze removed atoms vs missing atoms
    removed_map = removed_info.get('removed_map', {})
    matched_removed = []
    unmatched_missing = []
    
    for key in missing_in_modern:
        # Try exact match first
        if key in removed_map:
            matched_removed.append((key, removed_map[key], 'exact'))
        else:
            # Try without alt_loc (case where alt_loc differs)
            key_no_alt = (key[0], key[1], key[2], key[3], ' ')
            if key_no_alt in removed_map:
                matched_removed.append((key, removed_map[key_no_alt], 'no_alt'))
            else:
                unmatched_missing.append((key, legacy_dict[key]))
    
    # Check coordinate differences in common atoms
    coord_diffs = []
    for key in common:
        leg_xyz = legacy_dict[key].get('xyz', [])
        mod_xyz = modern_dict[key].get('xyz', [])
        if len(leg_xyz) == 3 and len(mod_xyz) == 3:
            max_diff = max(abs(leg_xyz[i] - mod_xyz[i]) for i in range(3))
            if max_diff > 0.001:  # More than 0.001 Angstrom difference
                coord_diffs.append({
                    'key': key,
                    'legacy_xyz': leg_xyz,
                    'modern_xyz': mod_xyz,
                    'max_diff': max_diff
                })
    
    # Analyze removed atoms by reason
    removed_by_reason = defaultdict(list)
    for atom_info in removed_info.get('removed_atoms', []):
        reason = atom_info.get('reason', 'unknown')
        removed_by_reason[reason].append(atom_info)
    
    # Print detailed analysis
    print(f"\nATOM COUNTS:")
    print(f"  Legacy: {len(legacy_atoms)}")
    print(f"  Modern: {len(modern_atoms)}")
    print(f"  Difference: {len(legacy_atoms) - len(modern_atoms)}")
    
    print(f"\nREMOVED ATOMS (from legacy JSON):")
    print(f"  Total removed: {removed_info.get('count', 0)}")
    print(f"  By reason:")
    for reason, count in removed_info.get('by_reason', {}).items():
        print(f"    {reason}: {count}")
    
    print(f"\nCOMPARISON:")
    print(f"  Missing in modern (should match removed): {len(missing_in_modern)}")
    print(f"    Matched with removed atoms: {len([m for m in matched_removed if m[2] == 'exact'])}")
    print(f"    Matched without alt_loc: {len([m for m in matched_removed if m[2] == 'no_alt'])}")
    print(f"    Unmatched: {len(unmatched_missing)}")
    print(f"  Extra in modern: {len(extra_in_modern)}")
    print(f"  Common: {len(common)}")
    print(f"  Common with coord diffs > 0.001Å: {len(coord_diffs)}")
    
    # Show examples of unmatched missing atoms
    if unmatched_missing:
        print(f"\nUNMATCHED MISSING ATOMS (first 10):")
        for key, atom in unmatched_missing[:10]:
            print(f"  {key}: {atom}")
    
    # Show examples of extra atoms in modern
    if extra_in_modern:
        print(f"\nEXTRA ATOMS IN MODERN (first 10):")
        for key in list(extra_in_modern)[:10]:
            atom = modern_dict[key]
            print(f"  {key}: {atom}")
    
    # Show examples of removed atoms by reason
    print(f"\nREMOVED ATOM SAMPLES BY REASON:")
    for reason, atoms in list(removed_by_reason.items())[:5]:
        print(f"\n  {reason} ({len(atoms)} total, showing first 3):")
        for atom in atoms[:3]:
            key = create_atom_key(atom) if 'chain_id' in atom else None
            print(f"    {key}: chain={atom.get('chain_id')}, "
                  f"res={atom.get('residue_name')}{atom.get('residue_seq')}, "
                  f"atom={atom.get('atom_name')}, "
                  f"alt_loc={atom.get('alt_loc', ' ')}")
    
    # Show coordinate differences if any
    if coord_diffs:
        print(f"\nCOORDINATE DIFFERENCES (first 5):")
        for diff in coord_diffs[:5]:
            print(f"  {diff['key']}:")
            print(f"    Legacy: {diff['legacy_xyz']}")
            print(f"    Modern: {diff['modern_xyz']}")
            print(f"    Max diff: {diff['max_diff']:.6f} Å")
    
    return {
        'pdb_id': pdb_id,
        'legacy_count': len(legacy_atoms),
        'modern_count': len(modern_atoms),
        'removed_count': removed_info.get('count', 0),
        'removed_by_reason': dict(removed_info.get('by_reason', {})),
        'missing_count': len(missing_in_modern),
        'extra_count': len(extra_in_modern),
        'matched_removed': len(matched_removed),
        'unmatched_missing': len(unmatched_missing),
        'coord_diffs_count': len(coord_diffs)
    }

def load_problematic_pdbs(project_root: Path) -> List[str]:
    """Load problematic PDB list."""
    problem_file = project_root / 'docs' / 'problematic_pdbs.txt'
    if not problem_file.exists():
        return []
    
    pdbs = []
    with open(problem_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            pdb_id = line.split()[0] if line.split() else None
            if pdb_id:
                pdbs.append(pdb_id)
    return pdbs

def main():
    project_root = Path(__file__).parent.parent.absolute()
    
    # Get PDBs to analyze
    if len(sys.argv) > 1:
        pdb_ids = sys.argv[1:]
    else:
        # Default: analyze first 5 problematic PDBs
        problematic_pdbs = load_problematic_pdbs(project_root)
        pdb_ids = problematic_pdbs[:5]
        print(f"Analyzing first 5 problematic PDBs: {pdb_ids}")
    
    results = []
    for pdb_id in pdb_ids:
        try:
            result = analyze_single_pdb(pdb_id, project_root)
            if 'error' not in result:
                results.append(result)
        except Exception as e:
            print(f"\nERROR analyzing {pdb_id}: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    if results:
        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        print(f"\n{'PDB':<10} {'Legacy':<8} {'Modern':<8} {'Removed':<8} {'Missing':<8} {'Extra':<8} {'Matched':<8} {'Unmatched':<8}")
        print("-" * 80)
        for r in results:
            print(f"{r['pdb_id']:<10} {r['legacy_count']:<8} {r['modern_count']:<8} "
                  f"{r['removed_count']:<8} {r['missing_count']:<8} {r['extra_count']:<8} "
                  f"{r['matched_removed']:<8} {r['unmatched_missing']:<8}")
        
        print(f"\nRemoved atoms by reason (across all PDBs):")
        all_reasons = defaultdict(int)
        for r in results:
            for reason, count in r.get('removed_by_reason', {}).items():
                all_reasons[reason] += count
        for reason, count in sorted(all_reasons.items(), key=lambda x: -x[1]):
            print(f"  {reason}: {count}")

if __name__ == '__main__':
    main()

