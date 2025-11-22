#!/usr/bin/env python3
"""
Analyze problematic PDBs using removed atom information.
Creates a summary of which atoms are different and why.
"""

import json
import sys
import os
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional

def load_problematic_pdbs():
    """Load problematic PDB list."""
    problem_file = Path('docs/problematic_pdbs.txt')
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
        return None
    except Exception as e:
        print(f"Error loading {json_file}: {e}")
        return None

def load_removed_atoms(json_file: str) -> Optional[Dict]:
    """Load removed atoms from legacy JSON."""
    if not os.path.exists(json_file):
        return None
    
    try:
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
    except Exception as e:
        print(f"Error loading removed atoms from {json_file}: {e}")
        return None

def create_atom_key(atom: Dict) -> Tuple:
    """Create a unique key for an atom."""
    return (
        atom.get('chain_id', ''),
        atom.get('residue_seq', 0),
        atom.get('atom_name', ''),
        atom.get('alt_loc', ' ')
    )

def compare_pdb(pdb_id: str) -> Dict:
    """Compare a single PDB and analyze differences."""
    legacy_file = f'data/json_legacy/{pdb_id}.json'
    modern_file = f'data/json/{pdb_id}.json'
    
    legacy_atoms = load_json_atoms(legacy_file)
    modern_atoms = load_json_atoms(modern_file)
    removed_info = load_removed_atoms(legacy_file)
    
    if legacy_atoms is None:
        return {'pdb_id': pdb_id, 'status': 'error', 'message': 'Legacy JSON not found'}
    if modern_atoms is None:
        return {'pdb_id': pdb_id, 'status': 'error', 'message': 'Modern JSON not found'}
    
    legacy_dict = {create_atom_key(a): a for a in legacy_atoms}
    modern_dict = {create_atom_key(a): a for a in modern_atoms}
    
    legacy_keys = set(legacy_dict.keys())
    modern_keys = set(modern_dict.keys())
    
    missing = legacy_keys - modern_keys
    extra = modern_keys - legacy_keys
    common = legacy_keys & modern_keys
    
    # Check if missing atoms match removed atoms
    matched_removed = 0
    unmatched_missing = []
    
    if removed_info:
        removed_map = removed_info.get('removed_map', {})
        for key in missing:
            atom = legacy_dict[key]
            # Try matching with different alt_loc values
            found = False
            for alt_loc in [' ', 'A', '1', 'B', 'C', 'D']:
                test_key = (key[0], key[1], key[2], alt_loc)
                if test_key in removed_map:
                    matched_removed += 1
                    found = True
                    break
            if not found:
                unmatched_missing.append(atom)
    
    return {
        'pdb_id': pdb_id,
        'legacy_count': len(legacy_atoms),
        'modern_count': len(modern_atoms),
        'missing_count': len(missing),
        'extra_count': len(extra),
        'common_count': len(common),
        'removed_info': removed_info,
        'matched_removed': matched_removed,
        'unmatched_missing': len(unmatched_missing)
    }

def main():
    problematic_pdbs = load_problematic_pdbs()
    
    if not problematic_pdbs:
        print("No problematic PDBs found in docs/problematic_pdbs.txt")
        sys.exit(1)
    
    print("=" * 80)
    print(f"ANALYZING {len(problematic_pdbs)} PROBLEMATIC PDBs")
    print("=" * 80)
    print()
    
    results = []
    for i, pdb_id in enumerate(problematic_pdbs, 1):
        if i % 20 == 0:
            print(f"Progress: {i}/{len(problematic_pdbs)}...")
        result = compare_pdb(pdb_id)
        results.append(result)
    
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    
    # Group by status
    by_status = defaultdict(list)
    for r in results:
        status = 'error' if 'error' in r else 'ok'
        by_status[status].append(r)
    
    print(f"Total analyzed: {len(results)}")
    print(f"  Errors: {len(by_status.get('error', []))}")
    print(f"  Successful: {len(by_status.get('ok', []))}")
    print()
    
    # Analyze differences
    has_differences = [r for r in results if r.get('missing_count', 0) > 0 or r.get('extra_count', 0) > 0]
    no_differences = [r for r in results if r.get('missing_count', 0) == 0 and r.get('extra_count', 0) == 0 and 'error' not in r]
    
    print(f"PDBs with differences: {len(has_differences)}")
    print(f"PDBs with no differences: {len(no_differences)}")
    print()
    
    # Show PDBs with removed atoms that match differences
    matched_cases = []
    for r in has_differences:
        if r.get('removed_info') and r.get('matched_removed', 0) > 0:
            matched_cases.append(r)
    
    print(f"PDBs where missing atoms match removed atoms: {len(matched_cases)}")
    if matched_cases:
        print("\nFirst 10 cases where removed atoms explain differences:")
        for r in matched_cases[:10]:
            removed = r.get('removed_info', {})
            num_removed = removed.get('num_removed', 0)
            matched = r.get('matched_removed', 0)
            missing = r.get('missing_count', 0)
            print(f"  {r['pdb_id']}: {matched}/{missing} missing atoms match {num_removed} removed")
    
    # Show PDBs with unmatched missing atoms
    unmatched_cases = [r for r in has_differences if r.get('unmatched_missing', 0) > 0]
    print(f"\nPDBs with unmatched missing atoms: {len(unmatched_cases)}")
    if unmatched_cases:
        print("\nFirst 10 cases with unmatched missing atoms:")
        for r in unmatched_cases[:10]:
            print(f"  {r['pdb_id']}: {r.get('unmatched_missing', 0)} unmatched missing atoms")
    
    # Save detailed results
    output_file = 'docs/problematic_pdbs_analysis.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nDetailed results saved to: {output_file}")
    
    # Create text summary
    summary_file = 'docs/problematic_pdbs_analysis.txt'
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PROBLEMATIC PDBs ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total analyzed: {len(results)}\n")
        f.write(f"PDBs with differences: {len(has_differences)}\n")
        f.write(f"PDBs with no differences: {len(no_differences)}\n\n")
        
        f.write("PDBs WITH DIFFERENCES:\n")
        f.write("-" * 80 + "\n")
        for r in sorted(has_differences, key=lambda x: x.get('missing_count', 0) + x.get('extra_count', 0), reverse=True)[:50]:
            f.write(f"{r['pdb_id']}: missing={r.get('missing_count', 0)}, extra={r.get('extra_count', 0)}")
            if r.get('removed_info'):
                num_removed = r.get('removed_info', {}).get('num_removed', 0)
                matched = r.get('matched_removed', 0)
                f.write(f", removed={num_removed}, matched={matched}")
            f.write("\n")
    
    print(f"Text summary saved to: {summary_file}")

if __name__ == '__main__':
    main()

