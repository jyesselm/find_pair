#!/usr/bin/env python3
"""
Diagnose why modern code isn't calculating frames.

This script checks:
1. Are atoms being parsed correctly?
2. Are residue types being detected correctly?
3. Are templates being loaded?
4. Is atom matching working?
"""

import json
import sys
from pathlib import Path


def analyze_modern_json(pdb_id: str, project_root: Path):
    """Analyze modern JSON to see why frames aren't being calculated."""
    modern_json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    legacy_json_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    
    if not modern_json_file.exists():
        print(f"ERROR: Modern JSON not found: {modern_json_file}")
        return
    
    if not legacy_json_file.exists():
        print(f"ERROR: Legacy JSON not found: {legacy_json_file}")
        return
    
    with open(modern_json_file) as f:
        modern_json = json.load(f)
    
    with open(legacy_json_file) as f:
        legacy_json = json.load(f)
    
    # Count frame calculations
    modern_frame_calcs = [r for r in modern_json.get('calculations', [])
                         if r.get('type') in ['base_frame_calc', 'frame_calc', 'ls_fitting']]
    legacy_frame_calcs = [r for r in legacy_json.get('calculations', [])
                          if r.get('type') in ['base_frame_calc', 'frame_calc', 'ls_fitting']]
    
    print(f"Modern frame calculations: {len(modern_frame_calcs)}")
    print(f"Legacy frame calculations: {len(legacy_frame_calcs)}")
    print()
    
    # Check pdb_atoms record
    modern_atoms_rec = None
    legacy_atoms_rec = None
    
    for rec in modern_json.get('calculations', []):
        if rec.get('type') == 'pdb_atoms':
            modern_atoms_rec = rec
            break
    
    for rec in legacy_json.get('calculations', []):
        if rec.get('type') == 'pdb_atoms':
            legacy_atoms_rec = rec
            break
    
    if modern_atoms_rec:
        modern_atoms = modern_atoms_rec.get('atoms', [])
        print(f"Modern PDB atoms: {len(modern_atoms)}")
        
        # Count by residue
        residues = {}
        for atom in modern_atoms:
            key = (atom.get('chain_id', ''), atom.get('residue_seq', 0), atom.get('insertion', ' '))
            if key not in residues:
                residues[key] = {
                    'name': atom.get('residue_name', ''),
                    'count': 0,
                    'atoms': []
                }
            residues[key]['count'] += 1
            residues[key]['atoms'].append(atom.get('atom_name', ''))
        
        print(f"Modern residues: {len(residues)}")
        print()
        
        # Show first few residues
        print("First 10 residues in modern JSON:")
        for i, (key, data) in enumerate(list(residues.items())[:10]):
            chain, seq, ins = key
            print(f"  {data['name']} {chain}:{seq}{ins} - {data['count']} atoms")
            # Show ring atoms if present
            ring_atoms = [a for a in data['atoms'] if any(ring in a for ring in ['C4', 'N3', 'C2', 'N1', 'C6', 'C5', 'N7', 'C8', 'N9'])]
            if ring_atoms:
                print(f"    Ring atoms: {ring_atoms[:5]}")
        print()
    
    if legacy_atoms_rec:
        legacy_atoms = legacy_atoms_rec.get('atoms', [])
        print(f"Legacy PDB atoms: {len(legacy_atoms)}")
        print()
    
    # Compare frame calculations
    print("Frame Calculation Comparison:")
    print("-" * 80)
    
    # Build maps
    legacy_map = {}
    for rec in legacy_frame_calcs:
        key = (rec.get('chain_id', ''), rec.get('residue_seq', 0), rec.get('insertion', ' '))
        legacy_map[key] = rec
    
    modern_map = {}
    for rec in modern_frame_calcs:
        key = (rec.get('chain_id', ''), rec.get('residue_seq', 0), rec.get('insertion', ' '))
        modern_map[key] = rec
    
    # Find residues in legacy but not modern
    missing_in_modern = set(legacy_map.keys()) - set(modern_map.keys())
    if missing_in_modern:
        print(f"\nResidues in Legacy but NOT in Modern ({len(missing_in_modern)}):")
        for key in list(missing_in_modern)[:10]:
            chain, seq, ins = key
            leg_rec = legacy_map[key]
            name = leg_rec.get('residue_name', '')
            base = leg_rec.get('base_type', '?')
            print(f"  {name} {base} {chain}:{seq}{ins}")
            if 'matched_atoms' in leg_rec:
                print(f"    Legacy matched: {leg_rec['matched_atoms']}")
        print()
    
    # Find residues in modern but not legacy
    missing_in_legacy = set(modern_map.keys()) - set(legacy_map.keys())
    if missing_in_legacy:
        print(f"\nResidues in Modern but NOT in Legacy ({len(missing_in_legacy)}):")
        for key in list(missing_in_legacy)[:10]:
            chain, seq, ins = key
            mod_rec = modern_map[key]
            name = mod_rec.get('residue_name', '')
            base = mod_rec.get('base_type', '?')
            print(f"  {name} {base} {chain}:{seq}{ins}")
        print()
    
    # Show matching residues
    matching = set(legacy_map.keys()) & set(modern_map.keys())
    if matching:
        print(f"\nMatching residues ({len(matching)}):")
        for key in list(matching)[:5]:
            chain, seq, ins = key
            leg_rec = legacy_map[key]
            mod_rec = modern_map[key]
            leg_atoms = leg_rec.get('matched_atoms', [])
            mod_atoms = mod_rec.get('matched_atoms', [])
            if set(leg_atoms) != set(mod_atoms):
                print(f"  {leg_rec.get('residue_name', '')} {chain}:{seq}{ins}:")
                print(f"    Legacy: {leg_atoms}")
                print(f"    Modern: {mod_atoms}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 diagnose_modern_frame_calculation.py <PDB_ID>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    project_root = Path(__file__).parent.parent
    
    analyze_modern_json(pdb_id, project_root)

