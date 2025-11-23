#!/usr/bin/env python3
"""
Analyze why modern JSON is missing frame calculations.

This script investigates specific PDBs to understand why residues
that have frame calculations in legacy JSON don't appear in modern JSON.
"""

import sys
import json
from pathlib import Path
from typing import Dict, List, Optional

from x3dna_json_compare import PdbFileReader, JsonComparator


def analyze_pdb(pdb_id: str, project_root: Path) -> Dict:
    """Analyze a single PDB to understand missing frame calculations."""
    legacy_file = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    pdb_file = project_root / 'data' / 'pdb' / f'{pdb_id}.pdb'
    
    if not all(f.exists() for f in [legacy_file, modern_file, pdb_file]):
        return {'error': 'Files not found'}
    
    # Load JSON files
    with open(legacy_file, 'r') as f:
        legacy_json = json.load(f)
    with open(modern_file, 'r') as f:
        modern_json = json.load(f)
    
    # Extract frame calculation records
    legacy_frames = []
    modern_frames = []
    
    for calc in legacy_json.get('calculations', []):
        if calc.get('type') in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
            legacy_frames.append(calc)
    
    for calc in modern_json.get('calculations', []):
        if calc.get('type') in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
            modern_frames.append(calc)
    
    # Build maps by (chain_id, residue_seq, insertion)
    legacy_map = {}
    modern_map = {}
    
    for rec in legacy_frames:
        key = (
            rec.get('chain_id', ''),
            rec.get('residue_seq', 0),
            rec.get('insertion', ' ')
        )
        if key not in legacy_map:
            legacy_map[key] = []
        legacy_map[key].append(rec)
    
    for rec in modern_frames:
        key = (
            rec.get('chain_id', ''),
            rec.get('residue_seq', 0),
            rec.get('insertion', ' ')
        )
        if key not in modern_map:
            modern_map[key] = []
        modern_map[key].append(rec)
    
    # Find missing residues
    missing_keys = set(legacy_map.keys()) - set(modern_map.keys())
    
    # Analyze missing residues
    missing_analysis = []
    pdb_reader = PdbFileReader(pdb_file)
    
    for key in list(missing_keys)[:10]:  # Analyze first 10 missing
        chain_id, residue_seq, insertion = key
        legacy_recs = legacy_map[key]
        
        # Get residue info from legacy
        residue_name = None
        base_type = None
        rms_fit = None
        num_matched = None
        matched_atoms = []
        
        for rec in legacy_recs:
            if rec.get('type') == 'base_frame_calc':
                residue_name = rec.get('residue_name', '')
                base_type = rec.get('base_type', '?')
                matched_atoms = rec.get('matched_atoms', [])
            elif rec.get('type') == 'ls_fitting':
                rms_fit = rec.get('rms_fit', 0.0)
                num_matched = rec.get('num_matched_atoms', rec.get('num_points', 0))
        
        # Try to get PDB line for this residue
        try:
            # Find atom lines for this residue
            atom_lines = pdb_reader.get_atom_lines_by_names(
                chain_id, residue_seq, insertion, matched_atoms[:3] if matched_atoms else []
            )
            has_atoms = len(atom_lines) > 0
        except:
            has_atoms = False
            atom_lines = []
        
        missing_analysis.append({
            'chain_id': chain_id,
            'residue_seq': residue_seq,
            'insertion': insertion,
            'residue_name': residue_name,
            'base_type': base_type,
            'rms_fit': rms_fit,
            'num_matched_atoms': num_matched,
            'matched_atoms': matched_atoms,
            'has_atoms_in_pdb': has_atoms,
            'atom_count': len(atom_lines)
        })
    
    return {
        'pdb_id': pdb_id,
        'legacy_frame_count': len(legacy_map),
        'modern_frame_count': len(modern_map),
        'missing_count': len(missing_keys),
        'missing_samples': missing_analysis
    }


def main():
    project_root = Path(__file__).parent.parent.parent.absolute()
    
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_missing_frames.py <pdb_id> [pdb_id2 ...]")
        print("\nExample:")
        print("  python3 analyze_missing_frames.py 157D 161D 1BNA")
        return 1
    
    pdb_ids = sys.argv[1:]
    
    print("=" * 80)
    print("MISSING FRAME ANALYSIS")
    print("=" * 80)
    print()
    
    for pdb_id in pdb_ids:
        print(f"Analyzing {pdb_id}...")
        result = analyze_pdb(pdb_id, project_root)
        
        if 'error' in result:
            print(f"  Error: {result['error']}")
            continue
        
        print(f"  Legacy frames: {result['legacy_frame_count']}")
        print(f"  Modern frames: {result['modern_frame_count']}")
        print(f"  Missing: {result['missing_count']}")
        print()
        
        if result['missing_samples']:
            print("  Sample missing residues:")
            for sample in result['missing_samples'][:5]:
                print(f"    {sample['residue_name']} {sample['chain_id']}:{sample['residue_seq']}{sample['insertion']}")
                print(f"      Base: {sample['base_type']}, RMS: {sample['rms_fit']:.4f}, Matched: {sample['num_matched_atoms']}")
                print(f"      Atoms in PDB: {sample['has_atoms_in_pdb']} ({sample['atom_count']} found)")
                print(f"      Matched atoms: {sample['matched_atoms']}")
                print()
        print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

