#!/usr/bin/env python3
"""
Debug a specific failing residue to understand the atom matching difference.
"""

import json
import sys
import os
from pathlib import Path

from x3dna_json_compare import PdbFileReader

def main():
    if len(sys.argv) < 4:
        print("Usage: debug_specific_failure.py <pdb_name> <chain_id> <seq_num>")
        sys.exit(1)
    
    pdb_name = sys.argv[1]
    chain_id = sys.argv[2]
    seq_num = int(sys.argv[3])
    
    json_file = 'docs/frame_calculation_failures.json'
    with open(json_file) as f:
        data = json.load(f)
    
    # Find the failure
    failure = None
    for f in data['failures']:
        if (f['pdb_name'] == pdb_name and 
            f['chain_id'] == chain_id and 
            f['seq_num'] == seq_num):
            failure = f
            break
    
    if not failure:
        print(f"Failure not found for {pdb_name} {chain_id}:{seq_num}")
        sys.exit(1)
    
    print(f"=== Failure Analysis for {pdb_name} {chain_id}:{seq_num} ===\n")
    print(f"Residue name: {failure['residue_name']}")
    print(f"Base type: {failure['base_type']}")
    print(f"Failure reason: {failure['failure_reason']}\n")
    
    print(f"Our num_matched: {failure['our']['num_matched']}")
    print(f"Legacy num_matched: {failure['legacy']['num_matched']}\n")
    
    print("Our matched atoms:")
    for i, atom in enumerate(failure['our']['matched_atoms'], 1):
        print(f"  {i}. {atom}")
    
    print("\nLegacy matched atoms:")
    for i, atom in enumerate(failure['legacy']['matched_atoms'], 1):
        print(f"  {i}. {atom}")
    
    print("\nAtoms in our list but not in legacy:")
    our_set = set(failure['our']['matched_atoms'])
    leg_set = set(failure['legacy']['matched_atoms'])
    only_ours = our_set - leg_set
    for atom in sorted(only_ours):
        print(f"  - {atom}")
    
    print("\nAtoms in legacy list but not in ours:")
    only_legacy = leg_set - our_set
    for atom in sorted(only_legacy):
        print(f"  - {atom}")
    
    print(f"\nOur RMS: {failure['our']['rms']:.6f}")
    print(f"Legacy RMS: {failure['legacy']['rms']:.6f}")
    print(f"RMS diff: {failure['differences']['rms_diff']:.6f}")
    print(f"Max rot diff: {failure['differences']['max_rot_diff']:.6f}")
    print(f"Max trans diff: {failure['differences']['max_trans_diff']:.6f}")
    
    # Check if we can find the PDB and inspect the residue
    pdb_file = Path(f"data/pdb/{pdb_name}.pdb")
    if pdb_file.exists():
        print(f"\n=== Inspecting PDB file: {pdb_file} ===")
        print(f"Looking for residue {chain_id}:{seq_num} {failure['residue_name']}\n")
        
        # Use PdbFileReader from lib
        reader = PdbFileReader(pdb_file)
        residue_atoms = reader.get_residue_atoms(chain_id, seq_num, ' ')
        
        atom_lines = []
        for atom_name, line_num, line in residue_atoms:
            atom_lines.append((atom_name, line))
        
        print(f"Found {len(atom_lines)} atoms in residue:")
        # Sort by atom name
        for atom_name, line in sorted(atom_lines):
            print(f"  {atom_name}: {line}")
        
        # Check for ring atoms
        ring_atoms_purine = [" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "]
        ring_atoms_pyrimidine = [" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "]
        
    print("\nRing atoms present in PDB:")
    if failure['base_type'] in ['A', 'G']:
        ring_atoms = ring_atoms_purine
    else:
        ring_atoms = ring_atoms_pyrimidine
    
    found_ring_atoms = []
    for ring_atom in ring_atoms:
        for atom_name, line in atom_lines:
            if atom_name == ring_atom:
                found_ring_atoms.append(ring_atom)
                print(f"  ✓ {ring_atom}")
                break
        else:
            print(f"  ✗ {ring_atom} (missing)")
    
    # Print actual PDB lines for matched atoms
    if pdb_file.exists():
        print("\n=== PDB Lines for Matched Atoms ===")
        print("\nOur matched atoms (with PDB lines):")
        reader = PdbFileReader(pdb_file)
        for atom_name in failure['our']['matched_atoms']:
            # Get PDB line for this atom
            atom_lines_dict = reader.get_atom_lines_by_names(
                chain_id, seq_num, ' ', [atom_name]
            )
            if atom_name in atom_lines_dict:
                line_num, line = atom_lines_dict[atom_name]
                print(f"  {atom_name}: Line {line_num}")
                print(f"    {line}")
            else:
                print(f"  {atom_name}: NOT FOUND IN PDB")
        
        print("\nLegacy matched atoms (with PDB lines):")
        for atom_name in failure['legacy']['matched_atoms']:
            # Get PDB line for this atom
            atom_lines_dict = reader.get_atom_lines_by_names(
                chain_id, seq_num, ' ', [atom_name]
            )
            if atom_name in atom_lines_dict:
                line_num, line = atom_lines_dict[atom_name]
                print(f"  {atom_name}: Line {line_num}")
                print(f"    {line}")
            else:
                print(f"  {atom_name}: NOT FOUND IN PDB")

if __name__ == '__main__':
    main()

