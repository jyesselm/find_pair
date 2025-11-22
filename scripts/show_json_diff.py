#!/usr/bin/env python3
"""
Tool to show detailed differences between generated and legacy JSON files.

Usage:
    python3 scripts/show_json_diff.py <pdb_id>
    python3 scripts/show_json_diff.py <pdb_id> --verbose
    python3 scripts/show_json_diff.py <pdb_id> --missing-only
    python3 scripts/show_json_diff.py <pdb_id> --extra-only
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional


def get_atom_key(atom: dict) -> Tuple[str, int, str, str]:
    """Create a unique key for an atom."""
    chain_id = atom.get("chain_id", "")
    residue_seq = atom.get("residue_seq", 0)
    insertion = atom.get("insertion", "")
    atom_name = atom.get("atom_name", "")
    return (chain_id, residue_seq, insertion, atom_name)


def format_atom_info(atom: dict, show_debug: bool = False) -> str:
    """Format atom information for display."""
    # Build residue identifier: sequence number + insertion code if present
    res_seq = atom.get('residue_seq', 'N/A')
    insertion = atom.get('insertion', '')
    if insertion and insertion != ' ':
        res_id = f"{res_seq}{insertion}"
    else:
        res_id = str(res_seq)
    
    parts = [
        f"atom_name={atom.get('atom_name', 'N/A')}",
        f"residue={atom.get('residue_name', 'N/A')}",
        f"chain={atom.get('chain_id', 'N/A')}",
        f"residue_seq={res_id}",  # Shows sequence number (from PDB columns 23-26)
    ]
    
    if atom.get("alt_loc") and atom.get("alt_loc") != ' ':
        parts.append(f"alt_loc={atom.get('alt_loc')}")
    
    if show_debug:
        if atom.get("line_number"):
            parts.append(f"line={atom.get('line_number')}")
        if atom.get("atom_serial"):
            parts.append(f"serial={atom.get('atom_serial')}")
        if atom.get("model_number"):
            parts.append(f"model={atom.get('model_number')}")
        if atom.get("occupancy") is not None:
            parts.append(f"occ={atom.get('occupancy')}")
        if atom.get("original_atom_name") and atom.get("original_atom_name") != atom.get("atom_name"):
            parts.append(f"orig_name={atom.get('original_atom_name')}")
    
    return " | ".join(parts)


def compare_atoms(gen_atom: dict, leg_atom: dict) -> List[str]:
    """Compare two atoms and return list of differences."""
    differences = []
    
    # Compare all fields
    fields_to_check = [
        "atom_name", "residue_name", "chain_id", "residue_seq",
        "record_type", "alt_loc", "insertion", "xyz"
    ]
    
    for field in fields_to_check:
        gen_val = gen_atom.get(field)
        leg_val = leg_atom.get(field)
        
        if field == "xyz":
            # Compare coordinates with tolerance
            if gen_val and leg_val:
                if len(gen_val) == 3 and len(leg_val) == 3:
                    for i, (g, l) in enumerate(zip(gen_val, leg_val)):
                        if abs(g - l) > 0.001:  # 0.001 Angstrom tolerance
                            differences.append(f"xyz[{i}]: {g} vs {l} (diff: {abs(g-l):.6f})")
        else:
            if gen_val != leg_val:
                differences.append(f"{field}: '{gen_val}' vs '{leg_val}'")
    
    return differences


def show_differences(pdb_id: str, verbose: bool = False, missing_only: bool = False, extra_only: bool = False):
    """Show differences between generated and legacy JSON files."""
    gen_file = Path(f"data/json/{pdb_id}.json")
    leg_file = Path(f"data/json_legacy/{pdb_id}.json")
    
    if not gen_file.exists():
        print(f"Error: Generated JSON file not found: {gen_file}")
        return
    
    if not leg_file.exists():
        print(f"Error: Legacy JSON file not found: {leg_file}")
        return
    
    # Load JSON files
    with open(gen_file) as f:
        gen_json = json.load(f)
    with open(leg_file) as f:
        leg_json = json.load(f)
    
    # Extract pdb_atoms records
    gen_atoms_rec = next((c for c in gen_json.get("calculations", []) if c.get("type") == "pdb_atoms"), None)
    leg_atoms_rec = next((c for c in leg_json.get("calculations", []) if c.get("type") == "pdb_atoms"), None)
    
    if not gen_atoms_rec or not leg_atoms_rec:
        print("Error: Could not find pdb_atoms records")
        return
    
    gen_atoms = gen_atoms_rec.get("atoms", [])
    leg_atoms = leg_atoms_rec.get("atoms", [])
    
    gen_count = gen_atoms_rec.get("num_atoms", len(gen_atoms))
    leg_count = leg_atoms_rec.get("num_atoms", len(leg_atoms))
    
    # Create maps
    gen_map = {get_atom_key(a): a for a in gen_atoms}
    leg_map = {get_atom_key(a): a for a in leg_atoms}
    
    gen_keys = set(gen_map.keys())
    leg_keys = set(leg_map.keys())
    
    missing = leg_keys - gen_keys  # In legacy but not generated
    extra = gen_keys - leg_keys    # In generated but not legacy
    common = gen_keys & leg_keys   # In both
    
    # Print summary
    print("=" * 80)
    print(f"Comparison: {pdb_id}")
    print("=" * 80)
    print(f"Generated atoms: {gen_count}")
    print(f"Legacy atoms:     {leg_count}")
    print(f"Difference:       {gen_count - leg_count}")
    print()
    print(f"Missing atoms (in legacy but not generated): {len(missing)}")
    print(f"Extra atoms (in generated but not legacy):   {len(extra)}")
    print(f"Common atoms:                                 {len(common)}")
    print()
    print("Note: 'residue_seq' is the SEQUENCE NUMBER from PDB file (columns 23-26),")
    print("      NOT the position in the chain (e.g., '28' means sequence number 28,")
    print("      not the 28th residue). Insertion codes are shown if present (e.g., '28A').")
    print()
    
    # Show missing atoms
    if missing and not extra_only:
        print("=" * 80)
        print("MISSING ATOMS (in legacy but not generated):")
        print("=" * 80)
        for i, key in enumerate(sorted(missing)[:50], 1):  # Show first 50
            atom = leg_map[key]
            print(f"{i:3d}. {format_atom_info(atom, show_debug=verbose)}")
            if verbose:
                xyz = atom.get("xyz", [])
                if xyz:
                    print(f"     xyz: [{xyz[0]:.6f}, {xyz[1]:.6f}, {xyz[2]:.6f}]")
        if len(missing) > 50:
            print(f"... and {len(missing) - 50} more missing atoms")
        print()
    
    # Show extra atoms
    if extra and not missing_only:
        print("=" * 80)
        print("EXTRA ATOMS (in generated but not legacy):")
        print("=" * 80)
        for i, key in enumerate(sorted(extra)[:50], 1):  # Show first 50
            atom = gen_map[key]
            print(f"{i:3d}. {format_atom_info(atom, show_debug=verbose)}")
            if verbose:
                xyz = atom.get("xyz", [])
                if xyz:
                    print(f"     xyz: [{xyz[0]:.6f}, {xyz[1]:.6f}, {xyz[2]:.6f}]")
        if len(extra) > 50:
            print(f"... and {len(extra) - 50} more extra atoms")
        print()
    
    # Show field differences in common atoms
    if verbose and common:
        print("=" * 80)
        print("FIELD DIFFERENCES IN COMMON ATOMS (first 20):")
        print("=" * 80)
        diff_count = 0
        for key in sorted(common)[:100]:  # Check first 100
            gen_atom = gen_map[key]
            leg_atom = leg_map[key]
            differences = compare_atoms(gen_atom, leg_atom)
            if differences:
                diff_count += 1
                if diff_count <= 20:
                    print(f"\n{format_atom_info(gen_atom, show_debug=False)}")
                    for diff in differences:
                        print(f"  - {diff}")
        if diff_count == 0:
            print("No field differences found in common atoms")
        print()
    
    # Summary
    print("=" * 80)
    if gen_count == leg_count and len(missing) == 0 and len(extra) == 0:
        print("✓ PERFECT MATCH!")
    else:
        print("✗ MISMATCH DETECTED")
    print("=" * 80)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    pdb_id = sys.argv[1].upper()
    verbose = "--verbose" in sys.argv or "-v" in sys.argv
    missing_only = "--missing-only" in sys.argv or "-m" in sys.argv
    extra_only = "--extra-only" in sys.argv or "-e" in sys.argv
    
    show_differences(pdb_id, verbose, missing_only, extra_only)


if __name__ == "__main__":
    main()

