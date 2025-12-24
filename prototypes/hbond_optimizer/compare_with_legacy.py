#!/usr/bin/env python3
"""
Compare H-bond optimizer output with legacy X3DNA output.

Usage:
    python compare_with_legacy.py 1BNA
    python compare_with_legacy.py data/pdb/1GID.pdb
"""

import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from optimizer import HBondOptimizer, Residue, HBond, format_hbond
from geometry import DONOR_CAPACITY, ACCEPTOR_CAPACITY


def parse_pdb_residues(pdb_path: Path) -> Dict[int, Residue]:
    """
    Parse a PDB file and extract nucleotide residues with their atoms.

    Returns dict mapping residue index (1-based) to Residue objects.
    """
    residues = {}
    current_res_idx = None
    current_chain = None
    current_res_name = None

    # Map 3-letter codes to single letters
    RES_MAP = {
        'ADE': 'A', 'A': 'A', 'DA': 'A', 'RA': 'A',
        'GUA': 'G', 'G': 'G', 'DG': 'G', 'RG': 'G',
        'CYT': 'C', 'C': 'C', 'DC': 'C', 'RC': 'C',
        'URA': 'U', 'U': 'U', 'RU': 'U',
        'THY': 'T', 'T': 'T', 'DT': 'T',
    }

    # Base and ribose atoms - N and O only (no C-H donors)
    BASE_ATOMS = {
        # Base N atoms
        'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9',
        # Base O atoms
        'O2', 'O4', 'O6',
        # Ribose O atoms
        "O2'", "O4'",
        # Ring carbons needed for geometry computation only
        'C2', 'C4', 'C5', 'C6', 'C8',
    }

    res_counter = 0

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_seq = line[22:26].strip()
            ins_code = line[26].strip()

            # Build residue key
            res_key = (chain, res_seq, ins_code)

            # Check if this is a new residue
            if res_key != (current_chain, current_res_idx, ins_code):
                current_chain = chain
                current_res_idx = res_seq
                res_counter += 1

                base_type = RES_MAP.get(res_name.upper())
                if base_type:
                    res_id = f"{chain}-{base_type}-{res_seq}{ins_code}".rstrip()
                    residues[res_counter] = Residue(
                        res_id=res_id,
                        base_type=base_type,
                        atoms={}
                    )

            # Add atom if it's a base atom
            if atom_name in BASE_ATOMS and res_counter in residues:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                residues[res_counter].atoms[atom_name] = np.array([x, y, z])

    return residues


def load_legacy_hbonds(json_path: Path) -> List[dict]:
    """Load legacy H-bond JSON file."""
    with open(json_path) as f:
        return json.load(f)


def compare_hbonds(optimizer: HBondOptimizer,
                   residues: Dict[int, Residue],
                   legacy_data: List[dict]) -> Tuple[int, int, int]:
    """
    Compare optimizer output with legacy H-bonds.

    Returns (matches, optimizer_only, legacy_only)
    """
    # Build a set of unique base pairs from legacy data
    processed_pairs = set()

    total_legacy_good = 0
    total_optimizer = 0
    total_matches = 0

    print("\n" + "="*80)
    print("H-BOND COMPARISON")
    print("="*80)

    for entry in legacy_data:
        base_i = entry['base_i']
        base_j = entry['base_j']

        # Skip if we've already processed this pair (legacy has both directions)
        pair_key = tuple(sorted([base_i, base_j]))
        if pair_key in processed_pairs:
            continue
        processed_pairs.add(pair_key)

        if base_i not in residues or base_j not in residues:
            continue

        res_i = residues[base_i]
        res_j = residues[base_j]

        # Get legacy "good" H-bonds (type == "-")
        legacy_good = []
        for hb in entry['hbonds']:
            if hb['type'].strip() == '-':
                legacy_good.append((hb['donor_atom'].strip(),
                                   hb['acceptor_atom'].strip(),
                                   hb['distance']))

        # Run optimizer
        optimizer.add_residue(res_i)
        optimizer.add_residue(res_j)
        opt_hbonds = optimizer.optimize_pair(res_i.res_id, res_j.res_id)

        # Convert optimizer output to comparable format
        opt_set = set()
        for hb in opt_hbonds:
            # Figure out which residue is donor
            if hb.donor_res_id == res_i.res_id:
                opt_set.add((hb.donor_atom, hb.acceptor_atom))
            else:
                opt_set.add((hb.acceptor_atom, hb.donor_atom))

        legacy_set = set((d, a) for d, a, _ in legacy_good)

        matches = opt_set & legacy_set
        opt_only = opt_set - legacy_set
        leg_only = legacy_set - opt_set

        total_legacy_good += len(legacy_good)
        total_optimizer += len(opt_hbonds)
        total_matches += len(matches)

        # Print comparison
        if legacy_good or opt_hbonds:
            print(f"\n{res_i.res_id} <-> {res_j.res_id}")
            print(f"  Legacy ({len(legacy_good)} good): ", end="")
            for d, a, dist in legacy_good:
                print(f"{d}-{a}({dist:.2f}) ", end="")
            print()

            print(f"  Optimizer ({len(opt_hbonds)}):     ", end="")
            for hb in opt_hbonds:
                d = hb.donor_atom if hb.donor_res_id == res_i.res_id else hb.acceptor_atom
                a = hb.acceptor_atom if hb.donor_res_id == res_i.res_id else hb.donor_atom
                print(f"{d}-{a}({hb.distance:.2f}) ", end="")
            print()

            if matches:
                print(f"  ✓ Matches: {matches}")
            if opt_only:
                print(f"  + Optimizer only: {opt_only}")
            if leg_only:
                print(f"  - Legacy only: {leg_only}")

    return total_matches, total_optimizer, total_legacy_good


def main():
    if len(sys.argv) < 2:
        print("Usage: python compare_with_legacy.py <PDB_ID or path>")
        sys.exit(1)

    pdb_arg = sys.argv[1]

    # Find the PDB and JSON files
    base_dir = Path(__file__).parent.parent.parent
    data_dir = base_dir / "data"

    if Path(pdb_arg).exists():
        pdb_path = Path(pdb_arg)
        pdb_id = pdb_path.stem
    else:
        pdb_id = pdb_arg.upper()
        pdb_path = data_dir / "pdb" / f"{pdb_id}.pdb"

    json_path = data_dir / "json_legacy" / "hbond_list" / f"{pdb_id}.json"

    if not pdb_path.exists():
        print(f"PDB file not found: {pdb_path}")
        sys.exit(1)

    if not json_path.exists():
        print(f"Legacy JSON not found: {json_path}")
        sys.exit(1)

    print(f"PDB: {pdb_path}")
    print(f"Legacy JSON: {json_path}")

    # Parse PDB
    residues = parse_pdb_residues(pdb_path)
    print(f"Loaded {len(residues)} nucleotide residues")

    # Load legacy data
    legacy_data = load_legacy_hbonds(json_path)
    print(f"Loaded {len(legacy_data)} legacy H-bond entries")

    # Create optimizer
    optimizer = HBondOptimizer(max_distance=3.5, min_alignment=0.3)

    # Compare
    matches, opt_total, leg_total = compare_hbonds(optimizer, residues, legacy_data)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Legacy 'good' H-bonds:    {leg_total}")
    print(f"Optimizer selected:       {opt_total}")
    print(f"Matching:                 {matches}")
    if leg_total > 0:
        print(f"Recall (legacy captured): {matches/leg_total*100:.1f}%")
    if opt_total > 0:
        print(f"Precision (opt correct):  {matches/opt_total*100:.1f}%")


if __name__ == '__main__':
    main()
