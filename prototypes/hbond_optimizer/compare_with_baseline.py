#!/usr/bin/env python3
"""
Compare H-bond optimizer output with C++ baseline output.

Usage:
    python compare_with_baseline.py 1GID
    python compare_with_baseline.py --test-set 100
"""

import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))

from optimizer import HBondOptimizer, Residue, HBond


def parse_pdb_residues(pdb_path: Path) -> Tuple[Dict[str, Residue], Dict[int, str]]:
    """
    Parse a PDB file and extract nucleotide residues.

    Returns:
        - Dict mapping res_id to Residue objects
        - Dict mapping 1-based index to res_id
    """
    residues = {}
    idx_to_res_id = {}

    # Map 3-letter codes to single letters
    RES_MAP = {
        'ADE': 'A', 'A': 'A', 'DA': 'A', 'RA': 'A',
        'GUA': 'G', 'G': 'G', 'DG': 'G', 'RG': 'G',
        'CYT': 'C', 'C': 'C', 'DC': 'C', 'RC': 'C',
        'URA': 'U', 'U': 'U', 'RU': 'U',
        'THY': 'T', 'T': 'T', 'DT': 'T',
    }

    # Base and backbone atoms for H-bonding
    BASE_ATOMS = {
        # Base N atoms
        'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9',
        # Base O atoms
        'O2', 'O4', 'O6',
        # Ribose atoms
        "O2'", "O3'", "O4'", "O5'",
        # Phosphate atoms
        'O1P', 'O2P', 'OP1', 'OP2',
        # Ring carbons needed for geometry
        'C2', 'C4', 'C5', 'C6', 'C8',
        # Sugar carbons for geometry
        "C1'", "C3'", "C4'", "C5'",
    }

    current_res_key = None
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

            res_key = (chain, res_seq, ins_code)

            if res_key != current_res_key:
                current_res_key = res_key
                res_counter += 1

                base_type = RES_MAP.get(res_name.upper())
                if base_type:
                    # Build res_id like "A-G-103"
                    res_id = f"{chain}-{base_type}-{res_seq}{ins_code}".rstrip()
                    residues[res_id] = Residue(
                        res_id=res_id,
                        base_type=base_type,
                        atoms={},
                        residue_code=res_name
                    )
                    idx_to_res_id[res_counter] = res_id

            # Get current res_id
            if res_counter in idx_to_res_id:
                res_id = idx_to_res_id[res_counter]
                if atom_name in BASE_ATOMS:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    residues[res_id].atoms[atom_name] = np.array([x, y, z])

    return residues, idx_to_res_id


def load_baseline_hbonds(json_path: Path) -> List[dict]:
    """Load baseline H-bond JSON file."""
    with open(json_path) as f:
        return json.load(f)


def normalize_atom_name(atom: str) -> str:
    """Normalize atom names for comparison."""
    atom = atom.strip()
    # Normalize phosphate oxygens
    if atom == 'OP1':
        return 'O1P'
    if atom == 'OP2':
        return 'O2P'
    return atom


def compare_single_pdb(pdb_id: str, base_dir: Path, verbose: bool = False) -> Tuple[int, int, int, List[str]]:
    """
    Compare optimizer output with baseline for a single PDB.

    Returns: (matches, baseline_only, optimizer_only, diff_details)
    """
    pdb_path = base_dir / "pdb" / f"{pdb_id}.pdb"
    json_path = base_dir / "json_baseline" / "hbond_list" / f"{pdb_id}.json"

    if not pdb_path.exists() or not json_path.exists():
        return 0, 0, 0, []

    # Parse PDB
    residues, idx_to_res_id = parse_pdb_residues(pdb_path)

    # Load baseline
    baseline_data = load_baseline_hbonds(json_path)

    # Create optimizer in baseline mode
    optimizer = HBondOptimizer(max_distance=4.0, min_alignment=0.3, baseline_mode=True)
    for res in residues.values():
        optimizer.add_residue(res)

    total_matches = 0
    total_baseline = 0
    total_optimizer = 0
    diff_details = []

    # Process each pair in baseline
    processed_pairs = set()

    for entry in baseline_data:
        res_id_i = entry.get('res_id_i')
        res_id_j = entry.get('res_id_j')

        if not res_id_i or not res_id_j:
            continue

        # Skip if already processed (baseline has both directions)
        pair_key = tuple(sorted([res_id_i, res_id_j]))
        if pair_key in processed_pairs:
            continue
        processed_pairs.add(pair_key)

        if res_id_i not in residues or res_id_j not in residues:
            continue

        # Get baseline "good" H-bonds (type == "-")
        baseline_good = set()
        for hb in entry['hbonds']:
            if hb['type'].strip() == '-':
                # Store as frozenset for direction-agnostic matching (normalize names)
                baseline_good.add(frozenset([
                    (res_id_i, normalize_atom_name(hb['donor_atom'])),
                    (res_id_j, normalize_atom_name(hb['acceptor_atom']))
                ]))

        # Run optimizer
        opt_hbonds = optimizer.optimize_pair(res_id_i, res_id_j)

        # Convert optimizer output (normalize names)
        opt_set = set()
        for hb in opt_hbonds:
            opt_set.add(frozenset([
                (hb.donor_res_id, normalize_atom_name(hb.donor_atom)),
                (hb.acceptor_res_id, normalize_atom_name(hb.acceptor_atom))
            ]))

        matches = baseline_good & opt_set
        baseline_only = baseline_good - opt_set
        opt_only = opt_set - baseline_good

        total_matches += len(matches)
        total_baseline += len(baseline_good)
        total_optimizer += len(opt_hbonds)

        if verbose and (baseline_only or opt_only):
            diff_details.append(f"\n{res_id_i} <-> {res_id_j}")
            if baseline_only:
                for hb in baseline_only:
                    atoms = list(hb)
                    diff_details.append(f"  - Baseline only: {atoms[0][1]}-{atoms[1][1]}")
            if opt_only:
                for hb in opt_only:
                    atoms = list(hb)
                    diff_details.append(f"  + Optimizer only: {atoms[0][1]}-{atoms[1][1]}")

    return total_matches, total_baseline - total_matches, total_optimizer - total_matches, diff_details


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_id', nargs='?', help='PDB ID to compare')
    parser.add_argument('--all', action='store_true', help='Run on all baseline PDBs')
    parser.add_argument('-n', type=int, default=100, help='Limit to first N PDBs (with --all)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show differences')
    args = parser.parse_args()

    base_dir = Path(__file__).parent.parent.parent / "data"

    if args.all:
        # Get all PDBs with baseline hbond files
        hbond_dir = base_dir / "json_baseline" / "hbond_list"
        pdb_ids = [f.stem for f in sorted(hbond_dir.glob("*.json"))]
        if args.n:
            pdb_ids = pdb_ids[:args.n]

        total_matches = 0
        total_baseline = 0
        total_opt = 0

        for pdb_id in pdb_ids:
            matches, bl_only, opt_only, details = compare_single_pdb(
                pdb_id, base_dir, verbose=args.verbose
            )
            total_matches += matches
            total_baseline += matches + bl_only
            total_opt += matches + opt_only

            if args.verbose and details:
                print(f"\n=== {pdb_id} ===")
                for d in details:
                    print(d)

        print(f"\n{'='*60}")
        print(f"SUMMARY ({len(pdb_ids)} PDBs)")
        print(f"{'='*60}")
        print(f"Baseline H-bonds:    {total_baseline}")
        print(f"Optimizer H-bonds:   {total_opt}")
        print(f"Matches:             {total_matches}")
        if total_baseline > 0:
            print(f"Recall:              {total_matches/total_baseline*100:.2f}%")
        if total_opt > 0:
            print(f"Precision:           {total_matches/total_opt*100:.2f}%")

    elif args.pdb_id:
        matches, bl_only, opt_only, details = compare_single_pdb(
            args.pdb_id.upper(), base_dir, verbose=True
        )

        for d in details:
            print(d)

        print(f"\n{'='*60}")
        print(f"SUMMARY")
        print(f"{'='*60}")
        total_bl = matches + bl_only
        total_opt = matches + opt_only
        print(f"Baseline H-bonds:    {total_bl}")
        print(f"Optimizer H-bonds:   {total_opt}")
        print(f"Matches:             {matches}")
        if total_bl > 0:
            print(f"Recall:              {matches/total_bl*100:.1f}%")
        if total_opt > 0:
            print(f"Precision:           {matches/total_opt*100:.1f}%")

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
