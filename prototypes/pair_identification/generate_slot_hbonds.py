#!/usr/bin/env python3
"""Generate H-bond JSON files using slot-based optimizer."""

import json
import sys
import numpy as np
from pathlib import Path
from multiprocessing import Pool
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent / "hbond_optimizer"))
sys.path.insert(0, str(Path(__file__).parent))

from optimizer import HBondOptimizer, Residue as OptResidue
from template_aligner import parse_pdb_residues


def process_pdb(args) -> Tuple[str, bool]:
    """Process a single PDB and generate slot-based H-bonds."""
    pdb_id, pdb_dir, output_dir = args

    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        return pdb_id, False

    try:
        # Parse residues
        residues = parse_pdb_residues(pdb_path)

        # Create optimizer and add all residues
        optimizer = HBondOptimizer(max_distance=4.0)

        opt_residues = {}
        for res_id, res in residues.items():
            opt_res = OptResidue(
                res_id=res_id,
                base_type=res.base_type,
                atoms={name: atom.coords for name, atom in res.atoms.items()}
            )
            optimizer.add_residue(opt_res)
            opt_residues[res_id] = opt_res

        # Find H-bonds for all residue pairs
        res_ids = list(residues.keys())
        all_hbonds = []

        for i, res_id1 in enumerate(res_ids):
            for res_id2 in res_ids[i+1:]:
                # Quick distance check - skip if residues too far apart
                res1 = residues[res_id1]
                res2 = residues[res_id2]

                # Use C1' or first available atom for distance check
                pos1 = None
                pos2 = None
                for atom_name in ["C1'", "N1", "N9", "C2"]:
                    if atom_name in res1.atoms and pos1 is None:
                        pos1 = res1.atoms[atom_name].coords
                    if atom_name in res2.atoms and pos2 is None:
                        pos2 = res2.atoms[atom_name].coords

                if pos1 is None or pos2 is None:
                    continue

                if np.linalg.norm(pos1 - pos2) > 15.0:  # Skip distant pairs
                    continue

                # Find H-bonds using slot optimizer
                hbonds = optimizer.optimize_pair(res_id1, res_id2)

                if hbonds:
                    # Group by base-base vs other contexts
                    base_atoms = {'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6'}

                    hbond_list = []
                    for hb in hbonds:
                        # Determine context
                        donor_is_base = hb.donor_atom in base_atoms
                        acceptor_is_base = hb.acceptor_atom in base_atoms

                        if donor_is_base and acceptor_is_base:
                            context = "base_base"
                        elif donor_is_base or acceptor_is_base:
                            context = "base_sugar"
                        else:
                            context = "sugar_sugar"

                        hbond_list.append({
                            "donor_atom": hb.donor_atom,
                            "acceptor_atom": hb.acceptor_atom,
                            "distance": round(hb.distance, 2),
                            "context": context,
                            "h_slot": hb.h_slot_idx,
                            "lp_slot": hb.lp_slot_idx,
                            "alignment": round(hb.alignment_score, 2)
                        })

                    all_hbonds.append({
                        "res_id_i": res_id1,
                        "res_id_j": res_id2,
                        "hbonds": hbond_list
                    })

        # Write output
        output_path = output_dir / f"{pdb_id}.json"
        with open(output_path, 'w') as f:
            json.dump(all_hbonds, f, indent=2)

        return pdb_id, True

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        return pdb_id, False


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate slot-based H-bond JSON files")
    parser.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/json/slot_hbonds"))
    parser.add_argument("--pdb-list", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-pdbs", type=int, default=None)

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Get PDB list
    if args.pdb_list:
        with open(args.pdb_list) as f:
            pdb_data = json.load(f)
        if isinstance(pdb_data, list):
            pdb_ids = pdb_data
        else:
            pdb_ids = list(next(iter(pdb_data.values())))
    else:
        pdb_ids = sorted([f.stem for f in args.pdb_dir.glob("*.pdb")])

    if args.max_pdbs:
        pdb_ids = pdb_ids[:args.max_pdbs]

    print(f"Generating slot-based H-bonds for {len(pdb_ids)} PDBs with {args.workers} workers...")

    work_args = [(pdb_id, args.pdb_dir, args.output_dir) for pdb_id in pdb_ids]

    with Pool(args.workers) as pool:
        results = pool.map(process_pdb, work_args)

    success = sum(1 for _, ok in results if ok)
    print(f"\nCompleted: {success}/{len(pdb_ids)} PDBs")


if __name__ == "__main__":
    main()
