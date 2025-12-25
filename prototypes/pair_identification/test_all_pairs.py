#!/usr/bin/env python3
"""Test cWW classification on ALL potential base pairs, not just ones already found.

This tests both precision (are we finding things that shouldn't be cWW?)
and recall (are we missing cWW pairs?).
"""

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set
from collections import defaultdict
import itertools

from cww_validator import (
    Residue, Atom, HBond, ValidationResult,
    parse_pdb_atoms, load_modern_hbonds, validate_cww_pair,
    normalize_bp_type, WC_HBOND_PATTERNS
)


def normalize_dssr_nt(nt: str) -> str:
    """Convert DSSR nt format (A.G1) to res_id format (A-G-1)."""
    parts = nt.split(".")
    if len(parts) != 2:
        return nt
    chain = parts[0]
    base_num = parts[1]

    i = len(base_num) - 1
    if i >= 0 and base_num[i].isalpha() and i > 0 and base_num[i-1].isdigit():
        i -= 1
    while i >= 0 and base_num[i].isdigit():
        i -= 1
    i += 1

    if i <= 0:
        return f"{chain}-{base_num}-0"

    base = base_num[:i]
    num = base_num[i:]
    return f"{chain}-{base}-{num}"


def load_dssr_cww_pairs(dssr_path: Path) -> Set[Tuple[str, str]]:
    """Load cWW pairs from DSSR as set of (res_id1, res_id2) tuples."""
    if not dssr_path.exists():
        return set()

    with open(dssr_path) as f:
        data = json.load(f)

    pairs = set()
    for p in data.get("pairs", []):
        if p.get("LW") == "cWW":
            res_id1 = normalize_dssr_nt(p.get("nt1", ""))
            res_id2 = normalize_dssr_nt(p.get("nt2", ""))
            pairs.add(tuple(sorted([res_id1, res_id2])))

    return pairs


def find_potential_pairs(
    residues: Dict[str, Residue],
    max_distance: float = 15.0,
) -> List[Tuple[str, str, float]]:
    """Find all residue pairs within distance cutoff.

    Returns list of (res_id1, res_id2, distance).
    """
    pairs = []
    res_list = list(residues.values())

    for i, res1 in enumerate(res_list):
        for res2 in res_list[i+1:]:
            # Check N1/N9 distance
            atom1 = res1.n1n9_atom
            atom2 = res2.n1n9_atom

            if atom1 is None or atom2 is None:
                continue

            dist = np.linalg.norm(atom2.coords - atom1.coords)
            if dist <= max_distance:
                pairs.append((res1.res_id, res2.res_id, dist))

    return pairs


def get_bp_type(res1: Residue, res2: Residue) -> str:
    """Get base pair type from two residues."""
    return res1.base_type + res2.base_type


def test_all_pairs(
    pdb_id: str,
    pdb_dir: Path,
    hbond_dir: Path,
    dssr_dir: Path,
    max_distance: float = 12.0,
    verbose: bool = False,
) -> Dict:
    """Test cWW classification on all potential pairs.

    Returns dict with precision, recall, and details.
    """
    # Load PDB
    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        return {"error": "PDB not found"}

    residues = parse_pdb_atoms(pdb_path)
    hbonds = load_modern_hbonds(hbond_dir / f"{pdb_id}.json")
    dssr_cww = load_dssr_cww_pairs(dssr_dir / f"{pdb_id}.json")

    # Find all potential pairs
    potential_pairs = find_potential_pairs(residues, max_distance)

    if verbose:
        print(f"Found {len(potential_pairs)} potential pairs within {max_distance}Å")
        print(f"DSSR cWW pairs: {len(dssr_cww)}")

    # Classify each pair
    true_positives = []
    false_positives = []
    false_negatives = []
    true_negatives = []

    for res_id1, res_id2, dist in potential_pairs:
        res1 = residues[res_id1]
        res2 = residues[res_id2]
        bp_type = get_bp_type(res1, res2)

        # Get H-bonds for this pair
        pair_key = tuple(sorted([res_id1, res_id2]))
        pair_hbonds = hbonds.get(pair_key, [])

        # Validate as cWW
        result = validate_cww_pair(res1, res2, pair_hbonds, bp_type)

        is_dssr_cww = pair_key in dssr_cww
        our_prediction = result.is_valid_cww

        if our_prediction and is_dssr_cww:
            true_positives.append((res_id1, res_id2, bp_type, result))
        elif our_prediction and not is_dssr_cww:
            false_positives.append((res_id1, res_id2, bp_type, result))
        elif not our_prediction and is_dssr_cww:
            false_negatives.append((res_id1, res_id2, bp_type, result))
        else:
            true_negatives.append((res_id1, res_id2, bp_type, result))

    # Calculate metrics
    tp = len(true_positives)
    fp = len(false_positives)
    fn = len(false_negatives)
    tn = len(true_negatives)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    result = {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "true_negatives": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }

    if verbose:
        print(f"\nResults:")
        print(f"  True Positives:  {tp} (correctly identified cWW)")
        print(f"  False Positives: {fp} (incorrectly called cWW)")
        print(f"  False Negatives: {fn} (missed cWW)")
        print(f"  True Negatives:  {tn} (correctly rejected)")
        print(f"\n  Precision: {precision:.1%}")
        print(f"  Recall:    {recall:.1%}")
        print(f"  F1 Score:  {f1:.1%}")

        if false_positives:
            print(f"\nFalse Positives (we called cWW, DSSR didn't):")
            for res1, res2, bp, res in false_positives[:10]:
                print(f"  {res1} - {res2} ({bp}) N1N9={res.n1n9_dist:.2f}Å hbonds={res.found_hbonds}")

        if false_negatives:
            print(f"\nFalse Negatives (DSSR called cWW, we didn't):")
            for res1, res2, bp, res in false_negatives[:10]:
                n1n9_str = f"{res.n1n9_dist:.2f}" if res.n1n9_dist else "N/A"
                print(f"  {res1} - {res2} ({bp}) N1N9={n1n9_str}Å hbonds={res.found_hbonds}/{res.expected_hbonds}")

    return result


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Test cWW classification on all potential pairs"
    )
    parser.add_argument(
        "--pdb-dir", type=Path,
        default=Path("data/pdb"),
        help="Directory with PDB files"
    )
    parser.add_argument(
        "--hbond-dir", type=Path,
        default=Path("data/json/all_hbond_list"),
        help="Directory with H-bond JSON files"
    )
    parser.add_argument(
        "--dssr-dir", type=Path,
        default=Path("data/json_dssr"),
        help="Directory with DSSR JSON files"
    )
    parser.add_argument(
        "--pdb", type=str, default=None,
        help="Specific PDB to analyze"
    )
    parser.add_argument(
        "--max-pdbs", type=int, default=50,
        help="Maximum PDBs to analyze"
    )
    parser.add_argument(
        "--max-distance", type=float, default=12.0,
        help="Maximum N1N9 distance to consider"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    if args.pdb:
        test_all_pairs(
            args.pdb,
            args.pdb_dir,
            args.hbond_dir,
            args.dssr_dir,
            args.max_distance,
            verbose=True
        )
        return

    # Test multiple PDBs
    dssr_files = sorted(args.dssr_dir.glob("*.json"))[:args.max_pdbs]

    total_tp = 0
    total_fp = 0
    total_fn = 0
    total_tn = 0

    for dssr_path in dssr_files:
        pdb_id = dssr_path.stem
        result = test_all_pairs(
            pdb_id,
            args.pdb_dir,
            args.hbond_dir,
            args.dssr_dir,
            args.max_distance,
            verbose=args.verbose
        )

        if "error" not in result:
            total_tp += result["true_positives"]
            total_fp += result["false_positives"]
            total_fn += result["false_negatives"]
            total_tn += result["true_negatives"]

    # Overall metrics
    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n{'='*60}")
    print(f"OVERALL RESULTS ({len(dssr_files)} PDBs)")
    print(f"{'='*60}")
    print(f"\nConfusion Matrix:")
    print(f"  True Positives:  {total_tp}")
    print(f"  False Positives: {total_fp}")
    print(f"  False Negatives: {total_fn}")
    print(f"  True Negatives:  {total_tn}")
    print(f"\nMetrics:")
    print(f"  Precision: {precision:.1%} (of pairs we call cWW, how many are correct?)")
    print(f"  Recall:    {recall:.1%} (of DSSR cWW pairs, how many do we find?)")
    print(f"  F1 Score:  {f1:.1%}")


if __name__ == "__main__":
    main()
