#!/usr/bin/env python3
"""Compare our C++ validation metrics against DSSR for standard WC pairs."""

import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import sys

sys.path.insert(0, str(Path(__file__).parent))
from lw_classifier import normalize_dssr_nt


def load_modern_validation(json_path: Path) -> dict:
    """Load modern pair_validation JSON."""
    if not json_path.exists():
        return {}
    with open(json_path) as f:
        data = json.load(f)
    result = {}
    for p in data:
        key = tuple(sorted([p["res_id_i"], p["res_id_j"]]))
        result[key] = p
    return result


def load_dssr_pairs(json_path: Path) -> dict:
    """Load DSSR JSON pairs."""
    if not json_path.exists():
        return {}
    with open(json_path) as f:
        data = json.load(f)
    result = {}
    for p in data.get("pairs", []):
        res_id1 = normalize_dssr_nt(p.get("nt1", ""))
        res_id2 = normalize_dssr_nt(p.get("nt2", ""))
        key = tuple(sorted([res_id1, res_id2]))
        result[key] = p
    return result


def process_pdb(args):
    """Process a single PDB and compare validation metrics."""
    pdb_id, modern_dir, dssr_dir = args

    validation_path = modern_dir / "pair_validation" / f"{pdb_id}.json"
    dssr_path = dssr_dir / f"{pdb_id}.json"

    if not validation_path.exists() or not dssr_path.exists():
        return None

    try:
        modern_val = load_modern_validation(validation_path)
        dssr_pairs = load_dssr_pairs(dssr_path)

        results = {
            "pdb_id": pdb_id,
            "comparisons": [],
        }

        for key, dssr_p in dssr_pairs.items():
            if dssr_p.get("LW") != "cWW":
                continue

            bp = dssr_p.get("bp", "")
            if not bp or bp[0] not in "GCAU" or bp[-1] not in "GCAU":
                continue

            seq = bp[0] + bp[-1]
            if seq not in ["GC", "CG", "AU", "UA"]:
                continue

            if key not in modern_val:
                continue

            modern_p = modern_val[key]
            calc = modern_p.get("calculated_values", {})

            # Compare metrics
            modern_dNN = calc.get("dNN", 0)
            dssr_n1n9 = dssr_p.get("N1N9_dist", 0)

            modern_angle = calc.get("plane_angle", 0)
            dssr_angle = dssr_p.get("interBase_angle", 0)

            results["comparisons"].append({
                "key": key,
                "seq": seq,
                "modern_dNN": modern_dNN,
                "dssr_n1n9": dssr_n1n9,
                "dNN_diff": abs(modern_dNN - dssr_n1n9) if dssr_n1n9 > 0 else None,
                "modern_angle": modern_angle,
                "dssr_angle": dssr_angle,
                "angle_diff": abs(modern_angle - dssr_angle) if dssr_angle > 0 else None,
            })

        return results

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        return None


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Compare validation metrics with DSSR")
    parser.add_argument("--modern-dir", type=Path, default=Path("data/json"))
    parser.add_argument("--dssr-dir", type=Path, default=Path("data/json_dssr"))
    parser.add_argument("--pdb-list", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-pdbs", type=int, default=None)

    args = parser.parse_args()

    # Get PDB list
    modern_pdbs = {f.stem for f in (args.modern_dir / "pair_validation").glob("*.json")}
    dssr_pdbs = {f.stem for f in args.dssr_dir.glob("*.json")}
    common_pdbs = sorted(modern_pdbs & dssr_pdbs)

    if args.pdb_list:
        with open(args.pdb_list) as f:
            pdb_list_data = json.load(f)
        if isinstance(pdb_list_data, list):
            allowed = set(pdb_list_data)
        else:
            allowed = set(next(iter(pdb_list_data.values())))
        common_pdbs = [p for p in common_pdbs if p in allowed]

    if args.max_pdbs:
        common_pdbs = common_pdbs[:args.max_pdbs]

    print(f"Comparing {len(common_pdbs)} PDBs...")

    work_args = [(pdb_id, args.modern_dir, args.dssr_dir) for pdb_id in common_pdbs]

    with Pool(args.workers) as pool:
        results = pool.map(process_pdb, work_args)

    # Aggregate
    all_comparisons = []
    for r in results:
        if r is None:
            continue
        all_comparisons.extend(r["comparisons"])

    print(f"\n{'='*70}")
    print(f"Validation Metrics Comparison (Modern C++ vs DSSR)")
    print(f"{'='*70}")
    print(f"Total comparisons: {len(all_comparisons)}")

    if all_comparisons:
        # N1N9 / dNN comparison
        dNN_diffs = [c["dNN_diff"] for c in all_comparisons if c["dNN_diff"] is not None]
        if dNN_diffs:
            print(f"\ndNN (modern) vs N1N9_dist (DSSR):")
            print(f"  Mean diff: {np.mean(dNN_diffs):.4f} Å")
            print(f"  Median diff: {np.median(dNN_diffs):.4f} Å")
            print(f"  Max diff: {np.max(dNN_diffs):.4f} Å")
            print(f"  <0.01Å: {sum(1 for d in dNN_diffs if d < 0.01)} ({100*sum(1 for d in dNN_diffs if d < 0.01)/len(dNN_diffs):.1f}%)")
            print(f"  <0.1Å: {sum(1 for d in dNN_diffs if d < 0.1)} ({100*sum(1 for d in dNN_diffs if d < 0.1)/len(dNN_diffs):.1f}%)")

        # Angle comparison
        angle_diffs = [c["angle_diff"] for c in all_comparisons if c["angle_diff"] is not None]
        if angle_diffs:
            print(f"\nplane_angle (modern) vs interBase_angle (DSSR):")
            print(f"  Mean diff: {np.mean(angle_diffs):.4f}°")
            print(f"  Median diff: {np.median(angle_diffs):.4f}°")
            print(f"  Max diff: {np.max(angle_diffs):.4f}°")
            print(f"  <0.1°: {sum(1 for d in angle_diffs if d < 0.1)} ({100*sum(1 for d in angle_diffs if d < 0.1)/len(angle_diffs):.1f}%)")
            print(f"  <1.0°: {sum(1 for d in angle_diffs if d < 1.0)} ({100*sum(1 for d in angle_diffs if d < 1.0)/len(angle_diffs):.1f}%)")

        # Sample comparisons
        print(f"\nSample comparisons (first 10):")
        for c in all_comparisons[:10]:
            print(f"  {c['key'][0]} - {c['key'][1]} ({c['seq']})")
            print(f"    dNN: modern={c['modern_dNN']:.3f} DSSR={c['dssr_n1n9']:.3f} diff={c['dNN_diff']:.4f} Å")
            print(f"    angle: modern={c['modern_angle']:.2f}° DSSR={c['dssr_angle']:.2f}° diff={c['angle_diff']:.4f}°")


if __name__ == "__main__":
    main()
