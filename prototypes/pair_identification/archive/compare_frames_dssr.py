#!/usr/bin/env python3
"""Compare our C++ frame calculation against DSSR for standard WC pairs."""

import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import sys

sys.path.insert(0, str(Path(__file__).parent))
from lw_classifier import normalize_dssr_nt


def load_modern_pairs(json_path: Path) -> dict:
    """Load modern base_pair JSON."""
    if not json_path.exists():
        return {}
    with open(json_path) as f:
        data = json.load(f)
    # Key by sorted res_id pair
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


def compute_n1n9_distance(res1_atoms: dict, res2_atoms: dict, base1: str, base2: str) -> float:
    """Compute N1-N9 or N9-N1 distance depending on base types."""
    # Purines (A, G) have N9, Pyrimidines (C, U, T) have N1 as glycosidic
    purine = {"A", "G"}

    if base1 in purine:
        atom1 = "N9"
    else:
        atom1 = "N1"

    if base2 in purine:
        atom2 = "N9"
    else:
        atom2 = "N1"

    if atom1 not in res1_atoms or atom2 not in res2_atoms:
        return float('nan')

    c1 = np.array(res1_atoms[atom1])
    c2 = np.array(res2_atoms[atom2])
    return float(np.linalg.norm(c1 - c2))


def process_pdb(args):
    """Process a single PDB and compare frames."""
    pdb_id, modern_dir, dssr_dir = args

    modern_path = modern_dir / "base_pair" / f"{pdb_id}.json"
    dssr_path = dssr_dir / f"{pdb_id}.json"

    if not modern_path.exists() or not dssr_path.exists():
        return None

    try:
        modern_pairs = load_modern_pairs(modern_path)
        dssr_pairs = load_dssr_pairs(dssr_path)

        results = {
            "pdb_id": pdb_id,
            "comparisons": [],
            "missing_in_modern": [],
            "extra_in_modern": [],
        }

        # Compare standard WC cWW pairs
        for key, dssr_p in dssr_pairs.items():
            if dssr_p.get("LW") != "cWW":
                continue

            bp = dssr_p.get("bp", "")
            if not bp or bp[0] not in "GCAU" or bp[-1] not in "GCAU":
                continue

            seq = bp[0] + bp[-1]
            if seq not in ["GC", "CG", "AU", "UA"]:
                continue

            if key not in modern_pairs:
                results["missing_in_modern"].append({
                    "key": key,
                    "dssr": dssr_p,
                })
                continue

            modern_p = modern_pairs[key]

            # Compare frame origins
            dssr_frame = dssr_p.get("frame", {})
            dssr_origin = dssr_frame.get("origin", [0, 0, 0])

            # Modern stores separate origins for each residue
            # Compute midpoint as the pair origin
            modern_org_i = np.array(modern_p["org_i"])
            modern_org_j = np.array(modern_p["org_j"])
            modern_origin = (modern_org_i + modern_org_j) / 2

            origin_diff = np.linalg.norm(np.array(dssr_origin) - modern_origin)

            # Compare N1N9 distance
            dssr_n1n9 = dssr_p.get("N1N9_dist", 0)

            # Compare inter-base angle
            dssr_angle = dssr_p.get("interBase_angle", 0)

            results["comparisons"].append({
                "key": key,
                "seq": seq,
                "dssr_n1n9": dssr_n1n9,
                "dssr_angle": dssr_angle,
                "dssr_origin": dssr_origin,
                "modern_origin": modern_origin.tolist(),
                "origin_diff": origin_diff,
                "dssr_bp_params": dssr_p.get("bp_params", []),
            })

        return results

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        return None


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Compare frames with DSSR")
    parser.add_argument("--modern-dir", type=Path, default=Path("data/json"))
    parser.add_argument("--dssr-dir", type=Path, default=Path("data/json_dssr"))
    parser.add_argument("--pdb-list", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-pdbs", type=int, default=None)

    args = parser.parse_args()

    # Get PDB list
    modern_pdbs = {f.stem for f in (args.modern_dir / "base_pair").glob("*.json")}
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
    all_missing = []

    for r in results:
        if r is None:
            continue
        all_comparisons.extend(r["comparisons"])
        all_missing.extend(r["missing_in_modern"])

    print(f"\n{'='*70}")
    print(f"Frame Comparison Results")
    print(f"{'='*70}")
    print(f"Total comparisons: {len(all_comparisons)}")
    print(f"Missing in modern: {len(all_missing)}")

    if all_comparisons:
        # Analyze origin differences
        origin_diffs = [c["origin_diff"] for c in all_comparisons]
        print(f"\nOrigin difference (modern vs DSSR midpoint):")
        print(f"  Mean: {np.mean(origin_diffs):.3f} Å")
        print(f"  Median: {np.median(origin_diffs):.3f} Å")
        print(f"  Max: {np.max(origin_diffs):.3f} Å")
        print(f"  <0.5Å: {sum(1 for d in origin_diffs if d < 0.5)} ({100*sum(1 for d in origin_diffs if d < 0.5)/len(origin_diffs):.1f}%)")
        print(f"  <1.0Å: {sum(1 for d in origin_diffs if d < 1.0)} ({100*sum(1 for d in origin_diffs if d < 1.0)/len(origin_diffs):.1f}%)")

        # N1N9 distances
        n1n9_values = [c["dssr_n1n9"] for c in all_comparisons if c["dssr_n1n9"] > 0]
        if n1n9_values:
            print(f"\nDSSR N1N9 distances:")
            print(f"  Mean: {np.mean(n1n9_values):.3f} Å")
            print(f"  Range: {np.min(n1n9_values):.3f} - {np.max(n1n9_values):.3f} Å")

        # Inter-base angles
        angles = [c["dssr_angle"] for c in all_comparisons if c["dssr_angle"] > 0]
        if angles:
            print(f"\nDSSR inter-base angles:")
            print(f"  Mean: {np.mean(angles):.1f}°")
            print(f"  Range: {np.min(angles):.1f}° - {np.max(angles):.1f}°")

        # Sample comparisons
        print(f"\nSample comparisons (first 10):")
        for c in all_comparisons[:10]:
            print(f"  {c['key'][0]} - {c['key'][1]} ({c['seq']})")
            print(f"    DSSR origin: {c['dssr_origin']}")
            print(f"    Modern origin: [{c['modern_origin'][0]:.3f}, {c['modern_origin'][1]:.3f}, {c['modern_origin'][2]:.3f}]")
            print(f"    Origin diff: {c['origin_diff']:.3f} Å")
            print(f"    N1N9: {c['dssr_n1n9']:.3f} Å, angle: {c['dssr_angle']:.1f}°")

    # Show missing pairs
    if all_missing:
        print(f"\nSample missing in modern (first 10):")
        for m in all_missing[:10]:
            print(f"  {m['key']}: {m['dssr'].get('bp', '')} {m['dssr'].get('LW', '')}")


if __name__ == "__main__":
    main()
