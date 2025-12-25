#!/usr/bin/env python3
"""Parallel cWW classification validation against DSSR."""

import json
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import sys

sys.path.insert(0, str(Path(__file__).parent))
from lw_classifier import LWClassifier, parse_pdb_residues, normalize_dssr_nt
from hbond_scorer import load_modern_hbonds


def process_pdb(args):
    """Process a single PDB and return stats."""
    pdb_id, pdb_dir, hbond_dir, dssr_dir, idealized_dir, exemplar_dir = args

    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        return None

    dssr_path = dssr_dir / f"{pdb_id}.json"
    if not dssr_path.exists():
        return None

    hbond_path = hbond_dir / f"{pdb_id}.json"
    if not hbond_path.exists():
        return None

    try:
        classifier = LWClassifier(idealized_dir, exemplar_dir)
        residues = parse_pdb_residues(pdb_path)
        hbonds = load_modern_hbonds(hbond_path)

        with open(dssr_path) as f:
            dssr = json.load(f)

        stats = {
            "correct": 0,
            "total": 0,
            "by_seq": defaultdict(lambda: {"c": 0, "t": 0}),
            "miss_reasons": defaultdict(int),
        }

        for pair in dssr.get("pairs", []):
            if pair.get("LW") != "cWW":
                continue

            res_id1 = normalize_dssr_nt(pair["nt1"])
            res_id2 = normalize_dssr_nt(pair["nt2"])

            if res_id1 not in residues or res_id2 not in residues:
                continue

            res1 = residues[res_id1]
            res2 = residues[res_id2]
            seq = res1.base_type + res2.base_type
            pair_key = tuple(sorted([res_id1, res_id2]))
            pair_hbonds = hbonds.get(pair_key, [])

            result = classifier.classify(res1, res2, pair_hbonds)

            stats["total"] += 1
            stats["by_seq"][seq]["t"] += 1

            if result.best_lw == "cWW":
                stats["correct"] += 1
                stats["by_seq"][seq]["c"] += 1
            else:
                if result.best_hbond_score == 0:
                    if len(pair_hbonds) == 0:
                        stats["miss_reasons"]["no_hbonds"] += 1
                    else:
                        stats["miss_reasons"]["pattern_mismatch"] += 1
                elif result.best_hbond_score == 1.0:
                    stats["miss_reasons"]["rmsd_prefers_other"] += 1
                else:
                    stats["miss_reasons"]["partial_match"] += 1

        # Convert defaultdicts to regular dicts for pickling
        stats["by_seq"] = {k: dict(v) for k, v in stats["by_seq"].items()}
        stats["miss_reasons"] = dict(stats["miss_reasons"])
        return stats

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        return None


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Parallel cWW validation")
    parser.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    parser.add_argument("--hbond-dir", type=Path, default=Path("data/json/all_hbond_list"))
    parser.add_argument("--dssr-dir", type=Path, default=Path("data/json_dssr"))
    parser.add_argument("--idealized-dir", type=Path, default=Path("basepair-idealized"))
    parser.add_argument("--exemplar-dir", type=Path, default=Path("basepair-exemplars"))
    parser.add_argument("--workers", type=int, default=10, help="Number of parallel workers")
    parser.add_argument("--max-pdbs", type=int, default=None, help="Max PDBs to process")
    parser.add_argument("--pdb-list", type=Path, default=None, help="JSON file with PDB list to use")

    args = parser.parse_args()

    # Find PDBs with both DSSR and H-bond data
    dssr_ids = {f.stem for f in args.dssr_dir.glob("*.json")}
    hbond_ids = {f.stem for f in args.hbond_dir.glob("*.json")}
    common_ids = sorted(dssr_ids & hbond_ids)

    # Filter to specific PDB list if provided
    if args.pdb_list:
        with open(args.pdb_list) as f:
            pdb_list_data = json.load(f)
        # Handle both list and dict formats
        if isinstance(pdb_list_data, list):
            allowed_ids = set(pdb_list_data)
        else:
            # Assume it's a dict with a list value
            allowed_ids = set(next(iter(pdb_list_data.values())))
        common_ids = [p for p in common_ids if p in allowed_ids]

    if args.max_pdbs:
        common_ids = common_ids[:args.max_pdbs]

    print(f"Processing {len(common_ids)} PDBs with {args.workers} workers...")

    # Prepare arguments for parallel processing
    work_args = [
        (pdb_id, args.pdb_dir, args.hbond_dir, args.dssr_dir,
         args.idealized_dir, args.exemplar_dir)
        for pdb_id in common_ids
    ]

    # Process in parallel
    with Pool(args.workers) as pool:
        results = pool.map(process_pdb, work_args)

    # Aggregate results
    total_stats = {
        "correct": 0,
        "total": 0,
        "by_seq": defaultdict(lambda: {"c": 0, "t": 0}),
        "miss_reasons": defaultdict(int),
    }

    for stats in results:
        if stats is None:
            continue
        total_stats["correct"] += stats["correct"]
        total_stats["total"] += stats["total"]
        for seq, counts in stats["by_seq"].items():
            total_stats["by_seq"][seq]["c"] += counts["c"]
            total_stats["by_seq"][seq]["t"] += counts["t"]
        for reason, count in stats["miss_reasons"].items():
            total_stats["miss_reasons"][reason] += count

    # Print results
    acc = 100 * total_stats["correct"] / total_stats["total"] if total_stats["total"] > 0 else 0
    print(f"\ncWW Accuracy: {total_stats['correct']}/{total_stats['total']} ({acc:.1f}%)")

    print(f"\nMiss reasons:")
    for reason, count in sorted(total_stats["miss_reasons"].items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count}")

    print(f"\nBy sequence (standard WC):")
    for seq in ["GC", "CG", "AU", "UA", "GU", "UG"]:
        d = total_stats["by_seq"].get(seq, {"c": 0, "t": 0})
        if d["t"] > 0:
            acc_seq = 100 * d["c"] / d["t"]
            miss = d["t"] - d["c"]
            print(f"  {seq}: {d['c']}/{d['t']} ({acc_seq:.0f}%) - {miss} misses")

    print(f"\nBy sequence (non-standard):")
    for seq in sorted(total_stats["by_seq"].keys()):
        if seq not in ["GC", "CG", "AU", "UA", "GU", "UG"]:
            d = total_stats["by_seq"][seq]
            if d["t"] >= 10:
                acc_seq = 100 * d["c"] / d["t"]
                miss = d["t"] - d["c"]
                print(f"  {seq}: {d['c']}/{d['t']} ({acc_seq:.0f}%) - {miss} misses")


if __name__ == "__main__":
    main()
