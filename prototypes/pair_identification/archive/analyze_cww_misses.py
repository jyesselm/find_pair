#!/usr/bin/env python3
"""Analyze cWW classification misses in detail."""

import json
from pathlib import Path
from collections import defaultdict
from typing import Dict, List
import argparse

from lw_classifier import (
    LWClassifier, parse_pdb_residues, normalize_dssr_nt
)
from hbond_scorer import load_modern_hbonds


def analyze_misses(
    pdb_ids: List[str],
    classifier: LWClassifier,
    pdb_dir: Path,
    hbond_dir: Path,
    dssr_dir: Path,
) -> Dict:
    """Analyze cWW misses and categorize reasons."""

    stats = {
        "total": 0,
        "correct": 0,
        "misses": [],
        "by_sequence": defaultdict(lambda: {"correct": 0, "total": 0}),
        "by_predicted": defaultdict(int),
        "by_reason": defaultdict(int),
    }

    for pdb_id in pdb_ids:
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
        if not pdb_path.exists():
            continue

        dssr_path = dssr_dir / f"{pdb_id}.json"
        if not dssr_path.exists():
            continue

        residues = parse_pdb_residues(pdb_path)
        hbonds = load_modern_hbonds(hbond_dir / f"{pdb_id}.json")

        with open(dssr_path) as f:
            dssr_data = json.load(f)

        for pair in dssr_data.get("pairs", []):
            if pair.get("LW") != "cWW":
                continue

            res_id1 = normalize_dssr_nt(pair.get("nt1", ""))
            res_id2 = normalize_dssr_nt(pair.get("nt2", ""))

            if res_id1 not in residues or res_id2 not in residues:
                continue

            res1 = residues[res_id1]
            res2 = residues[res_id2]
            sequence = res1.base_type + res2.base_type

            pair_key = tuple(sorted([res_id1, res_id2]))
            pair_hbonds = hbonds.get(pair_key, [])

            result = classifier.classify(res1, res2, pair_hbonds)

            stats["total"] += 1
            stats["by_sequence"][sequence]["total"] += 1

            if result.best_lw == "cWW":
                stats["correct"] += 1
                stats["by_sequence"][sequence]["correct"] += 1
            else:
                stats["by_predicted"][result.best_lw] += 1

                # Categorize reason for miss
                reason = "unknown"
                if result.best_hbond_score == 0:
                    if len(pair_hbonds) == 0:
                        reason = "no_hbonds_detected"
                    else:
                        reason = "hbond_pattern_mismatch"
                elif result.best_hbond_score == 1.0:
                    reason = "rmsd_preference"
                else:
                    reason = "partial_hbond_match"

                stats["by_reason"][reason] += 1

                stats["misses"].append({
                    "pdb": pdb_id,
                    "res1": res_id1,
                    "res2": res_id2,
                    "sequence": sequence,
                    "predicted": result.best_lw,
                    "pred_rmsd": result.best_rmsd,
                    "pred_hbond": result.best_hbond_score,
                    "second": result.second_lw,
                    "hbonds_found": len(pair_hbonds),
                    "hbond_atoms": [(h.donor_atom, h.acceptor_atom) for h in pair_hbonds],
                    "dssr_hbonds": pair.get("hbonds_desc", ""),
                    "reason": reason,
                })

    return stats


def main():
    parser = argparse.ArgumentParser(description="Analyze cWW classification misses")
    parser.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    parser.add_argument("--hbond-dir", type=Path, default=Path("data/json/all_hbond_list"))
    parser.add_argument("--dssr-dir", type=Path, default=Path("data/json_dssr"))
    parser.add_argument("--idealized-dir", type=Path, default=Path("basepair-idealized"))
    parser.add_argument("--exemplar-dir", type=Path, default=Path("basepair-exemplars"))
    parser.add_argument("--max-pdbs", type=int, default=500)
    parser.add_argument("--output", type=Path, default=None)

    args = parser.parse_args()

    # Get PDB IDs
    dssr_files = sorted(args.dssr_dir.glob("*.json"))[:args.max_pdbs]
    pdb_ids = [f.stem for f in dssr_files]

    print(f"Analyzing {len(pdb_ids)} PDBs...")

    classifier = LWClassifier(args.idealized_dir, args.exemplar_dir)
    stats = analyze_misses(pdb_ids, classifier, args.pdb_dir, args.hbond_dir, args.dssr_dir)

    # Summary
    acc = stats["correct"] / stats["total"] if stats["total"] > 0 else 0
    print(f"\n{'='*60}")
    print(f"cWW Classification Analysis ({len(pdb_ids)} PDBs)")
    print(f"{'='*60}")
    print(f"\nOverall: {stats['correct']}/{stats['total']} ({acc:.1%})")

    print(f"\nMisses by predicted class:")
    for pred, count in sorted(stats["by_predicted"].items(), key=lambda x: -x[1]):
        print(f"  {pred:6s}: {count:4d}")

    print(f"\nMisses by reason:")
    for reason, count in sorted(stats["by_reason"].items(), key=lambda x: -x[1]):
        print(f"  {reason:25s}: {count:4d}")

    print(f"\nBy sequence (showing problem sequences):")
    for seq, data in sorted(stats["by_sequence"].items(), key=lambda x: x[1]["total"], reverse=True):
        if data["total"] >= 10:
            acc = data["correct"] / data["total"]
            misses = data["total"] - data["correct"]
            if misses > 0:
                print(f"  {seq}: {data['correct']}/{data['total']} ({acc:.0%}) - {misses} misses")

    print(f"\nSample misses (first 30):")
    for m in stats["misses"][:30]:
        print(f"  {m['pdb']:6s} {m['res1']}-{m['res2']} ({m['sequence']}): "
              f"pred={m['predicted']} rmsd={m['pred_rmsd']:.2f} hbond={m['pred_hbond']:.0%} "
              f"reason={m['reason']}")
        if m["hbond_atoms"]:
            print(f"         hbonds: {m['hbond_atoms']}")

    if args.output:
        with open(args.output, 'w') as f:
            json.dump(stats, f, indent=2, default=str)
        print(f"\nFull results saved to {args.output}")


if __name__ == "__main__":
    main()
