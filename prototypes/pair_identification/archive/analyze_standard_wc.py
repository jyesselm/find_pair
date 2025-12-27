#!/usr/bin/env python3
"""Analyze standard WC (GC/CG/AU/UA) classification and parameters."""

import json
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import sys

sys.path.insert(0, str(Path(__file__).parent))
from lw_classifier import LWClassifier, parse_pdb_residues, normalize_dssr_nt
from hbond_scorer import load_modern_hbonds, HBond


def process_pdb(args):
    """Process a single PDB and return detailed analysis."""
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

        results = {
            "pdb_id": pdb_id,
            "misses": [],
            "correct": [],
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

            # Only standard WC pairs
            if seq not in ["GC", "CG", "AU", "UA"]:
                continue

            pair_key = tuple(sorted([res_id1, res_id2]))
            pair_hbonds = hbonds.get(pair_key, [])

            result = classifier.classify(res1, res2, pair_hbonds)

            # DSSR parameters
            dssr_params = {
                "bp": pair.get("bp", ""),
                "name": pair.get("name", ""),
                "Saenger": pair.get("Saenger", ""),
                "DSSR": pair.get("DSSR", ""),
                "hbonds_desc": pair.get("hbonds_desc", ""),
            }

            if result.best_lw == "cWW":
                # Correct - store for parameter comparison
                results["correct"].append({
                    "res1": res_id1,
                    "res2": res_id2,
                    "seq": seq,
                    "rmsd": result.best_rmsd,
                    "hbond_score": result.best_hbond_score,
                    "dssr": dssr_params,
                })
            else:
                # Miss - detailed analysis
                results["misses"].append({
                    "res1": res_id1,
                    "res2": res_id2,
                    "seq": seq,
                    "predicted": result.best_lw,
                    "pred_rmsd": result.best_rmsd,
                    "pred_hbond": result.best_hbond_score,
                    "second": result.second_lw,
                    "second_rmsd": result.second_rmsd,
                    "hbonds_found": [(h.donor_atom, h.acceptor_atom) for h in pair_hbonds],
                    "dssr": dssr_params,
                    "all_scores": result.all_results[:5],
                })

        return results

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        return None


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Analyze standard WC pairs")
    parser.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    parser.add_argument("--hbond-dir", type=Path, default=Path("data/json/all_hbond_list"))
    parser.add_argument("--dssr-dir", type=Path, default=Path("data/json_dssr"))
    parser.add_argument("--idealized-dir", type=Path, default=Path("basepair-idealized"))
    parser.add_argument("--exemplar-dir", type=Path, default=Path("basepair-exemplars"))
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--pdb-list", type=Path, default=None)
    parser.add_argument("--max-pdbs", type=int, default=None)
    parser.add_argument("--show-misses", type=int, default=50, help="Number of misses to show")

    args = parser.parse_args()

    # Find PDBs
    dssr_ids = {f.stem for f in args.dssr_dir.glob("*.json")}
    hbond_ids = {f.stem for f in args.hbond_dir.glob("*.json")}
    common_ids = sorted(dssr_ids & hbond_ids)

    if args.pdb_list:
        with open(args.pdb_list) as f:
            pdb_list_data = json.load(f)
        if isinstance(pdb_list_data, list):
            allowed_ids = set(pdb_list_data)
        else:
            allowed_ids = set(next(iter(pdb_list_data.values())))
        common_ids = [p for p in common_ids if p in allowed_ids]

    if args.max_pdbs:
        common_ids = common_ids[:args.max_pdbs]

    print(f"Processing {len(common_ids)} PDBs with {args.workers} workers...")

    work_args = [
        (pdb_id, args.pdb_dir, args.hbond_dir, args.dssr_dir,
         args.idealized_dir, args.exemplar_dir)
        for pdb_id in common_ids
    ]

    with Pool(args.workers) as pool:
        results = pool.map(process_pdb, work_args)

    # Aggregate
    all_misses = []
    all_correct = []
    miss_by_seq = defaultdict(list)
    miss_by_predicted = defaultdict(list)
    miss_by_reason = defaultdict(list)

    for r in results:
        if r is None:
            continue
        all_misses.extend(r["misses"])
        all_correct.extend(r["correct"])
        for m in r["misses"]:
            miss_by_seq[m["seq"]].append(m)
            miss_by_predicted[m["predicted"]].append(m)
            # Categorize reason
            if m["pred_hbond"] == 0:
                if len(m["hbonds_found"]) == 0:
                    miss_by_reason["no_hbonds"].append(m)
                else:
                    miss_by_reason["pattern_mismatch"].append(m)
            elif m["pred_hbond"] == 1.0:
                miss_by_reason["rmsd_prefers_other"].append(m)
            else:
                miss_by_reason["partial_hbond"].append(m)

    total = len(all_misses) + len(all_correct)
    acc = 100 * len(all_correct) / total if total > 0 else 0

    print(f"\n{'='*70}")
    print(f"Standard WC (GC/CG/AU/UA) Analysis")
    print(f"{'='*70}")
    print(f"\nOverall: {len(all_correct)}/{total} ({acc:.1f}%)")
    print(f"Misses: {len(all_misses)}")

    print(f"\nMisses by sequence:")
    for seq in ["GC", "CG", "AU", "UA"]:
        misses = miss_by_seq.get(seq, [])
        print(f"  {seq}: {len(misses)} misses")

    print(f"\nMisses by predicted class:")
    for pred, misses in sorted(miss_by_predicted.items(), key=lambda x: -len(x[1])):
        print(f"  {pred}: {len(misses)}")

    print(f"\nMisses by reason:")
    for reason, misses in sorted(miss_by_reason.items(), key=lambda x: -len(x[1])):
        print(f"  {reason}: {len(misses)}")

    # Show sample misses
    print(f"\n{'='*70}")
    print(f"Sample Misses (first {args.show_misses}):")
    print(f"{'='*70}")

    for i, m in enumerate(all_misses[:args.show_misses]):
        print(f"\n{i+1}. {m['seq']} predicted as {m['predicted']}")
        print(f"   Residues: {m['res1']} - {m['res2']}")
        print(f"   RMSD: {m['pred_rmsd']:.3f}Å, H-bond score: {m['pred_hbond']:.0%}")
        print(f"   H-bonds found: {m['hbonds_found']}")
        print(f"   DSSR hbonds: {m['dssr']['hbonds_desc']}")
        print(f"   Second best: {m['second']} (RMSD={m['second_rmsd']:.3f}Å)" if m['second_rmsd'] else "")
        if m['all_scores']:
            print(f"   Top scores:")
            for s in m['all_scores'][:3]:
                print(f"     {s['lw_class']:6s} RMSD={s['rmsd']:.3f} hbond={s['hbond_score']:.0%} combined={s['combined_score']:.3f}")

    # Analyze H-bond patterns in misses
    print(f"\n{'='*70}")
    print(f"H-bond Pattern Analysis in Misses:")
    print(f"{'='*70}")

    hbond_patterns = defaultdict(int)
    for m in all_misses:
        pattern = tuple(sorted(m["hbonds_found"]))
        hbond_patterns[pattern] += 1

    print(f"\nMost common H-bond patterns in misses:")
    for pattern, count in sorted(hbond_patterns.items(), key=lambda x: -x[1])[:20]:
        print(f"  {count:4d}: {pattern}")


if __name__ == "__main__":
    main()
