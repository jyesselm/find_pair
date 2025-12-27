#!/usr/bin/env python3
"""Calculate adjusted metrics excluding non-canonical and DSSR-questionable pairs.

This script:
1. Only considers canonical WC sequences (GC, CG, AU, UA, GU, UG)
2. Excludes pairs marked as "dssr_questionable" in our annotations
3. Provides breakdown of remaining false negatives by reason
"""

import json
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple

sys.path.insert(0, str(Path(__file__).parent))

from find_cww_pairs import (
    find_cww_pairs_single, load_dssr_cww_pairs, normalize_res_id_base,
    normalize_dssr_nt, CANONICAL_WC_SEQUENCES
)

# Modified residue registry
sys.path.insert(0, str(Path(__file__).parent.parent / "hbond_optimizer"))
from modified_registry import get_parent_base


def get_dssr_canonical_pairs(dssr_path: Path) -> Tuple[Set[Tuple[str, str]], int]:
    """Load only canonical WC cWW pairs from DSSR.

    Returns (set of canonical pairs, count of non-canonical pairs skipped)
    """
    if not dssr_path.exists():
        return set(), 0

    with open(dssr_path) as f:
        data = json.load(f)

    canonical = set()
    non_canonical_count = 0

    for p in data.get("pairs", []):
        if p.get("LW") != "cWW":
            continue

        nt1 = p.get("nt1", "")
        nt2 = p.get("nt2", "")

        res_id1 = normalize_dssr_nt(nt1)
        res_id2 = normalize_dssr_nt(nt2)

        # Normalize to parent bases
        norm1 = normalize_res_id_base(res_id1)
        norm2 = normalize_res_id_base(res_id2)

        # Get sequence
        parts1 = norm1.split("-")
        parts2 = norm2.split("-")
        base1 = parts1[1] if len(parts1) >= 2 else "?"
        base2 = parts2[1] if len(parts2) >= 2 else "?"
        sequence = base1 + base2

        if sequence in CANONICAL_WC_SEQUENCES:
            canonical.add(tuple(sorted([norm1, norm2])))
        else:
            non_canonical_count += 1

    return canonical, non_canonical_count


def load_annotations(annotation_dir: Path) -> Dict[str, Dict]:
    """Load all annotation files from a directory."""
    annotations = {}
    for f in annotation_dir.glob("*.json"):
        if f.name == "aggregate.json":
            continue
        with open(f) as fp:
            annotations[f.stem] = json.load(fp)
    return annotations


def is_dssr_questionable(fn_entry: dict) -> bool:
    """Check if a false negative is marked as dssr_questionable."""
    reasons = fn_entry.get("reasons", [])
    return "dssr_questionable" in reasons


def main():
    import argparse

    parser = argparse.ArgumentParser()

    project_root = Path(__file__).parent.parent.parent

    parser.add_argument("--pdb-dir", type=Path, default=project_root / "data" / "pdb")
    parser.add_argument("--dssr-dir", type=Path, default=project_root / "data" / "json_dssr")
    parser.add_argument("--hbond-dir", type=Path, default=project_root / "data" / "json" / "slot_hbonds")
    parser.add_argument("--template-dir", type=Path, default=project_root / "basepair-idealized")
    parser.add_argument("--annotation-dir", type=Path,
                        default=project_root / "prototypes" / "pair_identification" / "analysis_results" / "cww_analysis_100")
    parser.add_argument("--sample", type=int, help="Limit to N PDBs")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--output", "-o", type=Path)

    args = parser.parse_args()

    # Get PDBs with slot H-bond files
    hbond_files = list(args.hbond_dir.glob("*.json"))
    pdb_ids = [f.stem for f in hbond_files]

    if args.sample:
        pdb_ids = pdb_ids[:args.sample]

    # Load annotations if available
    annotations = {}
    if args.annotation_dir.exists():
        annotations = load_annotations(args.annotation_dir)
        print(f"Loaded annotations for {len(annotations)} PDBs")

    print(f"Processing {len(pdb_ids)} PDBs...")

    # Aggregated stats
    total_dssr_canonical = 0
    total_dssr_non_canonical = 0
    total_dssr_questionable = 0
    total_our_pairs = 0
    total_tp = 0
    total_fp = 0
    total_fn_real = 0  # FN after excluding questionable

    # Reason breakdown for real FN
    fn_reasons = defaultdict(int)
    fn_examples = defaultdict(list)

    from tqdm import tqdm

    for pdb_id in tqdm(pdb_ids, desc="Processing"):
        # Get canonical DSSR pairs
        dssr_path = args.dssr_dir / f"{pdb_id}.json"
        dssr_canonical, non_canonical_count = get_dssr_canonical_pairs(dssr_path)

        total_dssr_canonical += len(dssr_canonical)
        total_dssr_non_canonical += non_canonical_count

        if not dssr_canonical:
            continue

        # Run our finder
        result = find_cww_pairs_single(
            pdb_id, args.pdb_dir, args.dssr_dir, args.hbond_dir, args.template_dir
        )

        our_pairs = {
            tuple(sorted([normalize_res_id_base(p.res_id1), normalize_res_id_base(p.res_id2)]))
            for p in result.pairs
        }

        total_our_pairs += len(our_pairs)

        # Calculate TP, FP, FN
        tp = our_pairs & dssr_canonical
        fp = our_pairs - dssr_canonical
        fn_all = dssr_canonical - our_pairs

        total_tp += len(tp)
        total_fp += len(fp)

        # Check annotations for questionable pairs
        if pdb_id in annotations:
            ann = annotations[pdb_id]
            for fn_entry in ann.get("false_negatives", []):
                key = tuple(sorted([
                    normalize_res_id_base(fn_entry["res_id1"]),
                    normalize_res_id_base(fn_entry["res_id2"])
                ]))

                if key not in fn_all:
                    continue  # Already found or not in canonical set

                if is_dssr_questionable(fn_entry):
                    total_dssr_questionable += 1
                else:
                    # This is a real miss
                    total_fn_real += 1

                    # Categorize by reason
                    reasons = fn_entry.get("reasons", [])

                    # Primary reason categorization
                    if "no_hbonds" in reasons:
                        fn_reasons["no_hbonds"] += 1
                    elif "missing_hbonds" in reasons and "wrong_atoms" in reasons:
                        fn_reasons["hbond_pattern_mismatch"] += 1
                    elif "distance_issues" in reasons:
                        fn_reasons["distance_issues"] += 1
                    elif "geometric_outlier" in reasons:
                        fn_reasons["geometric_outlier"] += 1
                    elif "non_canonical" in reasons:
                        fn_reasons["non_canonical_saenger"] += 1
                    else:
                        fn_reasons["other"] += 1
                        if args.verbose and len(fn_examples["other"]) < 5:
                            fn_examples["other"].append({
                                "pdb": pdb_id,
                                "pair": key,
                                "reasons": reasons
                            })
        else:
            # No annotation - count all FN as real
            total_fn_real += len(fn_all)

    # Calculate metrics
    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0

    # Adjusted recall: exclude questionable pairs from denominator
    adjusted_dssr = total_dssr_canonical - total_dssr_questionable
    recall = total_tp / adjusted_dssr if adjusted_dssr > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n{'='*70}")
    print("ADJUSTED METRICS (canonical only, excluding DSSR-questionable)")
    print(f"{'='*70}")

    print(f"\n--- DSSR PAIRS ---")
    print(f"Total DSSR cWW pairs: {total_dssr_canonical + total_dssr_non_canonical}")
    print(f"  Non-canonical (AG, UU, etc.): {total_dssr_non_canonical} (excluded)")
    print(f"  Canonical (GC, AU, GU, etc.): {total_dssr_canonical}")
    print(f"  DSSR-questionable: {total_dssr_questionable} (excluded)")
    print(f"  Adjusted canonical: {adjusted_dssr}")

    print(f"\n--- OUR DETECTION ---")
    print(f"Our pairs found: {total_our_pairs}")
    print(f"True positives: {total_tp}")
    print(f"False positives: {total_fp}")
    print(f"False negatives (real): {total_fn_real}")

    print(f"\n--- ADJUSTED METRICS ---")
    print(f"Precision: {precision:.4f} ({precision:.2%})")
    print(f"Recall: {recall:.4f} ({recall:.2%})")
    print(f"F1 Score: {f1:.4f} ({f1:.2%})")

    if fn_reasons:
        print(f"\n--- FALSE NEGATIVE BREAKDOWN ---")
        for reason, count in sorted(fn_reasons.items(), key=lambda x: -x[1]):
            pct = 100 * count / total_fn_real if total_fn_real > 0 else 0
            print(f"  {reason}: {count} ({pct:.1f}%)")

    if args.verbose and fn_examples:
        print(f"\n--- EXAMPLES ---")
        for reason, examples in fn_examples.items():
            print(f"\n{reason}:")
            for ex in examples:
                print(f"  {ex['pdb']}: {ex['pair']} - {ex['reasons']}")

    if args.output:
        output_data = {
            "summary": {
                "pdb_count": len(pdb_ids),
                "dssr_total": total_dssr_canonical + total_dssr_non_canonical,
                "dssr_non_canonical_excluded": total_dssr_non_canonical,
                "dssr_canonical": total_dssr_canonical,
                "dssr_questionable_excluded": total_dssr_questionable,
                "dssr_adjusted": adjusted_dssr,
                "our_pairs": total_our_pairs,
                "true_positives": total_tp,
                "false_positives": total_fp,
                "false_negatives_real": total_fn_real,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            },
            "fn_breakdown": dict(fn_reasons),
        }
        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
