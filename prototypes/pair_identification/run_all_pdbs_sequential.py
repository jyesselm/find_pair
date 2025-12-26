#!/usr/bin/env python3
"""Run cWW pair finder on all PDBs sequentially (parallel crashes)."""

import json
import sys
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from find_cww_pairs import find_cww_pairs_single


def main():
    project_root = Path(__file__).parent.parent.parent

    pdb_dir = project_root / "data" / "pdb"
    dssr_dir = project_root / "data" / "json_dssr"
    hbond_dir = project_root / "data" / "json" / "slot_hbonds"
    template_dir = project_root / "basepair-idealized"

    # Get all PDBs with slot H-bond files
    hbond_files = list(hbond_dir.glob("*.json"))
    pdb_ids = [f.stem for f in hbond_files]

    print(f"Running cWW finder on {len(pdb_ids)} PDBs...")

    results = []
    errors = []

    for pdb_id in tqdm(pdb_ids, desc="Processing"):
        try:
            result = find_cww_pairs_single(
                pdb_id, pdb_dir, dssr_dir, hbond_dir, template_dir
            )
            results.append(result)
        except Exception as e:
            errors.append({"pdb_id": pdb_id, "error": str(e)})

    # Aggregate statistics
    total_dssr = sum(r.dssr_pairs for r in results)
    total_tp = sum(r.true_positives for r in results)
    total_fp = sum(r.false_positives for r in results)
    total_fn = sum(r.false_negatives for r in results)
    total_found = sum(r.final_pairs for r in results)

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / total_dssr if total_dssr > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n{'='*60}")
    print(f"AGGREGATE RESULTS ({len(results)} PDBs, {len(errors)} errors)")
    print(f"{'='*60}")
    print(f"DSSR cWW pairs: {total_dssr}")
    print(f"Our pairs found: {total_found}")
    print(f"True positives: {total_tp}")
    print(f"False positives: {total_fp}")
    print(f"False negatives: {total_fn}")
    print(f"\nPrecision: {precision:.4f} ({precision:.2%})")
    print(f"Recall: {recall:.4f} ({recall:.2%})")
    print(f"F1 Score: {f1:.4f} ({f1:.2%})")

    # Save results
    output_data = {
        "summary": {
            "pdb_count": len(results),
            "dssr_pairs": total_dssr,
            "found_pairs": total_found,
            "true_positives": total_tp,
            "false_positives": total_fp,
            "false_negatives": total_fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        },
        "errors": errors,
    }

    output_path = Path("/tmp/find_cww_results_tuned.json")
    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
