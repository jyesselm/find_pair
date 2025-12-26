#!/usr/bin/env python3
"""Analyze validation failures to understand threshold tuning needs.

Breakdown validation_failed and n1n9_out_of_range cases to see:
1. What are the actual N1N9 distances for out-of-range pairs?
2. What are the scores/hbond counts for validation failures?
"""

import json
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from find_cww_pairs import (
    parse_pdb, find_initial_candidates, filter_by_n1n9,
    load_slot_hbonds, load_dssr_cww_pairs,
    normalize_res_id_base, TemplateAligner, MAX_C1_DISTANCE, N1N9_RANGE_CWW,
    CANONICAL_WC_SEQUENCES, CandidatePair, get_n1n9_atom, compute_interbase_angle,
    count_wc_hbonds, normalize_bp_type
)


def analyze_validation_details(
    pdb_id: str,
    pdb_dir: Path,
    dssr_dir: Path,
    hbond_dir: Path,
    template_dir: Path,
) -> dict:
    """Get detailed info on why DSSR pairs failed validation."""
    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"

    if not pdb_path.exists():
        return {"error": f"PDB not found: {pdb_id}"}

    # Load DSSR pairs
    dssr_path = dssr_dir / f"{pdb_id}.json"
    dssr_pairs_raw = load_dssr_cww_pairs(dssr_path)

    if not dssr_pairs_raw:
        return {"error": "No DSSR pairs found"}

    # Parse PDB
    residues = parse_pdb(pdb_path)

    # Load H-bonds
    hbond_path = hbond_dir / f"{pdb_id}.json"
    slot_hbonds = load_slot_hbonds(hbond_path)

    # Create aligner
    aligner = TemplateAligner(idealized_dir=template_dir)

    # Find all candidates (non-filtered by canonical)
    initial = find_initial_candidates(residues, MAX_C1_DISTANCE)

    # Build lookup of N1N9 distances for all candidates
    n1n9_distances = {}
    for res_id1, res_id2, c1_dist in initial:
        res1 = residues.get(res_id1)
        res2 = residues.get(res_id2)

        if res1 is None or res2 is None:
            continue

        n1 = get_n1n9_atom(res1)
        n2 = get_n1n9_atom(res2)

        if n1 is None or n2 is None:
            continue

        n1n9_dist = np.linalg.norm(n2 - n1)
        key = tuple(sorted([normalize_res_id_base(res_id1), normalize_res_id_base(res_id2)]))
        n1n9_distances[key] = {
            "n1n9": n1n9_dist,
            "res_id1": res_id1,
            "res_id2": res_id2,
            "base1": res1.base_type,
            "base2": res2.base_type,
        }

    results = {
        "pdb_id": pdb_id,
        "n1n9_short": [],  # N1N9 < 8.0
        "n1n9_long": [],   # N1N9 > 10.0
        "validation_low_score": [],  # Score < 0.6
        "validation_no_hbonds": [],  # No H-bonds
        "ok": [],  # Should be captured
        "missing_atoms": [],  # Can't compute N1N9
    }

    for r1, r2 in dssr_pairs_raw:
        key = tuple(sorted([normalize_res_id_base(r1), normalize_res_id_base(r2)]))

        # Check sequence
        parts1 = key[0].split("-")
        parts2 = key[1].split("-")
        base1 = parts1[1] if len(parts1) >= 2 else "?"
        base2 = parts2[1] if len(parts2) >= 2 else "?"
        sequence = base1 + base2

        if sequence not in CANONICAL_WC_SEQUENCES:
            continue  # Skip non-canonical

        if key not in n1n9_distances:
            results["missing_atoms"].append({
                "key": key,
                "sequence": sequence,
            })
            continue

        info = n1n9_distances[key]
        n1n9 = info["n1n9"]

        # Check N1N9 range
        if n1n9 < N1N9_RANGE_CWW[0]:
            results["n1n9_short"].append({
                "key": key,
                "n1n9": n1n9,
                "sequence": sequence,
                "deficit": N1N9_RANGE_CWW[0] - n1n9,
            })
            continue

        if n1n9 > N1N9_RANGE_CWW[1]:
            results["n1n9_long"].append({
                "key": key,
                "n1n9": n1n9,
                "sequence": sequence,
                "excess": n1n9 - N1N9_RANGE_CWW[1],
            })
            continue

        # In N1N9 range - check validation
        # Get H-bonds
        hb_key = (info["res_id1"], info["res_id2"])
        hbonds = slot_hbonds.get(hb_key, [])
        if not hbonds:
            # Try reverse
            hb_key = (info["res_id2"], info["res_id1"])
            hbonds = slot_hbonds.get(hb_key, [])

        expected, found = count_wc_hbonds(sequence, hbonds)

        # Get RMSD
        res1 = residues.get(info["res_id1"])
        res2 = residues.get(info["res_id2"])

        rmsd = None
        if res1 and res2:
            result = aligner.classify_pair(res1, res2, lw_classes=["cWW"])
            if result.best_lw == "cWW":
                rmsd = result.best_rmsd

        # Compute angle
        angle = compute_interbase_angle(res1, res2) if res1 and res2 else 90.0

        # Compute score
        hbond_score = min(found / max(expected, 1), 1.0)
        if rmsd is not None:
            if rmsd <= 0.5:
                rmsd_score = 1.0
            elif rmsd >= 1.5:
                rmsd_score = 0.0
            else:
                rmsd_score = 1.0 - (rmsd - 0.5) / 1.0
        else:
            rmsd_score = 0.5

        if angle <= 15:
            angle_score = 1.0
        elif angle >= 30:
            angle_score = 0.0
        else:
            angle_score = 1.0 - (angle - 15) / 15

        score = 0.4 * hbond_score + 0.3 * rmsd_score + 0.3 * angle_score

        entry = {
            "key": key,
            "sequence": sequence,
            "n1n9": n1n9,
            "hbonds_found": found,
            "hbonds_expected": expected,
            "rmsd": rmsd,
            "angle": angle,
            "score": score,
        }

        # Categorize
        if found == 0:
            results["validation_no_hbonds"].append(entry)
        elif score < 0.6:
            results["validation_low_score"].append(entry)
        else:
            results["ok"].append(entry)

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser()

    project_root = Path(__file__).parent.parent.parent

    parser.add_argument("--pdb-dir", type=Path, default=project_root / "data" / "pdb")
    parser.add_argument("--dssr-dir", type=Path, default=project_root / "data" / "json_dssr")
    parser.add_argument("--hbond-dir", type=Path, default=project_root / "data" / "json" / "slot_hbonds")
    parser.add_argument("--template-dir", type=Path, default=project_root / "basepair-idealized")
    parser.add_argument("--sample", type=int, default=500)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--output", "-o", type=Path)

    args = parser.parse_args()

    # Get PDBs
    hbond_files = list(args.hbond_dir.glob("*.json"))
    pdb_ids = [f.stem for f in hbond_files[:args.sample]]

    print(f"Analyzing {len(pdb_ids)} PDBs...")

    all_n1n9_short = []
    all_n1n9_long = []
    all_low_score = []
    all_no_hbonds = []
    all_ok = []
    all_missing = []

    for i, pdb_id in enumerate(pdb_ids):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(pdb_ids)}...")

        results = analyze_validation_details(
            pdb_id, args.pdb_dir, args.dssr_dir, args.hbond_dir, args.template_dir
        )

        if "error" in results:
            continue

        all_n1n9_short.extend(results["n1n9_short"])
        all_n1n9_long.extend(results["n1n9_long"])
        all_low_score.extend(results["validation_low_score"])
        all_no_hbonds.extend(results["validation_no_hbonds"])
        all_ok.extend(results["ok"])
        all_missing.extend(results["missing_atoms"])

    print(f"\n{'='*60}")
    print("DETAILED VALIDATION ANALYSIS (canonical WC only)")
    print(f"{'='*60}")

    total_canonical = len(all_n1n9_short) + len(all_n1n9_long) + len(all_low_score) + len(all_no_hbonds) + len(all_ok) + len(all_missing)
    print(f"\nTotal canonical DSSR cWW pairs: {total_canonical}")

    print(f"\n--- N1N9 DISTANCE ISSUES ---")
    print(f"N1N9 too short (<8.0Å): {len(all_n1n9_short)}")
    if all_n1n9_short:
        distances = [x["n1n9"] for x in all_n1n9_short]
        print(f"  Range: {min(distances):.2f} - {max(distances):.2f}Å")
        print(f"  Mean: {np.mean(distances):.2f}Å")
        # How many would be captured at different thresholds?
        for thresh in [7.5, 7.0, 6.5, 6.0]:
            captured = sum(1 for d in distances if d >= thresh)
            print(f"  Would capture at {thresh}Å threshold: {captured} ({100*captured/len(distances):.1f}%)")

    print(f"\nN1N9 too long (>10.0Å): {len(all_n1n9_long)}")
    if all_n1n9_long:
        distances = [x["n1n9"] for x in all_n1n9_long]
        print(f"  Range: {min(distances):.2f} - {max(distances):.2f}Å")
        print(f"  Mean: {np.mean(distances):.2f}Å")
        for thresh in [10.5, 11.0, 11.5, 12.0]:
            captured = sum(1 for d in distances if d <= thresh)
            print(f"  Would capture at {thresh}Å threshold: {captured} ({100*captured/len(distances):.1f}%)")

    print(f"\n--- VALIDATION ISSUES ---")
    print(f"No H-bonds found: {len(all_no_hbonds)}")
    if all_no_hbonds and args.verbose:
        print("  Sample:")
        for x in all_no_hbonds[:5]:
            print(f"    {x['key']} ({x['sequence']}) n1n9={x['n1n9']:.2f} rmsd={x['rmsd']}")

    print(f"\nLow score (<0.6): {len(all_low_score)}")
    if all_low_score:
        scores = [x["score"] for x in all_low_score]
        hbonds = [x["hbonds_found"] for x in all_low_score]
        print(f"  Score range: {min(scores):.3f} - {max(scores):.3f}")
        print(f"  Mean score: {np.mean(scores):.3f}")
        print(f"  H-bonds: {min(hbonds)} - {max(hbonds)}, mean={np.mean(hbonds):.1f}")

        # What's causing low scores?
        low_hbond = sum(1 for x in all_low_score if x["hbonds_found"] < 2)
        low_rmsd = sum(1 for x in all_low_score if x["rmsd"] and x["rmsd"] > 1.0)
        high_angle = sum(1 for x in all_low_score if x["angle"] > 25)
        print(f"  Low H-bonds (<2): {low_hbond}")
        print(f"  High RMSD (>1.0): {low_rmsd}")
        print(f"  High angle (>25°): {high_angle}")

        # How many would pass with lower threshold?
        for thresh in [0.5, 0.4, 0.3]:
            passing = sum(1 for s in scores if s >= thresh)
            print(f"  Would pass at score≥{thresh}: {passing}")

    print(f"\nMissing atoms (can't compute N1N9): {len(all_missing)}")

    print(f"\n--- CORRECTLY CAPTURED ---")
    print(f"OK (validated & in range): {len(all_ok)}")

    # Summary
    total_issues = len(all_n1n9_short) + len(all_n1n9_long) + len(all_low_score) + len(all_no_hbonds) + len(all_missing)
    print(f"\n--- SUMMARY ---")
    print(f"Total issues: {total_issues} ({100*total_issues/total_canonical:.1f}% of canonical)")
    print(f"  N1N9 range: {len(all_n1n9_short) + len(all_n1n9_long)} ({100*(len(all_n1n9_short) + len(all_n1n9_long))/total_canonical:.1f}%)")
    print(f"  Validation: {len(all_low_score) + len(all_no_hbonds)} ({100*(len(all_low_score) + len(all_no_hbonds))/total_canonical:.1f}%)")
    print(f"  Missing atoms: {len(all_missing)} ({100*len(all_missing)/total_canonical:.1f}%)")

    if args.output:
        output_data = {
            "n1n9_short": all_n1n9_short[:100],  # Sample
            "n1n9_long": all_n1n9_long[:100],
            "low_score": all_low_score[:100],
            "no_hbonds": all_no_hbonds[:100],
            "summary": {
                "total_canonical": total_canonical,
                "n1n9_short": len(all_n1n9_short),
                "n1n9_long": len(all_n1n9_long),
                "low_score": len(all_low_score),
                "no_hbonds": len(all_no_hbonds),
                "ok": len(all_ok),
                "missing": len(all_missing),
            }
        }
        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
