#!/usr/bin/env python3
"""Analyze false negatives to determine greedy selection vs validation failures.

Key question: How many DSSR pairs that we miss are due to:
1. Greedy selection conflict - pair passes validation but loses to higher-scoring competitor
2. Validation failure - pair doesn't pass validation at all
"""

import json
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent))

from find_cww_pairs import (
    parse_pdb, find_initial_candidates, filter_by_n1n9,
    load_slot_hbonds, validate_candidates, load_dssr_cww_pairs,
    normalize_res_id_base, TemplateAligner, MAX_C1_DISTANCE, N1N9_RANGE_CWW,
    CANONICAL_WC_SEQUENCES, CandidatePair
)


@dataclass
class FNAnalysis:
    """Analysis of a false negative."""
    res_id1: str
    res_id2: str
    sequence: str
    reason: str  # 'greedy_conflict', 'validation_failed', 'no_n1n9', 'not_canonical', etc.
    details: dict = None

    def __post_init__(self):
        if self.details is None:
            self.details = {}


def analyze_single_pdb(
    pdb_id: str,
    pdb_dir: Path,
    dssr_dir: Path,
    hbond_dir: Path,
    template_dir: Path,
) -> Tuple[List[FNAnalysis], dict]:
    """Analyze all false negatives for a single PDB.

    Returns list of FNAnalysis and summary stats.
    """
    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"

    if not pdb_path.exists():
        return [], {"error": f"PDB not found: {pdb_id}"}

    # Load DSSR pairs
    dssr_path = dssr_dir / f"{pdb_id}.json"
    dssr_pairs_raw = load_dssr_cww_pairs(dssr_path)

    if not dssr_pairs_raw:
        return [], {"error": "No DSSR pairs found"}

    # Parse PDB
    residues = parse_pdb(pdb_path)

    # Find all candidates
    initial = find_initial_candidates(residues, MAX_C1_DISTANCE)

    # Filter by N1N9 (canonical only)
    n1n9_candidates = filter_by_n1n9(residues, initial, N1N9_RANGE_CWW, canonical_only=True)

    # Also get non-canonical filtered for checking
    n1n9_all = filter_by_n1n9(residues, initial, N1N9_RANGE_CWW, canonical_only=False)

    # Load H-bonds
    hbond_path = hbond_dir / f"{pdb_id}.json"
    slot_hbonds = load_slot_hbonds(hbond_path)

    # Create aligner
    aligner = TemplateAligner(idealized_dir=template_dir)

    # Validate WITHOUT greedy selection first - get ALL pairs that pass validation
    all_validated = validate_candidates(
        residues, n1n9_candidates, slot_hbonds, aligner, greedy_selection=False
    )
    validated_set = {
        tuple(sorted([normalize_res_id_base(p.res_id1), normalize_res_id_base(p.res_id2)]))
        for p in all_validated
    }

    # Get pairs with greedy selection
    greedy_validated = validate_candidates(
        residues, n1n9_candidates, slot_hbonds, aligner, greedy_selection=True
    )
    greedy_set = {
        tuple(sorted([normalize_res_id_base(p.res_id1), normalize_res_id_base(p.res_id2)]))
        for p in greedy_validated
    }

    # Build lookup for validated pairs
    validated_by_key = {
        tuple(sorted([normalize_res_id_base(p.res_id1), normalize_res_id_base(p.res_id2)])): p
        for p in all_validated
    }

    # Build lookup for n1n9 candidates
    n1n9_by_key = {
        tuple(sorted([normalize_res_id_base(p.res_id1), normalize_res_id_base(p.res_id2)])): p
        for p in n1n9_all
    }

    # Normalize DSSR pairs
    dssr_pairs_norm = {
        tuple(sorted([normalize_res_id_base(r1), normalize_res_id_base(r2)]))
        for r1, r2 in dssr_pairs_raw
    }

    # Find false negatives (DSSR pairs we don't find with greedy)
    false_negatives = dssr_pairs_norm - greedy_set

    analyses = []
    reason_counts = defaultdict(int)

    for fn_key in false_negatives:
        res_id1, res_id2 = fn_key

        # Determine sequence from residue IDs
        parts1 = res_id1.split("-")
        parts2 = res_id2.split("-")
        base1 = parts1[1] if len(parts1) >= 2 else "?"
        base2 = parts2[1] if len(parts2) >= 2 else "?"
        sequence = base1 + base2

        # Check if this pair passes validation but loses to greedy
        if fn_key in validated_set:
            # This pair PASSES validation but was eliminated by greedy selection
            pair_info = validated_by_key[fn_key]

            # Find competing pairs - which residue is shared?
            competing = []
            for gp in greedy_validated:
                gp_norm = tuple(sorted([
                    normalize_res_id_base(gp.res_id1),
                    normalize_res_id_base(gp.res_id2)
                ]))
                # Check if shares a residue
                if res_id1 in gp_norm or res_id2 in gp_norm:
                    competing.append({
                        "res_id1": gp.res_id1,
                        "res_id2": gp.res_id2,
                        "score": gp.score,
                        "hbond_count": gp.hbond_count,
                        "sequence": gp.sequence,
                    })

            analysis = FNAnalysis(
                res_id1=res_id1,
                res_id2=res_id2,
                sequence=sequence,
                reason="greedy_conflict",
                details={
                    "our_score": pair_info.score,
                    "our_hbonds": pair_info.hbond_count,
                    "our_rmsd": pair_info.rmsd_cww,
                    "competing_pairs": competing,
                }
            )
            reason_counts["greedy_conflict"] += 1

        elif fn_key in n1n9_by_key:
            # In N1N9 range but failed validation
            cand = n1n9_by_key[fn_key]

            # Check why it failed - need to revalidate to get details
            res1_found = res_id1 in residues or any(normalize_res_id_base(r) == res_id1 for r in residues)
            res2_found = res_id2 in residues or any(normalize_res_id_base(r) == res_id2 for r in residues)

            analysis = FNAnalysis(
                res_id1=res_id1,
                res_id2=res_id2,
                sequence=sequence,
                reason="validation_failed",
                details={
                    "n1n9_distance": cand.n1n9_distance,
                    "score": cand.score,
                    "hbond_count": cand.hbond_count,
                    "rmsd": cand.rmsd_cww,
                }
            )
            reason_counts["validation_failed"] += 1

        else:
            # Not in N1N9 range or other issue
            # Check if non-canonical sequence
            if sequence not in CANONICAL_WC_SEQUENCES:
                analysis = FNAnalysis(
                    res_id1=res_id1,
                    res_id2=res_id2,
                    sequence=sequence,
                    reason="non_canonical_sequence",
                )
                reason_counts["non_canonical_sequence"] += 1
            else:
                # Check if N1N9 out of range or missing atoms
                analysis = FNAnalysis(
                    res_id1=res_id1,
                    res_id2=res_id2,
                    sequence=sequence,
                    reason="n1n9_out_of_range",
                )
                reason_counts["n1n9_out_of_range"] += 1

        analyses.append(analysis)

    summary = {
        "pdb_id": pdb_id,
        "dssr_pairs": len(dssr_pairs_norm),
        "validated_no_greedy": len(validated_set),
        "validated_with_greedy": len(greedy_set),
        "false_negatives": len(false_negatives),
        "reason_counts": dict(reason_counts),
    }

    return analyses, summary


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Analyze false negative reasons")

    project_root = Path(__file__).parent.parent.parent

    parser.add_argument("--pdb-dir", type=Path, default=project_root / "data" / "pdb")
    parser.add_argument("--dssr-dir", type=Path, default=project_root / "data" / "json_dssr")
    parser.add_argument("--hbond-dir", type=Path, default=project_root / "data" / "json" / "slot_hbonds")
    parser.add_argument("--template-dir", type=Path, default=project_root / "basepair-idealized")
    parser.add_argument("--pdb", type=str, help="Single PDB to analyze")
    parser.add_argument("--sample", type=int, default=100, help="Number of PDBs to sample")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--output", "-o", type=Path, help="Output JSON file")

    args = parser.parse_args()

    if args.pdb:
        pdb_ids = [args.pdb]
    else:
        # Get list of PDBs with slot H-bond files
        hbond_files = list(args.hbond_dir.glob("*.json"))
        pdb_ids = [f.stem for f in hbond_files[:args.sample]]

    print(f"Analyzing {len(pdb_ids)} PDBs...")

    all_analyses = []
    all_summaries = []
    total_reasons = defaultdict(int)

    for i, pdb_id in enumerate(pdb_ids):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(pdb_ids)}...")

        analyses, summary = analyze_single_pdb(
            pdb_id, args.pdb_dir, args.dssr_dir, args.hbond_dir, args.template_dir
        )

        all_analyses.extend(analyses)
        all_summaries.append(summary)

        for reason, count in summary.get("reason_counts", {}).items():
            total_reasons[reason] += count

    # Print summary
    print(f"\n{'='*60}")
    print("FALSE NEGATIVE ANALYSIS")
    print(f"{'='*60}")

    total_fn = sum(s.get("false_negatives", 0) for s in all_summaries)
    total_dssr = sum(s.get("dssr_pairs", 0) for s in all_summaries)
    total_validated_no_greedy = sum(s.get("validated_no_greedy", 0) for s in all_summaries)
    total_validated_with_greedy = sum(s.get("validated_with_greedy", 0) for s in all_summaries)

    print(f"\nTotal DSSR pairs: {total_dssr}")
    print(f"Validated (no greedy): {total_validated_no_greedy}")
    print(f"Validated (with greedy): {total_validated_with_greedy}")
    print(f"Lost to greedy: {total_validated_no_greedy - total_validated_with_greedy}")
    print(f"Total false negatives: {total_fn}")

    print(f"\nReason breakdown:")
    for reason, count in sorted(total_reasons.items(), key=lambda x: -x[1]):
        pct = 100 * count / total_fn if total_fn > 0 else 0
        print(f"  {reason}: {count} ({pct:.1f}%)")

    # Show examples of greedy conflicts
    greedy_conflicts = [a for a in all_analyses if a.reason == "greedy_conflict"]
    if greedy_conflicts and args.verbose:
        print(f"\n\nExample greedy conflicts (showing first 10):")
        for gc in greedy_conflicts[:10]:
            print(f"\n  {gc.res_id1} - {gc.res_id2} ({gc.sequence})")
            print(f"    Our score: {gc.details['our_score']:.3f}, H-bonds: {gc.details['our_hbonds']}")
            print(f"    Competing pairs:")
            for cp in gc.details.get("competing_pairs", []):
                print(f"      {cp['res_id1']} - {cp['res_id2']} ({cp['sequence']}) "
                      f"score={cp['score']:.3f} hbonds={cp['hbond_count']}")

    # Save results
    if args.output:
        output_data = {
            "summary": {
                "pdb_count": len(pdb_ids),
                "total_dssr_pairs": total_dssr,
                "total_false_negatives": total_fn,
                "validated_without_greedy": total_validated_no_greedy,
                "validated_with_greedy": total_validated_with_greedy,
                "lost_to_greedy": total_validated_no_greedy - total_validated_with_greedy,
                "reason_counts": dict(total_reasons),
            },
            "per_pdb": all_summaries,
            "greedy_conflicts": [
                {
                    "res_id1": a.res_id1,
                    "res_id2": a.res_id2,
                    "sequence": a.sequence,
                    "details": a.details,
                }
                for a in greedy_conflicts
            ]
        }

        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults saved to {args.output}")


if __name__ == "__main__":
    main()
