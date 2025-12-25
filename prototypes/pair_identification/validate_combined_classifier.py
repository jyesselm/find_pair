#!/usr/bin/env python3
"""Validate combined LW classifier (RMSD + H-bond) against DSSR ground truth.

Tests classification accuracy on all pairs from multiple PDBs.
"""

import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
from collections import defaultdict
import argparse

from lw_classifier import (
    LWClassifier, LWClassificationResult,
    parse_pdb_residues, normalize_dssr_nt
)
from hbond_scorer import load_modern_hbonds, HBond


@dataclass
class ValidationStats:
    """Statistics for validation run."""
    total: int = 0
    correct: int = 0
    by_class: Dict[str, Dict] = field(default_factory=lambda: defaultdict(lambda: {"correct": 0, "total": 0}))
    mismatches: List[dict] = field(default_factory=list)


def load_dssr_pairs(dssr_path: Path) -> List[dict]:
    """Load base pairs from DSSR JSON."""
    if not dssr_path.exists():
        return []

    with open(dssr_path) as f:
        data = json.load(f)

    return data.get("pairs", [])


def validate_pdb(
    pdb_id: str,
    classifier: LWClassifier,
    pdb_dir: Path,
    hbond_dir: Path,
    dssr_dir: Path,
    stats: ValidationStats,
    lw_filter: str = None,
) -> None:
    """Validate classifier on a single PDB."""
    pdb_path = pdb_dir / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        return

    residues = parse_pdb_residues(pdb_path)
    hbonds = load_modern_hbonds(hbond_dir / f"{pdb_id}.json")
    dssr_pairs = load_dssr_pairs(dssr_dir / f"{pdb_id}.json")

    for pair in dssr_pairs:
        dssr_lw = pair.get("LW", "")
        if not dssr_lw:
            continue

        # Skip non-standard LW classes (e.g., "c.W", "--")
        if "." in dssr_lw or "-" in dssr_lw:
            continue

        # Filter by LW class if specified
        if lw_filter and dssr_lw != lw_filter:
            continue

        res_id1 = normalize_dssr_nt(pair.get("nt1", ""))
        res_id2 = normalize_dssr_nt(pair.get("nt2", ""))

        if res_id1 not in residues or res_id2 not in residues:
            continue

        res1 = residues[res_id1]
        res2 = residues[res_id2]

        # Get H-bonds for this pair
        pair_key = tuple(sorted([res_id1, res_id2]))
        pair_hbonds = hbonds.get(pair_key, [])

        result = classifier.classify(res1, res2, pair_hbonds)

        stats.total += 1
        stats.by_class[dssr_lw]["total"] += 1

        if result.best_lw == dssr_lw:
            stats.correct += 1
            stats.by_class[dssr_lw]["correct"] += 1
        else:
            stats.mismatches.append({
                "pdb": pdb_id,
                "res1": res_id1,
                "res2": res_id2,
                "sequence": result.sequence,
                "dssr": dssr_lw,
                "predicted": result.best_lw,
                "rmsd": result.best_rmsd,
                "hbond": result.best_hbond_score,
                "confidence": result.confidence,
                "second": result.second_lw,
            })


def main():
    parser = argparse.ArgumentParser(
        description="Validate combined LW classifier against DSSR"
    )
    parser.add_argument(
        "--pdb-dir", type=Path, default=Path("data/pdb"),
        help="Directory with PDB files"
    )
    parser.add_argument(
        "--hbond-dir", type=Path, default=Path("data/json/all_hbond_list"),
        help="Directory with H-bond JSON files"
    )
    parser.add_argument(
        "--dssr-dir", type=Path, default=Path("data/json_dssr"),
        help="Directory with DSSR JSON files"
    )
    parser.add_argument(
        "--idealized-dir", type=Path, default=Path("basepair-idealized"),
        help="Directory with idealized templates"
    )
    parser.add_argument(
        "--exemplar-dir", type=Path, default=Path("basepair-exemplars"),
        help="Directory with exemplar templates"
    )
    parser.add_argument(
        "--pdb", type=str, default=None,
        help="Specific PDB to validate"
    )
    parser.add_argument(
        "--max-pdbs", type=int, default=50,
        help="Maximum number of PDBs to validate"
    )
    parser.add_argument(
        "--lw-filter", type=str, default=None,
        help="Only validate specific LW class (e.g., cWW)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    classifier = LWClassifier(args.idealized_dir, args.exemplar_dir)
    stats = ValidationStats()

    if args.pdb:
        pdb_ids = [args.pdb]
    else:
        # Get PDB IDs from DSSR directory
        dssr_files = sorted(args.dssr_dir.glob("*.json"))[:args.max_pdbs]
        pdb_ids = [f.stem for f in dssr_files]

    print(f"Validating {len(pdb_ids)} PDBs...")
    if args.lw_filter:
        print(f"Filtering for LW class: {args.lw_filter}")

    for i, pdb_id in enumerate(pdb_ids):
        if args.verbose and i % 10 == 0:
            print(f"  Processing {i+1}/{len(pdb_ids)}: {pdb_id}")

        validate_pdb(
            pdb_id, classifier,
            args.pdb_dir, args.hbond_dir, args.dssr_dir,
            stats, args.lw_filter
        )

    # Results
    accuracy = stats.correct / stats.total if stats.total > 0 else 0
    print(f"\n{'='*60}")
    print(f"VALIDATION RESULTS ({len(pdb_ids)} PDBs)")
    print(f"{'='*60}")
    print(f"\nOverall Accuracy: {stats.correct}/{stats.total} ({accuracy:.1%})")

    print(f"\nBy LW class:")
    for lw in sorted(stats.by_class.keys()):
        c = stats.by_class[lw]
        acc = c["correct"] / c["total"] if c["total"] > 0 else 0
        print(f"  {lw:6s}: {c['correct']:4d}/{c['total']:<4d} ({acc:5.1%})")

    if args.verbose and stats.mismatches:
        print(f"\nSample mismatches (first 20):")
        for m in stats.mismatches[:20]:
            print(f"  {m['pdb']}: {m['res1']}-{m['res2']} ({m['sequence']}) "
                  f"DSSR={m['dssr']} Pred={m['predicted']} "
                  f"RMSD={m['rmsd']:.2f} hbond={m['hbond']:.0%}")


if __name__ == "__main__":
    main()
