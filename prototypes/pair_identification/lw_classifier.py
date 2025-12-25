#!/usr/bin/env python3
"""Combined LW classifier using template RMSD + H-bond pattern matching.

For each candidate pair:
1. Align to ALL LW class templates (cWW, tSH, tWH, cSW, etc.)
2. Compute RMSD for each template alignment
3. Compute H-bond score based on expected H-bond patterns
4. Combine scores to determine best classification + confidence
"""

import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

from core.pdb_parser import parse_pdb
from core.residue import normalize_base_type
from template_aligner import (
    TemplateAligner, AlignmentResult, ClassificationResult as AlignmentClassification,
    Residue, Atom
)
from hbond_scorer import HBondScorer, HBond, ScoringResult, load_modern_hbonds


@dataclass
class LWClassificationResult:
    """Result of classifying a pair using RMSD + H-bond scoring."""
    res_id1: str
    res_id2: str
    sequence: str

    # Best match
    best_lw: str
    best_rmsd: float
    best_hbond_score: float
    best_combined_score: float

    # Second best match
    second_lw: Optional[str] = None
    second_rmsd: Optional[float] = None
    second_hbond_score: Optional[float] = None
    second_combined_score: Optional[float] = None

    # All scored results
    all_results: List[dict] = field(default_factory=list)

    @property
    def confidence(self) -> float:
        """Confidence based on gap between best and second best combined score."""
        if self.second_combined_score is None:
            return 1.0
        gap = self.best_combined_score - self.second_combined_score
        # Normalize: 0.3 gap = high confidence
        return min(1.0, gap / 0.3)

    @property
    def rmsd_gap(self) -> Optional[float]:
        """Gap between best and second best RMSD."""
        if self.second_rmsd is None:
            return None
        return self.second_rmsd - self.best_rmsd

    def __repr__(self):
        conf_str = f"{self.confidence:.0%}"
        return (f"LWClassificationResult({self.res_id1}-{self.res_id2} {self.sequence}: "
                f"best={self.best_lw} RMSD={self.best_rmsd:.3f}Å hbond={self.best_hbond_score:.0%} "
                f"score={self.best_combined_score:.3f} conf={conf_str})")


class LWClassifier:
    """Classify base pairs using template RMSD + H-bond scoring."""

    # LW classes to try
    DEFAULT_LW_CLASSES = ["cWW", "tWW", "cWH", "tWH", "cWS", "tWS",
                          "cHH", "tHH", "cHS", "tHS", "cSS", "tSS"]

    def __init__(
        self,
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
        rmsd_weight: float = 0.4,
        hbond_weight: float = 0.6,
        rmsd_cutoff: float = 2.0,  # RMSD above this gets 0 score
    ):
        self.aligner = TemplateAligner(idealized_dir, exemplar_dir)
        self.hbond_scorer = HBondScorer()
        self.rmsd_weight = rmsd_weight
        self.hbond_weight = hbond_weight
        self.rmsd_cutoff = rmsd_cutoff

    def _rmsd_to_score(self, rmsd: float) -> float:
        """Convert RMSD to a 0-1 score (lower RMSD = higher score)."""
        if rmsd >= self.rmsd_cutoff:
            return 0.0
        return 1.0 - (rmsd / self.rmsd_cutoff)

    def classify(
        self,
        res1: Residue,
        res2: Residue,
        hbonds: List[HBond],
        lw_classes: Optional[List[str]] = None,
    ) -> LWClassificationResult:
        """Classify a pair by trying all LW class templates.

        Args:
            res1: First residue
            res2: Second residue
            hbonds: Observed H-bonds between the residues
            lw_classes: LW classes to try (default: all 12)

        Returns:
            LWClassificationResult with best match and confidence
        """
        if lw_classes is None:
            lw_classes = self.DEFAULT_LW_CLASSES

        sequence = res1.base_type + res2.base_type

        # Get template alignment results
        alignment_result = self.aligner.classify_pair(res1, res2, lw_classes)

        # Score each LW class using combined RMSD + H-bond
        scored_results = []

        for align in alignment_result.all_results:
            # Skip if RMSD is infinite (no template found)
            if align.rmsd == float('inf'):
                continue

            # Get RMSD score (0-1, higher is better)
            rmsd_score = self._rmsd_to_score(align.rmsd)

            # Get H-bond score for this LW class
            # Use the sequence from alignment (handles reversed sequences)
            hbond_seq = align.sequence.replace("(rev)", "")
            hbond_result = self.hbond_scorer.score_pattern(hbonds, align.lw_class, hbond_seq)

            # Combined score with bonus for more H-bonds matched
            # This helps discriminate cWW (3 H-bonds) from tWW (2 H-bonds)
            # when both have 100% match
            hbond_count_bonus = hbond_result.matched_hbonds * 0.05  # Larger bonus
            combined = (self.hbond_weight * hbond_result.score +
                       self.rmsd_weight * rmsd_score +
                       hbond_count_bonus)

            scored_results.append({
                "lw_class": align.lw_class,
                "sequence": align.sequence,
                "rmsd": align.rmsd,
                "rmsd_score": rmsd_score,
                "hbond_score": hbond_result.score,
                "hbond_matched": hbond_result.matched_hbonds,
                "hbond_expected": hbond_result.expected_hbonds,
                "combined_score": combined,
                "num_atoms": align.num_atoms_aligned,
            })

        # Sort by combined score (best first)
        scored_results.sort(key=lambda x: -x["combined_score"])

        if not scored_results:
            return LWClassificationResult(
                res_id1=res1.res_id,
                res_id2=res2.res_id,
                sequence=sequence,
                best_lw="unknown",
                best_rmsd=float('inf'),
                best_hbond_score=0.0,
                best_combined_score=0.0,
            )

        best = scored_results[0]
        second = scored_results[1] if len(scored_results) > 1 else None

        return LWClassificationResult(
            res_id1=res1.res_id,
            res_id2=res2.res_id,
            sequence=sequence,
            best_lw=best["lw_class"],
            best_rmsd=best["rmsd"],
            best_hbond_score=best["hbond_score"],
            best_combined_score=best["combined_score"],
            second_lw=second["lw_class"] if second else None,
            second_rmsd=second["rmsd"] if second else None,
            second_hbond_score=second["hbond_score"] if second else None,
            second_combined_score=second["combined_score"] if second else None,
            all_results=scored_results,
        )


def load_dssr_pairs(dssr_path: Path) -> List[dict]:
    """Load base pairs from DSSR JSON."""
    if not dssr_path.exists():
        return []

    with open(dssr_path) as f:
        data = json.load(f)

    return data.get("pairs", [])


def normalize_dssr_nt(nt: str) -> str:
    """Convert DSSR nt format (A.G1) to res_id format (A-G-1)."""
    parts = nt.split(".")
    if len(parts) != 2:
        return nt
    chain = parts[0]
    base_num = parts[1]

    # Find where the residue number starts - work backwards
    i = len(base_num) - 1
    if i >= 0 and base_num[i].isalpha() and i > 0 and base_num[i-1].isdigit():
        i -= 1
    while i >= 0 and base_num[i].isdigit():
        i -= 1
    i += 1

    if i <= 0:
        return f"{chain}-{base_num}-0"

    base = base_num[:i]
    num = base_num[i:]
    return f"{chain}-{base}-{num}"


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Classify base pairs using template RMSD + H-bond scoring"
    )
    parser.add_argument(
        "--pdb", type=str, required=True,
        help="PDB ID or path to PDB file"
    )
    parser.add_argument(
        "--res1", type=str, default=None,
        help="First residue ID (e.g., A-G-1)"
    )
    parser.add_argument(
        "--res2", type=str, default=None,
        help="Second residue ID (e.g., A-C-72)"
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
        "--all-dssr", action="store_true",
        help="Classify all DSSR pairs and compare"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    # Load PDB
    if Path(args.pdb).exists():
        pdb_path = Path(args.pdb)
        pdb_id = pdb_path.stem
    else:
        pdb_id = args.pdb
        pdb_path = args.pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            pdb_path = args.pdb_dir / f"{pdb_id.lower()}.pdb"

    if not pdb_path.exists():
        print(f"PDB file not found: {pdb_path}")
        return

    residues = parse_pdb(pdb_path)
    hbonds = load_modern_hbonds(args.hbond_dir / f"{pdb_id}.json")

    classifier = LWClassifier(args.idealized_dir, args.exemplar_dir)

    if args.all_dssr:
        # Classify all DSSR pairs and compare
        dssr_pairs = load_dssr_pairs(args.dssr_dir / f"{pdb_id}.json")

        if not dssr_pairs:
            print(f"No DSSR pairs found for {pdb_id}")
            return

        correct = 0
        total = 0
        mismatches = []

        for pair in dssr_pairs:
            dssr_lw = pair.get("LW", "")
            if not dssr_lw:
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

            total += 1
            if result.best_lw == dssr_lw:
                correct += 1
            else:
                mismatches.append({
                    "res1": res_id1,
                    "res2": res_id2,
                    "dssr": dssr_lw,
                    "predicted": result.best_lw,
                    "rmsd": result.best_rmsd,
                    "hbond": result.best_hbond_score,
                    "confidence": result.confidence,
                    "second": result.second_lw,
                })

        accuracy = correct / total if total > 0 else 0
        print(f"\nClassification Results for {pdb_id}:")
        print(f"  Accuracy: {correct}/{total} ({accuracy:.1%})")

        if mismatches and args.verbose:
            print(f"\nMismatches ({len(mismatches)}):")
            for m in mismatches[:20]:
                print(f"  {m['res1']}-{m['res2']}: DSSR={m['dssr']} Pred={m['predicted']} "
                      f"(RMSD={m['rmsd']:.3f} hbond={m['hbond']:.0%} conf={m['confidence']:.0%} 2nd={m['second']})")

        # Summary by LW class
        lw_counts = defaultdict(lambda: {"correct": 0, "total": 0})
        for pair in dssr_pairs:
            dssr_lw = pair.get("LW", "")
            if dssr_lw:
                lw_counts[dssr_lw]["total"] += 1

        for m in mismatches:
            lw_counts[m["dssr"]]["total"] -= 1  # We counted it above
        for pair in dssr_pairs:
            dssr_lw = pair.get("LW", "")
            res_id1 = normalize_dssr_nt(pair.get("nt1", ""))
            res_id2 = normalize_dssr_nt(pair.get("nt2", ""))
            if res_id1 in residues and res_id2 in residues:
                if dssr_lw not in [m["dssr"] for m in mismatches if m["res1"] == res_id1 and m["res2"] == res_id2]:
                    lw_counts[dssr_lw]["correct"] += 1

        print(f"\nBy LW class:")
        for lw in sorted(lw_counts.keys()):
            c = lw_counts[lw]
            acc = c["correct"] / c["total"] if c["total"] > 0 else 0
            print(f"  {lw:6s}: {c['correct']}/{c['total']} ({acc:.0%})")

        return

    # Single pair classification
    if args.res1 is None or args.res2 is None:
        print("Specify --res1 and --res2 for single pair, or use --all-dssr")
        return

    if args.res1 not in residues:
        print(f"Residue {args.res1} not found in PDB")
        return
    if args.res2 not in residues:
        print(f"Residue {args.res2} not found in PDB")
        return

    res1 = residues[args.res1]
    res2 = residues[args.res2]

    # Get H-bonds for this pair
    pair_key = tuple(sorted([args.res1, args.res2]))
    pair_hbonds = hbonds.get(pair_key, [])

    print(f"\nObserved H-bonds for {args.res1} - {args.res2}:")
    for hb in pair_hbonds:
        print(f"  {hb.donor_atom} -> {hb.acceptor_atom} ({hb.distance:.2f}Å)")

    result = classifier.classify(res1, res2, pair_hbonds)

    print(f"\n{result}")
    print(f"\nTop 5 classifications:")
    for r in result.all_results[:5]:
        print(f"  {r['lw_class']:6s} {r['sequence']:8s} "
              f"RMSD={r['rmsd']:.3f}Å hbond={r['hbond_matched']}/{r['hbond_expected']} "
              f"({r['hbond_score']:.0%}) combined={r['combined_score']:.3f}")


if __name__ == "__main__":
    main()
