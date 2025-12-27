#!/usr/bin/env python3
"""Score H-bond pattern match for LW class classification.

Each LW class has expected H-bond patterns. This module scores how well
observed H-bonds match those patterns.
"""

import json
import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict

# Import centralized H-bond patterns
from prototypes.pair_identification.hbond.patterns import HBOND_PATTERNS, HBondPattern as HBondPatternType


@dataclass
class HBond:
    """A hydrogen bond."""
    donor_atom: str
    acceptor_atom: str
    distance: float = 0.0

    def __repr__(self):
        return f"{self.donor_atom}-{self.acceptor_atom}"

    def matches(self, other: "HBond") -> bool:
        """Check if this H-bond matches another (same atoms)."""
        return (self.donor_atom == other.donor_atom and
                self.acceptor_atom == other.acceptor_atom)


@dataclass
class ScoringResult:
    """Result of scoring H-bonds against a pattern."""
    lw_class: str
    sequence: str
    expected_hbonds: int
    matched_hbonds: int
    score: float  # 0.0 to 1.0
    matches: List[Tuple[str, str]]  # Which H-bonds matched
    missing: List[Tuple[str, str]]  # Which expected H-bonds were missing


class HBondScorer:
    """Score H-bond patterns for LW classification."""

    def __init__(self, patterns: Dict = None):
        self.patterns = patterns or HBOND_PATTERNS

    def score_pattern(
        self,
        observed_hbonds: List[HBond],
        lw_class: str,
        sequence: str,
    ) -> ScoringResult:
        """Score how well observed H-bonds match expected pattern.

        Args:
            observed_hbonds: List of observed H-bonds
            lw_class: LW class to check (e.g., "cWW")
            sequence: Base pair sequence (e.g., "GC")

        Returns:
            ScoringResult with score 0.0 to 1.0
        """
        if lw_class not in self.patterns:
            return ScoringResult(
                lw_class=lw_class,
                sequence=sequence,
                expected_hbonds=0,
                matched_hbonds=0,
                score=0.0,
                matches=[],
                missing=[],
            )

        pattern = self.patterns[lw_class].get(sequence, [])
        if not pattern:
            # Try reversed sequence
            rev_seq = sequence[1] + sequence[0]
            pattern = self.patterns[lw_class].get(rev_seq, [])
            if pattern:
                # Need to swap donor/acceptor for reversed sequence
                pattern = [(a, d) for d, a in pattern]

        if not pattern:
            return ScoringResult(
                lw_class=lw_class,
                sequence=sequence,
                expected_hbonds=0,
                matched_hbonds=0,
                score=0.0,
                matches=[],
                missing=[],
            )

        # Check which expected H-bonds are present
        matches = []
        missing = []

        for donor, acceptor in pattern:
            found = False
            for hb in observed_hbonds:
                # Check both directions
                if ((hb.donor_atom == donor and hb.acceptor_atom == acceptor) or
                    (hb.donor_atom == acceptor and hb.acceptor_atom == donor)):
                    found = True
                    matches.append((donor, acceptor))
                    break

            if not found:
                missing.append((donor, acceptor))

        expected = len(pattern)
        matched = len(matches)
        score = matched / expected if expected > 0 else 0.0

        return ScoringResult(
            lw_class=lw_class,
            sequence=sequence,
            expected_hbonds=expected,
            matched_hbonds=matched,
            score=score,
            matches=matches,
            missing=missing,
        )

    def score_all_classes(
        self,
        observed_hbonds: List[HBond],
        sequence: str,
        lw_classes: List[str] = None,
    ) -> List[ScoringResult]:
        """Score H-bonds against all LW class patterns.

        Args:
            observed_hbonds: List of observed H-bonds
            sequence: Base pair sequence (e.g., "GC")
            lw_classes: LW classes to try (default: all)

        Returns:
            List of ScoringResult sorted by score (best first)
        """
        if lw_classes is None:
            lw_classes = list(self.patterns.keys())

        results = []
        for lw in lw_classes:
            result = self.score_pattern(observed_hbonds, lw, sequence)
            results.append(result)

        # Sort by score (best first)
        results.sort(key=lambda r: (-r.score, -r.matched_hbonds))
        return results


def load_modern_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[HBond]]:
    """Load H-bonds from modern JSON output."""
    if not hbond_path.exists():
        return {}

    with open(hbond_path) as f:
        data = json.load(f)

    result = {}
    for entry in data:
        res_id1 = entry.get("res_id_i", "")
        res_id2 = entry.get("res_id_j", "")

        hbonds = []
        for h in entry.get("hbonds", []):
            # Only base-base H-bonds
            if h.get("context") != "base_base":
                continue
            hbonds.append(HBond(
                donor_atom=h.get("donor_atom", ""),
                acceptor_atom=h.get("acceptor_atom", ""),
                distance=h.get("distance", 0.0),
            ))

        if hbonds:
            key = tuple(sorted([res_id1, res_id2]))
            result[key] = hbonds

    return result


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Score H-bond patterns for LW classification"
    )
    parser.add_argument(
        "--pdb", type=str, required=True,
        help="PDB ID"
    )
    parser.add_argument(
        "--res1", type=str, required=True,
        help="First residue ID"
    )
    parser.add_argument(
        "--res2", type=str, required=True,
        help="Second residue ID"
    )
    parser.add_argument(
        "--hbond-dir", type=Path, default=Path("data/json/all_hbond_list"),
        help="Directory with H-bond JSON files"
    )

    args = parser.parse_args()

    # Load H-bonds
    hbonds = load_modern_hbonds(args.hbond_dir / f"{args.pdb}.json")

    pair_key = tuple(sorted([args.res1, args.res2]))
    pair_hbonds = hbonds.get(pair_key, [])

    if not pair_hbonds:
        print(f"No base-base H-bonds found for {args.res1} - {args.res2}")
        return

    print(f"Observed H-bonds for {args.res1} - {args.res2}:")
    for hb in pair_hbonds:
        print(f"  {hb.donor_atom} -> {hb.acceptor_atom} ({hb.distance:.2f}Å)")

    # Determine sequence (need to load from PDB)
    # For now, extract from res_id
    base1 = args.res1.split("-")[1][-1].upper()
    base2 = args.res2.split("-")[1][-1].upper()
    sequence = base1 + base2

    print(f"\nSequence: {sequence}")

    # Score all patterns
    scorer = HBondScorer()
    results = scorer.score_all_classes(pair_hbonds, sequence)

    print(f"\nH-bond pattern scores:")
    for r in results[:10]:
        if r.expected_hbonds > 0:
            print(f"  {r.lw_class:6s} {r.sequence:4s}: {r.matched_hbonds}/{r.expected_hbonds} "
                  f"({r.score:.0%}) matched={r.matches}")


if __name__ == "__main__":
    main()
