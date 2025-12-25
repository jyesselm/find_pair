#!/usr/bin/env python3
"""Load and analyze H-bond data from modern X3DNA and DSSR.

This module provides utilities for loading H-bond information from both
modern X3DNA JSON output and DSSR JSON output, and comparing them.
"""

import json
import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict


@dataclass
class HBond:
    """A single hydrogen bond."""
    donor_atom: str
    acceptor_atom: str
    distance: float
    donor_annotation: str = ""  # e.g., "(amino)", "(imino)"
    acceptor_annotation: str = ""  # e.g., "(carbonyl)"

    def __repr__(self):
        return f"{self.donor_atom}-{self.acceptor_atom}[{self.distance:.2f}]"

    def matches(self, other: "HBond", dist_tolerance: float = 0.3) -> bool:
        """Check if this H-bond matches another (same atoms, similar distance)."""
        return (
            self.donor_atom == other.donor_atom and
            self.acceptor_atom == other.acceptor_atom and
            abs(self.distance - other.distance) < dist_tolerance
        )


@dataclass
class PairHBonds:
    """H-bonds between a pair of residues."""
    res_id1: str
    res_id2: str
    hbonds: List[HBond] = field(default_factory=list)
    source: str = ""  # "modern" or "dssr"

    def pair_key(self) -> Tuple[str, str]:
        """Return sorted tuple for order-independent comparison."""
        return tuple(sorted([self.res_id1, self.res_id2]))

    @property
    def base_base_hbonds(self) -> List[HBond]:
        """Return only base-base H-bonds."""
        base_atoms = {"N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6"}
        return [h for h in self.hbonds
                if h.donor_atom in base_atoms and h.acceptor_atom in base_atoms]


def parse_dssr_hbonds_desc(hbonds_desc: str) -> List[HBond]:
    """Parse DSSR hbonds_desc string into HBond objects.

    Format: "O6(carbonyl)-N4(amino)[2.83],N1(imino)-N3[2.88],N2(amino)-O2(carbonyl)[2.84]"
    """
    if not hbonds_desc:
        return []

    hbonds = []
    # Pattern: atom1(annotation)?-atom2(annotation)?[distance]
    # or atom1(annotation)?*atom2(annotation)?[distance] for non-WC
    pattern = r"(\w+)(?:\(([^)]+)\))?[*-](\w+)(?:\(([^)]+)\))?\[(\d+\.?\d*)\]"

    for match in re.finditer(pattern, hbonds_desc):
        donor = match.group(1)
        donor_ann = match.group(2) or ""
        acceptor = match.group(3)
        acceptor_ann = match.group(4) or ""
        distance = float(match.group(5))

        hbonds.append(HBond(
            donor_atom=donor,
            acceptor_atom=acceptor,
            distance=distance,
            donor_annotation=donor_ann,
            acceptor_annotation=acceptor_ann,
        ))

    return hbonds


def normalize_dssr_nt(nt: str) -> str:
    """Convert DSSR nt format (A.G1) to res_id format (A-G-1)."""
    parts = nt.split(".")
    if len(parts) != 2:
        return nt
    chain = parts[0]
    base_num = parts[1]

    # Find where the residue number starts
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


def load_modern_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], PairHBonds]:
    """Load H-bonds from modern X3DNA JSON output."""
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
            # Only include base-base H-bonds for comparison
            context = h.get("context", "")
            if context != "base_base":
                continue

            hbonds.append(HBond(
                donor_atom=h.get("donor_atom", ""),
                acceptor_atom=h.get("acceptor_atom", ""),
                distance=h.get("distance", 0.0),
            ))

        if hbonds:
            pair = PairHBonds(
                res_id1=res_id1,
                res_id2=res_id2,
                hbonds=hbonds,
                source="modern"
            )
            result[pair.pair_key()] = pair

    return result


def load_dssr_hbonds(dssr_path: Path, lw_filter: str = None) -> Dict[Tuple[str, str], PairHBonds]:
    """Load H-bonds from DSSR JSON output."""
    if not dssr_path.exists():
        return {}

    with open(dssr_path) as f:
        data = json.load(f)

    result = {}
    for p in data.get("pairs", []):
        if lw_filter and p.get("LW") != lw_filter:
            continue

        res_id1 = normalize_dssr_nt(p.get("nt1", ""))
        res_id2 = normalize_dssr_nt(p.get("nt2", ""))
        hbonds_desc = p.get("hbonds_desc", "")

        hbonds = parse_dssr_hbonds_desc(hbonds_desc)

        if hbonds:
            pair = PairHBonds(
                res_id1=res_id1,
                res_id2=res_id2,
                hbonds=hbonds,
                source="dssr"
            )
            result[pair.pair_key()] = pair

    return result


def compare_hbonds(
    modern_hbonds: Dict[Tuple[str, str], PairHBonds],
    dssr_hbonds: Dict[Tuple[str, str], PairHBonds],
    dist_tolerance: float = 0.3,
) -> Dict:
    """Compare H-bonds between modern and DSSR output.

    Returns statistics about matching and differences.
    """
    modern_keys = set(modern_hbonds.keys())
    dssr_keys = set(dssr_hbonds.keys())

    common_pairs = modern_keys & dssr_keys
    modern_only = modern_keys - dssr_keys
    dssr_only = dssr_keys - modern_keys

    # For common pairs, compare individual H-bonds
    matching_hbonds = 0
    total_dssr_hbonds = 0
    total_modern_hbonds = 0
    mismatched_pairs = []

    for key in common_pairs:
        modern_pair = modern_hbonds[key]
        dssr_pair = dssr_hbonds[key]

        total_dssr_hbonds += len(dssr_pair.hbonds)
        total_modern_hbonds += len(modern_pair.hbonds)

        # Count matching H-bonds
        matched = 0
        for dh in dssr_pair.hbonds:
            for mh in modern_pair.hbonds:
                if dh.matches(mh, dist_tolerance):
                    matched += 1
                    break

        matching_hbonds += matched

        if matched < len(dssr_pair.hbonds):
            mismatched_pairs.append({
                "pair": key,
                "dssr_hbonds": dssr_pair.hbonds,
                "modern_hbonds": modern_pair.hbonds,
                "matched": matched,
            })

    return {
        "common_pairs": len(common_pairs),
        "modern_only_pairs": len(modern_only),
        "dssr_only_pairs": len(dssr_only),
        "total_dssr_hbonds": total_dssr_hbonds,
        "total_modern_hbonds": total_modern_hbonds,
        "matching_hbonds": matching_hbonds,
        "mismatched_pairs": mismatched_pairs[:20],  # Limit for readability
    }


def analyze_pdb_hbonds(
    pdb_id: str,
    modern_dir: Path,
    dssr_dir: Path,
    lw_filter: str = "cWW",
    verbose: bool = False,
) -> Dict:
    """Analyze H-bond matching for a single PDB."""
    modern_hbonds = load_modern_hbonds(modern_dir / f"{pdb_id}.json")
    dssr_hbonds = load_dssr_hbonds(dssr_dir / f"{pdb_id}.json", lw_filter)

    comparison = compare_hbonds(modern_hbonds, dssr_hbonds)

    if verbose:
        print(f"\n{pdb_id}:")
        print(f"  Common pairs: {comparison['common_pairs']}")
        print(f"  DSSR-only pairs: {comparison['dssr_only_pairs']}")
        print(f"  Modern-only pairs: {comparison['modern_only_pairs']}")
        if comparison['total_dssr_hbonds'] > 0:
            match_pct = 100 * comparison['matching_hbonds'] / comparison['total_dssr_hbonds']
            print(f"  H-bond match: {comparison['matching_hbonds']}/{comparison['total_dssr_hbonds']} ({match_pct:.1f}%)")

        for mismatch in comparison['mismatched_pairs'][:5]:
            print(f"\n  Mismatch: {mismatch['pair']}")
            print(f"    DSSR:   {mismatch['dssr_hbonds']}")
            print(f"    Modern: {mismatch['modern_hbonds']}")

    return comparison


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Compare H-bonds between modern X3DNA and DSSR")
    parser.add_argument(
        "--modern-dir", type=Path,
        default=Path("data/json/all_hbond_list"),
        help="Directory with modern H-bond JSON files"
    )
    parser.add_argument(
        "--dssr-dir", type=Path,
        default=Path("data/json_dssr"),
        help="Directory with DSSR JSON files"
    )
    parser.add_argument(
        "--pdb", type=str, default=None,
        help="Specific PDB to analyze"
    )
    parser.add_argument(
        "--max-pdbs", type=int, default=50,
        help="Maximum PDBs to analyze"
    )
    parser.add_argument(
        "--lw-class", type=str, default="cWW",
        help="LW class to filter from DSSR"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    if args.pdb:
        analyze_pdb_hbonds(
            args.pdb, args.modern_dir, args.dssr_dir,
            args.lw_class, verbose=True
        )
        return

    # Analyze multiple PDBs
    dssr_files = sorted(args.dssr_dir.glob("*.json"))[:args.max_pdbs]

    total_stats = defaultdict(int)

    for dssr_path in dssr_files:
        pdb_id = dssr_path.stem
        result = analyze_pdb_hbonds(
            pdb_id, args.modern_dir, args.dssr_dir,
            args.lw_class, verbose=args.verbose
        )

        total_stats["common_pairs"] += result["common_pairs"]
        total_stats["dssr_only_pairs"] += result["dssr_only_pairs"]
        total_stats["total_dssr_hbonds"] += result["total_dssr_hbonds"]
        total_stats["matching_hbonds"] += result["matching_hbonds"]

    print(f"\n{'='*60}")
    print(f"H-BOND COMPARISON SUMMARY ({len(dssr_files)} PDBs)")
    print(f"{'='*60}")
    print(f"\nPair coverage:")
    print(f"  Common pairs: {total_stats['common_pairs']}")
    print(f"  DSSR-only pairs: {total_stats['dssr_only_pairs']}")

    if total_stats['total_dssr_hbonds'] > 0:
        match_pct = 100 * total_stats['matching_hbonds'] / total_stats['total_dssr_hbonds']
        print(f"\nH-bond matching:")
        print(f"  Total DSSR H-bonds: {total_stats['total_dssr_hbonds']}")
        print(f"  Matched by modern: {total_stats['matching_hbonds']} ({match_pct:.1f}%)")


if __name__ == "__main__":
    main()
