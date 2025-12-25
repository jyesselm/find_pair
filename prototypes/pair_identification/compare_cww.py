#!/usr/bin/env python3
"""Compare cWW base pairs across DSSR, legacy X3DNA, and modern X3DNA.

This script identifies which cWW pairs are found by each system and
reports coverage gaps.
"""

import json
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple
import argparse


@dataclass
class PairInfo:
    """Information about a base pair."""
    res_id1: str  # Normalized format: "chain-base-num"
    res_id2: str
    bp_type: str  # e.g., "GC", "AU", "GU"
    lw_class: Optional[str] = None  # e.g., "cWW"
    n1n9_dist: Optional[float] = None
    hbonds_desc: Optional[str] = None
    source: str = ""  # "dssr", "legacy", "modern"

    def pair_key(self) -> Tuple[str, str]:
        """Return sorted tuple for comparison (order-independent)."""
        return tuple(sorted([self.res_id1, self.res_id2]))

    def __hash__(self):
        return hash(self.pair_key())

    def __eq__(self, other):
        return self.pair_key() == other.pair_key()


def normalize_dssr_nt(nt: str) -> str:
    """Convert DSSR nt format (A.G1) to res_id format (A-G-1).

    DSSR format: chain.name+resnum
    Examples: A.G1, A.C72, A.2MG10, A.5MC40, A.PSU39
    The name can be 1-3 characters and may start with a digit.
    """
    parts = nt.split(".")
    if len(parts) != 2:
        return nt
    chain = parts[0]
    base_num = parts[1]

    # Find where the residue number starts
    # Look for the last stretch of digits (possibly with insertion code)
    # Work backwards from end to find where number starts
    i = len(base_num) - 1
    # Skip any trailing insertion code (letter at end)
    if i >= 0 and base_num[i].isalpha() and i > 0 and base_num[i-1].isdigit():
        i -= 1
    # Find start of number
    while i >= 0 and base_num[i].isdigit():
        i -= 1
    i += 1  # Move back to start of number

    if i <= 0:
        # No clear separation found
        return f"{chain}-{base_num}-0"

    base = base_num[:i]
    num = base_num[i:]

    return f"{chain}-{base}-{num}"


@dataclass
class DSSRPairInfo(PairInfo):
    """Extended pair info with DSSR-specific fields."""
    interBase_angle: Optional[float] = None
    planarity: Optional[float] = None
    saenger: Optional[str] = None
    name: Optional[str] = None


def load_dssr_pairs(dssr_path: Path, lw_filter: str = "cWW") -> List[PairInfo]:
    """Load pairs from DSSR JSON file, filtered by LW class."""
    if not dssr_path.exists():
        return []

    with open(dssr_path) as f:
        data = json.load(f)

    pairs = []
    for p in data.get("pairs", []):
        lw = p.get("LW", "")
        if lw_filter and lw != lw_filter:
            continue

        res_id1 = normalize_dssr_nt(p.get("nt1", ""))
        res_id2 = normalize_dssr_nt(p.get("nt2", ""))
        bp = p.get("bp", "-")
        bp_type = bp.replace("-", "") if bp else ""

        pairs.append(DSSRPairInfo(
            res_id1=res_id1,
            res_id2=res_id2,
            bp_type=bp_type,
            lw_class=lw,
            n1n9_dist=p.get("N1N9_dist"),
            hbonds_desc=p.get("hbonds_desc"),
            source="dssr",
            interBase_angle=p.get("interBase_angle"),
            planarity=p.get("planarity"),
            saenger=p.get("Saenger"),
            name=p.get("name"),
        ))

    return pairs


def load_legacy_pairs(legacy_path: Path, frame_calc_path: Path) -> List[PairInfo]:
    """Load pairs from legacy JSON file.

    Legacy uses base_i/base_j indices, so we need the frame_calc file to
    map indices to res_ids.
    """
    if not legacy_path.exists():
        return []

    # Build index to res_id mapping from frame_calc file
    idx_to_res_id = {}
    if frame_calc_path.exists():
        with open(frame_calc_path) as f:
            frames = json.load(f)
        for frame in frames:
            idx = frame.get("residue_idx")
            if idx is not None:
                chain = frame.get("chain_id", "?")
                base = frame.get("residue_name", "?").strip()
                num = frame.get("residue_seq", 0)
                res_id = f"{chain}-{base}-{num}"
                idx_to_res_id[idx] = res_id

    with open(legacy_path) as f:
        data = json.load(f)

    pairs = []
    for p in data:
        base_i = p.get("base_i")
        base_j = p.get("base_j")
        bp_type = p.get("bp_type", "")

        res_id1 = idx_to_res_id.get(base_i, f"?-?-{base_i}")
        res_id2 = idx_to_res_id.get(base_j, f"?-?-{base_j}")

        pairs.append(PairInfo(
            res_id1=res_id1,
            res_id2=res_id2,
            bp_type=bp_type,
            source="legacy"
        ))

    return pairs


def load_modern_pairs(modern_path: Path) -> List[PairInfo]:
    """Load pairs from modern JSON file."""
    if not modern_path.exists():
        return []

    with open(modern_path) as f:
        data = json.load(f)

    pairs = []
    for p in data:
        res_id1 = p.get("res_id_i", "")
        res_id2 = p.get("res_id_j", "")
        bp_type = p.get("bp_type", "")

        pairs.append(PairInfo(
            res_id1=res_id1,
            res_id2=res_id2,
            bp_type=bp_type,
            source="modern"
        ))

    return pairs


def load_modern_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[dict]]:
    """Load H-bonds from modern JSON, keyed by (res_id_i, res_id_j)."""
    if not hbond_path.exists():
        return {}

    with open(hbond_path) as f:
        data = json.load(f)

    hbonds = {}
    for entry in data:
        key = (entry.get("res_id_i", ""), entry.get("res_id_j", ""))
        hbonds[key] = entry.get("hbonds", [])

    return hbonds


def compare_pairs(
    dssr_pairs: List[PairInfo],
    legacy_pairs: List[PairInfo],
    modern_pairs: List[PairInfo],
) -> Dict[str, Set[Tuple[str, str]]]:
    """Compare pairs across the three sources.

    Returns dict with keys:
        - "dssr_only": pairs in DSSR but not in legacy or modern
        - "legacy_only": pairs in legacy but not in DSSR
        - "modern_only": pairs in modern but not in DSSR
        - "dssr_and_legacy": pairs in both DSSR and legacy
        - "dssr_and_modern": pairs in both DSSR and modern
        - "all_three": pairs in all three
        - "legacy_not_modern": pairs in legacy but not modern
        - "modern_not_legacy": pairs in modern but not legacy
    """
    dssr_keys = {p.pair_key() for p in dssr_pairs}
    legacy_keys = {p.pair_key() for p in legacy_pairs}
    modern_keys = {p.pair_key() for p in modern_pairs}

    return {
        "dssr_only": dssr_keys - legacy_keys - modern_keys,
        "legacy_only": legacy_keys - dssr_keys,
        "modern_only": modern_keys - dssr_keys,
        "dssr_and_legacy": dssr_keys & legacy_keys,
        "dssr_and_modern": dssr_keys & modern_keys,
        "all_three": dssr_keys & legacy_keys & modern_keys,
        "legacy_not_modern": legacy_keys - modern_keys,
        "modern_not_legacy": modern_keys - legacy_keys,
        "dssr_total": dssr_keys,
        "legacy_total": legacy_keys,
        "modern_total": modern_keys,
    }


def analyze_pdb(
    pdb_id: str,
    dssr_dir: Path,
    legacy_dir: Path,
    modern_dir: Path,
    frame_calc_dir: Path,
    lw_filter: str = "cWW",
    verbose: bool = False,
) -> Dict:
    """Analyze a single PDB for cWW pair coverage."""
    # Load pairs from each source
    dssr_pairs = load_dssr_pairs(dssr_dir / f"{pdb_id}.json", lw_filter)
    legacy_pairs = load_legacy_pairs(
        legacy_dir / f"{pdb_id}.json",
        frame_calc_dir / f"{pdb_id}.json"
    )
    modern_pairs = load_modern_pairs(modern_dir / f"{pdb_id}.json")

    # Compare
    comparison = compare_pairs(dssr_pairs, legacy_pairs, modern_pairs)

    result = {
        "pdb_id": pdb_id,
        "dssr_cww_count": len(dssr_pairs),
        "legacy_pair_count": len(legacy_pairs),
        "modern_pair_count": len(modern_pairs),
        "dssr_found_in_legacy": len(comparison["dssr_and_legacy"]),
        "dssr_found_in_modern": len(comparison["dssr_and_modern"]),
        "dssr_only": len(comparison["dssr_only"]),
        "legacy_only": len(comparison["legacy_only"]),
        "all_three": len(comparison["all_three"]),
    }

    if verbose and comparison["dssr_only"]:
        result["dssr_only_details"] = []
        for key in sorted(comparison["dssr_only"]):
            # Find the original pair info
            for p in dssr_pairs:
                if p.pair_key() == key:
                    detail = {
                        "res_id1": p.res_id1,
                        "res_id2": p.res_id2,
                        "bp_type": p.bp_type,
                        "n1n9_dist": p.n1n9_dist,
                        "hbonds_desc": p.hbonds_desc,
                    }
                    # Add DSSR-specific fields if available
                    if hasattr(p, "interBase_angle"):
                        detail["interBase_angle"] = p.interBase_angle
                    if hasattr(p, "planarity"):
                        detail["planarity"] = p.planarity
                    if hasattr(p, "saenger"):
                        detail["saenger"] = p.saenger
                    if hasattr(p, "name"):
                        detail["name"] = p.name
                    result["dssr_only_details"].append(detail)
                    break

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Compare cWW pairs across DSSR, legacy, and modern X3DNA"
    )
    parser.add_argument(
        "--dssr-dir", type=Path,
        default=Path("data/json_dssr"),
        help="Directory with DSSR JSON files"
    )
    parser.add_argument(
        "--legacy-dir", type=Path,
        default=Path("data/json_legacy/base_pair"),
        help="Directory with legacy base_pair JSON files"
    )
    parser.add_argument(
        "--modern-dir", type=Path,
        default=Path("data/json/base_pair"),
        help="Directory with modern base_pair JSON files"
    )
    parser.add_argument(
        "--frame-calc-dir", type=Path,
        default=Path("data/json_legacy/frame_calc"),
        help="Directory with frame_calc JSON files (for res_id mapping)"
    )
    parser.add_argument(
        "--lw-class", type=str, default="cWW",
        help="LW class to filter from DSSR (default: cWW)"
    )
    parser.add_argument(
        "--pdb", type=str, default=None,
        help="Specific PDB to analyze"
    )
    parser.add_argument(
        "--max-pdbs", type=int, default=100,
        help="Maximum PDBs to analyze"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )

    args = parser.parse_args()

    if args.pdb:
        # Analyze single PDB
        result = analyze_pdb(
            args.pdb,
            args.dssr_dir,
            args.legacy_dir,
            args.modern_dir,
            args.frame_calc_dir,
            args.lw_class,
            args.verbose
        )
        print(f"\n{args.pdb}:")
        print(f"  DSSR {args.lw_class}: {result['dssr_cww_count']}")
        print(f"  Legacy pairs: {result['legacy_pair_count']}")
        print(f"  Modern pairs: {result['modern_pair_count']}")
        print(f"  DSSR found in legacy: {result['dssr_found_in_legacy']}")
        print(f"  DSSR found in modern: {result['dssr_found_in_modern']}")
        print(f"  DSSR only: {result['dssr_only']}")
        return

    # Analyze multiple PDBs
    dssr_files = sorted(args.dssr_dir.glob("*.json"))[:args.max_pdbs]

    total_stats = defaultdict(int)
    pdbs_with_gaps = []

    for dssr_path in dssr_files:
        pdb_id = dssr_path.stem
        result = analyze_pdb(
            pdb_id,
            args.dssr_dir,
            args.legacy_dir,
            args.modern_dir,
            args.frame_calc_dir,
            args.lw_class,
            verbose=False
        )

        total_stats["dssr_total"] += result["dssr_cww_count"]
        total_stats["legacy_total"] += result["legacy_pair_count"]
        total_stats["modern_total"] += result["modern_pair_count"]
        total_stats["dssr_in_legacy"] += result["dssr_found_in_legacy"]
        total_stats["dssr_in_modern"] += result["dssr_found_in_modern"]
        total_stats["dssr_only"] += result["dssr_only"]
        total_stats["all_three"] += result["all_three"]

        if result["dssr_only"] > 0:
            pdbs_with_gaps.append((pdb_id, result["dssr_only"], result["dssr_cww_count"]))

    # Print summary
    print(f"\n{'='*60}")
    print(f"cWW PAIR COVERAGE ANALYSIS ({len(dssr_files)} PDBs)")
    print(f"{'='*60}")

    print(f"\nTotal pairs:")
    print(f"  DSSR cWW:     {total_stats['dssr_total']:>6}")
    print(f"  Legacy:       {total_stats['legacy_total']:>6}")
    print(f"  Modern:       {total_stats['modern_total']:>6}")

    print(f"\nDSSR coverage:")
    if total_stats['dssr_total'] > 0:
        legacy_pct = 100 * total_stats['dssr_in_legacy'] / total_stats['dssr_total']
        modern_pct = 100 * total_stats['dssr_in_modern'] / total_stats['dssr_total']
        print(f"  Found in legacy: {total_stats['dssr_in_legacy']:>6} ({legacy_pct:.1f}%)")
        print(f"  Found in modern: {total_stats['dssr_in_modern']:>6} ({modern_pct:.1f}%)")
        print(f"  DSSR only:       {total_stats['dssr_only']:>6} ({100-legacy_pct:.1f}%)")

    print(f"\nFound in all three: {total_stats['all_three']}")

    if pdbs_with_gaps and args.verbose:
        print(f"\n{'='*60}")
        print("PDBs with DSSR-only pairs (not in legacy/modern):")
        print(f"{'='*60}")
        for pdb_id, gap, total in sorted(pdbs_with_gaps, key=lambda x: -x[1])[:20]:
            print(f"  {pdb_id}: {gap}/{total} cWW pairs only in DSSR")


def analyze_dssr_only_pairs(
    dssr_dir: Path,
    legacy_dir: Path,
    modern_dir: Path,
    frame_calc_dir: Path,
    max_pdbs: int = 100,
    lw_filter: str = "cWW",
) -> None:
    """Analyze DSSR-only pairs in detail to understand why they're missed."""
    dssr_files = sorted(dssr_dir.glob("*.json"))[:max_pdbs]

    all_dssr_only = []

    for dssr_path in dssr_files:
        pdb_id = dssr_path.stem
        result = analyze_pdb(
            pdb_id, dssr_dir, legacy_dir, modern_dir, frame_calc_dir,
            lw_filter, verbose=True
        )
        if "dssr_only_details" in result:
            for detail in result["dssr_only_details"]:
                detail["pdb_id"] = pdb_id
                all_dssr_only.append(detail)

    if not all_dssr_only:
        print("No DSSR-only pairs found!")
        return

    # Analyze by bp_type
    bp_type_counts = defaultdict(int)
    for d in all_dssr_only:
        bp_type_counts[d["bp_type"]] += 1

    print(f"\n{'='*60}")
    print(f"DSSR-ONLY PAIRS ANALYSIS ({len(all_dssr_only)} pairs)")
    print(f"{'='*60}")

    print("\nBy base pair type:")
    for bp_type, count in sorted(bp_type_counts.items(), key=lambda x: -x[1]):
        print(f"  {bp_type}: {count}")

    # Show N1N9 distance distribution
    n1n9_dists = [d["n1n9_dist"] for d in all_dssr_only if d.get("n1n9_dist")]
    if n1n9_dists:
        import numpy as np
        print(f"\nN1N9 distance statistics:")
        print(f"  Mean: {np.mean(n1n9_dists):.2f} Å")
        print(f"  Min:  {np.min(n1n9_dists):.2f} Å")
        print(f"  Max:  {np.max(n1n9_dists):.2f} Å")

    # Analyze by interBase_angle
    high_angle = [d for d in all_dssr_only if d.get("interBase_angle", 0) > 30]
    low_angle = [d for d in all_dssr_only if d.get("interBase_angle", 0) <= 30]

    print(f"\nBy interBase_angle:")
    print(f"  High (>30°): {len(high_angle)} pairs - likely rejected by plane angle check")
    print(f"  Low (≤30°): {len(low_angle)} pairs - need investigation")

    # Show low angle pairs (these are the ones we should be finding)
    if low_angle:
        print(f"\nLow-angle DSSR-only pairs (should match):")
        for d in low_angle[:15]:
            n1n9 = d.get("n1n9_dist", 0)
            angle = d.get("interBase_angle", 0)
            saenger = d.get("saenger", "?")
            name = d.get("name", "?")
            print(f"  {d['pdb_id']}: {d['res_id1']} - {d['res_id2']} ({d['bp_type']}) "
                  f"N1N9={n1n9:.2f}Å angle={angle:.1f}° {name} {saenger}")

    # Show high angle pairs
    if high_angle:
        print(f"\nHigh-angle DSSR-only pairs (correctly rejected):")
        for d in high_angle[:10]:
            n1n9 = d.get("n1n9_dist", 0)
            angle = d.get("interBase_angle", 0)
            name = d.get("name", "?")
            print(f"  {d['pdb_id']}: {d['res_id1']} - {d['res_id2']} ({d['bp_type']}) "
                  f"N1N9={n1n9:.2f}Å angle={angle:.1f}° {name}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--analyze-dssr-only":
        # Run detailed analysis
        analyze_dssr_only_pairs(
            Path("data/json_dssr"),
            Path("data/json_legacy/base_pair"),
            Path("data/json/base_pair"),
            Path("data/json_legacy/frame_calc"),
            max_pdbs=100,
        )
    else:
        main()
