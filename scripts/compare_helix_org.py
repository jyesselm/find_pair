#!/usr/bin/env python3
"""
Compare helix organization between legacy and modern JSON outputs.

This script identifies differences in:
1. Helix structure (which pairs are grouped together)
2. Pair ordering within helices
3. Strand swap decisions
"""

import json
import sys
from pathlib import Path
from collections import defaultdict


def load_helix_org(json_path: Path) -> list:
    """Load helix organization JSON."""
    with open(json_path) as f:
        return json.load(f)


def build_pair_map(helices: list) -> dict:
    """Build a map from bp_idx to (helix_num, helix_pos, strand_swapped, base_i, base_j)."""
    pair_map = {}
    for helix in helices:
        helix_num = helix["helix_num"]
        for pair in helix["pairs"]:
            bp_idx = pair["bp_idx"]
            pair_map[bp_idx] = {
                "helix_num": helix_num,
                "helix_pos": pair["helix_pos"],
                "strand_swapped": pair["strand_swapped"],
                "base_i": pair["base_i"],
                "base_j": pair["base_j"]
            }
    return pair_map


def compare_helix_org(legacy_path: Path, modern_path: Path, verbose: bool = True):
    """Compare helix organization between legacy and modern."""
    legacy = load_helix_org(legacy_path)
    modern = load_helix_org(modern_path)

    legacy_map = build_pair_map(legacy)
    modern_map = build_pair_map(modern)

    # Statistics
    stats = {
        "total_pairs": len(legacy_map),
        "helix_match": 0,
        "helix_mismatch": 0,
        "swap_match": 0,
        "swap_mismatch": 0,
        "swap_mismatches": []
    }

    # Compare each pair
    for bp_idx in sorted(legacy_map.keys()):
        leg = legacy_map[bp_idx]
        mod = modern_map.get(bp_idx)

        if not mod:
            if verbose:
                print(f"  bp_idx {bp_idx}: Missing in modern")
            continue

        # Check if in same helix
        if leg["helix_num"] == mod["helix_num"]:
            stats["helix_match"] += 1
        else:
            stats["helix_mismatch"] += 1
            if verbose:
                print(f"  bp_idx {bp_idx}: Helix mismatch - legacy:{leg['helix_num']} vs modern:{mod['helix_num']}")

        # Check strand swap
        if leg["strand_swapped"] == mod["strand_swapped"]:
            stats["swap_match"] += 1
        else:
            stats["swap_mismatch"] += 1
            stats["swap_mismatches"].append({
                "bp_idx": bp_idx,
                "base_i": leg["base_i"],
                "base_j": leg["base_j"],
                "legacy_swap": leg["strand_swapped"],
                "modern_swap": mod["strand_swapped"],
                "legacy_helix": leg["helix_num"],
                "modern_helix": mod["helix_num"],
                "legacy_helix_pos": leg["helix_pos"],
                "modern_helix_pos": mod["helix_pos"]
            })

    return stats


def main():
    if len(sys.argv) < 2:
        print("Usage: python compare_helix_org.py <PDB_ID> [--verbose]")
        sys.exit(1)

    pdb_id = sys.argv[1]
    verbose = "--verbose" in sys.argv or "-v" in sys.argv

    base_dir = Path(__file__).parent.parent
    legacy_path = base_dir / "data" / "json_legacy" / "helix_organization" / f"{pdb_id}.json"
    modern_path = base_dir / "data" / "json" / "helix_organization" / f"{pdb_id}.json"

    if not legacy_path.exists():
        print(f"Error: Legacy file not found: {legacy_path}")
        sys.exit(1)
    if not modern_path.exists():
        print(f"Error: Modern file not found: {modern_path}")
        sys.exit(1)

    print(f"\nComparing helix organization for {pdb_id}")
    print("=" * 60)

    stats = compare_helix_org(legacy_path, modern_path, verbose=False)

    print(f"\nSummary:")
    print(f"  Total pairs: {stats['total_pairs']}")
    print(f"  Helix match: {stats['helix_match']} ({100*stats['helix_match']/stats['total_pairs']:.1f}%)")
    print(f"  Helix mismatch: {stats['helix_mismatch']}")
    print(f"  Swap match: {stats['swap_match']} ({100*stats['swap_match']/stats['total_pairs']:.1f}%)")
    print(f"  Swap mismatch: {stats['swap_mismatch']}")

    if stats["swap_mismatches"] and verbose:
        print(f"\nSwap mismatches ({len(stats['swap_mismatches'])} pairs):")
        for m in stats["swap_mismatches"][:20]:  # Show first 20
            print(f"  bp_idx={m['bp_idx']}: base=({m['base_i']},{m['base_j']}) "
                  f"legacy_swap={m['legacy_swap']} modern_swap={m['modern_swap']} "
                  f"helix={m['legacy_helix']}->{m['modern_helix']} "
                  f"pos={m['legacy_helix_pos']}->{m['modern_helix_pos']}")
        if len(stats["swap_mismatches"]) > 20:
            print(f"  ... and {len(stats['swap_mismatches']) - 20} more")

    return stats


if __name__ == "__main__":
    main()
