#!/usr/bin/env python3
"""
Compare helix organization across all test PDBs.
"""

import json
import subprocess
import sys
from pathlib import Path


def run_comparison(pdb_id: str, base_dir: Path) -> dict:
    """Run comparison for a single PDB."""
    legacy_path = base_dir / "data" / "json_legacy" / "helix_organization" / f"{pdb_id}.json"
    modern_path = base_dir / "data" / "json" / "helix_organization" / f"{pdb_id}.json"

    if not legacy_path.exists() or not modern_path.exists():
        return None

    with open(legacy_path) as f:
        legacy = json.load(f)
    with open(modern_path) as f:
        modern = json.load(f)

    # Build pair maps
    legacy_map = {}
    for helix in legacy:
        for pair in helix["pairs"]:
            bp_idx = pair["bp_idx"]
            legacy_map[bp_idx] = {
                "helix_num": helix["helix_num"],
                "helix_pos": pair["helix_pos"],
                "strand_swapped": pair["strand_swapped"]
            }

    modern_map = {}
    for helix in modern:
        for pair in helix["pairs"]:
            bp_idx = pair["bp_idx"]
            modern_map[bp_idx] = {
                "helix_num": helix["helix_num"],
                "helix_pos": pair["helix_pos"],
                "strand_swapped": pair["strand_swapped"]
            }

    total = len(legacy_map)
    helix_match = 0
    swap_match = 0

    for bp_idx in legacy_map:
        if bp_idx not in modern_map:
            continue
        if legacy_map[bp_idx]["helix_num"] == modern_map[bp_idx]["helix_num"]:
            helix_match += 1
        if legacy_map[bp_idx]["strand_swapped"] == modern_map[bp_idx]["strand_swapped"]:
            swap_match += 1

    return {
        "total": total,
        "helix_match": helix_match,
        "helix_pct": 100 * helix_match / total if total > 0 else 0,
        "swap_match": swap_match,
        "swap_pct": 100 * swap_match / total if total > 0 else 0
    }


def main():
    base_dir = Path(__file__).parent.parent

    # Get all legacy helix_organization files
    legacy_dir = base_dir / "data" / "json_legacy" / "helix_organization"
    if not legacy_dir.exists():
        print("No legacy helix_organization files found")
        return

    results = []
    for json_file in sorted(legacy_dir.glob("*.json")):
        pdb_id = json_file.stem
        result = run_comparison(pdb_id, base_dir)
        if result:
            results.append((pdb_id, result))

    print(f"\nHelix Organization Comparison ({len(results)} PDBs)")
    print("=" * 70)
    print(f"{'PDB':<8} {'Total':>6} {'Helix%':>8} {'Swap%':>8} {'Status':<20}")
    print("-" * 70)

    for pdb_id, r in sorted(results, key=lambda x: x[1]["helix_pct"]):
        status = ""
        if r["helix_pct"] < 50:
            status = "HELIX_DIFF"
        elif r["swap_pct"] < 90:
            status = "SWAP_DIFF"
        else:
            status = "OK"
        print(f"{pdb_id:<8} {r['total']:>6} {r['helix_pct']:>7.1f}% {r['swap_pct']:>7.1f}%  {status}")

    # Summary
    perfect_helix = sum(1 for _, r in results if r["helix_pct"] == 100)
    high_helix = sum(1 for _, r in results if r["helix_pct"] >= 90)
    perfect_swap = sum(1 for _, r in results if r["swap_pct"] == 100)
    high_swap = sum(1 for _, r in results if r["swap_pct"] >= 90)

    print("-" * 70)
    print(f"Helix match: {perfect_helix} perfect (100%), {high_helix} high (>=90%)")
    print(f"Swap match: {perfect_swap} perfect (100%), {high_swap} high (>=90%)")


if __name__ == "__main__":
    main()
