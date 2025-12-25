#!/usr/bin/env python3
"""
Compare base pairs found by legacy X3DNA with DSSR output.

Usage:
    python scripts/compare_legacy_dssr_pairs.py [--pdb PDB_ID] [--verbose] [--summary]
    python scripts/compare_legacy_dssr_pairs.py --wc-only   # Only compare Watson-Crick pairs
    python scripts/compare_legacy_dssr_pairs.py --stats     # Show pair type statistics
    python scripts/compare_legacy_dssr_pairs.py --frames    # Compare frames for matching pairs
"""

import argparse
import json
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
from collections import defaultdict, Counter
import numpy as np


@dataclass
class ResidueInfo:
    """Info about a residue."""
    chain: str
    name: str
    seq: int

    def to_dssr_format(self) -> str:
        """Convert to DSSR-style residue ID (e.g., 'A.G10')."""
        # DSSR uses format: chain.residue_name + seq
        # Need to handle DNA prefixes (DC, DG, DA, DT)
        return f"{self.chain}.{self.name.strip()}{self.seq}"


@dataclass
class BasePair:
    """A base pair with residue identifiers."""
    res1: str  # Residue ID in standardized format
    res2: str  # Residue ID in standardized format
    bp_type: str  # Pair type (e.g., "CG", "WC")
    source: str  # "legacy" or "dssr"
    # Frame data (optional)
    origin: Optional[np.ndarray] = None
    x_axis: Optional[np.ndarray] = None
    y_axis: Optional[np.ndarray] = None
    z_axis: Optional[np.ndarray] = None
    # Legacy-specific: individual base frames
    org_i: Optional[np.ndarray] = None
    org_j: Optional[np.ndarray] = None
    orien_i: Optional[np.ndarray] = None
    orien_j: Optional[np.ndarray] = None

    def __hash__(self):
        # Normalize order for comparison
        sorted_pair = tuple(sorted([self.res1, self.res2]))
        return hash(sorted_pair)

    def __eq__(self, other):
        if not isinstance(other, BasePair):
            return False
        self_sorted = tuple(sorted([self.res1, self.res2]))
        other_sorted = tuple(sorted([other.res1, other.res2]))
        return self_sorted == other_sorted

    def key(self) -> tuple:
        """Return a normalized key for comparison."""
        return tuple(sorted([self.res1, self.res2]))

    def is_wc(self) -> bool:
        """Check if this is a Watson-Crick pair."""
        if self.source == "dssr":
            return self.bp_type in ("WC", "cWW")
        else:
            # Legacy uses 2-letter codes
            return self.bp_type.upper() in ("CG", "GC", "AU", "UA", "AT", "TA")

    def is_wobble(self) -> bool:
        """Check if this is a Wobble pair."""
        if self.source == "dssr":
            return self.bp_type == "Wobble"
        else:
            return self.bp_type.upper() in ("GU", "UG", "GT", "TG")


class LegacyDSSRComparator:
    """Compare legacy X3DNA pairs with DSSR pairs."""

    def __init__(self, legacy_dir: Path, dssr_dir: Path):
        self.legacy_dir = legacy_dir
        self.dssr_dir = dssr_dir
        self.base_pair_dir = legacy_dir / "base_pair"
        self.pdb_atoms_dir = legacy_dir / "pdb_atoms"

    def get_common_pdbs(self) -> list[str]:
        """Get list of PDBs that have both legacy and DSSR data."""
        legacy_pdbs = set(p.stem for p in self.base_pair_dir.glob("*.json")
                         if not p.stem.endswith("_minimal_pairs_1_2"))
        dssr_pdbs = set(p.stem for p in self.dssr_dir.glob("*.json"))
        return sorted(legacy_pdbs & dssr_pdbs)

    def build_residue_map(self, pdb_id: str) -> dict[int, ResidueInfo]:
        """Build mapping from residue index to residue info."""
        pdb_atoms_file = self.pdb_atoms_dir / f"{pdb_id}.json"
        if not pdb_atoms_file.exists():
            return {}

        try:
            with open(pdb_atoms_file) as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"Warning: Failed to parse {pdb_atoms_file}: {e}")
            return {}

        # Get atoms list
        atoms = data[0]["atoms"] if isinstance(data, list) else data["atoms"]

        # Group atoms by residue (using first atom of each residue)
        residue_map = {}
        current_res_key = None
        current_res_idx = 0

        for atom in atoms:
            res_key = (atom["chain_id"], atom["residue_seq"])
            if res_key != current_res_key:
                current_res_key = res_key
                current_res_idx += 1
                residue_map[current_res_idx] = ResidueInfo(
                    chain=atom["chain_id"],
                    name=atom["residue_name"].strip(),
                    seq=atom["residue_seq"]
                )

        return residue_map

    def load_legacy_pairs(self, pdb_id: str) -> list[BasePair]:
        """Load base pairs from legacy JSON."""
        base_pair_file = self.base_pair_dir / f"{pdb_id}.json"
        if not base_pair_file.exists():
            return []

        try:
            with open(base_pair_file) as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"Warning: Failed to parse {base_pair_file}: {e}")
            return []

        # Build residue map
        res_map = self.build_residue_map(pdb_id)
        if not res_map:
            return []

        pairs = []
        for item in data:
            if item.get("type") != "base_pair":
                continue

            base_i = item["base_i"]
            base_j = item["base_j"]
            bp_type = item.get("bp_type", "")

            if base_i not in res_map or base_j not in res_map:
                continue

            res1 = res_map[base_i].to_dssr_format()
            res2 = res_map[base_j].to_dssr_format()

            # Extract frame data
            org_i = np.array(item.get("org_i")) if "org_i" in item else None
            org_j = np.array(item.get("org_j")) if "org_j" in item else None
            orien_i = np.array(item.get("orien_i")) if "orien_i" in item else None
            orien_j = np.array(item.get("orien_j")) if "orien_j" in item else None

            # Compute midpoint origin
            origin = None
            if org_i is not None and org_j is not None:
                origin = (org_i + org_j) / 2

            pairs.append(BasePair(
                res1, res2, bp_type, "legacy",
                origin=origin,
                org_i=org_i, org_j=org_j,
                orien_i=orien_i, orien_j=orien_j
            ))

        return pairs

    def load_dssr_pairs(self, pdb_id: str) -> list[BasePair]:
        """Load base pairs from DSSR JSON."""
        dssr_file = self.dssr_dir / f"{pdb_id}.json"
        if not dssr_file.exists():
            return []

        with open(dssr_file) as f:
            content = f.read().strip()
            if not content:
                return []
            data = json.loads(content)

        if "pairs" not in data:
            return []

        pairs = []
        for item in data["pairs"]:
            nt1 = item["nt1"]
            nt2 = item["nt2"]
            name = item.get("name", "")

            # Extract frame data
            origin = None
            x_axis = None
            y_axis = None
            z_axis = None

            if "frame" in item and item["frame"]:
                frame = item["frame"]
                if frame.get("origin"):
                    origin = np.array(frame["origin"])
                if frame.get("x_axis"):
                    x_axis = np.array(frame["x_axis"])
                if frame.get("y_axis"):
                    y_axis = np.array(frame["y_axis"])
                if frame.get("z_axis"):
                    z_axis = np.array(frame["z_axis"])

            pairs.append(BasePair(
                nt1, nt2, name, "dssr",
                origin=origin,
                x_axis=x_axis, y_axis=y_axis, z_axis=z_axis
            ))

        return pairs

    def compare_pdb(self, pdb_id: str, verbose: bool = False,
                    wc_only: bool = False, include_wobble: bool = False) -> dict:
        """Compare pairs for a single PDB."""
        legacy_pairs = self.load_legacy_pairs(pdb_id)
        dssr_pairs = self.load_dssr_pairs(pdb_id)

        # Filter for WC pairs if requested
        if wc_only:
            if include_wobble:
                legacy_pairs = [p for p in legacy_pairs if p.is_wc() or p.is_wobble()]
                dssr_pairs = [p for p in dssr_pairs if p.is_wc() or p.is_wobble()]
            else:
                legacy_pairs = [p for p in legacy_pairs if p.is_wc()]
                dssr_pairs = [p for p in dssr_pairs if p.is_wc()]

        # Convert to sets of keys for comparison
        legacy_keys = {p.key() for p in legacy_pairs}
        dssr_keys = {p.key() for p in dssr_pairs}

        common = legacy_keys & dssr_keys
        legacy_only = legacy_keys - dssr_keys
        dssr_only = dssr_keys - legacy_keys

        result = {
            "pdb_id": pdb_id,
            "legacy_count": len(legacy_pairs),
            "dssr_count": len(dssr_pairs),
            "common": len(common),
            "legacy_only": len(legacy_only),
            "dssr_only": len(dssr_only),
            "legacy_only_pairs": list(legacy_only),
            "dssr_only_pairs": list(dssr_only),
        }

        if verbose:
            print(f"\n{'='*60}")
            print(f"PDB: {pdb_id}")
            print(f"{'='*60}")
            print(f"Legacy pairs: {len(legacy_pairs)}")
            print(f"DSSR pairs:   {len(dssr_pairs)}")
            print(f"Common:       {len(common)}")
            print(f"Legacy only:  {len(legacy_only)}")
            print(f"DSSR only:    {len(dssr_only)}")

            if legacy_only:
                print(f"\nLegacy-only pairs:")
                # Find the BasePair objects for legacy-only pairs
                for pair in legacy_pairs:
                    if pair.key() in legacy_only:
                        print(f"  {pair.res1} - {pair.res2} ({pair.bp_type})")

            if dssr_only:
                print(f"\nDSSR-only pairs:")
                for pair in dssr_pairs:
                    if pair.key() in dssr_only:
                        print(f"  {pair.res1} - {pair.res2} ({pair.bp_type})")

        return result

    def compare_all(self, verbose: bool = False, max_pdbs: Optional[int] = None,
                    wc_only: bool = False, include_wobble: bool = False,
                    show_progress: bool = True) -> dict:
        """Compare all common PDBs."""
        pdbs = self.get_common_pdbs()
        if max_pdbs:
            pdbs = pdbs[:max_pdbs]

        results = []
        totals = {
            "legacy_count": 0,
            "dssr_count": 0,
            "common": 0,
            "legacy_only": 0,
            "dssr_only": 0,
        }

        for i, pdb_id in enumerate(pdbs):
            if show_progress and not verbose and (i + 1) % 500 == 0:
                print(f"  Processed {i + 1}/{len(pdbs)} PDBs...", file=sys.stderr)

            result = self.compare_pdb(pdb_id, verbose=verbose,
                                      wc_only=wc_only, include_wobble=include_wobble)
            results.append(result)

            for key in totals:
                totals[key] += result[key]

        return {
            "num_pdbs": len(pdbs),
            "totals": totals,
            "results": results,
        }

    def compare_frames(self, pdb_id: str, verbose: bool = False) -> dict:
        """Compare frames for matching pairs between legacy and DSSR."""
        legacy_pairs = self.load_legacy_pairs(pdb_id)
        dssr_pairs = self.load_dssr_pairs(pdb_id)

        # Build lookup by key
        legacy_by_key = {p.key(): p for p in legacy_pairs}
        dssr_by_key = {p.key(): p for p in dssr_pairs}

        # Find matching pairs
        common_keys = set(legacy_by_key.keys()) & set(dssr_by_key.keys())

        origin_diffs = []
        pairs_compared = 0

        for key in common_keys:
            leg = legacy_by_key[key]
            dssr = dssr_by_key[key]

            if leg.origin is not None and dssr.origin is not None:
                diff = np.linalg.norm(leg.origin - dssr.origin)
                origin_diffs.append(diff)
                pairs_compared += 1

                if verbose and diff > 0.1:  # Show large differences
                    print(f"  {leg.res1} - {leg.res2}: origin diff = {diff:.4f} Å")

        result = {
            "pdb_id": pdb_id,
            "pairs_compared": pairs_compared,
            "origin_diffs": origin_diffs,
        }

        if origin_diffs:
            result["mean_origin_diff"] = np.mean(origin_diffs)
            result["max_origin_diff"] = np.max(origin_diffs)
            result["min_origin_diff"] = np.min(origin_diffs)

        return result

    def compare_frames_all(self, max_pdbs: Optional[int] = None,
                           verbose: bool = False) -> dict:
        """Compare frames for all common PDBs."""
        pdbs = self.get_common_pdbs()
        if max_pdbs:
            pdbs = pdbs[:max_pdbs]

        all_origin_diffs = []
        results = []

        for i, pdb_id in enumerate(pdbs):
            if not verbose and (i + 1) % 500 == 0:
                print(f"  Processed {i + 1}/{len(pdbs)} PDBs...", file=sys.stderr)

            result = self.compare_frames(pdb_id, verbose=verbose)
            results.append(result)
            all_origin_diffs.extend(result.get("origin_diffs", []))

        summary = {
            "num_pdbs": len(pdbs),
            "total_pairs_compared": len(all_origin_diffs),
        }

        if all_origin_diffs:
            summary["mean_origin_diff"] = float(np.mean(all_origin_diffs))
            summary["median_origin_diff"] = float(np.median(all_origin_diffs))
            summary["max_origin_diff"] = float(np.max(all_origin_diffs))
            summary["std_origin_diff"] = float(np.std(all_origin_diffs))

            # Count pairs within thresholds
            diffs = np.array(all_origin_diffs)
            summary["within_0.01A"] = int(np.sum(diffs < 0.01))
            summary["within_0.1A"] = int(np.sum(diffs < 0.1))
            summary["within_0.5A"] = int(np.sum(diffs < 0.5))
            summary["within_1.0A"] = int(np.sum(diffs < 1.0))

        return {
            "summary": summary,
            "results": results,
        }

    def get_pair_type_stats(self, max_pdbs: Optional[int] = None) -> dict:
        """Get statistics on pair types."""
        pdbs = self.get_common_pdbs()
        if max_pdbs:
            pdbs = pdbs[:max_pdbs]

        legacy_types = Counter()
        dssr_types = Counter()

        for pdb_id in pdbs:
            legacy_pairs = self.load_legacy_pairs(pdb_id)
            dssr_pairs = self.load_dssr_pairs(pdb_id)

            for p in legacy_pairs:
                legacy_types[p.bp_type] += 1
            for p in dssr_pairs:
                dssr_types[p.bp_type] += 1

        return {
            "legacy": dict(legacy_types.most_common(30)),
            "dssr": dict(dssr_types.most_common(30)),
        }


def print_summary(comparison: dict):
    """Print a summary of comparison results."""
    totals = comparison["totals"]

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"PDBs compared:        {comparison['num_pdbs']}")
    print(f"Total legacy pairs:   {totals['legacy_count']}")
    print(f"Total DSSR pairs:     {totals['dssr_count']}")
    print(f"Total common:         {totals['common']}")
    print(f"Total legacy-only:    {totals['legacy_only']}")
    print(f"Total DSSR-only:      {totals['dssr_only']}")

    if totals['legacy_count'] > 0:
        legacy_match_rate = totals['common'] / totals['legacy_count'] * 100
        print(f"\nLegacy match rate:    {legacy_match_rate:.1f}% of legacy pairs found in DSSR")

    if totals['dssr_count'] > 0:
        dssr_match_rate = totals['common'] / totals['dssr_count'] * 100
        print(f"DSSR match rate:      {dssr_match_rate:.1f}% of DSSR pairs found in legacy")

    # Find PDbs with biggest differences
    results = comparison["results"]

    # PDbs where legacy has more pairs
    legacy_more = sorted(results, key=lambda x: x["legacy_only"], reverse=True)[:10]
    print(f"\nTop 10 PDBs where legacy has more pairs:")
    for r in legacy_more:
        if r["legacy_only"] > 0:
            print(f"  {r['pdb_id']}: +{r['legacy_only']} legacy-only (legacy={r['legacy_count']}, dssr={r['dssr_count']})")

    # PDbs where DSSR has more pairs
    dssr_more = sorted(results, key=lambda x: x["dssr_only"], reverse=True)[:10]
    print(f"\nTop 10 PDBs where DSSR has more pairs:")
    for r in dssr_more:
        if r["dssr_only"] > 0:
            print(f"  {r['pdb_id']}: +{r['dssr_only']} dssr-only (legacy={r['legacy_count']}, dssr={r['dssr_count']})")


def print_pair_type_stats(stats: dict):
    """Print pair type statistics."""
    print(f"\n{'='*60}")
    print("PAIR TYPE STATISTICS")
    print(f"{'='*60}")

    print("\nLegacy pair types (top 30):")
    for bp_type, count in stats["legacy"].items():
        print(f"  {bp_type:10s}: {count:6d}")

    print("\nDSSR pair types (top 30):")
    for bp_type, count in stats["dssr"].items():
        print(f"  {bp_type:15s}: {count:6d}")


def print_frame_comparison(comparison: dict):
    """Print frame comparison results."""
    summary = comparison["summary"]

    print(f"\n{'='*60}")
    print("FRAME COMPARISON (Origin Distance)")
    print(f"{'='*60}")
    print(f"PDBs compared:           {summary['num_pdbs']}")
    print(f"Total pairs compared:    {summary['total_pairs_compared']}")

    if summary['total_pairs_compared'] > 0:
        print(f"\nOrigin distance statistics:")
        print(f"  Mean:   {summary['mean_origin_diff']:.6f} Å")
        print(f"  Median: {summary['median_origin_diff']:.6f} Å")
        print(f"  Std:    {summary['std_origin_diff']:.6f} Å")
        print(f"  Max:    {summary['max_origin_diff']:.6f} Å")

        total = summary['total_pairs_compared']
        print(f"\nPairs within threshold:")
        print(f"  < 0.01 Å: {summary['within_0.01A']:6d} ({100*summary['within_0.01A']/total:.1f}%)")
        print(f"  < 0.10 Å: {summary['within_0.1A']:6d} ({100*summary['within_0.1A']/total:.1f}%)")
        print(f"  < 0.50 Å: {summary['within_0.5A']:6d} ({100*summary['within_0.5A']/total:.1f}%)")
        print(f"  < 1.00 Å: {summary['within_1.0A']:6d} ({100*summary['within_1.0A']/total:.1f}%)")


def main():
    parser = argparse.ArgumentParser(description="Compare legacy X3DNA pairs with DSSR pairs")
    parser.add_argument("--pdb", type=str, help="Compare a specific PDB")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show detailed output")
    parser.add_argument("--summary", "-s", action="store_true", help="Show summary only")
    parser.add_argument("--max-pdbs", type=int, help="Maximum number of PDBs to compare")
    parser.add_argument("--output", "-o", type=str, help="Save results to JSON file")
    parser.add_argument("--wc-only", action="store_true", help="Only compare Watson-Crick pairs")
    parser.add_argument("--include-wobble", action="store_true", help="Include wobble pairs with --wc-only")
    parser.add_argument("--stats", action="store_true", help="Show pair type statistics")
    parser.add_argument("--frames", action="store_true", help="Compare frames for matching pairs")

    args = parser.parse_args()

    # Find data directories
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
    legacy_dir = project_dir / "data" / "json_legacy"
    dssr_dir = project_dir / "data" / "json_dssr"

    comparator = LegacyDSSRComparator(legacy_dir, dssr_dir)

    if args.stats:
        stats = comparator.get_pair_type_stats(max_pdbs=args.max_pdbs)
        print_pair_type_stats(stats)
        return

    if args.frames:
        if args.pdb:
            result = comparator.compare_frames(args.pdb, verbose=True)
            print(f"\nPDB: {args.pdb}")
            print(f"Pairs compared: {result['pairs_compared']}")
            if result['pairs_compared'] > 0:
                print(f"Mean origin diff: {result['mean_origin_diff']:.6f} Å")
                print(f"Max origin diff:  {result['max_origin_diff']:.6f} Å")
        else:
            comparison = comparator.compare_frames_all(
                max_pdbs=args.max_pdbs,
                verbose=args.verbose
            )
            print_frame_comparison(comparison)
        return

    if args.pdb:
        result = comparator.compare_pdb(args.pdb, verbose=True,
                                        wc_only=args.wc_only,
                                        include_wobble=args.include_wobble)
        print(f"\nResult: {result['common']}/{result['legacy_count']} legacy pairs matched")
    else:
        if args.wc_only:
            filter_desc = "WC+Wobble" if args.include_wobble else "WC only"
            print(f"Filtering for {filter_desc} pairs...")

        comparison = comparator.compare_all(
            verbose=args.verbose,
            max_pdbs=args.max_pdbs,
            wc_only=args.wc_only,
            include_wobble=args.include_wobble
        )

        if args.output:
            # Remove the actual pair lists for the JSON output (too verbose)
            for r in comparison["results"]:
                r["legacy_only_pairs"] = len(r["legacy_only_pairs"])
                r["dssr_only_pairs"] = len(r["dssr_only_pairs"])

            with open(args.output, "w") as f:
                json.dump(comparison, f, indent=2)
            print(f"Results saved to {args.output}")

        if args.summary or not args.verbose:
            print_summary(comparison)


if __name__ == "__main__":
    main()
