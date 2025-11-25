#!/usr/bin/env python3
"""
Analyze differences between legacy and modern hbond_list records.

This script provides detailed analysis of H-bond differences to help
understand why legacy and modern differ.
"""

import json
import sys
from pathlib import Path
from collections import defaultdict

from x3dna_json_compare import JsonComparator


def analyze_hbond_differences(pdb_id: str, project_root: Path = Path(".")):
    """Analyze H-bond differences for a PDB file."""
    legacy_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
    modern_file = project_root / "data" / "json" / f"{pdb_id}.json"
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"

    comparator = JsonComparator(compare_pairs=True, force_recompute=True)
    result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)

    if not result.hbond_list_comparison:
        print(f"No hbond_list comparison found for {pdb_id}")
        return

    hb = result.hbond_list_comparison
    
    if not hb:
        print(f"No hbond_list comparison found for {pdb_id}")
        return
    
    # Handle both object and dict (from cache)
    if isinstance(hb, dict):
        # Reconstruct from dict
        from x3dna_json_compare.hbond_comparison import HBondListComparison
        hb_obj = HBondListComparison()
        hb_obj.total_legacy = hb.get('total_legacy', 0)
        hb_obj.total_modern = hb.get('total_modern', 0)
        hb_obj.common_count = hb.get('common_count', 0)
        hb_obj.missing_in_modern = hb.get('missing_in_modern', [])
        hb_obj.extra_in_modern = hb.get('extra_in_modern', [])
        hb_obj.mismatched_pairs = hb.get('mismatched_pairs', [])
        hb_obj.pair_comparisons = hb.get('pair_comparisons', {})
        hb = hb_obj
    else:
        # It's already an object, convert to dict for easier access
        hb_dict = {
            'total_legacy': hb.total_legacy,
            'total_modern': hb.total_modern,
            'common_count': hb.common_count,
            'missing_in_modern': hb.missing_in_modern,
            'extra_in_modern': hb.extra_in_modern,
            'mismatched_pairs': hb.mismatched_pairs,
            'pair_comparisons': hb.pair_comparisons
        }

    print(f"\n{'='*80}")
    print(f"H-bond Analysis for {pdb_id}")
    print(f"{'='*80}\n")

    print(f"Summary:")
    print(f"  Legacy pairs: {hb.total_legacy if hasattr(hb, 'total_legacy') else hb_dict['total_legacy']}")
    print(f"  Modern pairs: {hb.total_modern if hasattr(hb, 'total_modern') else hb_dict['total_modern']}")
    print(f"  Common pairs: {hb.common_count if hasattr(hb, 'common_count') else hb_dict['common_count']}")
    missing = hb.missing_in_modern if hasattr(hb, 'missing_in_modern') else hb_dict['missing_in_modern']
    extra = hb.extra_in_modern if hasattr(hb, 'extra_in_modern') else hb_dict['extra_in_modern']
    mismatched = hb.mismatched_pairs if hasattr(hb, 'mismatched_pairs') else hb_dict['mismatched_pairs']
    print(f"  Missing in modern: {len(missing)}")
    print(f"  Extra in modern: {len(extra)}")
    print(f"  Mismatched pairs: {len(mismatched)}\n")

    # Use dict for easier access
    missing = hb.missing_in_modern if hasattr(hb, 'missing_in_modern') else hb_dict['missing_in_modern']
    extra = hb.extra_in_modern if hasattr(hb, 'extra_in_modern') else hb_dict['extra_in_modern']
    mismatched = hb.mismatched_pairs if hasattr(hb, 'mismatched_pairs') else hb_dict['mismatched_pairs']
    pair_comps = hb.pair_comparisons if hasattr(hb, 'pair_comparisons') else hb_dict['pair_comparisons']
    
    # Analyze missing pairs
    if missing:
        print(f"{'='*80}")
        print(f"Missing Pairs in Modern ({len(missing)} total)")
        print(f"{'='*80}\n")

        for pair_info in missing:
            base_i = pair_info["base_i"]
            base_j = pair_info["base_j"]
            leg_rec = pair_info["legacy_record"]
            num_hbonds = pair_info["num_hbonds"]

            print(f"Pair ({base_i}, {base_j}): {num_hbonds} H-bonds in legacy")
            hbonds = leg_rec.get("hbonds", [])

            # Categorize H-bonds by type
            by_type = defaultdict(list)
            for hb in hbonds:
                hb_type = hb.get("type", " ")
                by_type[hb_type].append(hb)

            print(f"  H-bond types breakdown:")
            for hb_type, type_list in sorted(by_type.items()):
                type_name = {"-": "standard", "*": "non-standard", " ": "invalid"}.get(
                    hb_type, f"unknown ({hb_type})"
                )
                print(f"    {type_name} (type='{hb_type}'): {len(type_list)}")

            # Show sample H-bonds
            print(f"  Sample H-bonds (showing all):")
            for hb in hbonds:
                donor = hb.get("donor_atom", "").strip()
                acceptor = hb.get("acceptor_atom", "").strip()
                dist = hb.get("distance", 0.0)
                hb_type = hb.get("type", " ")
                print(
                    f"    {donor} -> {acceptor}: {dist:.3f} Å (type='{hb_type}')"
                )
            print()

    # Analyze mismatched pairs
    if mismatched:
        print(f"{'='*80}")
        print(f"Mismatched Pairs ({len(mismatched)} total)")
        print(f"{'='*80}\n")

        # Categorize mismatches
        missing_only = []  # Only missing H-bonds
        extra_only = []  # Only extra H-bonds
        count_mismatch = []  # Different counts

        for pair_info in mismatched:
            comp = pair_info["comparison"]
            if comp.missing_in_modern and not comp.extra_in_modern:
                missing_only.append(pair_info)
            elif comp.extra_in_modern and not comp.missing_in_modern:
                extra_only.append(pair_info)
            elif comp.legacy_num_hbonds != comp.modern_num_hbonds:
                count_mismatch.append(pair_info)

        print(f"Mismatch categories:")
        print(f"  Only missing H-bonds: {len(missing_only)}")
        print(f"  Only extra H-bonds: {len(extra_only)}")
        print(f"  Count differences: {len(count_mismatch)}\n")

        # Analyze missing H-bonds pattern
        if missing_only:
            print(f"{'-'*80}")
            print(f"Pairs with only missing H-bonds ({len(missing_only)} pairs):")
            print(f"{'-'*80}\n")

            missing_by_type = defaultdict(int)
            missing_count = 0

            for pair_info in missing_only[:5]:  # Show first 5
                base_i = pair_info["base_i"]
                base_j = pair_info["base_j"]
                comp = pair_info["comparison"]

                print(f"Pair ({base_i}, {base_j}):")
                print(
                    f"  Legacy: {comp.legacy_num_hbonds} H-bonds, Modern: {comp.modern_num_hbonds} H-bonds"
                )

                for hb_missing in comp.missing_in_modern:
                    missing_count += 1
                    hb_type = hb_missing["type"]
                    missing_by_type[hb_type] += 1
                    donor = hb_missing["donor_atom"].strip()
                    acceptor = hb_missing["acceptor_atom"].strip()
                    dist = hb_missing["distance"]
                    print(
                        f"    Missing: {donor} -> {acceptor}: {dist:.3f} Å (type='{hb_type}')"
                    )

            if len(missing_only) > 5:
                print(f"\n  ... and {len(missing_only) - 5} more pairs")

            print(f"\n  Missing H-bond type breakdown:")
            for hb_type, count in sorted(missing_by_type.items()):
                type_name = {
                    "-": "standard",
                    "*": "non-standard",
                    " ": "invalid",
                }.get(hb_type, f"unknown ({hb_type})")
                print(f"    {type_name} (type='{hb_type}'): {count}")

        # Show a sample of count mismatches
        if count_mismatch:
            print(f"\n{'-'*80}")
            print(f"Pairs with count differences (showing first 3):")
            print(f"{'-'*80}\n")

            for pair_info in count_mismatch[:3]:
                base_i = pair_info["base_i"]
                base_j = pair_info["base_j"]
                comp = pair_info["comparison"]
                leg_rec = pair_info["legacy_record"]
                mod_rec = pair_info["modern_record"]

                print(f"Pair ({base_i}, {base_j}):")
                print(
                    f"  Legacy: {comp.legacy_num_hbonds} H-bonds, Modern: {comp.modern_num_hbonds} H-bonds"
                )

                leg_hbonds = leg_rec.get("hbonds", [])
                mod_hbonds = mod_rec.get("hbonds", [])

                # Show legacy H-bond types
                leg_types = [h.get("type", " ") for h in leg_hbonds]
                mod_types = [h.get("type", " ") for h in mod_hbonds]

                print(f"  Legacy H-bond types: {leg_types}")
                print(f"  Modern H-bond types: {mod_types}")

                if comp.missing_in_modern:
                    print(f"  Missing in modern:")
                    for hb_missing in comp.missing_in_modern[:5]:
                        donor = hb_missing["donor_atom"].strip()
                        acceptor = hb_missing["acceptor_atom"].strip()
                        dist = hb_missing["distance"]
                        hb_type = hb_missing["type"]
                        print(
                            f"    {donor} -> {acceptor}: {dist:.3f} Å (type='{hb_type}')"
                        )
                print()


def main():
    if len(sys.argv) < 2:
        print("Usage: analyze_hbond_differences.py <pdb_id> [pdb_id2 ...]")
        print("Example: analyze_hbond_differences.py 1VBY")
        sys.exit(1)

    project_root = Path(".")
    for pdb_id in sys.argv[1:]:
        try:
            analyze_hbond_differences(pdb_id, project_root)
        except Exception as e:
            print(f"Error analyzing {pdb_id}: {e}")
            import traceback

            traceback.print_exc()


if __name__ == "__main__":
    main()

