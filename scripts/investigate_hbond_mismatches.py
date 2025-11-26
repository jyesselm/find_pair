#!/usr/bin/env python3
"""
Investigate H-bond mismatches between legacy and modern code.

This script provides detailed analysis of specific H-bond mismatches,
including the full record information to help identify root causes.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.compare_json import create_comparator, compare_single_pdb


def print_hbond_details(hbond: Dict, prefix: str = "  "):
    """Print detailed H-bond information."""
    donor = hbond.get("donor_atom", "").strip()
    acceptor = hbond.get("acceptor_atom", "").strip()
    distance = hbond.get("distance", 0.0)
    hb_type = hbond.get("type", " ")
    hbond_idx = hbond.get("hbond_idx", "?")
    
    print(f"{prefix}{hbond_idx}: {donor:>6} -> {acceptor:>6}  dist={distance:6.3f}  type='{hb_type}'")


def analyze_mismatch(pdb_id: str, base_i: int, base_j: int):
    """Analyze a specific mismatched H-bond pair in detail."""
    print(f"\n{'='*80}")
    print(f"Analyzing H-bond mismatch: {pdb_id} pair ({base_i}, {base_j})")
    print(f"{'='*80}\n")
    
    comparator = create_comparator(None, project_root)
    comparator.enable_cache = False
    comparator.force_recompute = True
    
    result = compare_single_pdb(pdb_id, project_root, False, comparator, False)
    
    if not result.hbond_list_comparison:
        print(f"ERROR: No hbond_list_comparison found for {pdb_id}")
        return
    
    hlc = result.hbond_list_comparison
    pair_key = (min(base_i, base_j), max(base_i, base_j))
    
    if pair_key not in hlc.pair_comparisons:
        print(f"ERROR: Pair ({base_i}, {base_j}) not found in comparisons")
        return
    
    comp = hlc.pair_comparisons[pair_key]
    
    # Find the full records
    legacy_record = None
    modern_record = None
    
    for mismatch in hlc.mismatched_pairs:
        if mismatch['base_i'] == pair_key[0] and mismatch['base_j'] == pair_key[1]:
            legacy_record = mismatch.get('legacy_record')
            modern_record = mismatch.get('modern_record')
            break
    
    if not legacy_record or not modern_record:
        print(f"ERROR: Could not find full records for pair ({base_i}, {base_j})")
        return
    
    legacy_hbonds = legacy_record.get("hbonds", [])
    modern_hbonds = modern_record.get("hbonds", [])
    
    print(f"Legacy: {len(legacy_hbonds)} H-bonds")
    print(f"Modern: {len(modern_hbonds)} H-bonds\n")
    
    print("LEGACY H-BONDS:")
    print("-" * 80)
    for hb in legacy_hbonds:
        print_hbond_details(hb)
    
    print("\nMODERN H-BONDS:")
    print("-" * 80)
    for hb in modern_hbonds:
        print_hbond_details(hb)
    
    print("\nCOMPARISON:")
    print("-" * 80)
    print(f"  Missing in modern: {len(comp.missing_in_modern)}")
    print(f"  Extra in modern: {len(comp.extra_in_modern)}")
    print(f"  Distance mismatches: {len(comp.mismatched_hbonds)}")
    
    if comp.missing_in_modern:
        print("\n  MISSING IN MODERN:")
        for hb in comp.missing_in_modern:
            print_hbond_details(hb, prefix="    ")
    
    if comp.extra_in_modern:
        print("\n  EXTRA IN MODERN:")
        for hb in comp.extra_in_modern:
            print_hbond_details(hb, prefix="    ")
    
    if comp.mismatched_hbonds:
        print("\n  DISTANCE MISMATCHES:")
        for hb in comp.mismatched_hbonds:
            donor = hb.get("donor_atom", "").strip()
            acceptor = hb.get("acceptor_atom", "").strip()
            leg_dist = hb.get("legacy_distance", 0.0)
            mod_dist = hb.get("modern_distance", 0.0)
            diff = hb.get("distance_diff", 0.0)
            hb_type = hb.get("type", " ")
            print(f"    {donor:>6} -> {acceptor:>6}  Legacy={leg_dist:6.3f}  Modern={mod_dist:6.3f}  Diff={diff:6.3f}  type='{hb_type}'")
    
    # Analyze pattern
    print("\n" + "="*80)
    print("PATTERN ANALYSIS:")
    print("="*80)
    
    # Check if it's a type mismatch
    if comp.missing_in_modern and comp.extra_in_modern:
        missing_hb = comp.missing_in_modern[0]
        extra_hb = comp.extra_in_modern[0]
        
        missing_donor = missing_hb.get("donor_atom", "").strip()
        missing_acceptor = missing_hb.get("acceptor_atom", "").strip()
        missing_dist = missing_hb.get("distance", 0.0)
        missing_type = missing_hb.get("type", " ")
        
        extra_donor = extra_hb.get("donor_atom", "").strip()
        extra_acceptor = extra_hb.get("acceptor_atom", "").strip()
        extra_dist = extra_hb.get("distance", 0.0)
        extra_type = extra_hb.get("type", " ")
        
        # Check if same atoms
        same_atoms = (missing_donor == extra_donor and missing_acceptor == extra_acceptor) or \
                     (missing_donor == extra_acceptor and missing_acceptor == extra_donor)
        
        if same_atoms and abs(missing_dist - extra_dist) < 0.01:
            print(f"\n🔍 TYPE MISMATCH DETECTED:")
            print(f"   Same atoms: {missing_donor} -> {missing_acceptor}")
            print(f"   Same distance: {missing_dist:.3f}")
            print(f"   Legacy type: '{missing_type}' (invalid)")
            print(f"   Modern type: '{extra_type}' (standard)")
            print(f"\n   ROOT CAUSE: Legacy marks this H-bond as type=' ' (invalid),")
            print(f"   but modern marks it as type='{extra_type}' (standard).")
            print(f"   This suggests a difference in:")
            print(f"   1. Conflict resolution (hb_atompair) - is it marked as a conflict?")
            print(f"   2. Validation (validate_hbonds) - is it processed correctly?")


def list_all_mismatches(pdb_id: str):
    """List all mismatched pairs for a PDB."""
    print(f"\n{'='*80}")
    print(f"All H-bond mismatches for {pdb_id}")
    print(f"{'='*80}\n")
    
    comparator = create_comparator(None, project_root)
    comparator.enable_cache = False
    comparator.force_recompute = True
    
    result = compare_single_pdb(pdb_id, project_root, False, comparator, False)
    
    if not result.hbond_list_comparison:
        print(f"No hbond_list_comparison found for {pdb_id}")
        return
    
    hlc = result.hbond_list_comparison
    
    print(f"Total pairs: {hlc.common_count}")
    print(f"Mismatched pairs: {len(hlc.mismatched_pairs)}\n")
    
    type_mismatches = []
    count_mismatches = []
    other_mismatches = []
    
    for mismatch in hlc.mismatched_pairs:
        base_i = mismatch['base_i']
        base_j = mismatch['base_j']
        comp = mismatch['comparison']
        
        # Categorize mismatch
        if comp.num_hbonds_match and len(comp.missing_in_modern) == 1 and len(comp.extra_in_modern) == 1:
            # Type mismatch pattern
            missing = comp.missing_in_modern[0]
            extra = comp.extra_in_modern[0]
            if (missing.get("donor_atom", "").strip() == extra.get("donor_atom", "").strip() and
                missing.get("acceptor_atom", "").strip() == extra.get("acceptor_atom", "").strip()):
                type_mismatches.append((base_i, base_j, missing, extra))
            else:
                other_mismatches.append((base_i, base_j, comp))
        elif not comp.num_hbonds_match:
            count_mismatches.append((base_i, base_j, comp))
        else:
            other_mismatches.append((base_i, base_j, comp))
    
    print(f"TYPE MISMATCHES (same atoms, different type): {len(type_mismatches)}")
    print(f"COUNT MISMATCHES (different number of H-bonds): {len(count_mismatches)}")
    print(f"OTHER MISMATCHES: {len(other_mismatches)}\n")
    
    if type_mismatches:
        print("TYPE MISMATCHES:")
        print("-" * 80)
        for base_i, base_j, missing, extra in type_mismatches[:10]:  # Show first 10
            missing_type = missing.get("type", " ")
            extra_type = extra.get("type", " ")
            donor = missing.get("donor_atom", "").strip()
            acceptor = missing.get("acceptor_atom", "").strip()
            dist = missing.get("distance", 0.0)
            print(f"  ({base_i:3d}, {base_j:3d}): {donor:>6} -> {acceptor:6}  dist={dist:6.3f}  Legacy type='{missing_type}'  Modern type='{extra_type}'")


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python3 scripts/investigate_hbond_mismatches.py list <PDB_ID>")
        print("  python3 scripts/investigate_hbond_mismatches.py analyze <PDB_ID> <base_i> <base_j>")
        print("\nExample:")
        print("  python3 scripts/investigate_hbond_mismatches.py list 1VBY")
        print("  python3 scripts/investigate_hbond_mismatches.py analyze 1VBY 45 62")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "list":
        if len(sys.argv) < 3:
            print("ERROR: PDB ID required for 'list' command")
            sys.exit(1)
        pdb_id = sys.argv[2]
        list_all_mismatches(pdb_id)
    
    elif command == "analyze":
        if len(sys.argv) < 5:
            print("ERROR: PDB ID, base_i, and base_j required for 'analyze' command")
            sys.exit(1)
        pdb_id = sys.argv[2]
        base_i = int(sys.argv[3])
        base_j = int(sys.argv[4])
        analyze_mismatch(pdb_id, base_i, base_j)
    
    else:
        print(f"ERROR: Unknown command '{command}'")
        sys.exit(1)


if __name__ == "__main__":
    main()

