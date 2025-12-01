#!/usr/bin/env python3
"""
Compare base pairs found by legacy and modern find_pair.

Usage:
    python3 scripts/compare_base_pairs.py [legacy_inp] [modern_inp]
"""

import sys
import re
from pathlib import Path
from dataclasses import dataclass
from typing import List, Set, Tuple, Optional
import argparse


@dataclass
class BasePair:
    """Represents a base pair from an .inp file."""
    bp_num: int
    res1: int
    res2: int
    flag: int
    bp_type: str = ""
    description: str = ""


def parse_legacy_inp(filepath: Path) -> List[BasePair]:
    """Parse legacy .inp file format."""
    pairs = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    in_pairs = False
    for line in lines:
        line = line.strip()
        
        # Skip header lines
        if line.startswith('#####'):
            break
        
        # Look for pair lines: "  947  1013   0 #    1 | ..."
        match = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+#\s*(\d+)\s*[|x+]?\s*(.*)', line)
        if match:
            res1 = int(match.group(1))
            res2 = int(match.group(2))
            flag = int(match.group(3))
            bp_num = int(match.group(4))
            description = match.group(5).strip()
            
            # Extract bp_type from description (e.g., "G-----C" -> "GC")
            bp_type_match = re.search(r'\[\.\.(\w)\](\w)[-*+]+(\w)\[\.\.(\w)\]', description)
            if bp_type_match:
                bp_type = bp_type_match.group(2) + bp_type_match.group(3)
            else:
                bp_type = ""
            
            pairs.append(BasePair(
                bp_num=bp_num,
                res1=res1,
                res2=res2,
                flag=flag,
                bp_type=bp_type,
                description=description
            ))
    
    return pairs


def parse_modern_inp(filepath: Path) -> List[BasePair]:
    """Parse modern .inp file format."""
    pairs = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        
        # Look for pair lines: "    1   947  1011     0 # GC"
        match = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*(?:#\s*(\w+))?', line)
        if match:
            bp_num = int(match.group(1))
            res1 = int(match.group(2))
            res2 = int(match.group(3))
            flag = int(match.group(4))
            bp_type = match.group(5) if match.group(5) else ""
            
            pairs.append(BasePair(
                bp_num=bp_num,
                res1=res1,
                res2=res2,
                flag=flag,
                bp_type=bp_type
            ))
    
    return pairs


def normalize_pair(res1: int, res2: int) -> Tuple[int, int]:
    """Normalize pair to (min, max) order."""
    return (min(res1, res2), max(res1, res2))


def compare_pairs(legacy: List[BasePair], modern: List[BasePair]) -> dict:
    """Compare legacy and modern base pairs."""
    
    # Create sets of normalized pairs
    legacy_set = {normalize_pair(p.res1, p.res2) for p in legacy}
    modern_set = {normalize_pair(p.res1, p.res2) for p in modern}
    
    # Find common and unique pairs
    common = legacy_set & modern_set
    legacy_only = legacy_set - modern_set
    modern_only = modern_set - legacy_set
    
    # Create lookup by residue
    legacy_by_res1 = {}
    for p in legacy:
        if p.res1 not in legacy_by_res1:
            legacy_by_res1[p.res1] = []
        legacy_by_res1[p.res1].append(p)
    
    modern_by_res1 = {}
    for p in modern:
        if p.res1 not in modern_by_res1:
            modern_by_res1[p.res1] = []
        modern_by_res1[p.res1].append(p)
    
    return {
        'legacy_count': len(legacy),
        'modern_count': len(modern),
        'common': sorted(common),
        'legacy_only': sorted(legacy_only),
        'modern_only': sorted(modern_only),
        'legacy_pairs': legacy,
        'modern_pairs': modern,
        'legacy_by_res1': legacy_by_res1,
        'modern_by_res1': modern_by_res1
    }


def print_comparison(results: dict):
    """Print detailed comparison."""
    print("=" * 70)
    print("Base Pair Comparison: Legacy vs Modern")
    print("=" * 70)
    print()
    print(f"Legacy pairs:  {results['legacy_count']}")
    print(f"Modern pairs:  {results['modern_count']}")
    print(f"Common pairs:  {len(results['common'])}")
    print(f"Legacy only:   {len(results['legacy_only'])}")
    print(f"Modern only:   {len(results['modern_only'])}")
    print()
    
    # Show common pairs
    if results['common']:
        print("-" * 70)
        print(f"COMMON PAIRS ({len(results['common'])}):")
        print("-" * 70)
        for pair in results['common'][:10]:
            print(f"  {pair[0]:5d} - {pair[1]:5d}")
        if len(results['common']) > 10:
            print(f"  ... and {len(results['common']) - 10} more")
    
    # Show legacy-only pairs
    if results['legacy_only']:
        print()
        print("-" * 70)
        print(f"LEGACY ONLY ({len(results['legacy_only'])}):")
        print("-" * 70)
        for pair in results['legacy_only']:
            # Find the original pair to get description
            for p in results['legacy_pairs']:
                if normalize_pair(p.res1, p.res2) == pair:
                    print(f"  {pair[0]:5d} - {pair[1]:5d}  ({p.bp_type})  {p.description[:50]}")
                    break
    
    # Show modern-only pairs
    if results['modern_only']:
        print()
        print("-" * 70)
        print(f"MODERN ONLY ({len(results['modern_only'])}):")
        print("-" * 70)
        for pair in results['modern_only']:
            # Find the original pair to get type
            for p in results['modern_pairs']:
                if normalize_pair(p.res1, p.res2) == pair:
                    print(f"  {pair[0]:5d} - {pair[1]:5d}  ({p.bp_type})")
                    break
    
    # Analyze pairing differences for same first residue
    print()
    print("-" * 70)
    print("PAIRING DIFFERENCES (same first residue, different partner):")
    print("-" * 70)
    
    # Find residues that appear in both but with different partners
    legacy_res1_set = set(results['legacy_by_res1'].keys())
    modern_res1_set = set(results['modern_by_res1'].keys())
    common_res1 = legacy_res1_set & modern_res1_set
    
    diff_count = 0
    for res1 in sorted(common_res1):
        leg_partners = [p.res2 for p in results['legacy_by_res1'][res1]]
        mod_partners = [p.res2 for p in results['modern_by_res1'][res1]]
        
        if set(leg_partners) != set(mod_partners):
            diff_count += 1
            leg_p = results['legacy_by_res1'][res1][0]
            mod_p = results['modern_by_res1'][res1][0]
            print(f"  Res {res1:5d}: Legacy pairs with {leg_partners} ({leg_p.bp_type})")
            print(f"            Modern pairs with {mod_partners} ({mod_p.bp_type})")
    
    if diff_count == 0:
        print("  (none)")
    
    print()
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description='Compare base pairs from legacy and modern .inp files')
    parser.add_argument('legacy_file', nargs='?', default='org/1H4S.inp',
                        help='Path to legacy .inp file')
    parser.add_argument('modern_file', nargs='?', default='modern_output.inp',
                        help='Path to modern .inp file')
    
    args = parser.parse_args()
    
    legacy_path = Path(args.legacy_file)
    modern_path = Path(args.modern_file)
    
    if not legacy_path.exists():
        print(f"Error: Legacy file not found: {legacy_path}")
        sys.exit(1)
    
    if not modern_path.exists():
        print(f"Error: Modern file not found: {modern_path}")
        sys.exit(1)
    
    print(f"Comparing:")
    print(f"  Legacy: {legacy_path}")
    print(f"  Modern: {modern_path}")
    print()
    
    # Parse files
    legacy_pairs = parse_legacy_inp(legacy_path)
    modern_pairs = parse_modern_inp(modern_path)
    
    # Compare
    results = compare_pairs(legacy_pairs, modern_pairs)
    
    # Print comparison
    print_comparison(results)


if __name__ == '__main__':
    main()

