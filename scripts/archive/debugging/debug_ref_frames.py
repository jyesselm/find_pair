#!/usr/bin/env python3
"""
Debug script for ref_frames comparison issues.

This script breaks down the ref_frames comparison into smaller, testable pieces:
1. Parse and compare individual base pairs
2. Check frame calculations step-by-step
3. Verify strand ordering logic
4. Compare with minimal test cases
"""

import sys
import json
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

@dataclass
class RefFrame:
    bp_num: int
    bp_type: str
    description: str
    origin: Tuple[float, float, float]
    x_axis: Tuple[float, float, float]
    y_axis: Tuple[float, float, float]
    z_axis: Tuple[float, float, float]
    res1: Optional[int] = None
    res2: Optional[int] = None

def parse_ref_frames(filepath: Path) -> List[RefFrame]:
    """Parse ref_frames.dat file."""
    frames = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    if not lines:
        return frames
    
    # First line: number of base pairs
    first_line = lines[0].strip()
    match = re.match(r'\s*(\d+)\s+base-pairs?', first_line)
    if not match:
        print(f"ERROR: Invalid first line: {first_line}")
        return frames
    
    num_bp = int(match.group(1))
    print(f"Found {num_bp} base pairs in {filepath.name}")
    
    # Parse each base pair (5 lines each)
    i = 1
    bp_count = 0
    while i < len(lines) and bp_count < num_bp:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        
        # Header line: ...     N bp_type   # description
        header_match = re.match(r'\.\.\.?\s*(\d+)\s+(\S+)\s+#\s*(.*)', line)
        if not header_match:
            i += 1
            continue
        
        bp_num = int(header_match.group(1))
        bp_type = header_match.group(2)
        description = header_match.group(3)
        
        # Extract residue numbers from description if possible
        res1, res2 = None, None
        # Format: "A:...1_:[..C]C - B:..20_:[..G]G"
        desc_match = re.search(r'A:.*?(\d+).*?B:.*?(\d+)', description)
        if desc_match:
            res1 = int(desc_match.group(1))
            res2 = int(desc_match.group(2))
        
        # Next 4 lines: origin, x-axis, y-axis, z-axis
        def parse_vector(line_str: str) -> Tuple[float, float, float]:
            parts = line_str.split('#')[0].strip().split()
            if len(parts) < 3:
                raise ValueError(f"Invalid vector line: {line_str}")
            return (float(parts[0]), float(parts[1]), float(parts[2]))
        
        try:
            origin = parse_vector(lines[i + 1])
            x_axis = parse_vector(lines[i + 2])
            y_axis = parse_vector(lines[i + 3])
            z_axis = parse_vector(lines[i + 4])
        except (IndexError, ValueError) as e:
            print(f"ERROR parsing frame {bp_num}: {e}")
            i += 1
            continue
        
        frames.append(RefFrame(
            bp_num=bp_num,
            bp_type=bp_type,
            description=description,
            origin=origin,
            x_axis=x_axis,
            y_axis=y_axis,
            z_axis=z_axis,
            res1=res1,
            res2=res2
        ))
        
        i += 5
        bp_count += 1
    
    return frames

def compare_vectors(v1: Tuple[float, float, float], 
                   v2: Tuple[float, float, float],
                   name: str,
                   tolerance: float = 0.01) -> Dict:
    """Compare two vectors and return detailed diff info."""
    diff = [abs(a - b) for a, b in zip(v1, v2)]
    max_diff = max(diff)
    is_opposite = all(abs(a + b) < tolerance for a, b in zip(v1, v2))
    is_same = max_diff < tolerance
    
    return {
        'name': name,
        'legacy': v1,
        'modern': v2,
        'diff': diff,
        'max_diff': max_diff,
        'is_same': is_same,
        'is_opposite': is_opposite,
        'matches': is_same or is_opposite
    }

def compare_frame_pair(legacy: RefFrame, modern: RefFrame, 
                      tolerance: float = 0.01) -> Dict:
    """Compare a single frame pair in detail."""
    results = {
        'bp_num': legacy.bp_num,
        'legacy_desc': legacy.description,
        'modern_desc': modern.description,
        'legacy_res': (legacy.res1, legacy.res2),
        'modern_res': (modern.res1, modern.res2),
        'residues_match': (legacy.res1, legacy.res2) == (modern.res1, modern.res2),
        'comparisons': {}
    }
    
    # Compare each component
    results['comparisons']['origin'] = compare_vectors(
        legacy.origin, modern.origin, 'origin', tolerance)
    results['comparisons']['x_axis'] = compare_vectors(
        legacy.x_axis, modern.x_axis, 'x-axis', tolerance)
    results['comparisons']['y_axis'] = compare_vectors(
        legacy.y_axis, modern.y_axis, 'y-axis', tolerance)
    results['comparisons']['z_axis'] = compare_vectors(
        legacy.z_axis, modern.z_axis, 'z-axis', tolerance)
    
    # Summary
    all_match = all(c['matches'] for c in results['comparisons'].values())
    results['all_match'] = all_match
    
    return results

def analyze_differences(legacy_frames: List[RefFrame], 
                       modern_frames: List[RefFrame],
                       tolerance: float = 0.01) -> Dict:
    """Analyze differences between legacy and modern frames."""
    results = {
        'total_legacy': len(legacy_frames),
        'total_modern': len(modern_frames),
        'by_position': [],
        'by_residue': {},
        'summary': {
            'perfect_matches': 0,
            'origin_matches': 0,
            'x_matches': 0,
            'y_matches': 0,
            'z_matches': 0,
            'y_inverted': 0,
            'z_inverted': 0,
            'different_pairs': 0
        }
    }
    
    # Compare by position
    min_len = min(len(legacy_frames), len(modern_frames))
    for i in range(min_len):
        leg = legacy_frames[i]
        mod = modern_frames[i]
        
        comparison = compare_frame_pair(leg, mod, tolerance)
        results['by_position'].append(comparison)
        
        # Update summary
        if comparison['all_match']:
            results['summary']['perfect_matches'] += 1
        if comparison['comparisons']['origin']['matches']:
            results['summary']['origin_matches'] += 1
        if comparison['comparisons']['x_axis']['matches']:
            results['summary']['x_matches'] += 1
        if comparison['comparisons']['y_axis']['matches']:
            results['summary']['y_matches'] += 1
        if comparison['comparisons']['y_axis']['is_opposite']:
            results['summary']['y_inverted'] += 1
        if comparison['comparisons']['z_axis']['matches']:
            results['summary']['z_matches'] += 1
        if comparison['comparisons']['z_axis']['is_opposite']:
            results['summary']['z_inverted'] += 1
        if not comparison['residues_match']:
            results['summary']['different_pairs'] += 1
    
    # Compare by residue pair
    legacy_by_res = {(f.res1, f.res2): f for f in legacy_frames 
                     if f.res1 is not None and f.res2 is not None}
    modern_by_res = {(f.res1, f.res2): f for f in modern_frames 
                     if f.res1 is not None and f.res2 is not None}
    
    common_pairs = set(legacy_by_res.keys()) & set(modern_by_res.keys())
    for res_pair in common_pairs:
        leg = legacy_by_res[res_pair]
        mod = modern_by_res[res_pair]
        comparison = compare_frame_pair(leg, mod, tolerance)
        results['by_residue'][res_pair] = comparison
    
    return results

def print_report(results: Dict, verbose: bool = False):
    """Print detailed analysis report."""
    print("=" * 80)
    print("REF_FRAMES COMPARISON ANALYSIS")
    print("=" * 80)
    print(f"\nTotal frames: Legacy={results['total_legacy']}, Modern={results['total_modern']}")
    
    summary = results['summary']
    print(f"\nSUMMARY:")
    print(f"  Perfect matches: {summary['perfect_matches']}/{results['total_legacy']}")
    print(f"  Origin matches: {summary['origin_matches']}/{results['total_legacy']}")
    print(f"  X-axis matches: {summary['x_matches']}/{results['total_legacy']}")
    print(f"  Y-axis matches: {summary['y_matches']}/{results['total_legacy']}")
    print(f"  Z-axis matches: {summary['z_matches']}/{results['total_legacy']}")
    print(f"  Y-axis inverted: {summary['y_inverted']}/{results['total_legacy']}")
    print(f"  Z-axis inverted: {summary['z_inverted']}/{results['total_legacy']}")
    print(f"  Different residue pairs: {summary['different_pairs']}/{results['total_legacy']}")
    
    if verbose:
        print(f"\nDETAILED COMPARISONS BY POSITION:")
        for comp in results['by_position']:
            print(f"\n  Base pair {comp['bp_num']}:")
            print(f"    Legacy: {comp['legacy_desc']}")
            print(f"    Modern: {comp['modern_desc']}")
            print(f"    Residues match: {comp['residues_match']}")
            for name, c in comp['comparisons'].items():
                status = "✓" if c['matches'] else "✗"
                if c['is_opposite']:
                    status += " (INVERTED)"
                print(f"    {name}: {status} (max_diff={c['max_diff']:.6f})")
    
    if results['by_residue']:
        print(f"\nCOMPARISONS BY RESIDUE PAIR ({len(results['by_residue'])} common pairs):")
        for res_pair, comp in sorted(results['by_residue'].items()):
            status = "✓" if comp['all_match'] else "✗"
            print(f"  {res_pair}: {status}")

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Debug ref_frames comparison issues')
    parser.add_argument('legacy_file', help='Path to legacy ref_frames.dat')
    parser.add_argument('modern_file', help='Path to modern ref_frames_modern.dat')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed comparisons')
    parser.add_argument('-t', '--tolerance', type=float, default=0.01,
                       help='Tolerance for comparisons')
    
    args = parser.parse_args()
    
    legacy_path = Path(args.legacy_file)
    modern_path = Path(args.modern_file)
    
    if not legacy_path.exists():
        print(f"ERROR: Legacy file not found: {legacy_path}")
        sys.exit(1)
    
    if not modern_path.exists():
        print(f"ERROR: Modern file not found: {modern_path}")
        sys.exit(1)
    
    print(f"Parsing legacy file: {legacy_path}")
    legacy_frames = parse_ref_frames(legacy_path)
    
    print(f"Parsing modern file: {modern_path}")
    modern_frames = parse_ref_frames(modern_path)
    
    if not legacy_frames or not modern_frames:
        print("ERROR: Could not parse frames")
        sys.exit(1)
    
    print("\nAnalyzing differences...")
    results = analyze_differences(legacy_frames, modern_frames, args.tolerance)
    
    print_report(results, args.verbose)
    
    # Exit with error if there are significant differences
    if results['summary']['perfect_matches'] < results['total_legacy']:
        sys.exit(1)

if __name__ == '__main__':
    main()

