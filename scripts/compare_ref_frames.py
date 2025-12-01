#!/usr/bin/env python3
"""
Compare ref_frames.dat (legacy) with ref_frames_modern.dat (modern) files.

Usage:
    python3 scripts/compare_ref_frames.py [legacy_file] [modern_file]
    python3 scripts/compare_ref_frames.py  # uses default files in current directory
    python3 scripts/compare_ref_frames.py --by-residue  # compare by residue pairs
"""

import sys
import re
from pathlib import Path
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict
import argparse


@dataclass
class RefFrame:
    """Reference frame for a base pair."""
    bp_num: int
    bp_type: str
    description: str
    origin: Tuple[float, float, float]
    x_axis: Tuple[float, float, float]
    y_axis: Tuple[float, float, float]
    z_axis: Tuple[float, float, float]
    res1: Optional[int] = None  # First residue index
    res2: Optional[int] = None  # Second residue index


def extract_residue_pair_from_inp(inp_file: Path, bp_num: int) -> Optional[Tuple[int, int]]:
    """Extract residue pair from .inp file for a given base pair number."""
    if not inp_file.exists():
        return None
    
    with open(inp_file, 'r') as f:
        for line in f:
            # Match lines like "  947  1013   0 #    1" or "    1   947  1013     0 # GC"
            # Legacy format: "  res1  res2   flag #    bp_num"
            legacy_match = re.match(r'\s*(\d+)\s+(\d+)\s+\d+\s+#\s*(\d+)', line)
            if legacy_match:
                res1, res2, num = int(legacy_match.group(1)), int(legacy_match.group(2)), int(legacy_match.group(3))
                if num == bp_num:
                    return (min(res1, res2), max(res1, res2))
            
            # Modern format: "    bp_num   res1  res2     flag # type"
            modern_match = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+#', line)
            if modern_match:
                num, res1, res2 = int(modern_match.group(1)), int(modern_match.group(2)), int(modern_match.group(3))
                if num == bp_num:
                    return (min(res1, res2), max(res1, res2))
    
    return None


def parse_ref_frames_file(filepath: Path, inp_file: Optional[Path] = None) -> List[RefFrame]:
    """Parse a ref_frames.dat file and return list of RefFrame objects."""
    frames = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    if not lines:
        return frames
    
    # First line: number of base pairs
    first_line = lines[0].strip()
    match = re.match(r'\s*(\d+)\s+base-pairs', first_line)
    if not match:
        raise ValueError(f"Invalid first line format: {first_line}")
    
    num_bp = int(match.group(1))
    
    # Parse each base pair (5 lines each)
    i = 1
    while i < len(lines):
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
        
        # Next 4 lines: origin, x-axis, y-axis, z-axis
        def parse_vector(line_str: str) -> Tuple[float, float, float]:
            """Parse a line with 3 floating point values."""
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
            print(f"Warning: Could not parse frame {bp_num}: {e}")
            i += 1
            continue
        
        # Try to get residue pair from .inp file
        res1, res2 = None, None
        if inp_file:
            pair = extract_residue_pair_from_inp(inp_file, bp_num)
            if pair:
                res1, res2 = pair
        
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
    
    return frames


def vector_diff(v1: Tuple[float, float, float], v2: Tuple[float, float, float]) -> float:
    """Calculate the maximum absolute difference between two vectors."""
    return max(abs(a - b) for a, b in zip(v1, v2))


def compare_frames_by_position(legacy: List[RefFrame], modern: List[RefFrame], 
                               tolerance: float = 0.01, verbose: bool = False) -> dict:
    """Compare two lists of reference frames by position number."""
    results = {
        'legacy_count': len(legacy),
        'modern_count': len(modern),
        'matched': 0,
        'origin_matches': 0,
        'x_axis_matches': 0,
        'y_axis_matches': 0,
        'z_axis_matches': 0,
        'differences': []
    }
    
    # Create lookup by bp_num
    legacy_by_num = {f.bp_num: f for f in legacy}
    modern_by_num = {f.bp_num: f for f in modern}
    
    # Compare common base pairs
    common_nums = set(legacy_by_num.keys()) & set(modern_by_num.keys())
    legacy_only = set(legacy_by_num.keys()) - set(modern_by_num.keys())
    modern_only = set(modern_by_num.keys()) - set(legacy_by_num.keys())
    
    results['common_count'] = len(common_nums)
    results['legacy_only'] = sorted(legacy_only)
    results['modern_only'] = sorted(modern_only)
    
    for bp_num in sorted(common_nums):
        leg = legacy_by_num[bp_num]
        mod = modern_by_num[bp_num]
        
        origin_diff = vector_diff(leg.origin, mod.origin)
        x_diff = vector_diff(leg.x_axis, mod.x_axis)
        y_diff = vector_diff(leg.y_axis, mod.y_axis)
        z_diff = vector_diff(leg.z_axis, mod.z_axis)
        
        origin_match = origin_diff <= tolerance
        x_match = x_diff <= tolerance
        y_match = y_diff <= tolerance
        z_match = z_diff <= tolerance
        
        if origin_match:
            results['origin_matches'] += 1
        if x_match:
            results['x_axis_matches'] += 1
        if y_match:
            results['y_axis_matches'] += 1
        if z_match:
            results['z_axis_matches'] += 1
        
        if origin_match and x_match and y_match and z_match:
            results['matched'] += 1
        else:
            diff_info = {
                'bp_num': bp_num,
                'legacy_type': leg.bp_type,
                'modern_type': mod.bp_type,
                'origin_diff': origin_diff,
                'x_axis_diff': x_diff,
                'y_axis_diff': y_diff,
                'z_axis_diff': z_diff
            }
            results['differences'].append(diff_info)
            
            if verbose:
                print(f"\nBase pair {bp_num}:")
                print(f"  Legacy:  {leg.bp_type} - {leg.description}")
                print(f"  Modern:  {mod.bp_type} - {mod.description}")
                print(f"  Origin diff: {origin_diff:.6f}")
                print(f"  X-axis diff: {x_diff:.6f}")
                print(f"  Y-axis diff: {y_diff:.6f}")
                print(f"  Z-axis diff: {z_diff:.6f}")
    
    return results


def compare_frames_by_residue(legacy: List[RefFrame], modern: List[RefFrame], 
                              tolerance: float = 0.01, verbose: bool = False) -> dict:
    """Compare two lists of reference frames by residue pair."""
    results = {
        'legacy_count': len(legacy),
        'modern_count': len(modern),
        'matched': 0,
        'origin_matches': 0,
        'x_axis_matches': 0,
        'y_axis_matches': 0,
        'z_axis_matches': 0,
        'differences': [],
        'comparison_type': 'by_residue_pair'
    }
    
    # Create lookup by residue pair (min, max)
    legacy_by_res = {}
    modern_by_res = {}
    
    for f in legacy:
        if f.res1 is not None and f.res2 is not None:
            key = (min(f.res1, f.res2), max(f.res1, f.res2))
            legacy_by_res[key] = f
    
    for f in modern:
        if f.res1 is not None and f.res2 is not None:
            key = (min(f.res1, f.res2), max(f.res1, f.res2))
            modern_by_res[key] = f
    
    # Compare common residue pairs
    common_pairs = set(legacy_by_res.keys()) & set(modern_by_res.keys())
    legacy_only = set(legacy_by_res.keys()) - set(modern_by_res.keys())
    modern_only = set(modern_by_res.keys()) - set(legacy_by_res.keys())
    
    results['common_count'] = len(common_pairs)
    results['legacy_only'] = sorted(legacy_only)
    results['modern_only'] = sorted(modern_only)
    results['legacy_with_residues'] = len(legacy_by_res)
    results['modern_with_residues'] = len(modern_by_res)
    
    for res_pair in sorted(common_pairs):
        leg = legacy_by_res[res_pair]
        mod = modern_by_res[res_pair]
        
        origin_diff = vector_diff(leg.origin, mod.origin)
        x_diff = vector_diff(leg.x_axis, mod.x_axis)
        y_diff = vector_diff(leg.y_axis, mod.y_axis)
        z_diff = vector_diff(leg.z_axis, mod.z_axis)
        
        origin_match = origin_diff <= tolerance
        x_match = x_diff <= tolerance
        y_match = y_diff <= tolerance
        z_match = z_diff <= tolerance
        
        if origin_match:
            results['origin_matches'] += 1
        if x_match:
            results['x_axis_matches'] += 1
        if y_match:
            results['y_axis_matches'] += 1
        if z_match:
            results['z_axis_matches'] += 1
        
        if origin_match and x_match and y_match and z_match:
            results['matched'] += 1
        else:
            diff_info = {
                'res_pair': res_pair,
                'legacy_pos': leg.bp_num,
                'modern_pos': mod.bp_num,
                'legacy_type': leg.bp_type,
                'modern_type': mod.bp_type,
                'origin_diff': origin_diff,
                'x_axis_diff': x_diff,
                'y_axis_diff': y_diff,
                'z_axis_diff': z_diff
            }
            results['differences'].append(diff_info)
            
            if verbose:
                print(f"\nResidue pair {res_pair}:")
                print(f"  Legacy:  pos={leg.bp_num}, {leg.bp_type}")
                print(f"  Modern:  pos={mod.bp_num}, {mod.bp_type}")
                print(f"  Origin diff: {origin_diff:.6f}")
                print(f"  X-axis diff: {x_diff:.6f}")
                print(f"  Y-axis diff: {y_diff:.6f}")
                print(f"  Z-axis diff: {z_diff:.6f}")
    
    return results


def print_summary(results: dict, tolerance: float):
    """Print comparison summary."""
    print("=" * 60)
    print("Reference Frames Comparison Summary")
    print("=" * 60)
    print(f"Tolerance: {tolerance}")
    
    if results.get('comparison_type') == 'by_residue_pair':
        print("Comparison mode: BY RESIDUE PAIR")
    else:
        print("Comparison mode: BY POSITION")
    
    print()
    print(f"Legacy file:  {results['legacy_count']} base pairs")
    print(f"Modern file:  {results['modern_count']} base pairs")
    
    if 'legacy_with_residues' in results:
        print(f"  Legacy with residue info: {results['legacy_with_residues']}")
        print(f"  Modern with residue info: {results['modern_with_residues']}")
    
    print(f"Common pairs: {results['common_count']}")
    print()
    
    if results['legacy_only']:
        if results.get('comparison_type') == 'by_residue_pair':
            print(f"Legacy only ({len(results['legacy_only'])}): {results['legacy_only'][:5]}...")
        else:
            print(f"Legacy only ({len(results['legacy_only'])}): {results['legacy_only']}")
    if results['modern_only']:
        if results.get('comparison_type') == 'by_residue_pair':
            print(f"Modern only ({len(results['modern_only'])}): {results['modern_only'][:5]}...")
        else:
            print(f"Modern only ({len(results['modern_only'])}): {results['modern_only']}")
    
    if results['common_count'] > 0:
        print()
        print("Match Statistics (for common base pairs):")
        print(f"  Perfect matches:  {results['matched']}/{results['common_count']} "
              f"({100*results['matched']/results['common_count']:.1f}%)")
        print(f"  Origin matches:   {results['origin_matches']}/{results['common_count']} "
              f"({100*results['origin_matches']/results['common_count']:.1f}%)")
        print(f"  X-axis matches:   {results['x_axis_matches']}/{results['common_count']} "
              f"({100*results['x_axis_matches']/results['common_count']:.1f}%)")
        print(f"  Y-axis matches:   {results['y_axis_matches']}/{results['common_count']} "
              f"({100*results['y_axis_matches']/results['common_count']:.1f}%)")
        print(f"  Z-axis matches:   {results['z_axis_matches']}/{results['common_count']} "
              f"({100*results['z_axis_matches']/results['common_count']:.1f}%)")
    
    if results['differences']:
        print()
        print(f"Differences ({len(results['differences'])}):")
        for diff in results['differences'][:10]:  # Show first 10
            if 'res_pair' in diff:
                print(f"  Pair {diff['res_pair']}: origin={diff['origin_diff']:.4f}, "
                      f"x={diff['x_axis_diff']:.4f}, y={diff['y_axis_diff']:.4f}, "
                      f"z={diff['z_axis_diff']:.4f}")
            else:
                print(f"  BP {diff['bp_num']}: origin={diff['origin_diff']:.4f}, "
                      f"x={diff['x_axis_diff']:.4f}, y={diff['y_axis_diff']:.4f}, "
                      f"z={diff['z_axis_diff']:.4f}")
        if len(results['differences']) > 10:
            print(f"  ... and {len(results['differences']) - 10} more")
    
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description='Compare ref_frames.dat (legacy) with ref_frames_modern.dat (modern)')
    parser.add_argument('legacy_file', nargs='?', default='ref_frames.dat',
                        help='Path to legacy ref_frames.dat (default: ref_frames.dat)')
    parser.add_argument('modern_file', nargs='?', default='ref_frames_modern.dat',
                        help='Path to modern ref_frames_modern.dat (default: ref_frames_modern.dat)')
    parser.add_argument('-t', '--tolerance', type=float, default=0.01,
                        help='Tolerance for floating point comparison (default: 0.01)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed differences')
    parser.add_argument('--by-residue', action='store_true',
                        help='Compare by residue pair instead of position')
    parser.add_argument('--legacy-inp', type=str, default=None,
                        help='Path to legacy .inp file for residue mapping')
    parser.add_argument('--modern-inp', type=str, default=None,
                        help='Path to modern .inp file for residue mapping')
    
    args = parser.parse_args()
    
    legacy_path = Path(args.legacy_file)
    modern_path = Path(args.modern_file)
    
    # Check files exist
    if not legacy_path.exists():
        print(f"Error: Legacy file not found: {legacy_path}")
        print("\nTo generate legacy ref_frames.dat:")
        print("  cd org && ./build/bin/find_pair_original ../data/pdb/1H4S.pdb")
        sys.exit(1)
    
    if not modern_path.exists():
        print(f"Error: Modern file not found: {modern_path}")
        print("\nTo generate modern ref_frames_modern.dat:")
        print("  ./build/find_pair_app data/pdb/1H4S.pdb output.inp")
        sys.exit(1)
    
    print(f"Comparing:")
    print(f"  Legacy: {legacy_path}")
    print(f"  Modern: {modern_path}")
    print()
    
    # Determine inp files for residue mapping
    legacy_inp = Path(args.legacy_inp) if args.legacy_inp else None
    modern_inp = Path(args.modern_inp) if args.modern_inp else None
    
    # Auto-detect inp files if --by-residue is used
    if args.by_residue:
        if legacy_inp is None:
            # Try common locations
            candidates = [
                legacy_path.parent / "org" / "1H4S.inp",
                Path("org/1H4S.inp"),
                legacy_path.with_suffix('.inp')
            ]
            for c in candidates:
                if c.exists():
                    legacy_inp = c
                    break
        
        if modern_inp is None:
            candidates = [
                modern_path.parent / "modern.inp",
                Path("modern.inp"),
                modern_path.with_suffix('.inp')
            ]
            for c in candidates:
                if c.exists():
                    modern_inp = c
                    break
        
        if legacy_inp:
            print(f"  Legacy .inp: {legacy_inp}")
        if modern_inp:
            print(f"  Modern .inp: {modern_inp}")
        print()
    
    # Parse files
    try:
        legacy_frames = parse_ref_frames_file(legacy_path, legacy_inp)
    except Exception as e:
        print(f"Error parsing legacy file: {e}")
        sys.exit(1)
    
    try:
        modern_frames = parse_ref_frames_file(modern_path, modern_inp)
    except Exception as e:
        print(f"Error parsing modern file: {e}")
        sys.exit(1)
    
    # Compare
    if args.by_residue:
        results = compare_frames_by_residue(legacy_frames, modern_frames, 
                                            tolerance=args.tolerance, verbose=args.verbose)
    else:
        results = compare_frames_by_position(legacy_frames, modern_frames, 
                                             tolerance=args.tolerance, verbose=args.verbose)
    
    # Print summary
    print_summary(results, args.tolerance)
    
    # Exit with error code if there are differences
    if results['matched'] < results['common_count']:
        sys.exit(1)
    
    return 0


if __name__ == '__main__':
    main()
