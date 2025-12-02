#!/usr/bin/env python3
"""
Compare reference frames between legacy and modern code.

This tool helps debug frame calculation and frame retrieval issues by comparing
frame origins and orientations from frame_calc JSON files.

Usage:
    python3 scripts/compare_frames.py <PDB_ID> [legacy_residue_idx ...]
    python3 scripts/compare_frames.py 9CF3 25 27  # Compare residues 25 and 27 (legacy indices)
    python3 scripts/compare_frames.py 9CF3 --all  # Compare all residues
    
Note: Residue indices must be legacy indices (1-based) from legacy JSON files.
"""

import sys
import json
import math
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional, List

def load_frame_data(json_file: Path) -> Dict[int, Dict]:
    """Load frame data from frame_calc JSON file."""
    frames = {}
    
    if not json_file.exists():
        return frames
    
    with open(json_file) as f:
        data = json.load(f)
    
    if not isinstance(data, list):
        return frames
    
    for record in data:
        if isinstance(record, dict):
            # CRITICAL: Always prefer legacy_residue_idx for comparisons
            # This ensures we're comparing the same residues between legacy and modern
            idx = record.get('legacy_residue_idx')
            if idx is None:
                idx = record.get('residue_idx')
                if idx is not None:
                    print(f"WARNING: Using residue_idx instead of legacy_residue_idx for record")
            
            translation = record.get('translation', [])
            rotation = record.get('rotation_matrix', [])
            
            if idx is not None and translation and len(translation) >= 3:
                frames[idx] = {
                    'origin': tuple(translation[:3]),
                    'rotation': rotation if rotation else None,
                    'residue_name': record.get('residue_name', '?'),
                    'chain_id': record.get('chain_id', '?'),
                    'rms_fit': record.get('rms_fit'),
                }
    
    return frames

def calculate_distance(p1: Tuple[float, float, float], 
                      p2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b)**2 for a, b in zip(p1, p2)))

def compare_frames(pdb_id: str, residue_indices: Optional[List[int]] = None):
    """Compare frames between legacy and modern."""
    legacy_file = Path(f"data/json_legacy/frame_calc/{pdb_id}.json")
    modern_file = Path(f"data/json/frame_calc/{pdb_id}.json")
    
    if not legacy_file.exists():
        print(f"ERROR: Legacy frame file not found: {legacy_file}")
        return
    
    if not modern_file.exists():
        print(f"ERROR: Modern frame file not found: {modern_file}")
        return
    
    legacy_frames = load_frame_data(legacy_file)
    modern_frames = load_frame_data(modern_file)
    
    if not legacy_frames:
        print(f"ERROR: No frame data in legacy file")
        return
    
    if not modern_frames:
        print(f"ERROR: No frame data in modern file")
        return
    
    # Determine which residues to compare
    if residue_indices:
        residues_to_check = residue_indices
    else:
        # Compare all common residues
        residues_to_check = sorted(set(legacy_frames.keys()) & set(modern_frames.keys()))
    
    print(f"Comparing frames for PDB: {pdb_id}")
    print(f"Residues to check: {len(residues_to_check)}")
    print("=" * 80)
    print()
    
    mismatches = []
    matches = []
    
    for idx in residues_to_check:
        if idx not in legacy_frames or idx not in modern_frames:
            continue
        
        leg_frame = legacy_frames[idx]
        mod_frame = modern_frames[idx]
        
        leg_origin = leg_frame['origin']
        mod_origin = mod_frame['origin']
        
        origin_distance = calculate_distance(leg_origin, mod_origin)
        
        if origin_distance > 0.001:  # Threshold for "different"
            mismatches.append((idx, origin_distance, leg_origin, mod_origin))
            status = "❌ MISMATCH"
        else:
            matches.append((idx, origin_distance))
            status = "✅ MATCH"
        
        print(f"Residue {idx:3d} ({leg_frame['residue_name']:>4s}, Chain {leg_frame['chain_id']}): {status}")
        print(f"  Legacy origin: ({leg_origin[0]:10.6f}, {leg_origin[1]:10.6f}, {leg_origin[2]:10.6f})")
        print(f"  Modern origin: ({mod_origin[0]:10.6f}, {mod_origin[1]:10.6f}, {mod_origin[2]:10.6f})")
        print(f"  Distance: {origin_distance:.6f} Å")
        
        if leg_frame.get('rms_fit') is not None:
            print(f"  RMS fit: legacy={leg_frame['rms_fit']:.6f}, modern={mod_frame.get('rms_fit', 'N/A')}")
        print()
    
    # Summary
    print("=" * 80)
    print(f"Summary:")
    print(f"  Total compared: {len(residues_to_check)}")
    print(f"  Matches: {len(matches)}")
    print(f"  Mismatches: {len(mismatches)}")
    
    if mismatches:
        print()
        print("Top mismatches:")
        mismatches.sort(key=lambda x: x[1], reverse=True)
        for idx, dist, leg_o, mod_o in mismatches[:10]:
            print(f"  Residue {idx}: {dist:.6f} Å difference")
            print(f"    Legacy: ({leg_o[0]:.3f}, {leg_o[1]:.3f}, {leg_o[2]:.3f})")
            print(f"    Modern: ({mod_o[0]:.3f}, {mod_o[1]:.3f}, {mod_o[2]:.3f})")

def main():
    parser = argparse.ArgumentParser(
        description='Compare reference frames between legacy and modern code')
    parser.add_argument('pdb_id', help='PDB ID (e.g., 9CF3)')
    parser.add_argument('residues', nargs='*', type=int,
                       help='Legacy residue indices (1-based) to compare (default: all)')
    parser.add_argument('--all', action='store_true',
                       help='Compare all residues (default if no residues specified)')
    
    args = parser.parse_args()
    
    residue_indices = args.residues if args.residues else None
    if not residue_indices and not args.all:
        # Default to all if nothing specified
        residue_indices = None
    
    compare_frames(args.pdb_id, residue_indices)

if __name__ == '__main__':
    main()

