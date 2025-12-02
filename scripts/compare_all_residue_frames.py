#!/usr/bin/env python3
"""
Compare residue reference frames (origin + orientation) between legacy and modern.

For each residue:
1. Extract frame from legacy (org, orien)
2. Extract frame from modern (org, orien)  
3. Compare the actual numerical values
4. Report differences

Usage:
    python3 scripts/compare_all_residue_frames.py <PDB_ID>
    python3 scripts/compare_all_residue_frames.py --all [--stop-on-fail]
"""

import sys
import json
import csv
import math
from pathlib import Path
from typing import Dict, Tuple, Optional
import argparse
from dataclasses import dataclass


@dataclass
class FrameDifference:
    """Difference between two frames."""
    residue_idx: int
    chain_id: str
    residue_seq: int
    residue_name: str
    origin_diff: float  # Euclidean distance in Angstroms
    orientation_diff: float  # Max absolute difference in rotation matrix


def extract_frames_from_base_pairs(pairs_json: list) -> Dict[int, Tuple]:
    """Extract frames by residue index from base_pair JSON.
    
    Returns: {legacy_residue_idx: (org, orien_flat)}
    """
    frames = {}
    
    for pair in pairs_json:
        # Residue i
        idx_i = pair.get('base_i')
        if idx_i and 'org_i' in pair and 'orien_i' in pair:
            org = pair['org_i']
            orien = pair['orien_i']
            
            # Both modern and legacy use nested 3x3 format: [[row1], [row2], [row3]]
            # Flatten to [a,b,c,d,e,f,g,h,i] for comparison
            if isinstance(orien, list) and len(orien) >= 3:
                if isinstance(orien[0], list) and len(orien[0]) >= 3:
                    # Nested 3x3: flatten row by row
                    orien_flat = []
                    for i in range(min(3, len(orien))):
                        for j in range(min(3, len(orien[i]))):
                            orien_flat.append(orien[i][j])
                else:
                    # Already flat (shouldn't happen)
                    orien_flat = orien[:9]
            else:
                orien_flat = orien
            
            frames[idx_i] = (org, orien_flat)
        
        # Residue j
        idx_j = pair.get('base_j')
        if idx_j and 'org_j' in pair and 'orien_j' in pair:
            org = pair['org_j']
            orien = pair['orien_j']
            
            # Flatten orien - both modern and legacy use nested 3x3 format
            if isinstance(orien, list) and len(orien) == 3:
                if isinstance(orien[0], list):
                    # Nested format: [[a,b,c], [d,e,f], [g,h,i]]
                    orien_flat = [orien[i][j] for i in range(3) for j in range(3)]
                else:
                    # Should not happen
                    orien_flat = orien
            else:
                orien_flat = orien
            
            frames[idx_j] = (org, orien_flat)
    
    return frames


def compare_frames(pdb_id: str, project_root: Path, 
                   tolerance_org: float = 0.01,
                   tolerance_orien: float = 0.001) -> dict:
    """Compare frames for one PDB."""
    result = {
        'pdb_id': pdb_id,
        'status': 'UNKNOWN',
        'residues_compared': 0,
        'perfect_matches': 0,
        'acceptable_matches': 0,
        'large_differences': 0,
        'max_origin_diff': 0.0,
        'max_orien_diff': 0.0,
        'differences': []
    }
    
    # Load modern base_pair
    modern_file = project_root / "data" / "json" / "base_pair" / f"{pdb_id}.json"
    if not modern_file.exists():
        result['status'] = 'SKIP'
        return result
    
    # Load legacy base_pair  
    legacy_file = project_root / "data" / "json_legacy" / "base_pair" / f"{pdb_id}.json"
    if not legacy_file.exists():
        result['status'] = 'SKIP'
        return result
    
    try:
        with open(modern_file) as f:
            modern_pairs = json.load(f)
        with open(legacy_file) as f:
            legacy_pairs = json.load(f)
    except Exception as e:
        result['status'] = 'ERROR'
        result['error'] = str(e)
        return result
    
    # Extract frames by residue index
    modern_frames = extract_frames_from_base_pairs(modern_pairs)
    legacy_frames = extract_frames_from_base_pairs(legacy_pairs)
    
    # Also load base_frame_calc for residue names
    residue_info = {}
    try:
        with open(project_root / "data" / "json" / "base_frame_calc" / f"{pdb_id}.json") as f:
            base_frames = json.load(f)
        for rec in base_frames:
            idx = rec.get('legacy_residue_idx', rec.get('residue_idx'))
            if idx:
                residue_info[idx] = {
                    'chain': rec.get('chain_id', '?'),
                    'seq': rec.get('residue_seq', 0),
                    'name': rec.get('residue_name', '???').strip()
                }
    except:
        pass
    
    # Compare common residues
    common_indices = sorted(set(modern_frames.keys()) & set(legacy_frames.keys()))
    result['residues_compared'] = len(common_indices)
    
    for idx in common_indices:
        mod_org, mod_orien = modern_frames[idx]
        leg_org, leg_orien = legacy_frames[idx]
        
        # Compare origin (Euclidean distance)
        if len(mod_org) == 3 and len(leg_org) == 3:
            org_diff = math.sqrt(sum((mod_org[i] - leg_org[i])**2 for i in range(3)))
        else:
            org_diff = float('inf')
        
        # Compare orientation (max absolute difference)
        if len(mod_orien) == 9 and len(leg_orien) == 9:
            orien_diff = max(abs(mod_orien[i] - leg_orien[i]) for i in range(9))
        else:
            orien_diff = float('inf')
        
        result['max_origin_diff'] = max(result['max_origin_diff'], org_diff)
        result['max_orien_diff'] = max(result['max_orien_diff'], orien_diff)
        
        # Categorize
        if org_diff < 0.001 and orien_diff < 0.0001:
            result['perfect_matches'] += 1
        elif org_diff < tolerance_org and orien_diff < tolerance_orien:
            result['acceptable_matches'] += 1
        else:
            result['large_differences'] += 1
            
            # Store details
            info = residue_info.get(idx, {'chain': '?', 'seq': 0, 'name': '???'})
            result['differences'].append(FrameDifference(
                residue_idx=idx,
                chain_id=info['chain'],
                residue_seq=info['seq'],
                residue_name=info['name'],
                origin_diff=org_diff,
                orientation_diff=orien_diff
            ))
    
    # Determine status
    if result['large_differences'] > 0:
        result['status'] = 'FAIL'
    elif result['acceptable_matches'] > 0:
        result['status'] = 'ACCEPTABLE'
    else:
        result['status'] = 'PERFECT'
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Compare residue reference frames between legacy and modern'
    )
    parser.add_argument('pdb_id', nargs='?', help='PDB ID to check')
    parser.add_argument('--all', action='store_true',
                       help='Check all validated PDBs')
    parser.add_argument('--stop-on-fail', action='store_true',
                       help='Stop on first PDB with large differences')
    parser.add_argument('--tolerance-org', type=float, default=0.01,
                       help='Tolerance for origin differences (Angstroms)')
    parser.add_argument('--tolerance-orien', type=float, default=0.001,
                       help='Tolerance for orientation differences')
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    
    if args.pdb_id:
        # Single PDB
        result = compare_frames(args.pdb_id, project_root, 
                               args.tolerance_org, args.tolerance_orien)
        
        print(f"\n{'='*70}")
        print(f"Frame Comparison: {result['pdb_id']}")
        print(f"{'='*70}")
        print(f"Status: {result['status']}")
        print(f"Residues compared: {result['residues_compared']}")
        print(f"Perfect matches: {result['perfect_matches']}")
        print(f"Acceptable matches: {result['acceptable_matches']}")
        print(f"Large differences: {result['large_differences']}")
        print(f"\nMax origin diff: {result['max_origin_diff']:.6f} Å")
        print(f"Max orientation diff: {result['max_orien_diff']:.8f}")
        
        if result['differences']:
            print(f"\nResidues with large differences:")
            for diff in result['differences'][:10]:
                print(f"  {diff.chain_id}:{diff.residue_seq} (idx={diff.residue_idx}): "
                      f"org={diff.origin_diff:.4f}Å, orien={diff.orientation_diff:.6f}")
        
        return 0 if result['status'] in ['PERFECT', 'ACCEPTABLE'] else 1
    
    elif args.all:
        # All validated PDBs
        status_csv = project_root / "data" / "index_validation_status.csv"
        
        pass_pdbs = []
        with open(status_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['match_status'] == 'PASS':
                    pass_pdbs.append(row['pdb_id'])
        
        print(f"Comparing frames for {len(pass_pdbs)} validated PDBs...")
        print(f"Tolerance: origin < {args.tolerance_org} Å, orientation < {args.tolerance_orien}")
        print()
        
        total_checked = 0
        total_perfect = 0
        total_acceptable = 0
        total_fail = 0
        total_skip = 0
        
        for i, pdb_id in enumerate(pass_pdbs):
            result = compare_frames(pdb_id, project_root,
                                   args.tolerance_org, args.tolerance_orien)
            total_checked += 1
            
            if result['status'] == 'PERFECT':
                total_perfect += 1
                if (i + 1) % 100 == 0:
                    print(f"  Progress: {i+1}/{len(pass_pdbs)} - {total_perfect} perfect, "
                          f"{total_acceptable} acceptable, {total_fail} issues")
            elif result['status'] == 'ACCEPTABLE':
                total_acceptable += 1
                print(f"  ⚠️  {pdb_id}: {result['large_differences']} residues with differences "
                      f"(max org={result['max_origin_diff']:.4f}Å)")
            elif result['status'] == 'FAIL':
                total_fail += 1
                print(f"  ❌ {pdb_id}: {result['large_differences']} residues with LARGE differences "
                      f"(max org={result['max_origin_diff']:.4f}Å, max orien={result['max_orien_diff']:.6f})")
                
                if args.stop_on_fail:
                    print(f"\n⚠️  Stopping on first failure: {pdb_id}")
                    print(f"   Large differences in {result['large_differences']} residues")
                    print(f"   Max origin diff: {result['max_origin_diff']:.6f} Å")
                    print(f"   Max orientation diff: {result['max_orien_diff']:.8f}")
                    return 1
            elif result['status'] == 'SKIP':
                total_skip += 1
        
        print(f"\n{'='*70}")
        print(f"FRAME COMPARISON COMPLETE")
        print(f"{'='*70}")
        print(f"Total checked: {total_checked}")
        print(f"✅ Perfect: {total_perfect}")
        print(f"⚠️  Acceptable: {total_acceptable}")
        print(f"❌ Large diffs: {total_fail}")
        print(f"⏭️  Skipped: {total_skip}")
        
        if total_fail == 0:
            print(f"\n✅ All frames match within tolerance!")
            return 0
        else:
            print(f"\n❌ Found {total_fail} PDBs with large frame differences")
            return 1
    
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())

