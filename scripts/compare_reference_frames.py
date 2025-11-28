#!/usr/bin/env python3
"""Compare reference frames between legacy and modern code for matching base pairs."""

import sys
import json
from pathlib import Path

def extract_base_pairs(inp_file, is_legacy=False):
    """Extract base pairs from .inp file."""
    if not Path(inp_file).exists():
        return []
    pairs = []
    with open(inp_file) as f:
        lines = f.readlines()
        for line in lines[5:]:
            if line.strip().startswith('#') or line.strip().startswith('#####'):
                continue
            parts = line.split()
            if is_legacy:
                if len(parts) >= 2:
                    try:
                        res1, res2 = int(parts[0]), int(parts[1])
                        if res1 > 10 and res2 > 10:
                            pairs.append((res1, res2))
                    except ValueError:
                        continue
            else:
                if len(parts) >= 3 and parts[0].isdigit():
                    try:
                        res1, res2 = int(parts[1]), int(parts[2])
                        pairs.append((res1, res2))
                    except ValueError:
                        continue
    return pairs

def normalize_pair(p1, p2):
    """Normalize a pair to always have smaller index first."""
    return (min(p1, p2), max(p1, p2))

def load_legacy_frames(pdb_name, project_root):
    """Load legacy reference frames from base_pair JSON (orien_i, orien_j, org_i, org_j)."""
    bp_file = project_root / f"data/json_legacy/base_pair/{pdb_name}.json"
    if not bp_file.exists():
        return {}
    
    frames = {}
    with open(bp_file) as f:
        data = json.load(f)
        if isinstance(data, list):
            for entry in data:
                # Legacy base_pair JSON has orien_i, orien_j, org_i, org_j for each pair
                # Extract frames for both residues in the pair
                if 'base_i' in entry and 'orien_i' in entry and 'org_i' in entry:
                    base_i = entry['base_i']
                    frames[base_i] = {
                        'rotation': entry['orien_i'],
                        'origin': entry['org_i']
                    }
                if 'base_j' in entry and 'orien_j' in entry and 'org_j' in entry:
                    base_j = entry['base_j']
                    frames[base_j] = {
                        'rotation': entry['orien_j'],
                        'origin': entry['org_j']
                    }
    return frames

def load_modern_frames(pdb_name, project_root):
    """Load modern reference frames from base_pair JSON or all_ref_frames."""
    # Try base_pair JSON first (has orien_i, orien_j, org_i, org_j)
    bp_file = project_root / f"data/json/base_pair/{pdb_name}.json"
    frames = {}
    
    if bp_file.exists():
        with open(bp_file) as f:
            data = json.load(f)
            if isinstance(data, list):
                for entry in data:
                    if 'base_i' in entry and 'orien_i' in entry and 'org_i' in entry:
                        base_i = entry['base_i']
                        frames[base_i] = {
                            'rotation': entry['orien_i'],
                            'origin': entry['org_i']
                        }
                    if 'base_j' in entry and 'orien_j' in entry and 'org_j' in entry:
                        base_j = entry['base_j']
                        frames[base_j] = {
                            'rotation': entry['orien_j'],
                            'origin': entry['org_j']
                        }
    
    # Also try all_ref_frames if available
    # (Would need to parse main JSON file for this)
    
    return frames

def compare_frames(legacy_frames, modern_frames, legacy_pairs, modern_pairs, project_root):
    """Compare reference frames for matching base pairs."""
    # Normalize pairs
    legacy_norm = {normalize_pair(p1, p2): (p1, p2) for p1, p2 in legacy_pairs}
    modern_norm = {normalize_pair(p1, p2): (p1, p2) for p1, p2 in modern_pairs}
    
    # Find matching pairs
    common_pairs = set(legacy_norm.keys()) & set(modern_norm.keys())
    
    differences = []
    tolerance = 0.01
    
    for pair_norm in common_pairs:
        leg_res1, leg_res2 = legacy_norm[pair_norm]
        mod_res1, mod_res2 = modern_norm[pair_norm]
        
        # Compare frames for residue 1
        if leg_res1 in legacy_frames and mod_res1 in modern_frames:
            leg_frame1 = legacy_frames[leg_res1]
            mod_frame1 = modern_frames[mod_res1]
            
            diff1 = compare_frame_data(leg_frame1, mod_frame1, leg_res1)
            if diff1:
                differences.append({
                    'residue': leg_res1,
                    'pair': pair_norm,
                    'differences': diff1
                })
        
        # Compare frames for residue 2
        if leg_res2 in legacy_frames and mod_res2 in modern_frames:
            leg_frame2 = legacy_frames[leg_res2]
            mod_frame2 = modern_frames[mod_res2]
            
            diff2 = compare_frame_data(leg_frame2, mod_frame2, leg_res2)
            if diff2:
                differences.append({
                    'residue': leg_res2,
                    'pair': pair_norm,
                    'differences': diff2
                })
    
    return differences, len(common_pairs)

def compare_frame_data(leg_frame, mod_frame, residue_idx):
    """Compare two frame data structures."""
    differences = []
    tolerance = 0.01
    
    # Compare origin
    if 'origin' in leg_frame and 'origin' in mod_frame:
        leg_org = leg_frame['origin']
        mod_org = mod_frame['origin']
        if isinstance(leg_org, list) and isinstance(mod_org, list) and len(leg_org) == 3 and len(mod_org) == 3:
            for i, (l, m) in enumerate(zip(leg_org, mod_org)):
                diff = abs(l - m)
                if diff > tolerance:
                    differences.append({
                        'component': f'origin[{i}]',
                        'legacy': l,
                        'modern': m,
                        'diff': diff
                    })
    
    # Compare rotation matrix
    if 'rotation' in leg_frame and 'rotation' in mod_frame:
        leg_rot = leg_frame['rotation']
        mod_rot = mod_frame['rotation']
        if isinstance(leg_rot, list) and isinstance(mod_rot, list):
            # Handle both 3x3 and 9-element formats
            leg_flat = flatten_matrix(leg_rot)
            mod_flat = flatten_matrix(mod_rot)
            
            if len(leg_flat) == 9 and len(mod_flat) == 9:
                for i, (l, m) in enumerate(zip(leg_flat, mod_flat)):
                    diff = abs(l - m)
                    if diff > tolerance:
                        differences.append({
                            'component': f'rotation[{i//3}][{i%3}]',
                            'legacy': l,
                            'modern': m,
                            'diff': diff
                        })
    
    return differences

def flatten_matrix(matrix):
    """Flatten a matrix (handle both 2D and 1D formats)."""
    if isinstance(matrix, list):
        if len(matrix) > 0 and isinstance(matrix[0], list):
            # 2D format: [[a,b,c], [d,e,f], [g,h,i]]
            return [item for row in matrix for item in row]
        else:
            # Already flat: [a,b,c,d,e,f,g,h,i]
            return matrix
    return []

def main():
    project_root = Path(__file__).parent.parent
    
    if len(sys.argv) < 2:
        test_pdbs = ['6V9Q', '7EH2', '1A34']
    else:
        test_pdbs = sys.argv[1:]
    
    print("Comparing reference frames for matching base pairs...\n")
    print("="*80)
    print()
    
    for pdb_name in test_pdbs:
        print(f"Processing {pdb_name}...")
        
        legacy_inp = project_root / f"{pdb_name}_legacy.inp"
        modern_inp = project_root / f"{pdb_name}_modern.inp"
        
        if not legacy_inp.exists() or not modern_inp.exists():
            print(f"  ⚠️  Skipped (missing .inp files)\n")
            continue
        
        legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True)
        modern_pairs = extract_base_pairs(modern_inp, is_legacy=False)
        
        legacy_frames = load_legacy_frames(pdb_name, project_root)
        modern_frames = load_modern_frames(pdb_name, project_root)
        
        if not legacy_frames or not modern_frames:
            print(f"  ⚠️  Skipped (missing frame JSON files)\n")
            continue
        
        differences, num_pairs = compare_frames(
            legacy_frames, modern_frames, legacy_pairs, modern_pairs, project_root
        )
        
        print(f"  Pairs: {num_pairs}")
        print(f"  Legacy frames: {len(legacy_frames)}, Modern frames: {len(modern_frames)}")
        
        if not differences:
            print(f"  ✅ ALL REFERENCE FRAMES MATCH (tolerance: 0.01)")
        else:
            print(f"  ❌ Reference frame differences: {len(differences)} residues")
            for diff in differences[:3]:
                print(f"    Residue {diff['residue']} (pair {diff['pair']}): {len(diff['differences'])} differences")
                for d in diff['differences'][:2]:
                    print(f"      {d['component']}: legacy={d['legacy']:.6f}, modern={d['modern']:.6f}, diff={d['diff']:.6f}")
        
        print()

if __name__ == "__main__":
    main()

