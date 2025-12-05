#!/usr/bin/env python3
"""
Compare legacy test_hbond_initial output with modern initial H-bonds
"""

import json
import subprocess
import sys
import os
import re

def get_legacy_initial(pdb_file, res1, res2):
    """Get legacy initial H-bonds from test_hbond_initial tool"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    tool_path = os.path.join(project_root, 'org', 'build', 'bin', 'test_hbond_initial')
    pdb_path = os.path.join(project_root, pdb_file)
    
    try:
        result = subprocess.run(
            [tool_path, pdb_path, str(res1), str(res2)],
            capture_output=True, text=True, cwd=project_root
        )
        
        if result.returncode != 0:
            print(f"Error running legacy tool: {result.stderr}", file=sys.stderr)
            return None
        
        # Extract JSON from output
        output = result.stdout
        json_start = output.find('{')
        if json_start == -1:
            print("Could not find JSON in legacy output", file=sys.stderr)
            return None
        
        json_str = output[json_start:]
        return json.loads(json_str)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return None

def get_modern_initial(pdb_file, res1, res2):
    """Get modern initial H-bonds from compare_hbond_stages"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    tool_path = os.path.join(project_root, 'build', 'compare_hbond_stages')
    pdb_path = os.path.join(project_root, pdb_file)
    
    try:
        result = subprocess.run(
            [tool_path, pdb_path, str(res1), str(res2)],
            capture_output=True, text=True, cwd=project_root
        )
        
        if result.returncode != 0:
            print(f"Error running modern tool: {result.stderr}", file=sys.stderr)
            return None
        
        # Parse output to extract Stage 1 H-bonds
        output = result.stdout
        stage1_start = output.find("Stage 1:")
        if stage1_start == -1:
            return None
        
        # Extract modern H-bonds from output - only get the numbered list
        modern_hbonds = []
        lines = output[stage1_start:].split('\n')
        seen = set()  # Track seen H-bonds to avoid duplicates
        
        for line in lines:
            # Look for numbered H-bonds like "  1.    N6  ->    N3 , dist=3.440"
            if re.match(r'\s+\d+\.\s+', line) and '->' in line and 'dist=' in line:
                # Parse line like "  1.    N6  ->    N3 , dist=3.440, type=-"
                parts = line.split('->')
                if len(parts) == 2:
                    # Get donor (everything after the number and dot)
                    donor_part = parts[0].strip()
                    donor = donor_part.split('.', 1)[1].strip() if '.' in donor_part else donor_part
                    
                    # Get acceptor and distance
                    acceptor_part = parts[1].split(',')[0].strip()
                    dist_part = [p for p in parts[1].split(',') if 'dist=' in p]
                    if dist_part:
                        dist_str = dist_part[0].split('=')[1].strip()
                        dist_str = dist_str.rstrip(')')
                        dist = float(dist_str)
                        
                        # Create unique key
                        hb_key = (donor, acceptor_part, dist)
                        if hb_key not in seen:
                            seen.add(hb_key)
                            modern_hbonds.append({
                                'donor_atom': donor,
                                'acceptor_atom': acceptor_part,
                                'distance': dist
                            })
        
        return {'hbonds': modern_hbonds}
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return None

def normalize_atom_name(name):
    """Normalize atom name for comparison"""
    return name.strip()

def compare_hbonds(legacy_data, modern_data):
    """Compare legacy and modern H-bonds"""
    if not legacy_data or not modern_data:
        return
    
    legacy_hbonds = legacy_data.get('hbonds', [])
    modern_hbonds = modern_data.get('hbonds', [])
    
    print("=" * 60)
    print("Initial H-bond Comparison")
    print("=" * 60)
    print(f"\nLegacy: {len(legacy_hbonds)} H-bonds")
    print(f"Modern: {len(modern_hbonds)} H-bonds")
    
    print("\nLegacy H-bonds:")
    for i, hb in enumerate(legacy_hbonds, 1):
        donor = normalize_atom_name(hb.get('donor_atom', ''))
        acceptor = normalize_atom_name(hb.get('acceptor_atom', ''))
        dist = hb.get('distance', 0)
        print(f"  {i}. {donor} -> {acceptor}, dist={dist:.3f}")
    
    print("\nModern H-bonds:")
    for i, hb in enumerate(modern_hbonds, 1):
        donor = normalize_atom_name(hb.get('donor_atom', ''))
        acceptor = normalize_atom_name(hb.get('acceptor_atom', ''))
        dist = hb.get('distance', 0)
        print(f"  {i}. {donor} -> {acceptor}, dist={dist:.3f}")
    
    # Find matches
    matches = []
    legacy_matched = [False] * len(legacy_hbonds)
    modern_matched = [False] * len(modern_hbonds)
    
    for i, lhb in enumerate(legacy_hbonds):
        l_donor = normalize_atom_name(lhb.get('donor_atom', ''))
        l_acceptor = normalize_atom_name(lhb.get('acceptor_atom', ''))
        l_dist = lhb.get('distance', 0)
        
        for j, mhb in enumerate(modern_hbonds):
            if modern_matched[j]:
                continue
            
            m_donor = normalize_atom_name(mhb.get('donor_atom', ''))
            m_acceptor = normalize_atom_name(mhb.get('acceptor_atom', ''))
            m_dist = mhb.get('distance', 0)
            
            # Check if same atoms (order matters for H-bonds)
            if l_donor == m_donor and l_acceptor == m_acceptor:
                matches.append((i, j, l_dist, m_dist))
                legacy_matched[i] = True
                modern_matched[j] = True
                break
    
    print(f"\nMatches: {len(matches)}")
    for leg_idx, mod_idx, l_dist, m_dist in matches:
        lhb = legacy_hbonds[leg_idx]
        mhb = modern_hbonds[mod_idx]
        print(f"  {lhb.get('donor_atom')} -> {lhb.get('acceptor_atom')}: "
              f"Legacy {l_dist:.3f}, Modern {m_dist:.3f}")
    
    print(f"\nMissing in modern ({sum(1 for m in legacy_matched if not m)}):")
    for i, matched in enumerate(legacy_matched):
        if not matched:
            hb = legacy_hbonds[i]
            print(f"  - {hb.get('donor_atom')} -> {hb.get('acceptor_atom')}, "
                  f"dist={hb.get('distance'):.3f}")
    
    print(f"\nExtra in modern ({sum(1 for m in modern_matched if not m)}):")
    for i, matched in enumerate(modern_matched):
        if not matched:
            hb = modern_hbonds[i]
            print(f"  + {hb.get('donor_atom')} -> {hb.get('acceptor_atom')}, "
                  f"dist={hb.get('distance'):.3f}")
    
    # Show seidx info
    if 'seidx_i' in legacy_data:
        print(f"\nLegacy seidx ranges:")
        print(f"  Residue i: {legacy_data['seidx_i']}")
        print(f"  Residue j: {legacy_data['seidx_j']}")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <pdb_file> <residue1> <residue2>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    res1 = int(sys.argv[2])
    res2 = int(sys.argv[3])
    
    print("Getting legacy initial H-bonds...")
    legacy_data = get_legacy_initial(pdb_file, res1, res2)
    
    print("Getting modern initial H-bonds...")
    modern_data = get_modern_initial(pdb_file, res1, res2)
    
    compare_hbonds(legacy_data, modern_data)

