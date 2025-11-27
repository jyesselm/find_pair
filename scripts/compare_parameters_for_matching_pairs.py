#!/usr/bin/env python3
"""Compare parameters for matching base pairs between legacy and modern code."""

import sys
import subprocess
from pathlib import Path
import re

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

def extract_step_parameters(analyze_output):
    """Extract step parameters from analyze output."""
    params = {}
    lines = analyze_output.split('\n')
    in_step = False
    
    step_num = 0
    for line in lines:
        if "Step Parameters" in line:
            in_step = True
            continue
        if "Helical Parameters" in line:
            break
        if in_step and line.strip() and not line.startswith('#') and not line.startswith('='):
            parts = line.split()
            if len(parts) >= 7:
                try:
                    step_num = int(parts[0])
                    values = [float(p) for p in parts[1:7]]
                    params[step_num] = {
                        'shift': values[0],
                        'slide': values[1],
                        'rise': values[2],
                        'tilt': values[3],
                        'roll': values[4],
                        'twist': values[5]
                    }
                except (ValueError, IndexError):
                    continue
    return params

def extract_helical_parameters(analyze_output):
    """Extract helical parameters from analyze output."""
    params = {}
    lines = analyze_output.split('\n')
    in_helix = False
    
    step_num = 0
    for line in lines:
        if "Helical Parameters" in line:
            in_helix = True
            continue
        if in_helix and line.strip() and not line.startswith('#') and not line.startswith('='):
            parts = line.split()
            if len(parts) >= 7:
                try:
                    step_num = int(parts[0])
                    values = [float(p) for p in parts[1:7]]
                    params[step_num] = {
                        'x_disp': values[0],
                        'y_disp': values[1],
                        'h_rise': values[2],
                        'inclination': values[3],
                        'tip': values[4],
                        'h_twist': values[5]
                    }
                except (ValueError, IndexError):
                    continue
    return params

def compare_parameters_for_pdb(pdb_name, project_root):
    """Compare parameters for a single PDB when pairs match."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_name}.pdb"
    legacy_inp = project_root / f"{pdb_name}_legacy.inp"
    modern_inp = project_root / f"{pdb_name}_modern.inp"
    
    if not pdb_file.exists():
        return None
    
    # Get base pairs
    legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True) if legacy_inp.exists() else []
    modern_pairs = extract_base_pairs(modern_inp, is_legacy=False) if modern_inp.exists() else []
    
    # Normalize pairs
    legacy_norm = set(normalize_pair(p1, p2) for p1, p2 in legacy_pairs)
    modern_norm = set(normalize_pair(p1, p2) for p1, p2 in modern_pairs)
    
    # Only proceed if pairs match
    if legacy_norm != modern_norm:
        return {
            'pdb': pdb_name,
            'match': False,
            'reason': f"Pairs don't match: legacy={len(legacy_norm)}, modern={len(modern_norm)}"
        }
    
    if len(legacy_norm) == 0:
        return {
            'pdb': pdb_name,
            'match': False,
            'reason': "No pairs found"
        }
    
    # Run modern analyze to get parameters
    if not modern_inp.exists():
        return {
            'pdb': pdb_name,
            'match': False,
            'reason': "Modern .inp file not found"
        }
    
    result = subprocess.run(
        f"./build/analyze_app {modern_inp}",
        shell=True,
        cwd=project_root,
        capture_output=True,
        text=True
    )
    
    modern_step_params = extract_step_parameters(result.stdout)
    modern_helix_params = extract_helical_parameters(result.stdout)
    
    # For legacy, we'd need to parse bp_step.par and bp_helical.par files
    # For now, we'll just verify modern has parameters for all steps
    num_steps = len(legacy_pairs) - 1
    step_missing = []
    helix_missing = []
    
    for i in range(1, num_steps + 1):
        if i not in modern_step_params:
            step_missing.append(i)
        if i not in modern_helix_params:
            helix_missing.append(i)
    
    return {
        'pdb': pdb_name,
        'match': True,
        'num_pairs': len(legacy_norm),
        'num_steps': num_steps,
        'step_params_found': len(modern_step_params),
        'helix_params_found': len(modern_helix_params),
        'step_missing': step_missing,
        'helix_missing': helix_missing,
        'step_params': modern_step_params,
        'helix_params': modern_helix_params
    }

def main():
    if len(sys.argv) < 2:
        # Get list of PDBs that have matching pairs from batch comparison
        project_root = Path(__file__).parent.parent
        
        # Test on known matching PDBs first
        test_pdbs = ['6V9Q', '7EH2', '1A34']
        print("Testing known perfect matches first...")
    else:
        test_pdbs = sys.argv[1:]
    
    project_root = Path(__file__).parent.parent
    
    results = []
    matching_results = []
    
    for pdb in test_pdbs:
        print(f"Analyzing {pdb}...", end=" ", flush=True)
        result = compare_parameters_for_pdb(pdb, project_root)
        if result:
            results.append(result)
            if result.get('match'):
                matching_results.append(result)
                print(f"✅ {result['num_pairs']} pairs, {result['num_steps']} steps")
                if result['step_missing'] or result['helix_missing']:
                    print(f"   ⚠️  Missing params: step={result['step_missing']}, helix={result['helix_missing']}")
            else:
                print(f"❌ {result.get('reason', 'Unknown error')}")
    
    # Summary
    print(f"\n{'='*80}")
    print("PARAMETER VERIFICATION SUMMARY")
    print(f"{'='*80}\n")
    
    print(f"Total tested: {len(results)}")
    print(f"Perfect matches (pairs match): {len(matching_results)}")
    
    if matching_results:
        print("\nParameter completeness for matching pairs:")
        for r in matching_results:
            print(f"  {r['pdb']}: {r['step_params_found']}/{r['num_steps']} step params, "
                  f"{r['helix_params_found']}/{r['num_steps']} helix params")
            if r['step_missing']:
                print(f"    Missing step params: {r['step_missing']}")
            if r['helix_missing']:
                print(f"    Missing helix params: {r['helix_missing']}")

if __name__ == "__main__":
    main()

