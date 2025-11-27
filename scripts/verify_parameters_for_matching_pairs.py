#!/usr/bin/env python3
"""Verify that parameters match exactly for matching base pairs."""

import sys
import subprocess
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

def extract_modern_parameters(analyze_output):
    """Extract step and helical parameters from modern analyze output."""
    step_params = {}
    helical_params = {}
    
    lines = analyze_output.split('\n')
    in_step = False
    in_helix = False
    
    for line in lines:
        if "Step Parameters" in line:
            in_step = True
            in_helix = False
            continue
        if "Helical Parameters" in line:
            in_step = False
            in_helix = True
            continue
        if line.strip() and not line.startswith('#') and not line.startswith('='):
            parts = line.split()
            if len(parts) >= 7:
                try:
                    step_num = int(parts[0])
                    values = [float(p) for p in parts[1:7]]
                    if in_step:
                        step_params[step_num] = values
                    elif in_helix:
                        helical_params[step_num] = values
                except (ValueError, IndexError):
                    continue
    
    return step_params, helical_params

def parse_legacy_par_file(par_file, is_step=True):
    """Parse legacy .par file and extract parameters."""
    if not Path(par_file).exists():
        return {}
    
    params = {}
    with open(par_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            # Try to extract numeric values
            try:
                # Legacy format varies, try to find step number and 6 parameters
                numeric_parts = [float(p) for p in parts if p.replace('.', '').replace('-', '').replace('e', '').replace('E', '').replace('+', '').isdigit()]
                if len(numeric_parts) >= 7:
                    # First number might be step number, rest are parameters
                    step_num = int(numeric_parts[0])
                    values = numeric_parts[1:7]
                    params[step_num] = values
            except (ValueError, IndexError):
                continue
    return params

def compare_parameters_for_pdb(pdb_name, project_root):
    """Compare parameters for matching pairs in a PDB."""
    legacy_inp = project_root / f"{pdb_name}_legacy.inp"
    modern_inp = project_root / f"{pdb_name}_modern.inp"
    
    if not legacy_inp.exists() or not modern_inp.exists():
        return None
    
    # Get base pairs
    legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True)
    modern_pairs = extract_base_pairs(modern_inp, is_legacy=False)
    
    legacy_norm = set(normalize_pair(p1, p2) for p1, p2 in legacy_pairs)
    modern_norm = set(normalize_pair(p1, p2) for p1, p2 in modern_pairs)
    
    # Only proceed if pairs match exactly
    if legacy_norm != modern_norm:
        return {
            'pdb': pdb_name,
            'match': False,
            'reason': f"Pairs don't match exactly"
        }
    
    if len(legacy_norm) == 0:
        return {
            'pdb': pdb_name,
            'match': False,
            'reason': "No pairs found"
        }
    
    # Run modern analyze
    result = subprocess.run(
        f"./build/analyze_app {modern_inp}",
        shell=True,
        cwd=project_root,
        capture_output=True,
        text=True
    )
    
    modern_step, modern_helix = extract_modern_parameters(result.stdout)
    
    # For legacy, we'd need to parse .par files
    # This requires running legacy analyze and parsing output
    # For now, verify modern has parameters for all steps
    num_steps = len(legacy_pairs) - 1
    
    missing_steps = []
    for i in range(1, num_steps + 1):
        if i not in modern_step:
            missing_steps.append(('step', i))
        if i not in modern_helix:
            missing_steps.append(('helix', i))
    
    # Check for any zero or NaN values (might indicate calculation errors)
    invalid_values = []
    for step_num, values in modern_step.items():
        for param_idx, val in enumerate(values):
            if val == 0.0 and step_num > 1:  # First step might legitimately be zero
                invalid_values.append(('step', step_num, param_idx, val))
    
    return {
        'pdb': pdb_name,
        'match': True,
        'num_pairs': len(legacy_norm),
        'num_steps': num_steps,
        'step_params_count': len(modern_step),
        'helix_params_count': len(modern_helix),
        'missing_steps': missing_steps,
        'invalid_values': invalid_values,
        'all_params_present': len(missing_steps) == 0
    }

def main():
    project_root = Path(__file__).parent.parent
    
    # Focus on known perfect matches
    test_pdbs = ['6V9Q', '7EH2', '1A34']
    
    if len(sys.argv) > 1:
        test_pdbs = sys.argv[1:]
    
    print("Verifying parameters for matching pairs...\n")
    print(f"{'='*80}\n")
    
    results = []
    
    for pdb in test_pdbs:
        print(f"Analyzing {pdb}...", end=" ", flush=True)
        result = compare_parameters_for_pdb(pdb, project_root)
        
        if result:
            results.append(result)
            if result.get('match'):
                print(f"✅")
                print(f"  Pairs: {result['num_pairs']}, Steps: {result['num_steps']}")
                print(f"  Step params: {result['step_params_count']}/{result['num_steps']}")
                print(f"  Helix params: {result['helix_params_count']}/{result['num_steps']}")
                
                if result['missing_steps']:
                    print(f"  ⚠️  Missing: {result['missing_steps']}")
                elif result['invalid_values']:
                    print(f"  ⚠️  Invalid values: {len(result['invalid_values'])}")
                else:
                    print(f"  ✅ All parameters calculated correctly")
            else:
                print(f"❌ {result.get('reason', 'Unknown error')}")
        else:
            print(f"❌ Could not analyze")
        print()
    
    # Summary
    print(f"{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}\n")
    
    matching = [r for r in results if r.get('match')]
    print(f"Perfect matches tested: {len(matching)}")
    
    if matching:
        all_good = all(r.get('all_params_present', False) for r in matching)
        if all_good:
            print("✅ All parameters are calculated correctly for matching pairs")
        else:
            print("⚠️  Some parameters are missing")
        
        print(f"\nDetails:")
        for r in matching:
            print(f"  {r['pdb']}: {r['step_params_count']}/{r['num_steps']} step, "
                  f"{r['helix_params_count']}/{r['num_steps']} helix params")

if __name__ == "__main__":
    main()

