#!/usr/bin/env python3
"""Compare step and helical parameters exactly for matching pairs."""

import sys
import subprocess
import re
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

def parse_legacy_step_par(par_file):
    """Parse legacy bp_step.par file."""
    if not Path(par_file).exists():
        return {}
    
    params = {}
    with open(par_file) as f:
        # Legacy format has:
        # Line 1: "   7 # base-pairs"
        # Line 2: "   0 # ***local base-pair & step parameters***"
        # Line 3: "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist"
        # Then data lines with base pair info + 13 values: base_pair_type + 6 bp params + 6 step params
        # BUT: First line has 0.000 for all step params (it's the first base pair, no step)
        # Subsequent lines represent steps between consecutive pairs
        
        line_num = 0
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse data line
            # Format: base_pair_type + 6 bp params (Shear Stretch Stagger Buckle Prop-Tw Opening) + 6 step params (Shift Slide Rise Tilt Roll Twist)
            parts = line.split()
            if len(parts) >= 13:
                try:
                    # Check if first token is a base pair type (e.g., "U-A", "C-G")
                    # If it contains '-', it's a base pair type, otherwise it's a number
                    if '-' in parts[0]:
                        # Format: "U-A" + 6 bp params + 6 step params
                        # Step number is line number (1-based), but first line is base pair 1, step 1 is between bp 1 and bp 2
                        step_num = line_num  # line_num starts at 0, becomes 1, 2, 3...
                        # Step parameters are at positions 7-12 (0-indexed: 6-11)
                        values = [float(p) for p in parts[7:13]]
                        if len(values) == 6:
                            # Skip first line if all step params are zero (it's the first base pair)
                            if line_num == 0 and all(abs(v) < 0.001 for v in values):
                                line_num += 1
                                continue
                            
                            params[step_num] = {
                                'shift': values[0],
                                'slide': values[1],
                                'rise': values[2],
                                'tilt': values[3],
                                'roll': values[4],
                                'twist': values[5]
                            }
                            line_num += 1
                except (ValueError, IndexError) as e:
                    continue
            elif len(parts) >= 7:  # Fallback: maybe just step params
                try:
                    step_num = int(parts[0])
                    values = [float(p) for p in parts[1:7]]
                    if len(values) == 6:
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

def parse_legacy_helix_par(par_file):
    """Parse legacy bp_helical.par file."""
    if not Path(par_file).exists():
        return {}
    
    params = {}
    with open(par_file) as f:
        # Legacy format has:
        # Line 1: "   7 # base-pairs"
        # Line 2: "   1 # ***local base-pair & helical parameters***"
        # Line 3: "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     X-disp    Y-disp    h-Rise    Incl.     Tip     h-Twist"
        # Then data lines with base pair info + 13 values: base_pair_type + 6 bp params + 6 helical params
        line_num = 0
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse data line
            # Format: base_pair_type + 6 bp params + 6 helical params (X-disp Y-disp h-Rise Incl. Tip h-Twist)
            parts = line.split()
            if len(parts) >= 13:
                try:
                    # Check if first token is a base pair type (e.g., "U-A", "C-G")
                    if '-' in parts[0]:
                        # Format: "U-A" + 6 bp params + 6 helical params
                        step_num = line_num
                        # Helical parameters are at positions 7-12 (0-indexed: 6-11)
                        values = [float(p) for p in parts[7:13]]
                        if len(values) == 6:
                            # Skip first line if all helical params are zero
                            if line_num == 0 and all(abs(v) < 0.001 for v in values):
                                line_num += 1
                                continue
                            
                            params[step_num] = {
                                'x_disp': values[0],
                                'y_disp': values[1],
                                'h_rise': values[2],
                                'inclination': values[3],
                                'tip': values[4],
                                'h_twist': values[5]
                            }
                            line_num += 1
                except (ValueError, IndexError):
                    continue
            elif len(parts) >= 7:  # Fallback
                try:
                    step_num = int(parts[0])
                    values = [float(p) for p in parts[1:7]]
                    if len(values) == 6:
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

def extract_modern_parameters(analyze_output):
    """Extract parameters from modern analyze output."""
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
                        step_params[step_num] = {
                            'shift': values[0],
                            'slide': values[1],
                            'rise': values[2],
                            'tilt': values[3],
                            'roll': values[4],
                            'twist': values[5]
                        }
                    elif in_helix:
                        helical_params[step_num] = {
                            'x_disp': values[0],
                            'y_disp': values[1],
                            'h_rise': values[2],
                            'inclination': values[3],
                            'tip': values[4],
                            'h_twist': values[5]
                        }
                except (ValueError, IndexError):
                    continue
    
    return step_params, helical_params

def compare_pdb_parameters(pdb_name, project_root):
    """Compare parameters for a PDB where pairs match."""
    legacy_inp = project_root / f"{pdb_name}_legacy.inp"
    modern_inp = project_root / f"{pdb_name}_modern.inp"
    
    if not legacy_inp.exists() or not modern_inp.exists():
        return None
    
    # Verify pairs match
    legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True)
    modern_pairs = extract_base_pairs(modern_inp, is_legacy=False)
    
    legacy_norm = set(normalize_pair(p1, p2) for p1, p2 in legacy_pairs)
    modern_norm = set(normalize_pair(p1, p2) for p1, p2 in modern_pairs)
    
    if legacy_norm != modern_norm or len(legacy_norm) == 0:
        return None
    
    # Run legacy analyze to get .par files
    # Legacy analyze needs to run from project root so it can find data/pdb/
    work_dir = project_root
    legacy_step_par = work_dir / "bp_step.par"
    legacy_helix_par = work_dir / "bp_helical.par"
    
    # Remove existing .par files
    if legacy_step_par.exists():
        legacy_step_par.unlink()
    if legacy_helix_par.exists():
        legacy_helix_par.unlink()
    
    # Run legacy analyze from project root
    legacy_inp_rel = legacy_inp.relative_to(project_root)
    result = subprocess.run(
        f"org/build/bin/analyze_original {legacy_inp_rel} --no-json",
        shell=True,
        cwd=project_root,
        capture_output=True,
        text=True,
        timeout=60
    )
    
    if result.returncode != 0:
        print(f"  ⚠️  Legacy analyze failed: {result.stderr[:200]}")
        return None
    
    # Run modern analyze
    result_modern = subprocess.run(
        f"./build/analyze_app {modern_inp}",
        shell=True,
        cwd=project_root,
        capture_output=True,
        text=True,
        timeout=60
    )
    
    modern_step, modern_helix = extract_modern_parameters(result_modern.stdout)
    
    # Parse legacy files
    legacy_step = parse_legacy_step_par(legacy_step_par)
    legacy_helix = parse_legacy_helix_par(legacy_helix_par)
    
    # Create mapping: find matching consecutive pair sequences
    # Order doesn't matter - we match by residue indices
    def find_matching_step_sequence(legacy_pairs, modern_pairs):
        """Find matching consecutive pair sequences and map their step parameters."""
        matched_steps = []  # List of (legacy_step_num, modern_step_num, bp_pair_seq)
        
        # Build maps: pair sequence -> step number
        # Legacy: step i is between bp i and bp i+1
        legacy_pair_to_step = {}
        for i in range(len(legacy_pairs) - 1):
            bp1 = normalize_pair(legacy_pairs[i][0], legacy_pairs[i][1])
            bp2 = normalize_pair(legacy_pairs[i+1][0], legacy_pairs[i+1][1])
            legacy_pair_to_step[(bp1, bp2)] = i + 1  # 1-based step numbering
        
        # Modern: step i is between bp i and bp i+1
        modern_pair_to_step = {}
        for i in range(len(modern_pairs) - 1):
            bp1 = normalize_pair(modern_pairs[i][0], modern_pairs[i][1])
            bp2 = normalize_pair(modern_pairs[i+1][0], modern_pairs[i+1][1])
            modern_pair_to_step[(bp1, bp2)] = i + 1  # 1-based step numbering
        
        # Find matching consecutive sequences
        for (leg_bp1, leg_bp2), leg_step in legacy_pair_to_step.items():
            if (leg_bp1, leg_bp2) in modern_pair_to_step:
                mod_step = modern_pair_to_step[(leg_bp1, leg_bp2)]
                matched_steps.append((leg_step, mod_step, (leg_bp1, leg_bp2)))
        
        return matched_steps
    
    # Find matching step sequences
    matched_steps = find_matching_step_sequence(legacy_pairs, modern_pairs)
    
    # Compare
    step_diffs = []
    helix_diffs = []
    tolerance = 0.01
    
    for leg_step, mod_step, (bp1, bp2) in matched_steps:
        if leg_step in legacy_step and mod_step in modern_step:
            leg = legacy_step[leg_step]
            mod = modern_step[mod_step]
            
            for param_name in ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']:
                leg_val = leg[param_name]
                mod_val = mod[param_name]
                diff = abs(leg_val - mod_val)
                if diff > tolerance:
                    step_diffs.append({
                        'legacy_step': leg_step,
                        'modern_step': mod_step,
                        'bp_pair': f"{bp1[0]}-{bp1[1]}/{bp2[0]}-{bp2[1]}",
                        'param': param_name,
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': diff
                    })
        elif leg_step not in legacy_step:
            step_diffs.append({
                'legacy_step': leg_step,
                'modern_step': mod_step,
                'bp_pair': f"{bp1[0]}-{bp1[1]}/{bp2[0]}-{bp2[1]}",
                'param': 'ALL',
                'error': 'Missing in legacy'
            })
        elif mod_step not in modern_step:
            step_diffs.append({
                'legacy_step': leg_step,
                'modern_step': mod_step,
                'bp_pair': f"{bp1[0]}-{bp1[1]}/{bp2[0]}-{bp2[1]}",
                'param': 'ALL',
                'error': 'Missing in modern'
            })
    
    # Use same matching for helical parameters
    matched_helix_steps = matched_steps  # Same matching logic
    
    for leg_step, mod_step, (bp1, bp2) in matched_helix_steps:
        if leg_step in legacy_helix and mod_step in modern_helix:
            leg = legacy_helix[leg_step]
            mod = modern_helix[mod_step]
            
            for param_name in ['x_disp', 'y_disp', 'h_rise', 'inclination', 'tip', 'h_twist']:
                leg_val = leg[param_name]
                mod_val = mod[param_name]
                diff = abs(leg_val - mod_val)
                if diff > tolerance:
                    helix_diffs.append({
                        'legacy_step': leg_step,
                        'modern_step': mod_step,
                        'bp_pair': f"{bp1[0]}-{bp1[1]}/{bp2[0]}-{bp2[1]}",
                        'param': param_name,
                        'legacy': leg_val,
                        'modern': mod_val,
                        'diff': diff
                    })
        elif leg_step not in legacy_helix:
            helix_diffs.append({
                'legacy_step': leg_step,
                'modern_step': mod_step,
                'bp_pair': f"{bp1[0]}-{bp1[1]}/{bp2[0]}-{bp2[1]}",
                'param': 'ALL',
                'error': 'Missing in legacy'
            })
        elif mod_step not in modern_helix:
            helix_diffs.append({
                'legacy_step': leg_step,
                'modern_step': mod_step,
                'bp_pair': f"{bp1[0]}-{bp1[1]}/{bp2[0]}-{bp2[1]}",
                'param': 'ALL',
                'error': 'Missing in modern'
            })
    
    return {
        'pdb': pdb_name,
        'num_pairs': len(legacy_pairs),
        'num_matched_steps': len(matched_steps),
        'legacy_step_count': len(legacy_step),
        'modern_step_count': len(modern_step),
        'legacy_helix_count': len(legacy_helix),
        'modern_helix_count': len(modern_helix),
        'step_diffs': step_diffs,
        'helix_diffs': helix_diffs,
        'step_match': len(step_diffs) == 0,
        'helix_match': len(helix_diffs) == 0,
        'all_match': len(step_diffs) == 0 and len(helix_diffs) == 0
    }

def main():
    project_root = Path(__file__).parent.parent
    
    if len(sys.argv) < 2:
        # Use known perfect matches
        test_pdbs = ['6V9Q', '7EH2', '1A34']
    else:
        test_pdbs = sys.argv[1:]
    
    print("Comparing parameters for matching pairs...\n")
    print("="*80)
    print()
    
    results = []
    
    for pdb in test_pdbs:
        print(f"Processing {pdb}...")
        result = compare_pdb_parameters(pdb, project_root)
        
        if result is None:
            print(f"  ⚠️  Skipped (pairs don't match or files missing)\n")
            continue
        
        results.append(result)
        
        print(f"  Pairs: {result['num_pairs']}, Matched step sequences: {result['num_matched_steps']}")
        print(f"  Legacy step params: {result['legacy_step_count']}, modern: {result['modern_step_count']}")
        print(f"  Legacy helix params: {result['legacy_helix_count']}, modern: {result['modern_helix_count']}")
        
        if result['all_match']:
            print(f"  ✅ ALL PARAMETERS MATCH (tolerance: 0.01)")
        else:
            if result['step_diffs']:
                print(f"  ❌ Step parameter differences: {len(result['step_diffs'])}")
                for diff in result['step_diffs'][:3]:
                    if 'error' in diff:
                        bp_info = diff.get('bp_pair', 'N/A')
                        print(f"    {bp_info} (L:{diff['legacy_step']}, M:{diff['modern_step']}): {diff['error']}")
                    else:
                        bp_info = diff.get('bp_pair', 'N/A')
                        print(f"    {bp_info} {diff['param']}: legacy={diff['legacy']:.3f}, "
                              f"modern={diff['modern']:.3f}, diff={diff['diff']:.3f}")
            
            if result['helix_diffs']:
                print(f"  ❌ Helical parameter differences: {len(result['helix_diffs'])}")
                for diff in result['helix_diffs'][:3]:
                    if 'error' in diff:
                        bp_info = diff.get('bp_pair', 'N/A')
                        print(f"    {bp_info} (L:{diff['legacy_step']}, M:{diff['modern_step']}): {diff['error']}")
                    else:
                        bp_info = diff.get('bp_pair', 'N/A')
                        print(f"    {bp_info} {diff['param']}: legacy={diff['legacy']:.3f}, "
                              f"modern={diff['modern']:.3f}, diff={diff['diff']:.3f}")
        
        print()
    
    # Summary
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print()
    
    if not results:
        print("No results to summarize.")
        return
    
    all_match = [r for r in results if r['all_match']]
    some_diffs = [r for r in results if not r['all_match']]
    
    print(f"Total tested: {len(results)}")
    print(f"Perfect parameter matches: {len(all_match)}")
    print(f"Some differences: {len(some_diffs)}")
    
    if all_match:
        print(f"\n✅ Perfect matches:")
        for r in all_match:
            print(f"  {r['pdb']}: {r['num_steps']} steps, all parameters match")
    
    if some_diffs:
        print(f"\n⚠️  Cases with differences:")
        for r in some_diffs:
            print(f"  {r['pdb']}: {len(r['step_diffs'])} step diffs, {len(r['helix_diffs'])} helix diffs")

if __name__ == "__main__":
    main()

