#!/usr/bin/env python3
"""Side-by-side comparison of legacy vs modern output."""

import sys
import subprocess
from pathlib import Path

def extract_base_pairs(inp_file, is_legacy=False):
    """Extract base pairs from .inp file."""
    pairs = []
    with open(inp_file) as f:
        lines = f.readlines()
        # Skip header: line 1 (PDB), line 2 (outfile), line 3 (duplex), line 4 (num bp), line 5 (flags)
        for line in lines[5:]:
            if line.strip().startswith('#') or line.strip().startswith('#####'):
                continue
            parts = line.split()
            if is_legacy:
                if len(parts) >= 2:
                    try:
                        res1, res2 = int(parts[0]), int(parts[1])
                        if res1 > 10 and res2 > 10:  # Skip header artifacts
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

def normalize_pairs(pairs):
    """Normalize pairs (sort and handle reverse)."""
    normalized = []
    for p1, p2 in pairs:
        # Always put smaller index first
        if p1 < p2:
            normalized.append((p1, p2))
        else:
            normalized.append((p2, p1))
    return sorted(set(normalized))

def run_analyze_get_params(inp_file):
    """Run analyze and extract parameters."""
    result = subprocess.run(
        f"./build/analyze_app {inp_file}",
        shell=True,
        capture_output=True,
        text=True,
        cwd=Path(__file__).parent.parent
    )
    
    step_params = []
    helical_params = []
    
    lines = result.stdout.split('\n')
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
                    values = [float(p) for p in parts[1:7]]
                    if in_step:
                        step_params.append(values)
                    elif in_helix:
                        helical_params.append(values)
                except ValueError:
                    continue
    
    return step_params, helical_params

def main():
    if len(sys.argv) < 2:
        print("Usage: compare_side_by_side.py <pdb1> [pdb2] ...")
        sys.exit(1)
    
    project_root = Path(__file__).parent.parent
    
    for pdb in sys.argv[1:]:
        print(f"\n{'='*80}")
        print(f"COMPARISON: {pdb}")
        print(f"{'='*80}\n")
        
        pdb_file = project_root / "data" / "pdb" / f"{pdb}.pdb"
        legacy_inp = project_root / f"{pdb}_legacy.inp"
        modern_inp = project_root / f"{pdb}_modern.inp"
        
        # Run find_pair if needed
        if not legacy_inp.exists():
            print("Running legacy find_pair...")
            subprocess.run(
                f"org/build/bin/find_pair_original data/pdb/{pdb}.pdb {pdb}_legacy.inp",
                shell=True,
                cwd=project_root,
                capture_output=True
            )
        
        if not modern_inp.exists():
            print("Running modern find_pair...")
            subprocess.run(
                f"./build/find_pair_app data/pdb/{pdb}.pdb {pdb}_modern.inp",
                shell=True,
                cwd=project_root,
                capture_output=True
            )
        
        # Extract base pairs
        if legacy_inp.exists():
            legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True)
            legacy_norm = normalize_pairs(legacy_pairs)
        else:
            legacy_pairs = []
            legacy_norm = []
        
        if modern_inp.exists():
            modern_pairs = extract_base_pairs(modern_inp, is_legacy=False)
            modern_norm = normalize_pairs(modern_pairs)
        else:
            modern_pairs = []
            modern_norm = []
        
        # Compare base pairs
        print("BASE PAIRS:")
        print(f"  Legacy: {len(legacy_pairs)} pairs -> {len(legacy_norm)} unique")
        print(f"  Modern: {len(modern_pairs)} pairs -> {len(modern_norm)} unique")
        
        if set(legacy_norm) == set(modern_norm):
            print("  ✅ Base pairs match!")
        else:
            only_legacy = set(legacy_norm) - set(modern_norm)
            only_modern = set(modern_norm) - set(legacy_norm)
            if only_legacy:
                print(f"  ⚠️  Only in legacy: {sorted(only_legacy)[:5]}")
            if only_modern:
                print(f"  ⚠️  Only in modern: {sorted(only_modern)[:5]}")
        
        # Get parameters from modern
        if modern_inp.exists():
            print("\nPARAMETERS (Modern Output):")
            step_params, helical_params = run_analyze_get_params(modern_inp)
            
            if step_params:
                print(f"\n  Step Parameters: {len(step_params)} steps")
                print("    #    Shift    Slide    Rise     Tilt     Roll     Twist")
                for i, params in enumerate(step_params[:5], 1):
                    print(f"    {i:2d}  {params[0]:8.2f} {params[1]:8.2f} {params[2]:8.2f} "
                          f"{params[3]:8.2f} {params[4]:8.2f} {params[5]:8.2f}")
                if len(step_params) > 5:
                    print(f"    ... ({len(step_params) - 5} more)")
            
            if helical_params:
                print(f"\n  Helical Parameters: {len(helical_params)} steps")
                print("    #  X-disp   Y-disp   h-Rise   Incl.    Tip      h-Twist")
                for i, params in enumerate(helical_params[:5], 1):
                    print(f"    {i:2d}  {params[0]:8.2f} {params[1]:8.2f} {params[2]:8.2f} "
                          f"{params[3]:8.2f} {params[4]:8.2f} {params[5]:8.2f}")
                if len(helical_params) > 5:
                    print(f"    ... ({len(helical_params) - 5} more)")
        
        print()

if __name__ == "__main__":
    main()

