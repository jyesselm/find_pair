#!/usr/bin/env python3
"""Quick comparison summary - compares base pairs and shows modern output."""

import sys
import subprocess
from pathlib import Path

def run_find_pair(pdb_name, legacy=False):
    """Run find_pair and extract base pairs."""
    pdb_file = f"data/pdb/{pdb_name}.pdb"
    inp_file = f"{pdb_name}_legacy.inp" if legacy else f"{pdb_name}_modern.inp"
    
    if legacy:
        cmd = f"org/build/bin/find_pair_original {pdb_file} {inp_file}"
    else:
        cmd = f"./build/find_pair_app {pdb_file} {inp_file}"
    
    subprocess.run(cmd, shell=True, capture_output=True)
    
    if not Path(inp_file).exists():
        return None
    
    pairs = []
    with open(inp_file) as f:
        lines = f.readlines()
        # Skip header lines (PDB path, output file, duplex number, num pairs, flags)
        start_idx = 4  # Start after flags line
        for line in lines[start_idx:]:
            parts = line.split()
            if legacy:
                # Legacy format: bp_num res1 res2 flag # comment
                # Skip lines that are criteria or helix info
                if line.strip().startswith('#'):
                    continue
                if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                    # Skip if it's the criteria line (starts with #####)
                    if line.strip().startswith('#####'):
                        continue
                    try:
                        res1 = int(parts[0])
                        res2 = int(parts[1])
                        # Validate reasonable residue numbers (skip flags line like "1 1")
                        if res1 > 10 and res2 > 10:  # Reasonable threshold
                            pairs.append((res1, res2))
                    except (ValueError, IndexError):
                        continue
            else:
                # Modern format: bp_num res1 res2 flag # comment
                if len(parts) >= 3 and parts[0].isdigit():
                    try:
                        bp_num = int(parts[0])
                        res1 = int(parts[1])
                        res2 = int(parts[2])
                        pairs.append((res1, res2))
                    except (ValueError, IndexError):
                        continue
    return pairs

def main():
    if len(sys.argv) < 2:
        print("Usage: compare_output_summary.py <pdb1> [pdb2] ...")
        sys.exit(1)
    
    pdbs = sys.argv[1:]
    
    for pdb in pdbs:
        print(f"\n{'='*70}")
        print(f"Comparison: {pdb}")
        print(f"{'='*70}\n")
        
        print("Running find_pair...")
        legacy_pairs = run_find_pair(pdb, legacy=True)
        modern_pairs = run_find_pair(pdb, legacy=False)
        
        if legacy_pairs:
            print(f"  Legacy: {len(legacy_pairs)} base pairs")
            print(f"    {legacy_pairs}")
        else:
            print("  Legacy: Failed or not found")
        
        if modern_pairs:
            print(f"  Modern: {len(modern_pairs)} base pairs")
            print(f"    {modern_pairs}")
        else:
            print("  Modern: Failed or not found")
        
        if legacy_pairs and modern_pairs:
            if set(legacy_pairs) == set(modern_pairs):
                print(f"\n✅ Base pairs match! ({len(legacy_pairs)} pairs)")
            else:
                print(f"\n⚠️  Base pairs differ:")
                only_legacy = set(legacy_pairs) - set(modern_pairs)
                only_modern = set(modern_pairs) - set(legacy_pairs)
                if only_legacy:
                    print(f"    Only in legacy: {only_legacy}")
                if only_modern:
                    print(f"    Only in modern: {only_modern}")
        
        # Run modern analyze
        print(f"\nModern analyze output:")
        inp_file = f"{pdb}_modern.inp"
        if Path(inp_file).exists():
            result = subprocess.run(
                f"./build/analyze_app {inp_file}",
                shell=True,
                capture_output=True,
                text=True
            )
            output = result.stdout
            # Extract key info
            for line in output.split('\n'):
                if 'Calculated' in line or 'Step Parameters' in line or 'Helical Parameters' in line:
                    print(f"  {line}")
                elif line.strip() and line[0] in '123456789':
                    # Parameter line
                    print(f"  {line}")
                    if output.split('\n').index(line) > 50:
                        break

if __name__ == "__main__":
    main()

