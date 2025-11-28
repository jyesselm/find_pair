#!/usr/bin/env python3
"""Compare legacy vs modern output without JSON generation."""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path

# Colors for output
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
BLUE = '\033[94m'
RESET = '\033[0m'

def run_command(cmd, cwd=None, capture_output=True):
    """Run a command and return result."""
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=cwd, capture_output=capture_output,
            text=True, check=False
        )
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def compare_inp_files(legacy_inp, modern_inp):
    """Compare two .inp files."""
    if not os.path.exists(legacy_inp):
        return False, "Legacy .inp file not found"
    if not os.path.exists(modern_inp):
        return False, "Modern .inp file not found"
    
    with open(legacy_inp) as f:
        legacy_lines = f.readlines()
    with open(modern_inp) as f:
        modern_lines = f.readlines()
    
    if len(legacy_lines) != len(modern_lines):
        return False, f"Line count differs: legacy={len(legacy_lines)}, modern={len(modern_lines)}"
    
    differences = []
    for i, (leg, mod) in enumerate(zip(legacy_lines, modern_lines)):
        if leg != mod:
            differences.append(f"Line {i+1}: legacy='{leg.strip()}' vs modern='{mod.strip()}'")
    
    if differences:
        return False, f"{len(differences)} differences found:\n" + "\n".join(differences[:5])
    return True, "Files are identical"

def parse_par_file(par_file):
    """Parse a .par file and return parameter values."""
    if not os.path.exists(par_file):
        return None
    
    params = []
    with open(par_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Try to parse numeric values
            parts = line.split()
            try:
                values = [float(p) for p in parts if p.replace('.', '').replace('-', '').isdigit()]
                if len(values) >= 6:  # At least 6 parameters
                    params.append(values[:6])  # Take first 6
            except:
                continue
    return params

def compare_par_files(legacy_par, modern_par, par_type="Step"):
    """Compare two .par files."""
    legacy_params = parse_par_file(legacy_par)
    modern_params = parse_par_file(modern_par)
    
    if legacy_params is None:
        return False, "Legacy .par file not found or could not be parsed"
    if modern_params is None:
        return False, "Modern .par file not found or could not be parsed"
    
    if len(legacy_params) != len(modern_params):
        return False, f"Parameter count differs: legacy={len(legacy_params)}, modern={len(modern_params)}"
    
    differences = []
    tolerance = 0.01
    for i, (leg, mod) in enumerate(zip(legacy_params, modern_params)):
        for j, (l, m) in enumerate(zip(leg, mod)):
            diff = abs(l - m)
            if diff > tolerance:
                differences.append(f"  Step {i+1}, param {j+1}: legacy={l:.3f}, modern={m:.3f}, diff={diff:.3f}")
    
    if differences:
        return False, f"{len(differences)} parameters differ (tolerance={tolerance}):\n" + "\n".join(differences[:10])
    return True, f"All {len(legacy_params)} steps match within tolerance"

def compare_pdb(pdb_name, project_root):
    """Compare outputs for a single PDB."""
    pdb_file = os.path.join(project_root, "data/pdb", f"{pdb_name}.pdb")
    
    if not os.path.exists(pdb_file):
        print(f"{RED}ERROR: PDB file not found: {pdb_file}{RESET}")
        return False
    
    print(f"\n{'='*70}")
    print(f"{BLUE}Comparing: {pdb_name}{RESET}")
    print(f"{'='*70}\n")
    
    # Create temp directory
    temp_dir = tempfile.mkdtemp(prefix=f"compare_{pdb_name}_")
    
    try:
        # Try multiple possible locations for legacy binaries
        legacy_bin = None
        for name in ["find_pair", "find_pair_original"]:
            for path in [
                os.path.join(project_root, "org", "build", "bin", name),
                os.path.join(project_root, "org", "bin", name),
                os.path.join(project_root, "org", "build", name),
            ]:
                if os.path.exists(path):
                    legacy_bin = path
                    break
            if legacy_bin:
                break
        
        modern_bin = os.path.join(project_root, "build", "find_pair_app")
        
        legacy_analyze = None
        for name in ["analyze", "analyze_original"]:
            for path in [
                os.path.join(project_root, "org", "build", "bin", name),
                os.path.join(project_root, "org", "bin", name),
                os.path.join(project_root, "org", "build", name),
            ]:
                if os.path.exists(path):
                    legacy_analyze = path
                    break
            if legacy_analyze:
                break
        
        modern_analyze = os.path.join(project_root, "build", "analyze_app")
        
        # Check binaries
        if legacy_bin is None:
            print(f"{YELLOW}WARNING: Legacy find_pair binary not found. Skipping legacy comparison.{RESET}")
        elif not os.path.exists(legacy_bin):
            print(f"{YELLOW}WARNING: Legacy binary not found: {legacy_bin}{RESET}")
            legacy_bin = None
        
        if not os.path.exists(modern_bin):
            print(f"{RED}ERROR: Modern binary not found: {modern_bin}{RESET}")
            print(f"  Please build the project first with 'make release'{RESET}")
            return False
        
        # Step 1: Run find_pair
        print(f"{BLUE}--- Step 1: Running find_pair ---{RESET}")
        
        legacy_inp = os.path.join(temp_dir, f"{pdb_name}_legacy.inp")
        modern_inp = os.path.join(temp_dir, f"{pdb_name}_modern.inp")
        
        if legacy_bin:
            print("  Running legacy find_pair...")
            # Use absolute path for binary
            abs_bin = os.path.abspath(legacy_bin)
            # Run from project root, use relative path for PDB
            rel_pdb = os.path.relpath(pdb_file, project_root)
            rel_inp = os.path.relpath(legacy_inp, project_root)
            
            cmd = f"{abs_bin} --no-json {rel_pdb} {rel_inp}"
            success, out, err = run_command(cmd, cwd=project_root, capture_output=False)
            
            # Check for .inp file in various possible locations
            possible_locations = [
                legacy_inp,
                os.path.join(os.path.dirname(pdb_file), f"{pdb_name}.inp"),
                os.path.join(project_root, f"{pdb_name}.inp"),
                os.path.join(project_root, os.path.basename(rel_inp)),
                rel_inp,  # Relative path from project root
            ]
            
            found = False
            for loc in possible_locations:
                if os.path.exists(loc) and loc != legacy_inp:
                    if not found:
                        shutil.move(loc, legacy_inp)
                        found = True
                    else:
                        # Remove duplicate
                        os.remove(loc)
            
            if not os.path.exists(legacy_inp):
                print(f"    {YELLOW}Warning: Legacy find_pair may have failed - .inp file not found{RESET}")
                print(f"      Checked locations: {possible_locations[:3]}")
        
        if os.path.exists(modern_bin):
            print("  Running modern find_pair...")
            success, out, err = run_command(
                f"{modern_bin} {pdb_file} {modern_inp}",
                capture_output=False
            )
            if not success:
                print(f"    {YELLOW}Warning: Modern find_pair may have failed{RESET}")
        
        # Step 2: Compare .inp files
        print(f"\n{BLUE}--- Step 2: Comparing .inp files ---{RESET}")
        match, msg = compare_inp_files(legacy_inp, modern_inp)
        if match:
            print(f"  {GREEN}✓{RESET} {msg}")
        else:
            print(f"  {RED}✗{RESET} {msg}")
        
        # Step 3: Run analyze
        print(f"\n{BLUE}--- Step 3: Running analyze ---{RESET}")
        
        # Create working directories for analyze
        legacy_work = os.path.join(temp_dir, "legacy_work")
        modern_work = os.path.join(temp_dir, "modern_work")
        os.makedirs(legacy_work, exist_ok=True)
        os.makedirs(modern_work, exist_ok=True)
        
        legacy_work_inp = os.path.join(legacy_work, f"{pdb_name}.inp")
        modern_work_inp = os.path.join(modern_work, f"{pdb_name}.inp")
        
        if os.path.exists(legacy_inp):
            shutil.copy(legacy_inp, legacy_work_inp)
        if os.path.exists(modern_inp):
            shutil.copy(modern_inp, modern_work_inp)
        
        legacy_step_par = os.path.join(legacy_work, "bp_step.par")
        legacy_helix_par = os.path.join(legacy_work, "bp_helical.par")
        modern_step_par = os.path.join(temp_dir, f"{pdb_name}_modern_step.txt")
        modern_helix_par = os.path.join(temp_dir, f"{pdb_name}_modern_helix.txt")
        
        if legacy_analyze and os.path.exists(legacy_analyze) and os.path.exists(legacy_work_inp):
            print("  Running legacy analyze...")
            success, out, err = run_command(
                f"../../org/build/bin/analyze --no-json {os.path.basename(legacy_work_inp)}",
                cwd=legacy_work,
                capture_output=False
            )
        
        if os.path.exists(modern_analyze) and os.path.exists(modern_work_inp):
            print("  Running modern analyze...")
            success, out, err = run_command(
                f"{modern_analyze} {modern_work_inp}",
                capture_output=True
            )
            if success:
                # Extract step parameters
                step_lines = []
                helix_lines = []
                in_step = False
                in_helix = False
                for line in out.split('\n'):
                    if "Step Parameters" in line:
                        in_step = True
                        in_helix = False
                        continue
                    if "Helical Parameters" in line:
                        in_step = False
                        in_helix = True
                        continue
                    if in_step and line.strip() and not line.startswith('='):
                        step_lines.append(line)
                    if in_helix and line.strip() and not line.startswith('='):
                        helix_lines.append(line)
                
                with open(modern_step_par, 'w') as f:
                    f.write('\n'.join(step_lines))
                with open(modern_helix_par, 'w') as f:
                    f.write('\n'.join(helix_lines))
        
        # Step 4: Compare parameters
        print(f"\n{BLUE}--- Step 4: Comparing Step Parameters ---{RESET}")
        match, msg = compare_par_files(legacy_step_par, modern_step_par, "Step")
        if match:
            print(f"  {GREEN}✓{RESET} {msg}")
        else:
            print(f"  {RED}✗{RESET} {msg}")
        
        print(f"\n{BLUE}--- Step 5: Comparing Helical Parameters ---{RESET}")
        match, msg = compare_par_files(legacy_helix_par, modern_helix_par, "Helical")
        if match:
            print(f"  {GREEN}✓{RESET} {msg}")
        else:
            print(f"  {RED}✗{RESET} {msg}")
        
        print(f"\n{YELLOW}Temp directory: {temp_dir}{RESET}")
        print(f"  (Files preserved for manual inspection)")
        
    except Exception as e:
        print(f"{RED}ERROR: {e}{RESET}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

def main():
    if len(sys.argv) < 2:
        print("Usage: compare_output_no_json.py <pdb1> [pdb2] [pdb3] ...")
        print("Example: compare_output_no_json.py 6V9Q 7EH2")
        sys.exit(1)
    
    project_root = Path(__file__).parent.parent
    pdbs = sys.argv[1:]
    
    for pdb in pdbs:
        compare_pdb(pdb, str(project_root))
        print()

if __name__ == "__main__":
    main()

