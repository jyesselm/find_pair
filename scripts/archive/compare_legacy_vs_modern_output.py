#!/usr/bin/env python3
"""
Compare legacy and modern code outputs for find_pair and analyze.
"""

import sys
import os
import subprocess
import tempfile
import json
from pathlib import Path

def run_legacy_find_pair(pdb_file, output_dir):
    """Run legacy find_pair and return .inp file path."""
    # Check if legacy binary exists
    legacy_bin = Path(__file__).parent.parent / "org" / "build" / "bin" / "find_pair"
    if not legacy_bin.exists():
        print(f"Warning: Legacy find_pair not found at {legacy_bin}")
        return None
    
    inp_file = output_dir / "legacy_output.inp"
    
    # Run legacy find_pair
    cmd = [str(legacy_bin), pdb_file]
    print(f"Running legacy find_pair: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(output_dir.parent),
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # Legacy find_pair creates .inp file based on PDB name
        pdb_stem = Path(pdb_file).stem
        default_inp = output_dir.parent / f"{pdb_stem}.inp"
        
        if default_inp.exists():
            default_inp.rename(inp_file)
            print(f"Legacy find_pair output: {inp_file}")
            return inp_file
        else:
            print(f"Warning: Legacy find_pair did not create {default_inp}")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return None
    except Exception as e:
        print(f"Error running legacy find_pair: {e}")
        return None

def run_modern_find_pair(pdb_file, output_dir):
    """Run modern find_pair_app and return .inp file path."""
    modern_bin = Path(__file__).parent.parent / "build" / "find_pair_app"
    if not modern_bin.exists():
        print(f"Error: Modern find_pair_app not found at {modern_bin}")
        return None
    
    inp_file = output_dir / "modern_output.inp"
    
    cmd = [str(modern_bin), pdb_file, str(inp_file)]
    print(f"Running modern find_pair_app: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0 and inp_file.exists():
            print(f"Modern find_pair_app output: {inp_file}")
            return inp_file
        else:
            print(f"Error: Modern find_pair_app failed")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return None
    except Exception as e:
        print(f"Error running modern find_pair_app: {e}")
        return None

def parse_inp_file(inp_file):
    """Parse .inp file and return structure."""
    if not inp_file or not inp_file.exists():
        return None
    
    data = {
        "pdb_file": None,
        "output_file": None,
        "duplex_number": None,
        "num_base_pairs": None,
        "flags": None,
        "base_pairs": []
    }
    
    with open(inp_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 5:
        return data
    
    # Line 1: PDB file
    data["pdb_file"] = lines[0].strip()
    
    # Line 2: Output file
    data["output_file"] = lines[1].strip()
    
    # Line 3: Duplex number
    try:
        data["duplex_number"] = int(lines[2].split()[0])
    except:
        pass
    
    # Line 4: Number of base pairs
    try:
        data["num_base_pairs"] = int(lines[3].split()[0])
    except:
        pass
    
    # Line 5: Flags
    try:
        parts = lines[4].split()
        data["flags"] = int(parts[0]) if parts else None
    except:
        pass
    
    # Remaining lines: Base pairs
    for line in lines[5:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) >= 3:
            try:
                bp_num = int(parts[0])
                res1 = int(parts[1])
                res2 = int(parts[2])
                flag = int(parts[3]) if len(parts) > 3 else 0
                comment = ' '.join(parts[4:]) if len(parts) > 4 else ""
                data["base_pairs"].append({
                    "bp_num": bp_num,
                    "residue1": res1,
                    "residue2": res2,
                    "flag": flag,
                    "comment": comment
                })
            except:
                pass
    
    return data

def compare_inp_files(legacy_inp, modern_inp):
    """Compare two .inp files."""
    print("\n" + "="*70)
    print("COMPARING .inp FILES")
    print("="*70)
    
    legacy_data = parse_inp_file(legacy_inp)
    modern_data = parse_inp_file(modern_inp)
    
    if not legacy_data or not modern_data:
        print("Error: Could not parse one or both .inp files")
        return False
    
    differences = []
    matches = []
    
    # Compare structure
    if legacy_data["num_base_pairs"] != modern_data["num_base_pairs"]:
        differences.append(
            f"Base pair count: Legacy={legacy_data['num_base_pairs']}, "
            f"Modern={modern_data['num_base_pairs']}"
        )
    else:
        matches.append(f"Base pair count: {legacy_data['num_base_pairs']}")
    
    if legacy_data["duplex_number"] != modern_data["duplex_number"]:
        differences.append(
            f"Duplex number: Legacy={legacy_data['duplex_number']}, "
            f"Modern={modern_data['duplex_number']}"
        )
    
    # Compare base pairs
    if len(legacy_data["base_pairs"]) != len(modern_data["base_pairs"]):
        differences.append(
            f"Number of base pair lines: Legacy={len(legacy_data['base_pairs'])}, "
            f"Modern={len(modern_data['base_pairs'])}"
        )
    else:
        matches.append(f"Number of base pair lines: {len(legacy_data['base_pairs'])}")
        
        # Compare each base pair
        for i, (leg_bp, mod_bp) in enumerate(zip(legacy_data["base_pairs"], modern_data["base_pairs"])):
            if leg_bp["residue1"] != mod_bp["residue1"] or leg_bp["residue2"] != mod_bp["residue2"]:
                differences.append(
                    f"Base pair {i+1}: Legacy=({leg_bp['residue1']}, {leg_bp['residue2']}), "
                    f"Modern=({mod_bp['residue1']}, {mod_bp['residue2']})"
                )
    
    # Print results
    if matches:
        print("\n✅ MATCHES:")
        for match in matches:
            print(f"  {match}")
    
    if differences:
        print("\n❌ DIFFERENCES:")
        for diff in differences:
            print(f"  {diff}")
    else:
        print("\n✅ All .inp file fields match!")
    
    return len(differences) == 0

def run_legacy_analyze(inp_file, output_dir):
    """Run legacy analyze and return output file."""
    legacy_bin = Path(__file__).parent.parent / "org" / "build" / "bin" / "analyze"
    if not legacy_bin.exists():
        print(f"Warning: Legacy analyze not found at {legacy_bin}")
        return None
    
    cmd = [str(legacy_bin), str(inp_file)]
    print(f"\nRunning legacy analyze: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(output_dir.parent),
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # Legacy analyze creates .outp file
        # Extract output filename from .inp file
        inp_data = parse_inp_file(inp_file)
        if inp_data and inp_data["output_file"]:
            outp_file = output_dir.parent / inp_data["output_file"]
            if outp_file.exists():
                print(f"Legacy analyze output: {outp_file}")
                return outp_file
        
        print(f"STDOUT: {result.stdout[:500]}")
        if result.stderr:
            print(f"STDERR: {result.stderr[:500]}")
        return None
    except Exception as e:
        print(f"Error running legacy analyze: {e}")
        return None

def run_modern_analyze(inp_file, output_dir):
    """Run modern analyze_app."""
    modern_bin = Path(__file__).parent.parent / "build" / "analyze_app"
    if not modern_bin.exists():
        print(f"Error: Modern analyze_app not found at {modern_bin}")
        return False
    
    cmd = [str(modern_bin), str(inp_file)]
    print(f"\nRunning modern analyze_app: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            print("Modern analyze_app completed successfully")
            print(f"STDOUT: {result.stdout}")
            return True
        else:
            print(f"Error: Modern analyze_app failed")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return False
    except Exception as e:
        print(f"Error running modern analyze_app: {e}")
        return False

def main():
    if len(sys.argv) < 2:
        print("Usage: compare_legacy_vs_modern_output.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = Path(sys.argv[1])
    if not pdb_file.exists():
        print(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)
    
    # Create temporary output directory
    output_dir = Path(tempfile.mkdtemp(prefix="x3dna_compare_"))
    print(f"Using output directory: {output_dir}")
    
    try:
        # Step 1: Compare find_pair outputs
        print("\n" + "="*70)
        print("STEP 1: COMPARING find_pair OUTPUTS")
        print("="*70)
        
        legacy_inp = run_legacy_find_pair(pdb_file, output_dir)
        modern_inp = run_modern_find_pair(pdb_file, output_dir)
        
        if legacy_inp and modern_inp:
            inp_match = compare_inp_files(legacy_inp, modern_inp)
        else:
            print("Skipping .inp comparison (one or both files missing)")
            inp_match = False
        
        # Step 2: Run analyze on both .inp files
        print("\n" + "="*70)
        print("STEP 2: RUNNING analyze")
        print("="*70)
        
        if legacy_inp:
            legacy_outp = run_legacy_analyze(legacy_inp, output_dir)
        
        if modern_inp:
            modern_analyze_success = run_modern_analyze(modern_inp, output_dir)
        
        # Summary
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        print(f".inp file comparison: {'✅ PASS' if inp_match else '❌ FAIL'}")
        
        print(f"\nOutput files saved to: {output_dir}")
        print("You can manually compare the outputs for detailed differences.")
        
    finally:
        # Keep output directory for manual inspection
        print(f"\nOutput directory preserved: {output_dir}")

if __name__ == "__main__":
    main()

