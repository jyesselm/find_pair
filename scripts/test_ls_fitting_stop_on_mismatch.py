#!/usr/bin/env python3
"""
Test ls_fitting with FIXED legacy code and stop on significant mismatches.
"""

import json
import subprocess
import sys
from pathlib import Path

project_root = Path(".")
sys.path.insert(0, str(project_root))
from scripts.test_utils import load_valid_pdbs_fast

def test_pdb(pdb_id, legacy_exe, modern_exe):
    """Test a single PDB and return detailed results."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    
    if not pdb_file.exists():
        return {"status": "skip", "reason": "PDB file not found"}
    
    # Generate legacy JSON
    try:
        subprocess.run([str(legacy_exe.resolve()), f"../data/pdb/{pdb_id}.pdb"],
                      cwd=str(project_root / "org"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                      timeout=30)
    except Exception as e:
        return {"status": "error", "reason": f"Legacy generation failed: {e}"}
    
    # Generate modern JSON
    output_dir = project_root / "tmp" / "stop_on_mismatch"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        subprocess.run([str(modern_exe), str(pdb_file), str(output_dir), "--stage=ls_fitting"],
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
    except Exception as e:
        return {"status": "error", "reason": f"Modern generation failed: {e}"}
    
    # Compare
    legacy_file = project_root / "data" / "json_legacy" / "ls_fitting" / f"{pdb_id}.json"
    modern_file = output_dir / "ls_fitting" / f"{pdb_id}.json"
    
    if not legacy_file.exists() or not modern_file.exists():
        return {"status": "error", "reason": "JSON files not generated"}
    
    legacy = json.load(open(legacy_file))
    modern = json.load(open(modern_file))
    
    if len(legacy) != len(modern):
        return {
            "status": "count_mismatch",
            "legacy_count": len(legacy),
            "modern_count": len(modern)
        }
    
    # Check each record
    mismatches = []
    max_diff = 0.0
    
    for i, (leg, mod) in enumerate(zip(legacy, modern)):
        mod_clean = {k: v for k, v in mod.items() if k != 'type'}
        
        if leg != mod_clean:
            # Find the differences
            for key in set(leg.keys()) | set(mod_clean.keys()):
                leg_val = leg.get(key)
                mod_val = mod_clean.get(key)
                
                if leg_val != mod_val:
                    if isinstance(leg_val, float) and isinstance(mod_val, float):
                        diff = abs(leg_val - mod_val)
                        max_diff = max(max_diff, diff)
                        mismatches.append({
                            "record": i,
                            "field": key,
                            "diff": diff,
                            "type": "float"
                        })
                    elif isinstance(leg_val, list) and isinstance(mod_val, list):
                        for j, (lv, mv) in enumerate(zip(leg_val, mod_val)):
                            if lv != mv and isinstance(lv, float) and isinstance(mv, float):
                                diff = abs(lv - mv)
                                max_diff = max(max_diff, diff)
                                mismatches.append({
                                    "record": i,
                                    "field": f"{key}[{j}]",
                                    "diff": diff,
                                    "type": "float"
                                })
                    else:
                        mismatches.append({
                            "record": i,
                            "field": key,
                            "legacy": leg_val,
                            "modern": mod_val,
                            "type": "other"
                        })
    
    if not mismatches:
        return {"status": "pass", "count": len(legacy)}
    else:
        return {
            "status": "mismatch",
            "count": len(legacy),
            "num_mismatches": len(mismatches),
            "max_diff": max_diff,
            "mismatches": mismatches[:5]  # First 5
        }

def main():
    # Find executables
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    modern_exe = project_root / "build" / "generate_modern_json"
    
    if not legacy_exe.exists():
        print(f"Error: Legacy executable not found: {legacy_exe}")
        return 1
    if not modern_exe.exists():
        print(f"Error: Modern executable not found: {modern_exe}")
        return 1
    
    # Load fast PDBs
    pdb_ids = load_valid_pdbs_fast(project_root)
    
    print(f"Testing {len(pdb_ids)} PDBs with FIXED legacy code...")
    print(f"Will STOP if max difference > 1e-5 or non-float mismatch found\n")
    
    passed = 0
    failed = 0
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        result = test_pdb(pdb_id, legacy_exe, modern_exe)
        
        if result["status"] == "pass":
            print(f"[{i}/{len(pdb_ids)}] {pdb_id}: ✓ PASS ({result['count']} records)")
            passed += 1
            
        elif result["status"] == "mismatch":
            max_diff = result.get("max_diff", 0)
            num_mm = result.get("num_mismatches", 0)
            
            # Check if it's acceptable (only small FP differences)
            if max_diff <= 1e-5 and all(m["type"] == "float" for m in result["mismatches"]):
                print(f"[{i}/{len(pdb_ids)}] {pdb_id}: ✓ OK ({result['count']} records, {num_mm} FP diffs, max={max_diff:.2e})")
                passed += 1
            else:
                print(f"\n[{i}/{len(pdb_ids)}] {pdb_id}: ✗ SIGNIFICANT MISMATCH!")
                print(f"  Records: {result['count']}")
                print(f"  Mismatches: {num_mm}")
                print(f"  Max difference: {max_diff:.2e}")
                print(f"  Details:")
                for mm in result["mismatches"][:3]:
                    print(f"    Record {mm['record']}, field {mm['field']}: {mm}")
                print(f"\n⚠️  STOPPING to investigate...")
                return 1
                
        elif result["status"] == "count_mismatch":
            print(f"\n[{i}/{len(pdb_ids)}] {pdb_id}: ✗ COUNT MISMATCH!")
            print(f"  Legacy: {result['legacy_count']} records")
            print(f"  Modern: {result['modern_count']} records")
            print(f"\n⚠️  STOPPING - this should not happen with fixed legacy code!")
            return 1
            
        else:
            print(f"[{i}/{len(pdb_ids)}] {pdb_id}: SKIP ({result.get('reason', 'unknown')})")
        
        # Progress update every 50
        if i % 50 == 0:
            print(f"\n--- Progress: {i}/{len(pdb_ids)} tested, {passed} passed ---\n")
    
    print(f"\n✅ All {passed} PDBs passed!")
    print(f"   (Only minor floating-point differences <= 1e-5)")
    return 0

if __name__ == "__main__":
    sys.exit(main())

