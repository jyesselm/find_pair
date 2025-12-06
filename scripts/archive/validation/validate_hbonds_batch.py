#!/usr/bin/env python3
"""Batch H-bond validation - 10 workers, uses valid_pdbs_fast.json."""

import json
import subprocess
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def validate_one(name: str) -> tuple[str, str]:
    """Validate one PDB, return (name, status)."""
    pdb_path = f"data/pdb/{name}.pdb"
    
    try:
        # Generate modern JSON (hbonds only)
        result = subprocess.run(
            ["./build/generate_modern_json", pdb_path, "data/json", "--stage=hbonds"],
            capture_output=True, text=True, timeout=60
        )
        
        if result.returncode != 0:
            return name, "GEN_FAIL"
        
        legacy_path = f"data/json_legacy/hbond_list/{name}.json"
        modern_path = f"data/json/hbond_list/{name}.json"
        
        if not os.path.exists(legacy_path):
            return name, "NO_LEGACY"
        if not os.path.exists(modern_path):
            return name, "NO_MODERN"
        
        # Compare
        with open(modern_path) as f:
            modern = json.load(f)
        with open(legacy_path) as f:
            legacy = json.load(f)
        
        def normalize(p):
            return (min(p), max(p))
        
        modern_pairs = set(normalize((e['base_i'], e['base_j'])) for e in modern)
        legacy_pairs = set(normalize((e['base_i']-1, e['base_j']-1)) for e in legacy)
        
        # Clean up immediately
        try:
            os.unlink(modern_path)
        except:
            pass
        
        if modern_pairs == legacy_pairs:
            return name, "PASS"
        return name, f"MISMATCH:m={len(modern_pairs)},l={len(legacy_pairs)}"
        
    except subprocess.TimeoutExpired:
        return name, "TIMEOUT"
    except Exception as e:
        return name, f"ERROR:{str(e)[:30]}"

def main():
    # Load fast PDBs from valid_pdbs_fast.json
    with open("data/valid_pdbs_fast.json") as f:
        data = json.load(f)
    
    pdbs = data["valid_pdbs_with_atoms_and_frames"]
    total = len(pdbs)
    
    print(f"Validating {total} fast PDBs with 10 workers...", flush=True)
    print(flush=True)
    
    passed = 0
    failures = []
    no_legacy = 0
    completed = 0
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(validate_one, name): name for name in pdbs}
        
        for future in as_completed(futures):
            name, status = future.result()
            completed += 1
            
            if status == "PASS":
                passed += 1
            elif status == "NO_LEGACY":
                no_legacy += 1
            else:
                failures.append((name, status))
            
            # Progress every 100
            if completed % 100 == 0:
                print(f"Progress: {completed}/{total} - PASS:{passed} NO_LEGACY:{no_legacy} FAIL:{len(failures)}", flush=True)
    
    print(f"Progress: {completed}/{total} - PASS:{passed} NO_LEGACY:{no_legacy} FAIL:{len(failures)}", flush=True)
    print(flush=True)
    print("=" * 60, flush=True)
    
    testable = passed + len(failures)
    if testable > 0:
        print(f"FINAL: {passed}/{testable} passed ({100*passed/testable:.1f}%)", flush=True)
        print(f"  No legacy data: {no_legacy}", flush=True)
    
    if failures:
        print(f"\nFailures ({len(failures)}):", flush=True)
        for name, status in failures[:20]:
            print(f"  {name}: {status}", flush=True)
        if len(failures) > 20:
            print(f"  ... and {len(failures)-20} more", flush=True)

if __name__ == "__main__":
    main()
