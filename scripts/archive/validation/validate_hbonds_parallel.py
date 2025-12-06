#!/usr/bin/env python3
"""Parallel H-bond validation script."""

import json
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def validate_pdb(pdb_path: str) -> tuple[str, str]:
    """Validate a single PDB's H-bond output."""
    pdb_path = Path(pdb_path)
    name = pdb_path.stem
    
    try:
        # Generate modern JSON
        result = subprocess.run(
            ["./build/generate_modern_json", str(pdb_path), "data/json", "--stage=hbonds"],
            capture_output=True, text=True, timeout=60
        )
        
        if result.returncode != 0:
            return name, "GEN_FAIL"
        
        # Check if legacy exists
        legacy_path = Path(f"data/json_legacy/hbond_list/{name}.json")
        modern_path = Path(f"data/json/hbond_list/{name}.json")
        
        if not legacy_path.exists():
            return name, "NO_LEGACY"
        
        if not modern_path.exists():
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
        
        if modern_pairs == legacy_pairs:
            return name, "PASS"
        else:
            return name, f"MISMATCH:m={len(modern_pairs)},l={len(legacy_pairs)}"
            
    except Exception as e:
        return name, f"ERROR:{str(e)[:30]}"

def main():
    # Get all PDBs
    pdb_dir = Path("data/pdb")
    pdbs = sorted(pdb_dir.glob("*.pdb"))
    total = len(pdbs)
    
    print(f"Validating {total} PDBs with 10 workers...")
    print()
    
    results = {"PASS": 0, "NO_LEGACY": 0, "MISMATCH": 0, "ERROR": 0, "OTHER": 0}
    failures = []
    
    with ProcessPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(validate_pdb, str(p)): p for p in pdbs}
        
        completed = 0
        for future in as_completed(futures):
            name, status = future.result()
            completed += 1
            
            if status == "PASS":
                results["PASS"] += 1
            elif status == "NO_LEGACY":
                results["NO_LEGACY"] += 1
            elif status.startswith("MISMATCH"):
                results["MISMATCH"] += 1
                failures.append((name, status))
            elif status.startswith("ERROR"):
                results["ERROR"] += 1
                failures.append((name, status))
            else:
                results["OTHER"] += 1
            
            # Progress every 100
            if completed % 100 == 0 or completed == total:
                pct = 100 * completed / total
                print(f"Progress: {completed}/{total} ({pct:.1f}%) - "
                      f"PASS:{results['PASS']} NO_LEGACY:{results['NO_LEGACY']} "
                      f"MISMATCH:{results['MISMATCH']} ERROR:{results['ERROR']}")
    
    print()
    print("=" * 60)
    print(f"FINAL: {results['PASS']}/{results['PASS'] + results['MISMATCH']} passed")
    print(f"  - Passed: {results['PASS']}")
    print(f"  - No legacy data: {results['NO_LEGACY']}")
    print(f"  - Mismatches: {results['MISMATCH']}")
    print(f"  - Errors: {results['ERROR']}")
    
    if failures:
        print()
        print("Failures (first 10):")
        for name, status in failures[:10]:
            print(f"  {name}: {status}")

if __name__ == "__main__":
    main()

