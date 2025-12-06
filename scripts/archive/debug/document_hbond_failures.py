#!/usr/bin/env python3
"""Document all H-bond validation failures."""

import json
import subprocess
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def validate_one(name: str) -> tuple[str, str, dict]:
    """Validate one PDB, return (name, status, details)."""
    pdb_path = f"data/pdb/{name}.pdb"
    details = {}
    
    try:
        result = subprocess.run(
            ["./build/generate_modern_json", pdb_path, "data/json", "--stage=hbonds"],
            capture_output=True, text=True, timeout=60
        )
        
        if result.returncode != 0:
            return name, "GEN_FAIL", {"error": result.stderr[:100]}
        
        legacy_path = f"data/json_legacy/hbond_list/{name}.json"
        modern_path = f"data/json/hbond_list/{name}.json"
        
        if not os.path.exists(legacy_path):
            return name, "NO_LEGACY", {}
        if not os.path.exists(modern_path):
            return name, "NO_MODERN", {}
        
        with open(modern_path) as f:
            modern = json.load(f)
        with open(legacy_path) as f:
            legacy = json.load(f)
        
        def normalize(p):
            return (min(p), max(p))
        
        modern_pairs = set(normalize((e['base_i'], e['base_j'])) for e in modern)
        legacy_pairs = set(normalize((e['base_i']-1, e['base_j']-1)) for e in legacy)
        
        # Clean up
        try:
            os.unlink(modern_path)
        except:
            pass
        
        if modern_pairs == legacy_pairs:
            return name, "PASS", {}
        
        # Document the differences
        only_modern = modern_pairs - legacy_pairs
        only_legacy = legacy_pairs - modern_pairs
        
        details = {
            "modern_count": len(modern_pairs),
            "legacy_count": len(legacy_pairs),
            "only_in_modern": sorted(list(only_modern))[:10],
            "only_in_legacy": sorted(list(only_legacy))[:10],
        }
        
        return name, "MISMATCH", details
        
    except subprocess.TimeoutExpired:
        return name, "TIMEOUT", {}
    except Exception as e:
        return name, f"ERROR", {"error": str(e)[:100]}

def main():
    # Load fast PDBs
    with open("data/valid_pdbs_fast.json") as f:
        data = json.load(f)
    pdbs = data["valid_pdbs_with_atoms_and_frames"]
    
    print(f"Checking {len(pdbs)} PDBs for failures...", flush=True)
    
    failures = []
    completed = 0
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(validate_one, name): name for name in pdbs}
        
        for future in as_completed(futures):
            name, status, details = future.result()
            completed += 1
            
            if status != "PASS":
                failures.append({"pdb": name, "status": status, **details})
            
            if completed % 500 == 0:
                print(f"Progress: {completed}/{len(pdbs)}", flush=True)
    
    # Sort failures by type
    failures.sort(key=lambda x: (x["status"], x["pdb"]))
    
    # Write to file
    with open("data/hbond_failures.json", "w") as f:
        json.dump(failures, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"Total failures: {len(failures)}")
    
    # Group by status
    by_status = {}
    for f in failures:
        status = f["status"]
        if status not in by_status:
            by_status[status] = []
        by_status[status].append(f)
    
    for status, items in sorted(by_status.items()):
        print(f"\n{status} ({len(items)}):")
        for item in items[:10]:
            if status == "MISMATCH":
                print(f"  {item['pdb']}: modern={item['modern_count']}, legacy={item['legacy_count']}")
                if item.get('only_in_modern'):
                    print(f"    Only in modern: {item['only_in_modern'][:3]}")
                if item.get('only_in_legacy'):
                    print(f"    Only in legacy: {item['only_in_legacy'][:3]}")
            else:
                print(f"  {item['pdb']}: {item.get('error', '')[:50]}")
        if len(items) > 10:
            print(f"  ... and {len(items)-10} more")
    
    print(f"\nResults saved to data/hbond_failures.json")

if __name__ == "__main__":
    main()

