#!/usr/bin/env python3
"""Find only actual MISMATCH failures in H-bond validation."""

import json
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

# Skip known problematic PDBs
SKIP = {"1F8V", "1FFZ", "9CJI"}  # From stage2_exclusions.json

def validate_one(name: str) -> tuple[str, str, dict]:
    """Validate one PDB."""
    if name in SKIP:
        return name, "SKIP", {}
    
    try:
        result = subprocess.run(
            ["./build/generate_modern_json", f"data/pdb/{name}.pdb", "data/json", "--stage=hbonds"],
            capture_output=True, text=True, timeout=60
        )
        
        if result.returncode != 0:
            return name, "GEN_FAIL", {}
        
        legacy_path = f"data/json_legacy/hbond_list/{name}.json"
        modern_path = f"data/json/hbond_list/{name}.json"
        
        if not os.path.exists(legacy_path) or not os.path.exists(modern_path):
            return name, "NO_DATA", {}
        
        with open(modern_path) as f:
            modern = json.load(f)
        with open(legacy_path) as f:
            legacy = json.load(f)
        
        def normalize(p): return (min(p), max(p))
        
        modern_pairs = set(normalize((e['base_i'], e['base_j'])) for e in modern)
        legacy_pairs = set(normalize((e['base_i']-1, e['base_j']-1)) for e in legacy)
        
        try:
            os.unlink(modern_path)
        except:
            pass
        
        if modern_pairs == legacy_pairs:
            return name, "PASS", {}
        
        return name, "MISMATCH", {
            "modern": len(modern_pairs),
            "legacy": len(legacy_pairs),
            "diff": len(modern_pairs) - len(legacy_pairs),
            "only_modern": sorted(list(modern_pairs - legacy_pairs)),
            "only_legacy": sorted(list(legacy_pairs - modern_pairs))
        }
    except subprocess.TimeoutExpired:
        return name, "TIMEOUT", {}
    except Exception as e:
        return name, "ERROR", {"msg": str(e)[:80]}

# Load all fast PDBs
with open("data/valid_pdbs_fast.json") as f:
    pdbs = json.load(f)["valid_pdbs_with_atoms_and_frames"]

print(f"Testing {len(pdbs)} PDBs (skipping {len(SKIP)} known issues)...", flush=True)

results = {"PASS": 0, "MISMATCH": [], "TIMEOUT": [], "ERROR": [], "SKIP": 0, "NO_DATA": 0}

with ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(validate_one, name): name for name in pdbs}
    done = 0
    for future in as_completed(futures):
        name, status, details = future.result()
        done += 1
        
        if status == "PASS":
            results["PASS"] += 1
        elif status == "MISMATCH":
            results["MISMATCH"].append({"pdb": name, **details})
        elif status == "TIMEOUT":
            results["TIMEOUT"].append(name)
        elif status == "ERROR":
            results["ERROR"].append({"pdb": name, **details})
        elif status == "SKIP":
            results["SKIP"] += 1
        else:
            results["NO_DATA"] += 1
        
        if done % 100 == 0:
            print(f"Progress: {done}/{len(pdbs)} - PASS:{results['PASS']} MISMATCH:{len(results['MISMATCH'])}", flush=True)

print(f"\n{'='*60}", flush=True)
print(f"PASS: {results['PASS']}", flush=True)
print(f"MISMATCH: {len(results['MISMATCH'])}", flush=True)
print(f"TIMEOUT: {len(results['TIMEOUT'])}", flush=True)
print(f"ERROR: {len(results['ERROR'])}", flush=True)
print(f"SKIP: {results['SKIP']}", flush=True)

if results["MISMATCH"]:
    print(f"\n=== MISMATCH DETAILS ===", flush=True)
    for m in results["MISMATCH"]:
        print(f"\n{m['pdb']}: modern={m['modern']}, legacy={m['legacy']} (diff={m['diff']})", flush=True)
        if m['only_modern']:
            print(f"  Only in modern: {m['only_modern'][:5]}", flush=True)
        if m['only_legacy']:
            print(f"  Only in legacy: {m['only_legacy'][:5]}", flush=True)

# Save results
with open("data/hbond_validation_results.json", "w") as f:
    json.dump({
        "pass_count": results["PASS"],
        "mismatches": results["MISMATCH"],
        "timeouts": results["TIMEOUT"],
        "errors": results["ERROR"]
    }, f, indent=2)

print(f"\nResults saved to data/hbond_validation_results.json", flush=True)

