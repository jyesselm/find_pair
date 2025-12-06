#!/usr/bin/env python3
"""Quick H-bond validation on first 200 PDBs with 10 workers."""

import json
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

def validate_one(name: str) -> tuple[str, str, dict]:
    """Validate one PDB."""
    try:
        result = subprocess.run(
            ["./build/generate_modern_json", f"data/pdb/{name}.pdb", "data/json", "--stage=hbonds"],
            capture_output=True, text=True, timeout=30
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
            "only_modern": sorted(list(modern_pairs - legacy_pairs))[:5],
            "only_legacy": sorted(list(legacy_pairs - modern_pairs))[:5]
        }
    except subprocess.TimeoutExpired:
        return name, "TIMEOUT", {}
    except Exception as e:
        return name, "ERROR", {"msg": str(e)[:50]}

# Load first 200 PDBs
with open("data/valid_pdbs_fast.json") as f:
    pdbs = json.load(f)["valid_pdbs_with_atoms_and_frames"][:200]

print(f"Testing {len(pdbs)} PDBs with 10 workers...", flush=True)

failures = []
with ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(validate_one, name): name for name in pdbs}
    done = 0
    for future in as_completed(futures):
        name, status, details = future.result()
        done += 1
        if status not in ["PASS", "NO_DATA"]:
            failures.append({"pdb": name, "status": status, **details})
        if done % 50 == 0:
            print(f"Progress: {done}/200", flush=True)

print(f"\n{'='*60}", flush=True)
print(f"Failures: {len(failures)}", flush=True)
for f in failures:
    print(json.dumps(f), flush=True)

