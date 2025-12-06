#!/usr/bin/env python3
"""Find first H-bond mismatch and stop for investigation."""

import json
import subprocess
import os

# Skip known problematic PDBs
SKIP = {"1F8V", "1FFZ", "9CJI"}

# Load all fast PDBs
with open("data/valid_pdbs_fast.json") as f:
    pdbs = json.load(f)["valid_pdbs_with_atoms_and_frames"]

print(f"Testing {len(pdbs)} PDBs (stopping on first mismatch)...", flush=True)

for i, name in enumerate(pdbs):
    if name in SKIP:
        continue
    
    if (i+1) % 100 == 0:
        print(f"Progress: {i+1}/{len(pdbs)}", flush=True)
    
    try:
        result = subprocess.run(
            ["./build/generate_modern_json", f"data/pdb/{name}.pdb", "data/json", "--stage=hbonds"],
            capture_output=True, text=True, timeout=60
        )
        
        if result.returncode != 0:
            continue
        
        legacy_path = f"data/json_legacy/hbond_list/{name}.json"
        modern_path = f"data/json/hbond_list/{name}.json"
        
        if not os.path.exists(legacy_path) or not os.path.exists(modern_path):
            continue
        
        with open(modern_path) as f:
            modern = json.load(f)
        with open(legacy_path) as f:
            legacy = json.load(f)
        
        def normalize(p): return (min(p), max(p))
        
        modern_pairs = set(normalize((e['base_i'], e['base_j'])) for e in modern)
        legacy_pairs = set(normalize((e['base_i']-1, e['base_j']-1)) for e in legacy)
        
        if modern_pairs != legacy_pairs:
            print(f"\n{'='*60}")
            print(f"MISMATCH FOUND: {name}")
            print(f"Modern pairs: {len(modern_pairs)}")
            print(f"Legacy pairs: {len(legacy_pairs)}")
            
            only_modern = sorted(list(modern_pairs - legacy_pairs))
            only_legacy = sorted(list(legacy_pairs - modern_pairs))
            
            print(f"\nOnly in modern ({len(only_modern)}):")
            for p in only_modern[:10]:
                print(f"  {p}")
            
            print(f"\nOnly in legacy ({len(only_legacy)}):")
            for p in only_legacy[:10]:
                print(f"  {p}")
            
            # Keep modern file for investigation
            print(f"\nModern JSON kept at: {modern_path}")
            print(f"Legacy JSON at: {legacy_path}")
            break
        else:
            os.unlink(modern_path)
            
    except subprocess.TimeoutExpired:
        continue
    except Exception as e:
        continue
else:
    print("No mismatches found!")

