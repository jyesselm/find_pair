#!/usr/bin/env python3
import json
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
legacy_file = project_root / "data/json_legacy/1H4S.json"
modern_file = project_root / "data/json/1H4S_legacy.json"

if not modern_file.exists():
    print(f"ERROR: {modern_file} not found")
    print("Run: ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
    sys.exit(1)

with open(legacy_file) as f:
    legacy = json.load(f)
with open(modern_file) as f:
    modern = json.load(f)

leg_calcs = {c.get('residue_seq', 0): c for c in legacy['calculations'] if c.get('type') == 'base_frame_calc'}
mod_calcs = {(c.get('chain_id', ''), c.get('residue_seq', 0)): c for c in modern['calculations'] if c.get('type') == 'base_frame_calc'}

print(f"Legacy calculations: {len(leg_calcs)}")
print(f"Modern calculations: {len(mod_calcs)}")

exact = 0
diffs = 0

for seq, leg_calc in list(leg_calcs.items())[:10]:
    chain = leg_calc.get('chain_id', '')
    mod_calc = mod_calcs.get((chain, seq))
    if not mod_calc:
        continue
    
    leg_atoms = leg_calc.get('matched_atoms', [])
    mod_atoms = mod_calc.get('matched_atoms', [])
    
    if leg_atoms == mod_atoms:
        exact += 1
    else:
        diffs += 1
        if diffs <= 3:
            print(f"\nDiff: {chain}:{seq} ({leg_calc.get('base_type')})")
            print(f"  Legacy: {leg_atoms}")
            print(f"  Modern: {mod_atoms}")

print(f"\nExact matches: {exact}")
print(f"Differences: {diffs}")

