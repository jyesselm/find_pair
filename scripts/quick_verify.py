#!/usr/bin/env python3
"""Quick verification script - no shell dependencies"""

import json
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent

legacy_file = project_root / "data" / "json_legacy" / "1H4S.json"
modern_file = project_root / "data" / "json" / "1H4S_legacy.json"

if not modern_file.exists():
    print(f"❌ Modern JSON file not found: {modern_file}")
    print("Please run: ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
    sys.exit(1)

with open(legacy_file) as f:
    legacy = json.load(f)
with open(modern_file) as f:
    modern = json.load(f)

leg_calcs = [c for c in legacy['calculations'] if c.get('type') == 'base_frame_calc']
mod_calcs = {f\"{c.get('chain_id', '')}:{c.get('residue_seq', 0)}\": c 
             for c in modern['calculations'] if c.get('type') == 'base_frame_calc'}

print(f"Legacy: {len(leg_calcs)} residues")
print(f"Modern: {len(mod_calcs)} residues")
print()

exact = 0
set_match = 0
real_diff = 0

for leg_calc in leg_calcs:
    chain = leg_calc.get('chain_id', '')
    seq = leg_calc.get('residue_seq', 0)
    key = f\"{chain}:{seq}\"
    if key not in mod_calcs:
        continue
    
    mod_calc = mod_calcs[key]
    leg_atoms = leg_calc.get('matched_atoms', [])
    mod_atoms = mod_calc.get('matched_atoms', [])
    
    if leg_atoms == mod_atoms:
        exact += 1
        set_match += 1
    elif set(leg_atoms) == set(mod_atoms):
        set_match += 1
    else:
        real_diff += 1
        if real_diff <= 3:  # Show first few differences
            print(f\"Difference {real_diff}: {leg_calc.get('base_type', '?')} {chain}:{seq}\")
            leg_only = sorted(set(leg_atoms) - set(mod_atoms))
            mod_only = sorted(set(mod_atoms) - set(leg_atoms))
            print(f\"  Legacy only: {leg_only}\")
            print(f\"  Modern only: {mod_only}\")
            print()

total = len([c for c in leg_calcs if f\"{c.get('chain_id', '')}:{c.get('residue_seq', 0)}\" in mod_calcs])

print(f"Exact matches: {exact}/{total} ({exact*100/total:.1f}%)")
print(f"Set matches: {set_match}/{total} ({set_match*100/total:.1f}%)")
print(f"Real differences: {real_diff}/{total} ({real_diff*100/total:.1f}%)")
print()

# Spot checks
leg_g = [c for c in leg_calcs if c.get('base_type')=='G' and c.get('chain_id')=='T' and c.get('residue_seq')==4]
if leg_g:
    leg_g = leg_g[0]
    key = f\"{leg_g.get('chain_id', '')}:{leg_g.get('residue_seq', 0)}\"
    if key in mod_calcs:
        mod_g = mod_calcs[key]
        print(f"Purine G T:4:")
        print(f"  Legacy: {leg_g['matched_atoms']}")
        print(f"  Modern: {mod_g['matched_atoms']}")
        print(f"  Match: {leg_g['matched_atoms'] == mod_g['matched_atoms']}")
        print()

if real_diff == 0 and exact == total:
    print("✅ PERFECT MATCH! All residues match exactly!")
    sys.exit(0)
elif real_diff == 0:
    print(f"✅ ATOM SET MATCH! ({exact} exact, {total-exact} order differences)")
    sys.exit(0)
else:
    print(f"❌ Still {real_diff} real differences")
    sys.exit(1)

