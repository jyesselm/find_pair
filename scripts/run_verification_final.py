#!/usr/bin/env python3
"""Final verification script"""

import subprocess
import sys
import json
from pathlib import Path

project_root = Path(__file__).parent.parent

print("=" * 80)
print("LEGACY MODE VERIFICATION")
print("=" * 80)
print()

# Step 1: Rebuild
print("Step 1: Rebuilding code...")
result = subprocess.run(["make", "-j8"], cwd=project_root / "build", 
                       capture_output=True, text=True)
if result.returncode != 0:
    print("❌ Build failed!")
    print(result.stderr)
    sys.exit(1)
print("✅ Build successful")
print()

# Step 2: Generate JSON
print("Step 2: Generating JSON with legacy mode...")
result = subprocess.run([
    str(project_root / "build" / "generate_modern_json"),
    str(project_root / "data" / "pdb" / "1H4S.pdb"),
    str(project_root / "data" / "json" / "1H4S_legacy.json"),
    "--legacy"
], capture_output=True, text=True)
if result.returncode != 0:
    print("❌ Generation failed!")
    print(result.stderr)
    sys.exit(1)
print(result.stdout.split('\n')[-10:])
print("✅ JSON generated")
print()

# Step 3: Load and compare
print("Step 3: Comparing JSON files...")
legacy_file = project_root / "data" / "json_legacy" / "1H4S.json"
modern_file = project_root / "data" / "json" / "1H4S_legacy.json"

with open(legacy_file) as f:
    legacy = json.load(f)
with open(modern_file) as f:
    modern = json.load(f)

leg_calcs = {c.get('residue_seq', 0): c for c in legacy['calculations'] 
             if c.get('type') == 'base_frame_calc'}
mod_calcs = {(c.get('chain_id', ''), c.get('residue_seq', 0)): c 
             for c in modern['calculations'] 
             if c.get('type') == 'base_frame_calc'}

print(f"Legacy: {len(leg_calcs)} residues")
print(f"Modern: {len(mod_calcs)} residues")
print()

exact = 0
set_match = 0
real_diff = 0

for seq, leg_calc in sorted(leg_calcs.items()):
    chain = leg_calc.get('chain_id', '')
    key = (chain, seq)
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

total = len(leg_calcs)
print(f"Exact matches: {exact}/{total} ({exact*100/total:.1f}%)")
print(f"Set matches: {set_match}/{total} ({set_match*100/total:.1f}%)")
print(f"Real differences: {real_diff}/{total} ({real_diff*100/total:.1f}%)")
print()

# Spot checks
print("=" * 80)
print("SPOT CHECKS")
print("=" * 80)
print()

# Purine G T:4
leg_g = [c for c in legacy['calculations'] if c.get('type')=='base_frame_calc' and c.get('base_type')=='G' and c.get('chain_id')=='T' and c.get('residue_seq')==4]
mod_g = [c for c in modern['calculations'] if c.get('type')=='base_frame_calc' and c.get('base_type')=='G' and c.get('chain_id')=='T' and c.get('residue_seq')==4]

if leg_g and mod_g:
    print(f"Purine G T:4:")
    print(f"  Legacy: {leg_g[0]['matched_atoms']}")
    print(f"  Modern: {mod_g[0]['matched_atoms']}")
    print(f"  Match: {leg_g[0]['matched_atoms'] == mod_g[0]['matched_atoms']}")
    print(f"  Legacy num: {leg_g[0].get('num_matched_atoms')}, Modern num: {mod_g[0].get('num_matched_atoms')}")
    print()

# Pyrimidine U
leg_u = [c for c in legacy['calculations'] if c.get('type')=='base_frame_calc' and c.get('base_type')=='U'][:1]
if leg_u:
    leg_u = leg_u[0]
    mod_u = [c for c in modern['calculations'] if c.get('type')=='base_frame_calc' and c.get('base_type')=='U' and c.get('chain_id')==leg_u.get('chain_id') and c.get('residue_seq')==leg_u.get('residue_seq')]
    if mod_u:
        mod_u = mod_u[0]
        print(f"Pyrimidine U {leg_u.get('chain_id')}:{leg_u.get('residue_seq')}:")
        print(f"  Legacy: {leg_u['matched_atoms']}")
        print(f"  Modern: {mod_u['matched_atoms']}")
        print(f"  Match: {leg_u['matched_atoms'] == mod_u['matched_atoms']}")
        print(f"  Legacy num: {leg_u.get('num_matched_atoms')}, Modern num: {mod_u.get('num_matched_atoms')}")
        print()

# Final verdict
print("=" * 80)
print("VERDICT")
print("=" * 80)
if real_diff == 0 and exact == total:
    print("✅ PERFECT MATCH! All residues match exactly!")
elif real_diff == 0:
    print(f"✅ ATOM SET MATCH! ({exact} exact, {total-exact} order differences)")
else:
    print(f"❌ Still {real_diff} real differences")
    sys.exit(1)

