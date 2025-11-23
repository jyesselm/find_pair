#!/usr/bin/env python3
"""Rerun verification - uses subprocess with clean environment"""

import subprocess
import json
import sys
import os
from pathlib import Path

project_root = Path(__file__).parent.parent
os.chdir(project_root)

print("=" * 80)
print("LEGACY MODE VERIFICATION - REBUILD AND RERUN")
print("=" * 80)
print()

# Clean environment
env = os.environ.copy()
env.pop('PS1', None)
env.pop('PROMPT', None)

# Step 1: Rebuild
print("Step 1: Rebuilding code...")
try:
    result = subprocess.run(
        ["make", "-j8"],
        cwd=project_root / "build",
        env=env,
        capture_output=True,
        text=True,
        timeout=300
    )
    if result.returncode != 0:
        print("❌ Build failed!")
        print(result.stderr[-500:])
        sys.exit(1)
    print("✅ Build successful")
    print(result.stdout.split('\n')[-5:])
except Exception as e:
    print(f"❌ Build error: {e}")
    sys.exit(1)
print()

# Step 2: Generate JSON
print("Step 2: Generating JSON with --legacy flag...")
modern_json = project_root / "data" / "json" / "1H4S_legacy.json"
try:
    result = subprocess.run(
        [str(project_root / "build" / "generate_modern_json"),
         str(project_root / "data" / "pdb" / "1H4S.pdb"),
         str(modern_json),
         "--legacy"],
        env=env,
        capture_output=True,
        text=True,
        timeout=120
    )
    if result.returncode != 0:
        print("❌ Generation failed!")
        print(result.stderr[-500:])
        sys.exit(1)
    print("✅ JSON generated")
    lines = [l for l in result.stdout.split('\n') if l.strip()]
    print('\n'.join(lines[-5:]))
except Exception as e:
    print(f"❌ Generation error: {e}")
    sys.exit(1)
print()

# Step 3: Verify
print("Step 3: Verifying JSON match...")
legacy_file = project_root / "data" / "json_legacy" / "1H4S.json"

if not modern_json.exists():
    print(f"❌ Modern JSON file not found: {modern_json}")
    sys.exit(1)

with open(legacy_file) as f:
    legacy = json.load(f)
with open(modern_json) as f:
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
first_diff = None

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
        if first_diff is None:
            first_diff = (leg_calc, mod_calc, leg_atoms, mod_atoms)

total = len([c for c in leg_calcs if f\"{c.get('chain_id', '')}:{c.get('residue_seq', 0)}\" in mod_calcs])

print(f"Exact matches: {exact}/{total} ({exact*100/total:.1f}%)")
print(f"Set matches: {set_match}/{total} ({set_match*100/total:.1f}%)")
print(f"Real differences: {real_diff}/{total} ({real_diff*100/total:.1f}%)")
print()

if first_diff:
    leg_calc, mod_calc, leg_atoms, mod_atoms = first_diff
    print(f"First difference: {leg_calc.get('base_type', '?')} {leg_calc.get('chain_id', '?')}:{leg_calc.get('residue_seq', 0)}")
    leg_only = sorted(set(leg_atoms) - set(mod_atoms))
    mod_only = sorted(set(mod_atoms) - set(leg_atoms))
    print(f"  Legacy only: {leg_only}")
    print(f"  Modern only: {mod_only}")
    print()

# Spot check: Purine G T:4
leg_g = [c for c in leg_calcs if c.get('base_type')=='G' and c.get('chain_id')=='T' and c.get('residue_seq')==4]
if leg_g:
    leg_g = leg_g[0]
    key = f\"{leg_g.get('chain_id', '')}:{leg_g.get('residue_seq', 0)}\"
    if key in mod_calcs:
        mod_g = mod_calcs[key]
        print(f"Spot check - Purine G T:4:")
        print(f"  Legacy: {leg_g['matched_atoms']}")
        print(f"  Modern: {mod_g['matched_atoms']}")
        print(f"  Match: {leg_g['matched_atoms'] == mod_g['matched_atoms']}")
        print()

print("=" * 80)
if real_diff == 0 and exact == total:
    print("✅ PERFECT MATCH! All residues match exactly!")
    sys.exit(0)
elif real_diff == 0:
    print(f"✅ ATOM SET MATCH! ({exact} exact, {total-exact} order differences)")
    sys.exit(0)
else:
    print(f"❌ Still {real_diff} real differences")
    sys.exit(1)

