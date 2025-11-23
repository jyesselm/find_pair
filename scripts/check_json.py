#!/usr/bin/env python3
"""Check JSON files to see if legacy mode was applied"""

import json
from pathlib import Path

project_root = Path(__file__).parent.parent

legacy_file = project_root / "data" / "json_legacy" / "1H4S.json"
modern_file = project_root / "data" / "json" / "1H4S_legacy.json"

print("=" * 80)
print("CHECKING JSON FILES")
print("=" * 80)
print()

if not modern_file.exists():
    print(f"❌ Modern JSON file not found: {modern_file}")
    print("Need to generate with: ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
    exit(1)

with open(legacy_file) as f:
    legacy = json.load(f)
with open(modern_file) as f:
    modern = json.load(f)

leg_calcs = [c for c in legacy['calculations'] if c.get('type') == 'base_frame_calc']
mod_calcs = [c for c in modern['calculations'] if c.get('type') == 'base_frame_calc']

print(f"Legacy: {len(leg_calcs)} base_frame_calc records")
print(f"Modern: {len(mod_calcs)} base_frame_calc records")
print()

# Check for G T:4 (purine)
leg_g = [c for c in leg_calcs if c.get('base_type') == 'G' and c.get('chain_id') == 'T' and c.get('residue_seq') == 4]
mod_g = [c for c in mod_calcs if c.get('base_type') == 'G' and c.get('chain_id') == 'T' and c.get('residue_seq') == 4]

print("=" * 80)
print("PURINE G T:4 CHECK")
print("=" * 80)

if leg_g and mod_g:
    leg = leg_g[0]
    mod = mod_g[0]
    
    leg_atoms = leg.get('matched_atoms', [])
    mod_atoms = mod.get('matched_atoms', [])
    
    print(f"Legacy matched_atoms: {leg_atoms}")
    print(f"Modern matched_atoms: {mod_atoms}")
    print()
    
    # Check for C1' (should NOT be in legacy mode)
    has_c1_legacy = any("C1'" in a for a in leg_atoms)
    has_c1_modern = any("C1'" in a for a in mod_atoms)
    
    # Check for C4 (should NOT be in legacy mode)
    has_c4_legacy = any(" C4 " == a for a in leg_atoms)
    has_c4_modern = any(" C4 " == a for a in mod_atoms)
    
    # Check for H (should be in purines in legacy mode)
    has_h_legacy = any(" H" in a or a == " H" for a in leg_atoms)
    has_h_modern = any(" H" in a or a == " H" for a in mod_atoms)
    
    print(f"Legacy - C1': {has_c1_legacy}, C4: {has_c4_legacy}, H: {has_h_legacy}")
    print(f"Modern - C1': {has_c1_modern}, C4: {has_c4_modern}, H: {has_h_modern}")
    print()
    
    exact_match = leg_atoms == mod_atoms
    set_match = set(leg_atoms) == set(mod_atoms)
    
    if exact_match:
        print("✅ EXACT MATCH!")
    elif set_match:
        print("⚠️  SET MATCH (order differs)")
    else:
        print("❌ DIFFERENT SETS")
        leg_only = sorted(set(leg_atoms) - set(mod_atoms))
        mod_only = sorted(set(mod_atoms) - set(leg_atoms))
        print(f"  Legacy only: {leg_only}")
        print(f"  Modern only: {mod_only}")
        
        # Check if modern needs regeneration
        if has_c1_modern or has_c4_modern or not has_h_modern:
            print()
            print("❌ MODERN FILE DOES NOT HAVE LEGACY MODE!")
            print("Need to regenerate with --legacy flag:")
            print("  ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
else:
    print(f"Legacy G T:4 found: {len(leg_g) > 0}")
    print(f"Modern G T:4 found: {len(mod_g) > 0}")

# Check pyrimidine
print()
print("=" * 80)
print("PYRIMIDINE CHECK (first U)")
print("=" * 80)

leg_u = [c for c in leg_calcs if c.get('base_type') == 'U']
mod_u = [c for c in mod_calcs if c.get('base_type') == 'U']

if leg_u and mod_u:
    # Match by chain and seq
    leg = leg_u[0]
    chain = leg.get('chain_id')
    seq = leg.get('residue_seq')
    mod = [c for c in mod_u if c.get('chain_id') == chain and c.get('residue_seq') == seq]
    
    if mod:
        mod = mod[0]
        leg_atoms = leg.get('matched_atoms', [])
        mod_atoms = mod.get('matched_atoms', [])
        
        print(f"U {chain}:{seq}")
        print(f"Legacy: {leg_atoms}")
        print(f"Modern: {mod_atoms}")
        print()
        
        has_c1_modern = any("C1'" in a for a in mod_atoms)
        has_c4_modern = any(" C4 " == a for a in mod_atoms)
        has_n7_modern = any(" N7 " == a for a in mod_atoms)
        
        has_c1_legacy = any("C1'" in a for a in leg_atoms)
        has_c4_legacy = any(" C4 " == a for a in leg_atoms)
        has_n7_legacy = any(" N7 " == a for a in leg_atoms)
        
        print(f"Legacy - C1': {has_c1_legacy}, C4: {has_c4_legacy}, N7: {has_n7_legacy}")
        print(f"Modern - C1': {has_c1_modern}, C4: {has_c4_modern}, N7: {has_n7_modern}")
        print()
        
        if leg_atoms == mod_atoms:
            print("✅ EXACT MATCH!")
        elif set(leg_atoms) == set(mod_atoms):
            print("⚠️  SET MATCH (order differs)")
        else:
            print("❌ DIFFERENT SETS")
            if has_c1_modern or has_c4_modern or not has_n7_modern:
                print()
                print("❌ MODERN FILE DOES NOT HAVE LEGACY MODE!")
                print("Need to regenerate with --legacy flag")

print()
print("=" * 80)
print("VERDICT")
print("=" * 80)

# Quick overall check
needs_regeneration = False
for mod_calc in mod_calcs[:20]:  # Check first 20
    atoms = mod_calc.get('matched_atoms', [])
    if any("C1'" in a for a in atoms):
        needs_regeneration = True
        print(f"Found C1' in {mod_calc.get('base_type')} {mod_calc.get('chain_id')}:{mod_calc.get('residue_seq')}")
        break
    if any(" C4 " == a for a in atoms):
        needs_regeneration = True
        print(f"Found C4 in {mod_calc.get('base_type')} {mod_calc.get('chain_id')}:{mod_calc.get('residue_seq')}")
        break

if needs_regeneration:
    print()
    print("❌ JSON FILE NEEDS REGENERATION WITH --legacy FLAG")
    print("Run: ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
else:
    print("✅ JSON file appears to be in legacy mode (no C1', no C4)")

