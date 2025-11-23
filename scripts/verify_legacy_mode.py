#!/usr/bin/env python3
"""Comprehensive verification that legacy mode matches legacy JSON output exactly"""

import json
import sys
from pathlib import Path
from collections import defaultdict

def main():
    project_root = Path(__file__).parent.parent
    legacy_file = project_root / "data/json_legacy/1H4S.json"
    modern_file = project_root / "data/json/1H4S_legacy.json"
    
    print("=" * 80)
    print("COMPREHENSIVE LEGACY MODE VERIFICATION")
    print("=" * 80)
    print()
    
    # Check files exist
    if not legacy_file.exists():
        print(f"ERROR: Legacy file not found: {legacy_file}")
        sys.exit(1)
    
    if not modern_file.exists():
        print(f"ERROR: Modern file not found: {modern_file}")
        print("Please run: ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
        sys.exit(1)
    
    # Load JSON files
    print("Loading JSON files...")
    with open(legacy_file) as f:
        legacy_data = json.load(f)
    with open(modern_file) as f:
        modern_data = json.load(f)
    
    # Extract base_frame_calc entries
    leg_calcs = {c.get('residue_seq', 0): c for c in legacy_data.get('calculations', []) 
                 if c.get('type') == 'base_frame_calc'}
    mod_calcs = {(c.get('chain_id', ''), c.get('residue_seq', 0)): c 
                 for c in modern_data.get('calculations', []) 
                 if c.get('type') == 'base_frame_calc'}
    
    print(f"Legacy base_frame_calc entries: {len(leg_calcs)}")
    print(f"Modern base_frame_calc entries: {len(mod_calcs)}")
    print()
    
    # Statistics
    total = 0
    exact_matches = 0
    atom_set_matches = 0
    atom_order_matches = 0
    real_diffs = 0
    
    # Pattern analysis
    patterns = defaultdict(int)
    
    # Detailed comparison
    print("=" * 80)
    print("DETAILED COMPARISON")
    print("=" * 80)
    print()
    
    examples = []
    
    for seq, leg_calc in sorted(leg_calcs.items()):
        chain = leg_calc.get('chain_id', '')
        base = leg_calc.get('base_type', '')
        key = (chain, seq)
        
        if key not in mod_calcs:
            print(f"WARNING: Modern missing {base} {chain}:{seq}")
            continue
        
        mod_calc = mod_calcs[key]
        total += 1
        
        leg_atoms = leg_calc.get('matched_atoms', [])
        mod_atoms = mod_calc.get('matched_atoms', [])
        
        # Exact match (same list)
        if leg_atoms == mod_atoms:
            exact_matches += 1
            atom_set_matches += 1
            atom_order_matches += 1
        else:
            # Check if same set (ignoring order)
            leg_set = set(leg_atoms)
            mod_set = set(mod_atoms)
            
            if leg_set == mod_set:
                atom_set_matches += 1
                # Order is different
                only_leg = leg_set - mod_set
                only_mod = mod_set - leg_set
                
                # Store example
                if len(examples) < 5:
                    examples.append({
                        'residue': f"{base} {chain}:{seq}",
                        'legacy': leg_atoms,
                        'modern': mod_atoms,
                        'issue': 'order_only'
                    })
            else:
                # Real differences
                real_diffs += 1
                only_leg = leg_set - mod_set
                only_mod = mod_set - leg_set
                
                # Track pattern
                pattern_key = f"Legacy_only: {sorted(only_leg)}, Modern_only: {sorted(only_mod)}"
                patterns[pattern_key] += 1
                
                # Store example
                if len(examples) < 10:
                    examples.append({
                        'residue': f"{base} {chain}:{seq}",
                        'legacy': leg_atoms,
                        'modern': mod_atoms,
                        'legacy_only': sorted(only_leg),
                        'modern_only': sorted(only_mod),
                        'issue': 'set_difference'
                    })
    
    # Print statistics
    print("=" * 80)
    print("STATISTICS")
    print("=" * 80)
    print(f"Total residues compared: {total}")
    print(f"Exact matches (same list): {exact_matches}/{total} ({exact_matches*100/total:.1f}%)")
    print(f"Atom set matches (same atoms, order may differ): {atom_set_matches}/{total} ({atom_set_matches*100/total:.1f}%)")
    print(f"Atom order matches: {atom_order_matches}/{total} ({atom_order_matches*100/total:.1f}%)")
    print(f"Real differences (different atom sets): {real_diffs}/{total} ({real_diffs*100/total:.1f}%)")
    print()
    
    # Print patterns
    if patterns:
        print("=" * 80)
        print("DIFFERENCE PATTERNS")
        print("=" * 80)
        for pattern, count in sorted(patterns.items(), key=lambda x: x[1], reverse=True):
            print(f"{count}x: {pattern}")
        print()
    
    # Print examples
    if examples:
        print("=" * 80)
        print("EXAMPLE DIFFERENCES")
        print("=" * 80)
        for i, ex in enumerate(examples[:10], 1):
            print(f"\n{i}. {ex['residue']}:")
            print(f"   Legacy: {ex['legacy']}")
            print(f"   Modern: {ex['modern']}")
            if ex['issue'] == 'set_difference':
                if ex.get('legacy_only'):
                    print(f"   Only in legacy: {ex['legacy_only']}")
                if ex.get('modern_only'):
                    print(f"   Only in modern: {ex['modern_only']}")
            elif ex['issue'] == 'order_only':
                print(f"   (Same atoms, different order)")
        print()
    
    # Check specific base types
    print("=" * 80)
    print("BASE TYPE ANALYSIS")
    print("=" * 80)
    print()
    
    for base_type in ['A', 'G', 'C', 'T', 'U']:
        base_leg = {seq: c for seq, c in leg_calcs.items() if c.get('base_type') == base_type}
        base_mod = {k: c for k, c in mod_calcs.items() if c.get('base_type') == base_type}
        
        if not base_leg:
            continue
        
        print(f"\n{base_type} ({len(base_leg)} residues):")
        
        base_exact = 0
        base_set_match = 0
        base_real_diff = 0
        
        for seq, leg_calc in base_leg.items():
            chain = leg_calc.get('chain_id', '')
            key = (chain, seq)
            
            if key not in base_mod:
                continue
            
            mod_calc = base_mod[key]
            leg_atoms = leg_calc.get('matched_atoms', [])
            mod_atoms = mod_calc.get('matched_atoms', [])
            
            if leg_atoms == mod_atoms:
                base_exact += 1
                base_set_match += 1
            elif set(leg_atoms) == set(mod_atoms):
                base_set_match += 1
            else:
                base_real_diff += 1
        
        print(f"  Exact matches: {base_exact}/{len(base_leg)} ({base_exact*100/len(base_leg):.1f}%)")
        print(f"  Set matches: {base_set_match}/{len(base_leg)} ({base_set_match*100/len(base_leg):.1f}%)")
        print(f"  Real differences: {base_real_diff}/{len(base_leg)} ({base_real_diff*100/len(base_leg):.1f}%)")
    
    print()
    
    # Final verdict
    print("=" * 80)
    print("VERDICT")
    print("=" * 80)
    
    if real_diffs == 0:
        if atom_order_matches == total:
            print("✅ PERFECT MATCH! All residues match exactly (same atoms, same order)")
        else:
            print(f"✅ ATOM SET MATCH! All residues have same atom sets (but {total - atom_order_matches} have different order)")
        print()
        print("Legacy mode is working correctly!")
    else:
        print(f"❌ MISMATCH: {real_diffs} residues have different atom sets")
        print()
        print("Legacy mode needs further investigation.")
        sys.exit(1)

if __name__ == "__main__":
    main()

