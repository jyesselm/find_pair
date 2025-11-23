#!/usr/bin/env python3
"""Test if legacy mode matches legacy JSON output"""

import json
from pathlib import Path
from x3dna_json_compare import JsonComparator

def main():
    project_root = Path(__file__).parent.parent
    legacy_file = project_root / "data/json_legacy/1H4S.json"
    modern_file = project_root / "data/json/1H4S_legacy.json"
    
    if not modern_file.exists():
        print(f"ERROR: Modern JSON file not found: {modern_file}")
        print("Please run: ./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S_legacy.json --legacy")
        return
    
    comparator = JsonComparator(enable_cache=False)
    result = comparator.compare_files(legacy_file, modern_file, None, '1H4S')
    
    if result.frame_comparison:
        fc = result.frame_comparison
        print(f"Legacy: {fc.total_legacy} residues")
        print(f"Modern: {fc.total_modern} residues")
        print(f"Missing: {len(fc.missing_residues)}")
        print(f"Mismatched: {len(fc.mismatched_calculations)}")
        print()
        
        # Check exact matches
        exact_matches = 0
        real_diffs = 0
        order_diffs = 0
        
        for mismatch in fc.mismatched_calculations:
            leg_atoms = set(mismatch.legacy_matched_atoms)
            mod_atoms = set(mismatch.modern_matched_atoms)
            
            if leg_atoms == mod_atoms:
                exact_matches += 1
                if list(mismatch.legacy_matched_atoms) != list(mismatch.modern_matched_atoms):
                    order_diffs += 1
            else:
                real_diffs += 1
                if real_diffs <= 5:
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    only_leg = sorted(leg_atoms - mod_atoms)
                    only_mod = sorted(mod_atoms - leg_atoms)
                    print(f"Diff {real_diffs}: {chain_id}:{residue_seq}{insertion}")
                    if only_leg:
                        print(f"  Only legacy: {only_leg}")
                    if only_mod:
                        print(f"  Only modern: {only_mod}")
        
        print()
        print(f"Exact atom set matches: {exact_matches}/{len(fc.mismatched_calculations)}")
        print(f"Atom set differences: {real_diffs}")
        print(f"Order-only differences: {order_diffs}")
        
        if real_diffs == 0:
            if order_diffs == 0:
                print("✅ PERFECT MATCH!")
            else:
                print(f"✅ Atom sets match (order differs in {order_diffs} cases)")
        else:
            print(f"❌ Still {real_diffs} atom set differences")
    
    # Also show example comparisons
    print("\n" + "=" * 80)
    print("EXAMPLE COMPARISONS")
    print("=" * 80)
    
    with open(legacy_file) as f:
        legacy_data = json.load(f)
    with open(modern_file) as f:
        modern_data = json.load(f)
    
    # Find examples
    examples = []
    for leg_calc in legacy_data.get('calculations', []):
        if leg_calc.get('type') == 'base_frame_calc':
            base = leg_calc.get('base_type', '')
            if base in ['G', 'U']:
                leg_atoms = leg_calc.get('matched_atoms', [])
                
                # Find matching modern
                for mod_calc in modern_data.get('calculations', []):
                    if (mod_calc.get('type') == 'base_frame_calc' and
                        mod_calc.get('base_type') == base and
                        mod_calc.get('chain_id') == leg_calc.get('chain_id') and
                        mod_calc.get('residue_seq') == leg_calc.get('residue_seq')):
                        mod_atoms = mod_calc.get('matched_atoms', [])
                        examples.append((base, leg_calc.get('chain_id'), 
                                       leg_calc.get('residue_seq'),
                                       leg_atoms, mod_atoms))
                        break
                
                if len(examples) >= 2:
                    break
    
    for base, chain, seq, leg_atoms, mod_atoms in examples:
        print(f"\n{base} {chain}:{seq}:")
        print(f"  Legacy: {leg_atoms}")
        print(f"  Modern: {mod_atoms}")
        print(f"  Match: {leg_atoms == mod_atoms}")

if __name__ == "__main__":
    main()

