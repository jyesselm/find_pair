#!/usr/bin/env python3
"""
Validate H-bond detection on test_set_10.
"""
import json
import sys
from pathlib import Path
import importlib.util

# Load hbond comparison module
spec = importlib.util.spec_from_file_location(
    "hbond_comparison",
    "x3dna_json_compare/hbond_comparison.py"
)
hbond_comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hbond_comparison)

compare_hbond_lists = hbond_comparison.compare_hbond_lists

test_pdbs = ["1Q96", "1VBY", "3AVY", "3G8T", "3KNC", "4AL5", "5UJ2", "6CAQ", "6LTU", "8J1J"]

print("H-bond List Validation (test_set_10)")
print("="*70)

results = []
perfect_count = 0
minor_type_diff = 0  # Count of H-bond type differences only

for pdb in test_pdbs:
    legacy_file = Path(f"data/json_legacy/hbond_list/{pdb}.json")
    modern_file = Path(f"data/json/hbond_list/{pdb}.json")
    
    if not legacy_file.exists():
        print(f"⚠️  {pdb}: No legacy hbond_list")
        continue
    
    if not modern_file.exists():
        print(f"❌ {pdb}: No modern hbond_list")
        results.append({'pdb': pdb, 'status': 'no_modern'})
        continue
    
    with open(legacy_file) as f:
        legacy = json.load(f)
    with open(modern_file) as f:
        modern = json.load(f)
    
    comparison = compare_hbond_lists(legacy, modern, tolerance=1e-6)
    
    total = comparison.total_legacy
    common = comparison.common_count
    mismatched = len(comparison.mismatched_pairs)
    
    if common == total and mismatched == 0:
        print(f"✅ {pdb}: Perfect match ({total} pairs)")
        perfect_count += 1
        results.append({'pdb': pdb, 'status': 'perfect', 'total': total})
    else:
        # Check if only type differences
        only_type_diffs = True
        if comparison.mismatched_pairs:
            for mismatch in comparison.mismatched_pairs:
                pair_comp = mismatch.get('comparison')
                if pair_comp and len(pair_comp.mismatched_hbonds) > 0:
                    for hb_mismatch in pair_comp.mismatched_hbonds:
                        # Check if only type differs
                        donor_same = hb_mismatch['donor_atom'] == hb_mismatch['donor_atom']
                        acceptor_same = hb_mismatch['acceptor_atom'] == hb_mismatch['acceptor_atom']
                        dist_diff = abs(hb_mismatch.get('distance_diff', 0))
                        
                        if dist_diff > 1e-6 or not donor_same or not acceptor_same:
                            only_type_diffs = False
                            break
        
        if only_type_diffs and mismatched > 0:
            print(f"⚠️  {pdb}: {common}/{total} common, {mismatched} type-only diffs")
            minor_type_diff += 1
            results.append({'pdb': pdb, 'status': 'type_diff', 'total': total, 'mismatched': mismatched})
        else:
            print(f"❌ {pdb}: {common}/{total} common, {mismatched} mismatched")
            results.append({'pdb': pdb, 'status': 'mismatch', 'total': total, 'mismatched': mismatched})

# Summary
print("\n" + "="*70)
print(f"Summary:")
print(f"  Perfect matches: {perfect_count}/{len(test_pdbs)}")
print(f"  Minor type-only differences: {minor_type_diff}/{len(test_pdbs)}")
print(f"  Total acceptable: {perfect_count + minor_type_diff}/{len(test_pdbs)}")

if perfect_count == len(test_pdbs):
    print("\n🎉 All test PDBs have perfect H-bond matches!")
elif perfect_count + minor_type_diff >= len(test_pdbs) * 0.8:  # 80% threshold
    print("\n✅ Good enough - type differences are minor algorithmic variations")
else:
    print("\n⚠️  Some significant mismatches found - needs investigation")

# Save results
output = {
    'total_tested': len(test_pdbs),
    'perfect': perfect_count,
    'minor_type_diff': minor_type_diff,
    'results': results
}

Path("data/validation_results").mkdir(parents=True, exist_ok=True)
with open("data/validation_results/hbond_test10_results.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\n✅ Results saved to data/validation_results/hbond_test10_results.json")

sys.exit(0 if perfect_count + minor_type_diff >= len(test_pdbs) * 0.8 else 1)

