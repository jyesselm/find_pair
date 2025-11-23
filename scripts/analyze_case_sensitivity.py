#!/usr/bin/env python3
"""
Analyze case sensitivity differences between legacy and modern code.

This script investigates why lowercase bases (especially 'a') have much worse
match rates than uppercase bases.

Usage:
    python3 scripts/analyze_case_sensitivity.py
    python3 scripts/analyze_case_sensitivity.py --verbose
    python3 scripts/analyze_case_sensitivity.py --base a
"""

import json
import sys
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List

from x3dna_json_compare import JsonComparator, ComparisonResult


def analyze_case_sensitivity(project_root: Path, base_type: str = None,
                            verbose: bool = False) -> None:
    """Analyze case sensitivity differences."""
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"
    
    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)
    
    comparator = JsonComparator(enable_cache=False)
    
    print("=" * 80)
    print("CASE SENSITIVITY ANALYSIS")
    print("=" * 80)
    print(f"\nProcessing {len(common)} files...\n")
    
    case_stats = defaultdict(lambda: {
        'total': 0,
        'exact_matches': 0,
        'set_matches': 0,
        'real_diffs': 0,
        'missing': 0,
        'examples': []
    })
    
    import multiprocessing
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    max_workers = multiprocessing.cpu_count()
    
    def compare_pdb(pdb_id):
        """Compare a single PDB and extract case-related differences."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        
        if not all([legacy_file.exists(), modern_file.exists()]):
            return []
        
        try:
            result = comparator.compare_files(legacy_file, modern_file, None, pdb_id)
            
            pdb_differences = []
            
            if result.frame_comparison:
                fc = result.frame_comparison
                
                # Check missing residues
                for missing in fc.missing_residues:
                    missing_base = missing.base_type
                    if base_type and missing_base != base_type:
                        continue
                    if missing_base.islower():
                        pdb_differences.append({
                            'type': 'missing',
                            'pdb_id': pdb_id,
                            'residue': f"{missing.chain_id}:{missing.residue_seq}{missing.insertion}",
                            'residue_name': missing.residue_name,
                            'base_type': missing_base,
                        })
                
                # Check mismatched calculations
                for mismatch in fc.mismatched_calculations:
                    mismatch_base = mismatch.legacy_record.get('base_type', '?')
                    if base_type and mismatch_base != base_type:
                        continue
                    
                    leg_atoms = set(mismatch.legacy_matched_atoms)
                    mod_atoms = set(mismatch.modern_matched_atoms)
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    
                    if leg_atoms != mod_atoms:
                        only_legacy = sorted(leg_atoms - mod_atoms)
                        only_modern = sorted(mod_atoms - leg_atoms)
                        
                        pdb_differences.append({
                            'type': 'mismatch',
                            'pdb_id': pdb_id,
                            'residue': f"{chain_id}:{residue_seq}{insertion}",
                            'residue_name': mismatch.legacy_record.get('residue_name', ''),
                            'base_type': mismatch_base,
                            'only_legacy': only_legacy,
                            'only_modern': only_modern,
                            'legacy_atoms': sorted(leg_atoms),
                            'modern_atoms': sorted(mod_atoms),
                        })
            
            return pdb_differences
        except Exception as e:
            return []
    
    # Process in parallel
    completed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }
        
        for future in as_completed(future_to_pdb):
            pdb_diffs = future.result()
            
            for diff in pdb_diffs:
                base = diff['base_type']
                stats = case_stats[base]
                stats['total'] += 1
                
                if diff['type'] == 'missing':
                    stats['missing'] += 1
                else:
                    if diff['legacy_atoms'] == diff['modern_atoms']:
                        stats['exact_matches'] += 1
                        stats['set_matches'] += 1
                    elif set(diff['legacy_atoms']) == set(diff['modern_atoms']):
                        stats['set_matches'] += 1
                    else:
                        stats['real_diffs'] += 1
                    
                    if len(stats['examples']) < 10:
                        stats['examples'].append(diff)
            
            completed += 1
            if completed % 100 == 0:
                print(f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)",
                      flush=True)
    
    print(f"\n{'=' * 80}")
    print("CASE SENSITIVITY STATISTICS")
    print(f"{'=' * 80}\n")
    
    # Compare uppercase vs lowercase
    uppercase_bases = ['A', 'G', 'C', 'T', 'U']
    lowercase_bases = ['a', 'g', 'c', 't', 'u']
    
    print("UPPERCASE vs LOWERCASE COMPARISON")
    print("-" * 80)
    
    for upper, lower in zip(uppercase_bases, lowercase_bases):
        upper_stats = case_stats.get(upper, {})
        lower_stats = case_stats.get(lower, {})
        
        upper_total = upper_stats.get('total', 0)
        lower_total = lower_stats.get('total', 0)
        
        if upper_total == 0 and lower_total == 0:
            continue
        
        print(f"\n{upper} (uppercase):")
        if upper_total > 0:
            exact_rate = upper_stats.get('exact_matches', 0) / upper_total * 100
            diff_rate = upper_stats.get('real_diffs', 0) / upper_total * 100
            print(f"  Total: {upper_total}")
            print(f"  Exact matches: {upper_stats.get('exact_matches', 0)} ({exact_rate:.1f}%)")
            print(f"  Real differences: {upper_stats.get('real_diffs', 0)} ({diff_rate:.1f}%)")
        else:
            print(f"  No residues found")
        
        print(f"{lower} (lowercase):")
        if lower_total > 0:
            exact_rate = lower_stats.get('exact_matches', 0) / lower_total * 100
            diff_rate = lower_stats.get('real_diffs', 0) / lower_total * 100
            print(f"  Total: {lower_total}")
            print(f"  Exact matches: {lower_stats.get('exact_matches', 0)} ({exact_rate:.1f}%)")
            print(f"  Real differences: {lower_stats.get('real_diffs', 0)} ({diff_rate:.1f}%)")
            print(f"  Missing: {lower_stats.get('missing', 0)}")
        else:
            print(f"  No residues found")
    
    # Show examples for lowercase bases
    print(f"\n{'=' * 80}")
    print("LOWERCASE BASE EXAMPLES")
    print(f"{'=' * 80}\n")
    
    for base in lowercase_bases:
        stats = case_stats.get(base, {})
        if stats.get('total', 0) == 0:
            continue
        
        print(f"{base} (lowercase) - {stats['total']} total, "
              f"{stats['real_diffs']} real differences:")
        
        for i, example in enumerate(stats['examples'][:5], 1):
            print(f"  {i}. {example['pdb_id']} {example['residue']}: "
                  f"{example['residue_name']}")
            if example['type'] == 'mismatch':
                if example['only_legacy']:
                    print(f"     Only in legacy: {example['only_legacy']}")
                if example['only_modern']:
                    print(f"     Only in modern: {example['only_modern']}")
            print()
    
    # Check residue name patterns
    print(f"{'=' * 80}")
    print("RESIDUE NAME ANALYSIS")
    print(f"{'=' * 80}\n")
    
    residue_names = defaultdict(int)
    for base in lowercase_bases:
        stats = case_stats.get(base, {})
        for example in stats['examples']:
            residue_names[example['residue_name']] += 1
    
    print("Most common residue names in lowercase base differences:")
    for residue_name, count in sorted(residue_names.items(), key=lambda x: x[1], reverse=True)[:20]:
        print(f"  {residue_name}: {count} occurrences")
    
    # Save detailed report
    output_file = project_root / "docs" / "CASE_SENSITIVITY_ANALYSIS.json"
    with open(output_file, 'w') as f:
        json.dump({
            'case_statistics': dict(case_stats),
            'residue_name_counts': dict(residue_names),
        }, f, indent=2)
    
    print(f"\n💾 Detailed report saved to: {output_file}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Analyze case sensitivity differences"
    )
    parser.add_argument(
        "--base",
        help="Focus on specific base type (a, g, c, t, u)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show verbose output"
    )
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    analyze_case_sensitivity(project_root, args.base, args.verbose)


if __name__ == "__main__":
    main()

