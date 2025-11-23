#!/usr/bin/env python3
"""
Analyze differences for a specific base type.

This script extracts and analyzes residues of a specific base type that have
differences between legacy and modern JSON files.

Usage:
    python3 scripts/analyze_base_type_differences.py A
    python3 scripts/analyze_base_type_differences.py A --limit 50
    python3 scripts/analyze_base_type_differences.py a --verbose
"""

import json
import sys
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

from x3dna_json_compare import JsonComparator, ComparisonResult


def analyze_base_type(base_type: str, project_root: Path, limit: int = 100,
                     verbose: bool = False) -> None:
    """Analyze differences for a specific base type."""
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"
    
    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)
    
    comparator = JsonComparator(enable_cache=False)
    
    print("=" * 80)
    print(f"ANALYZING BASE TYPE: {base_type}")
    print("=" * 80)
    print(f"\nProcessing {len(common)} files...\n")
    
    differences = []
    patterns = defaultdict(list)
    
    completed = 0
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import multiprocessing
    
    max_workers = multiprocessing.cpu_count()
    
    def compare_pdb(pdb_id):
        """Compare a single PDB and extract base type differences."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        
        if not all([legacy_file.exists(), modern_file.exists()]):
            return []
        
        try:
            result = comparator.compare_files(legacy_file, modern_file, None, pdb_id)
            
            pdb_differences = []
            
            if result.frame_comparison:
                fc = result.frame_comparison
                
                for mismatch in fc.mismatched_calculations:
                    mismatch_base = mismatch.legacy_record.get('base_type', '?')
                    if mismatch_base != base_type:
                        continue
                    
                    leg_atoms = set(mismatch.legacy_matched_atoms)
                    mod_atoms = set(mismatch.modern_matched_atoms)
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    
                    if leg_atoms != mod_atoms:
                        only_legacy = sorted(leg_atoms - mod_atoms)
                        only_modern = sorted(mod_atoms - leg_atoms)
                        
                        diff_info = {
                            'pdb_id': pdb_id,
                            'residue': f"{chain_id}:{residue_seq}{insertion}",
                            'residue_name': mismatch.legacy_record.get('residue_name', ''),
                            'base_type': mismatch_base,
                            'legacy_atoms': sorted(leg_atoms),
                            'modern_atoms': sorted(mod_atoms),
                            'only_legacy': only_legacy,
                            'only_modern': only_modern,
                            'legacy_count': len(leg_atoms),
                            'modern_count': len(mod_atoms),
                            'rms_legacy': mismatch.legacy_record.get('rms_fit', 0),
                            'rms_modern': mismatch.modern_record.get('rms_fit', 0),
                            'rms_diff': abs(mismatch.legacy_record.get('rms_fit', 0) - 
                                          mismatch.modern_record.get('rms_fit', 0)),
                        }
                        
                        pdb_differences.append(diff_info)
                        
                        # Track pattern
                        pattern_key = (tuple(only_legacy), tuple(only_modern))
                        patterns[pattern_key].append(diff_info)
            
            return pdb_differences
        except Exception as e:
            return []
    
    # Process in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }
        
        for future in as_completed(future_to_pdb):
            pdb_diffs = future.result()
            differences.extend(pdb_diffs)
            completed += 1
            
            if completed % 100 == 0:
                print(f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)", 
                      flush=True)
    
    print(f"\n{'=' * 80}")
    print(f"SUMMARY FOR BASE TYPE: {base_type}")
    print(f"{'=' * 80}\n")
    
    print(f"Total differences found: {len(differences)}")
    print(f"Unique patterns: {len(patterns)}\n")
    
    # Show top patterns
    sorted_patterns = sorted(patterns.items(), key=lambda x: len(x[1]), reverse=True)
    
    print("=" * 80)
    print("TOP 10 DIFFERENCE PATTERNS")
    print("=" * 80)
    print()
    
    for i, (pattern, examples) in enumerate(sorted_patterns[:10], 1):
        only_legacy, only_modern = pattern
        print(f"Pattern {i}: {len(examples)} occurrences ({len(examples)*100/len(differences):.1f}%)")
        if only_legacy:
            print(f"  Only in legacy: {list(only_legacy)}")
        if only_modern:
            print(f"  Only in modern: {list(only_modern)}")
        print(f"  Example: {examples[0]['pdb_id']} {examples[0]['residue']} "
              f"({examples[0]['residue_name']})")
        if verbose and examples[0]['rms_diff'] > 0.001:
            print(f"  RMS diff: {examples[0]['rms_diff']:.6f}")
        print()
    
    # Show examples
    print("=" * 80)
    print(f"EXAMPLE DIFFERENCES (First {min(limit, len(differences))})")
    print("=" * 80)
    print()
    
    for i, diff in enumerate(differences[:limit], 1):
        print(f"{i}. {diff['pdb_id']} {diff['residue']}: {diff['residue_name']} ({diff['base_type']})")
        print(f"   Legacy: {diff['legacy_count']} atoms - {diff['legacy_atoms']}")
        print(f"   Modern: {diff['modern_count']} atoms - {diff['modern_atoms']}")
        if diff['only_legacy']:
            print(f"   Only in legacy: {diff['only_legacy']}")
        if diff['only_modern']:
            print(f"   Only in modern: {diff['only_modern']}")
        if diff['rms_diff'] > 0.001:
            print(f"   RMS: Legacy={diff['rms_legacy']:.6f}, Modern={diff['rms_modern']:.6f}, "
                  f"Diff={diff['rms_diff']:.6f}")
        print()
    
    # Save detailed report
    output_file = project_root / "docs" / f"BASE_TYPE_{base_type}_ANALYSIS.json"
    with open(output_file, 'w') as f:
        json.dump({
            'base_type': base_type,
            'total_differences': len(differences),
            'unique_patterns': len(patterns),
            'patterns': {
                f"Pattern_{i}": {
                    'count': len(examples),
                    'only_legacy': list(pattern[0]),
                    'only_modern': list(pattern[1]),
                    'example': examples[0] if examples else None
                }
                for i, (pattern, examples) in enumerate(sorted_patterns[:20], 1)
            },
            'examples': differences[:limit],
        }, f, indent=2)
    
    print(f"\n💾 Detailed report saved to: {output_file}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Analyze differences for a specific base type"
    )
    parser.add_argument(
        "base_type",
        help="Base type to analyze (A, G, C, T, U, a, g, c, t, u, etc.)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=100,
        help="Limit number of examples shown (default: 100)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show verbose output including RMS differences"
    )
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    analyze_base_type(args.base_type, project_root, args.limit, args.verbose)


if __name__ == "__main__":
    main()

