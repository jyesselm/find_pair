#!/usr/bin/env python3
"""
Analyze residues with high RMS differences.

This script extracts and analyzes residues with RMS differences above a threshold
to understand if these are correlated with atom set differences.

Usage:
    python3 scripts/analyze_high_rms_differences.py
    python3 scripts/analyze_high_rms_differences.py --threshold 1.0
    python3 scripts/analyze_high_rms_differences.py --threshold 0.1 --limit 50
"""

import json
import sys
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List

from x3dna_json_compare import JsonComparator, ComparisonResult


def analyze_high_rms(project_root: Path, threshold: float = 0.1,
                     limit: int = 100) -> None:
    """Analyze residues with high RMS differences."""
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"
    
    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)
    
    comparator = JsonComparator(enable_cache=False)
    
    print("=" * 80)
    print(f"ANALYZING HIGH RMS DIFFERENCES (threshold: {threshold} Å)")
    print("=" * 80)
    print(f"\nProcessing {len(common)} files...\n")
    
    high_rms_residues = []
    base_type_stats = defaultdict(lambda: {'count': 0, 'total_rms': 0.0, 'max_rms': 0.0})
    
    import multiprocessing
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    max_workers = multiprocessing.cpu_count()
    
    def compare_pdb(pdb_id):
        """Compare a single PDB and extract high RMS differences."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        
        if not all([legacy_file.exists(), modern_file.exists()]):
            return []
        
        try:
            result = comparator.compare_files(legacy_file, modern_file, None, pdb_id)
            
            pdb_high_rms = []
            
            if result.frame_comparison:
                fc = result.frame_comparison
                
                for mismatch in fc.mismatched_calculations:
                    leg_rms = mismatch.legacy_record.get('rms_fit', 0.0)
                    mod_rms = mismatch.modern_record.get('rms_fit', 0.0)
                    rms_diff = abs(leg_rms - mod_rms)
                    
                    if rms_diff >= threshold:
                        chain_id, residue_seq, insertion = mismatch.residue_key
                        base_type = mismatch.legacy_record.get('base_type', '?')
                        leg_atoms = set(mismatch.legacy_matched_atoms)
                        mod_atoms = set(mismatch.modern_matched_atoms)
                        
                        has_atom_diff = leg_atoms != mod_atoms
                        
                        rms_info = {
                            'pdb_id': pdb_id,
                            'residue': f"{chain_id}:{residue_seq}{insertion}",
                            'residue_name': mismatch.legacy_record.get('residue_name', ''),
                            'base_type': base_type,
                            'rms_legacy': leg_rms,
                            'rms_modern': mod_rms,
                            'rms_diff': rms_diff,
                            'has_atom_difference': has_atom_diff,
                            'legacy_atoms': sorted(leg_atoms),
                            'modern_atoms': sorted(mod_atoms),
                            'only_legacy': sorted(leg_atoms - mod_atoms) if has_atom_diff else [],
                            'only_modern': sorted(mod_atoms - leg_atoms) if has_atom_diff else [],
                        }
                        
                        pdb_high_rms.append(rms_info)
            
            return pdb_high_rms
        except Exception as e:
            return []
    
    # Process in parallel
    completed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }
        
        for future in as_completed(future_to_pdb):
            pdb_high_rms = future.result()
            high_rms_residues.extend(pdb_high_rms)
            
            for rms_info in pdb_high_rms:
                base_type = rms_info['base_type']
                stats = base_type_stats[base_type]
                stats['count'] += 1
                stats['total_rms'] += rms_info['rms_diff']
                stats['max_rms'] = max(stats['max_rms'], rms_info['rms_diff'])
            
            completed += 1
            if completed % 100 == 0:
                print(f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)",
                      flush=True)
    
    print(f"\n{'=' * 80}")
    print(f"SUMMARY: HIGH RMS DIFFERENCES (≥{threshold} Å)")
    print(f"{'=' * 80}\n")
    
    print(f"Total residues with RMS diff ≥ {threshold} Å: {len(high_rms_residues)}")
    
    # Correlation with atom differences
    with_atom_diff = sum(1 for r in high_rms_residues if r['has_atom_difference'])
    without_atom_diff = len(high_rms_residues) - with_atom_diff
    
    print(f"  With atom set differences: {with_atom_diff} ({with_atom_diff*100/len(high_rms_residues):.1f}%)")
    print(f"  Without atom set differences: {without_atom_diff} ({without_atom_diff*100/len(high_rms_residues):.1f}%)")
    print()
    
    # Base type statistics
    print("BASE TYPE STATISTICS")
    print("-" * 80)
    for base_type in sorted(base_type_stats.keys()):
        stats = base_type_stats[base_type]
        avg_rms = stats['total_rms'] / stats['count'] if stats['count'] > 0 else 0
        print(f"{base_type}: {stats['count']} residues")
        print(f"  Average RMS diff: {avg_rms:.6f} Å")
        print(f"  Maximum RMS diff: {stats['max_rms']:.6f} Å")
        print()
    
    # Sort by RMS difference
    high_rms_residues.sort(key=lambda x: x['rms_diff'], reverse=True)
    
    print("=" * 80)
    print(f"TOP {min(limit, len(high_rms_residues))} HIGHEST RMS DIFFERENCES")
    print("=" * 80)
    print()
    
    for i, rms_info in enumerate(high_rms_residues[:limit], 1):
        print(f"{i}. {rms_info['pdb_id']} {rms_info['residue']}: "
              f"{rms_info['residue_name']} ({rms_info['base_type']})")
        print(f"   RMS: Legacy={rms_info['rms_legacy']:.6f}, "
              f"Modern={rms_info['rms_modern']:.6f}, "
              f"Diff={rms_info['rms_diff']:.6f} Å")
        if rms_info['has_atom_difference']:
            print(f"   ⚠️  Has atom set differences!")
            if rms_info['only_legacy']:
                print(f"      Only in legacy: {rms_info['only_legacy']}")
            if rms_info['only_modern']:
                print(f"      Only in modern: {rms_info['only_modern']}")
        else:
            print(f"   ✓ Same atom sets (RMS difference likely due to numerical precision)")
        print()
    
    # Distribution analysis
    print("=" * 80)
    print("RMS DIFFERENCE DISTRIBUTION")
    print("=" * 80)
    print()
    
    ranges = [
        (0.1, 0.5, "0.1 - 0.5 Å"),
        (0.5, 1.0, "0.5 - 1.0 Å"),
        (1.0, 2.0, "1.0 - 2.0 Å"),
        (2.0, float('inf'), "> 2.0 Å"),
    ]
    
    for min_rms, max_rms, label in ranges:
        count = sum(1 for r in high_rms_residues 
                   if min_rms <= r['rms_diff'] < max_rms)
        if count > 0:
            print(f"{label}: {count} residues ({count*100/len(high_rms_residues):.1f}%)")
    
    # Save detailed report
    output_file = project_root / "docs" / f"HIGH_RMS_ANALYSIS_{threshold}.json"
    with open(output_file, 'w') as f:
        json.dump({
            'threshold': threshold,
            'total_residues': len(high_rms_residues),
            'with_atom_differences': with_atom_diff,
            'without_atom_differences': without_atom_diff,
            'base_type_stats': dict(base_type_stats),
            'top_residues': high_rms_residues[:limit],
        }, f, indent=2)
    
    print(f"\n💾 Detailed report saved to: {output_file}")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Analyze residues with high RMS differences"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="RMS difference threshold in Angstroms (default: 0.1)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=100,
        help="Limit number of examples shown (default: 100)"
    )
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    analyze_high_rms(project_root, args.threshold, args.limit)


if __name__ == "__main__":
    main()

