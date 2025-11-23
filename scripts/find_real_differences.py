#!/usr/bin/env python3
"""
Find only REAL differences between legacy and modern code.

This script:
1. Compares JSON files
2. Identifies only residues with actual differences (not just atom ordering)
3. Reports missing residues, different atom sets, and significant RMS differences
4. Saves a focused report of real issues
"""

import json
import sys
from pathlib import Path
from collections import defaultdict
from x3dna_json_compare import JsonComparator, ComparisonResult
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing

def find_real_differences(project_root: Path):
    """Find residues with real differences (not just ordering)."""
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"
    
    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)
    
    comparator = JsonComparator(enable_cache=False)
    max_workers = multiprocessing.cpu_count()
    
    print("=" * 80)
    print("FINDING REAL DIFFERENCES (excluding atom ordering)")
    print("=" * 80)
    print(f"\nProcessing {len(common)} files...\n")
    
    real_differences = []
    missing_residues = []
    
    def compare_pdb(pdb_id):
        """Compare and return only real differences."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        
        if not all([legacy_file.exists(), modern_file.exists()]):
            return (pdb_id, [], [])
        
        try:
            result = comparator.compare_files(legacy_file, modern_file, None, pdb_id)
            
            pdb_differences = []
            pdb_missing = []
            
            if result.frame_comparison:
                fc = result.frame_comparison
                
                # Collect missing residues (these are real issues)
                for missing in fc.missing_residues:
                    pdb_missing.append({
                        'pdb_id': pdb_id,
                        'residue': f"{missing.chain_id}:{missing.residue_seq}{missing.insertion}",
                        'residue_name': missing.residue_name,
                        'base_type': missing.base_type,
                        'matched_atoms': missing.matched_atoms,
                        'num_matched_atoms': missing.num_matched_atoms,
                        'rms_fit': missing.rms_fit,
                    })
                
                # Collect mismatches with REAL differences
                for mismatch in fc.mismatched_calculations:
                    leg_atoms = set(mismatch.legacy_matched_atoms)
                    mod_atoms = set(mismatch.modern_matched_atoms)
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    
                    # Check for real differences
                    only_legacy = leg_atoms - mod_atoms
                    only_modern = mod_atoms - leg_atoms
                    leg_rms = mismatch.legacy_record.get('rms_fit', 0)
                    mod_rms = mismatch.modern_record.get('rms_fit', 0)
                    rms_diff = abs(leg_rms - mod_rms)
                    
                    # Only include if:
                    # 1. Different atom sets (not just ordering)
                    # 2. Significant RMS difference (> 0.1 Å)
                    has_different_atoms = bool(only_legacy or only_modern)
                    has_significant_rms = rms_diff > 0.1
                    
                    if has_different_atoms or has_significant_rms:
                        pdb_differences.append({
                            'pdb_id': pdb_id,
                            'residue': f"{chain_id}:{residue_seq}{insertion}",
                            'residue_name': mismatch.legacy_record.get('residue_name', ''),
                            'base_type': mismatch.legacy_record.get('base_type', '?'),
                            'legacy_atoms': sorted(leg_atoms),
                            'modern_atoms': sorted(mod_atoms),
                            'only_legacy': sorted(only_legacy),
                            'only_modern': sorted(only_modern),
                            'legacy_count': len(leg_atoms),
                            'modern_count': len(mod_atoms),
                            'rms_legacy': leg_rms,
                            'rms_modern': mod_rms,
                            'rms_diff': rms_diff,
                            'has_different_atoms': has_different_atoms,
                            'has_significant_rms': has_significant_rms,
                        })
            
            return (pdb_id, pdb_differences, pdb_missing)
        except Exception as e:
            return (pdb_id, [], [])
    
    # Process in parallel
    completed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }
        
        for future in as_completed(future_to_pdb):
            pdb_id, differences, missing = future.result()
            if differences or missing:
                real_differences.extend(differences)
                missing_residues.extend(missing)
            completed += 1
            
            if completed % 100 == 0 or completed == len(common):
                print(f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)", flush=True)
    
    print(f"\n{'=' * 80}")
    print("SUMMARY OF REAL DIFFERENCES")
    print(f"{'=' * 80}\n")
    
    print(f"Missing residues (in legacy but not modern): {len(missing_residues)}")
    print(f"Residues with different atom sets: {sum(1 for d in real_differences if d['has_different_atoms'])}")
    print(f"Residues with significant RMS differences (>0.1Å): {sum(1 for d in real_differences if d['has_significant_rms'])}")
    print(f"Total residues with real differences: {len(real_differences)}\n")
    
    # Group by type of difference
    atom_set_diffs = [d for d in real_differences if d['has_different_atoms']]
    rms_diffs = [d for d in real_differences if d['has_significant_rms']]
    
    # Analyze patterns
    print("=" * 80)
    print("ATOM SET DIFFERENCES")
    print("=" * 80)
    
    if atom_set_diffs:
        # Group by pattern
        atom_patterns = defaultdict(list)
        for diff in atom_set_diffs:
            pattern_key = (
                tuple(sorted(diff['only_legacy'])),
                tuple(sorted(diff['only_modern']))
            )
            atom_patterns[pattern_key].append(diff)
        
        print(f"\nUnique atom difference patterns: {len(atom_patterns)}\n")
        
        # Show top patterns
        sorted_patterns = sorted(atom_patterns.items(), key=lambda x: len(x[1]), reverse=True)
        for i, (pattern, diffs) in enumerate(sorted_patterns[:10], 1):
            only_legacy, only_modern = pattern
            print(f"Pattern {i}: {len(diffs)} occurrences")
            if only_legacy:
                print(f"  Only in legacy: {list(only_legacy)}")
            if only_modern:
                print(f"  Only in modern: {list(only_modern)}")
            print(f"  Example: {diffs[0]['pdb_id']} {diffs[0]['residue']} ({diffs[0]['base_type']})")
            print()
    else:
        print("\n✅ No atom set differences found!")
    
    # RMS differences
    print("=" * 80)
    print("SIGNIFICANT RMS DIFFERENCES (>0.1Å)")
    print("=" * 80)
    
    if rms_diffs:
        rms_diffs.sort(key=lambda x: x['rms_diff'], reverse=True)
        print(f"\nTop 20 RMS differences:\n")
        for i, diff in enumerate(rms_diffs[:20], 1):
            print(f"{i}. {diff['pdb_id']} {diff['residue']}: {diff['residue_name']} ({diff['base_type']})")
            print(f"   RMS: Legacy={diff['rms_legacy']:.6f}, Modern={diff['rms_modern']:.6f}, Diff={diff['rms_diff']:.6f}")
            if diff['has_different_atoms']:
                print(f"   Also has different atoms!")
            print()
    else:
        print("\n✅ No significant RMS differences found!")
    
    # Missing residues
    print("=" * 80)
    print("MISSING RESIDUES")
    print("=" * 80)
    
    if missing_residues:
        print(f"\nResidues in legacy but missing in modern:\n")
        for missing in missing_residues[:20]:
            print(f"  {missing['pdb_id']} {missing['residue']}: {missing['residue_name']} ({missing['base_type']})")
            print(f"    Atoms: {missing['matched_atoms']}")
            print()
        if len(missing_residues) > 20:
            print(f"  ... and {len(missing_residues) - 20} more\n")
    else:
        print("\n✅ No missing residues found!")
    
    # Save detailed report
    output_file = project_root / "docs" / "REAL_DIFFERENCES.json"
    with open(output_file, 'w') as f:
        json.dump({
            'summary': {
                'total_real_differences': len(real_differences),
                'missing_residues': len(missing_residues),
                'atom_set_differences': len(atom_set_diffs),
                'significant_rms_differences': len(rms_diffs),
            },
            'missing_residues': missing_residues,
            'atom_set_differences': atom_set_diffs,
            'rms_differences': rms_diffs,
            'all_real_differences': real_differences,
        }, f, indent=2)
    
    print(f"\n💾 Detailed report saved to: {output_file}")
    
    # Generate markdown report
    report_file = project_root / "docs" / "REAL_DIFFERENCES_REPORT.md"
    with open(report_file, 'w') as f:
        f.write("# Real Differences Report\n\n")
        f.write("This report contains only **real differences** - excluding atom ordering.\n\n")
        f.write("## Summary\n\n")
        f.write(f"- **Missing residues** (in legacy but not modern): {len(missing_residues)}\n")
        f.write(f"- **Residues with different atom sets**: {len(atom_set_diffs)}\n")
        f.write(f"- **Residues with significant RMS differences (>0.1Å)**: {len(rms_diffs)}\n")
        f.write(f"- **Total residues with real differences**: {len(real_differences)}\n\n")
        
        if missing_residues:
            f.write("## Missing Residues\n\n")
            f.write("Residues that appear in legacy JSON but are missing from modern JSON:\n\n")
            for missing in missing_residues:
                f.write(f"- **{missing['pdb_id']}** {missing['residue']}: {missing['residue_name']} ({missing['base_type']})\n")
                f.write(f"  - Atoms: {missing['matched_atoms']}\n")
                f.write(f"  - RMS fit: {missing['rms_fit']:.6f}\n\n")
        
        if atom_set_diffs:
            f.write("## Atom Set Differences\n\n")
            f.write("Residues where modern and legacy match different sets of atoms:\n\n")
            for diff in atom_set_diffs[:20]:
                f.write(f"- **{diff['pdb_id']}** {diff['residue']}: {diff['residue_name']} ({diff['base_type']})\n")
                if diff['only_legacy']:
                    f.write(f"  - Only in legacy: {diff['only_legacy']}\n")
                if diff['only_modern']:
                    f.write(f"  - Only in modern: {diff['only_modern']}\n")
                f.write(f"  - Legacy: {diff['legacy_count']} atoms, Modern: {diff['modern_count']} atoms\n")
                f.write(f"  - RMS diff: {diff['rms_diff']:.6f}\n\n")
        
        if rms_diffs:
            f.write("## Significant RMS Differences (>0.1Å)\n\n")
            f.write("Residues with significant differences in RMS fit:\n\n")
            for diff in rms_diffs[:20]:
                f.write(f"- **{diff['pdb_id']}** {diff['residue']}: {diff['residue_name']} ({diff['base_type']})\n")
                f.write(f"  - RMS: Legacy={diff['rms_legacy']:.6f}, Modern={diff['rms_modern']:.6f}, Diff={diff['rms_diff']:.6f}\n\n")
    
    print(f"📄 Markdown report saved to: {report_file}")

if __name__ == '__main__':
    project_root = Path(__file__).parent.parent
    find_real_differences(project_root)

