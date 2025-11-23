#!/usr/bin/env python3
"""
Analyze the 69 mismatched residues from comprehensive comparison.

This script:
1. Identifies all mismatched residues
2. Analyzes each one to understand the difference
3. Suggests fixes or legacy mode usage
"""

import json
import sys
from pathlib import Path
from collections import defaultdict
from x3dna_json_compare import JsonComparator, ComparisonResult

def analyze_mismatched_residues(project_root: Path):
    """Find and analyze all mismatched residues."""
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"
    
    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)
    
    comparator = JsonComparator(enable_cache=False)
    
    print("=" * 80)
    print("ANALYZING MISMATCHED RESIDUES")
    print("=" * 80)
    print(f"\nProcessing {len(common)} files to find mismatches...\n")
    
    all_mismatches = []
    files_with_mismatches = defaultdict(list)
    
    # Process files and collect mismatches
    completed = 0
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import multiprocessing
    
    max_workers = multiprocessing.cpu_count()
    
    def compare_pdb(pdb_id):
        """Compare a single PDB file and return mismatches."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        pdb_file = project_root / f"data/pdb/{pdb_id}.pdb"
        
        if not all([legacy_file.exists(), modern_file.exists(), pdb_file.exists()]):
            return (pdb_id, [])
        
        try:
            result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
            mismatches = []
            
            if result.frame_comparison:
                fc = result.frame_comparison
                
                # Collect missing residues
                for missing in fc.missing_residues:
                    mismatches.append({
                        'pdb_id': pdb_id,
                        'type': 'missing',
                        'residue': f"{missing.chain_id}:{missing.residue_seq}{missing.insertion}",
                        'residue_name': missing.residue_name,
                        'base_type': missing.base_type,
                        'matched_atoms': missing.matched_atoms,
                    })
                
                # Collect mismatched calculations
                for mismatch in fc.mismatched_calculations:
                    leg_atoms = set(mismatch.legacy_matched_atoms)
                    mod_atoms = set(mismatch.modern_matched_atoms)
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    
                    mismatches.append({
                        'pdb_id': pdb_id,
                        'type': 'mismatch',
                        'residue': f"{chain_id}:{residue_seq}{insertion}",
                        'residue_name': mismatch.legacy_record.get('residue_name', ''),
                        'base_type': mismatch.legacy_record.get('base_type', '?'),
                        'legacy_atoms': sorted(leg_atoms),
                        'modern_atoms': sorted(mod_atoms),
                        'only_legacy': sorted(leg_atoms - mod_atoms),
                        'only_modern': sorted(mod_atoms - leg_atoms),
                        'rms_legacy': mismatch.legacy_record.get('rms_fit', 0),
                        'rms_modern': mismatch.modern_record.get('rms_fit', 0),
                        'rms_diff': abs(mismatch.legacy_record.get('rms_fit', 0) - mismatch.modern_record.get('rms_fit', 0)),
                    })
            
            return (pdb_id, mismatches)
        except Exception as e:
            return (pdb_id, [{'type': 'error', 'error': str(e)}])
    
    # Process in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }
        
        for future in as_completed(future_to_pdb):
            pdb_id, mismatches = future.result()
            if mismatches:
                all_mismatches.extend(mismatches)
                files_with_mismatches[pdb_id].extend(mismatches)
            completed += 1
            
            if completed % 100 == 0:
                print(f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)", flush=True)
    
    print(f"\n{'=' * 80}")
    print(f"FOUND {len(all_mismatches)} MISMATCHED RESIDUES IN {len(files_with_mismatches)} FILES")
    print(f"{'=' * 80}\n")
    
    # Categorize mismatches
    missing_residues = [m for m in all_mismatches if m.get('type') == 'missing']
    atom_mismatches = [m for m in all_mismatches if m.get('type') == 'mismatch']
    
    print(f"Missing residues: {len(missing_residues)}")
    print(f"Atom mismatches: {len(atom_mismatches)}")
    
    # Analyze atom mismatches
    print(f"\n{'=' * 80}")
    print("ATOM MISMATCH ANALYSIS")
    print(f"{'=' * 80}\n")
    
    # Group by pattern
    patterns = defaultdict(list)
    
    for mismatch in atom_mismatches:
        if mismatch['only_legacy'] or mismatch['only_modern']:
            pattern_key = (
                tuple(sorted(mismatch['only_legacy'])),
                tuple(sorted(mismatch['only_modern']))
            )
            patterns[pattern_key].append(mismatch)
    
    print(f"Unique mismatch patterns: {len(patterns)}\n")
    
    # Print top patterns
    sorted_patterns = sorted(patterns.items(), key=lambda x: len(x[1]), reverse=True)
    
    for i, (pattern, mismatches) in enumerate(sorted_patterns[:10], 1):
        only_legacy, only_modern = pattern
        print(f"Pattern {i}: {len(mismatches)} occurrences")
        if only_legacy:
            print(f"  Only in legacy: {list(only_legacy)}")
        if only_modern:
            print(f"  Only in modern: {list(only_modern)}")
        print(f"  Example: {mismatches[0]['pdb_id']} {mismatches[0]['residue']}")
        print()
    
    # Save detailed report
    output_file = project_root / "docs" / "MISMATCHED_RESIDUES_DETAILED.json"
    with open(output_file, 'w') as f:
        json.dump({
            'total_mismatches': len(all_mismatches),
            'total_files': len(files_with_mismatches),
            'missing_residues': missing_residues,
            'atom_mismatches': atom_mismatches,
            'patterns': {
                str(k): len(v) for k, v in sorted_patterns
            }
        }, f, indent=2)
    
    print(f"\nDetailed report saved to: {output_file}")
    
    # Print all mismatched residues
    print(f"\n{'=' * 80}")
    print("ALL MISMATCHED RESIDUES")
    print(f"{'=' * 80}\n")
    
    for mismatch in sorted(atom_mismatches, key=lambda x: (x['pdb_id'], x['residue'])):
        print(f"{mismatch['pdb_id']} {mismatch['residue']}: {mismatch['residue_name']} ({mismatch['base_type']})")
        if mismatch['only_legacy']:
            print(f"  Only in legacy: {mismatch['only_legacy']}")
        if mismatch['only_modern']:
            print(f"  Only in modern: {mismatch['only_modern']}")
        if mismatch['rms_diff'] > 0.001:
            print(f"  RMS diff: {mismatch['rms_diff']:.6f}")
        print()
    
    return all_mismatches, files_with_mismatches

if __name__ == '__main__':
    project_root = Path(__file__).parent.parent
    analyze_mismatched_residues(project_root)

