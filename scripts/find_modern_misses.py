#!/usr/bin/env python3
"""Find cases where modern code missed pairs that legacy found."""

import sys
from pathlib import Path
import subprocess

def extract_base_pairs(inp_file, is_legacy=False):
    """Extract base pairs from .inp file."""
    if not Path(inp_file).exists():
        return []
    pairs = []
    with open(inp_file) as f:
        lines = f.readlines()
        for line in lines[5:]:
            if line.strip().startswith('#') or line.strip().startswith('#####'):
                continue
            parts = line.split()
            if is_legacy:
                if len(parts) >= 2:
                    try:
                        res1, res2 = int(parts[0]), int(parts[1])
                        if res1 > 10 and res2 > 10:
                            pairs.append((res1, res2))
                    except ValueError:
                        continue
            else:
                if len(parts) >= 3 and parts[0].isdigit():
                    try:
                        res1, res2 = int(parts[1]), int(parts[2])
                        pairs.append((res1, res2))
                    except ValueError:
                        continue
    return pairs

def normalize_pair(p1, p2):
    """Normalize a pair to always have smaller index first."""
    return (min(p1, p2), max(p1, p2))

def main():
    project_root = Path(__file__).parent.parent
    
    if len(sys.argv) > 1:
        pdbs = sys.argv[1:]
    else:
        # Test on all PDBs we have .inp files for
        pdbs = sorted([p.stem.replace('_legacy', '').replace('_modern', '') 
                      for p in project_root.glob("*_legacy.inp")])[:100]
    
    print("Finding cases where modern missed pairs that legacy found...\n")
    print(f"{'='*80}\n")
    
    modern_misses = []
    
    for pdb in pdbs:
        legacy_inp = project_root / f"{pdb}_legacy.inp"
        modern_inp = project_root / f"{pdb}_modern.inp"
        
        if not legacy_inp.exists() or not modern_inp.exists():
            continue
        
        # Run find_pair if needed
        if not legacy_inp.exists():
            subprocess.run(
                f"org/build/bin/find_pair_original data/pdb/{pdb}.pdb {pdb}_legacy.inp",
                shell=True,
                cwd=project_root,
                capture_output=True
            )
        
        if not modern_inp.exists():
            subprocess.run(
                f"./build/find_pair_app data/pdb/{pdb}.pdb {pdb}_modern.inp",
                shell=True,
                cwd=project_root,
                capture_output=True
            )
        
        legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True)
        modern_pairs = extract_base_pairs(modern_inp, is_legacy=False)
        
        legacy_norm = set(normalize_pair(p1, p2) for p1, p2 in legacy_pairs)
        modern_norm = set(normalize_pair(p1, p2) for p1, p2 in modern_pairs)
        
        only_legacy = legacy_norm - modern_norm
        only_modern = modern_norm - legacy_norm
        
        if only_legacy:
            modern_misses.append({
                'pdb': pdb,
                'legacy_count': len(legacy_pairs),
                'modern_count': len(modern_pairs),
                'legacy_only': len(only_legacy),
                'modern_only': len(only_modern),
                'legacy_only_pairs': sorted(only_legacy),
                'modern_only_pairs': sorted(only_modern)
            })
    
    if not modern_misses:
        print("✅ GOOD NEWS: Modern code never missed pairs that legacy found!")
        print("   (Based on tested PDBs)")
        return
    
    print(f"Found {len(modern_misses)} cases where modern missed pairs:\n")
    
    # Sort by number of missed pairs
    modern_misses.sort(key=lambda x: x['legacy_only'], reverse=True)
    
    total_missed = 0
    total_extra = 0
    
    for case in modern_misses:
        print(f"{case['pdb']}:")
        print(f"  Legacy: {case['legacy_count']} pairs")
        print(f"  Modern: {case['modern_count']} pairs")
        print(f"  ❌ Modern MISSED: {case['legacy_only']} pairs")
        if case['modern_only'] > 0:
            print(f"  ⚠️  Modern found EXTRA: {case['modern_only']} pairs")
        print(f"  Missed pairs: {case['legacy_only_pairs'][:5]}")
        if len(case['legacy_only_pairs']) > 5:
            print(f"    ... and {len(case['legacy_only_pairs']) - 5} more")
        print()
        
        total_missed += case['legacy_only']
        total_extra += case['modern_only']
    
    print(f"{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}\n")
    print(f"Total cases where modern missed pairs: {len(modern_misses)}")
    print(f"Total pairs missed: {total_missed}")
    print(f"Total pairs modern found extra: {total_extra}")
    print(f"\nAverage missed per case: {total_missed/len(modern_misses):.1f}")
    print(f"Average extra per case: {total_extra/len(modern_misses):.1f}")

if __name__ == "__main__":
    main()

