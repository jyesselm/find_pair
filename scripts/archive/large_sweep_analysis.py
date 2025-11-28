#!/usr/bin/env python3
"""Large-scale analysis across many PDBs."""

import sys
import subprocess
from pathlib import Path
from collections import defaultdict

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

def analyze_pdb(pdb_name, project_root, run_if_needed=True):
    """Analyze differences for a single PDB."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_name}.pdb"
    legacy_inp = project_root / f"{pdb_name}_legacy.inp"
    modern_inp = project_root / f"{pdb_name}_modern.inp"
    
    if not pdb_file.exists():
        return None
    
    # Run find_pair if needed and requested
    if run_if_needed:
        if not legacy_inp.exists():
            try:
                subprocess.run(
                    f"org/build/bin/find_pair_original data/pdb/{pdb_name}.pdb {pdb_name}_legacy.inp",
                    shell=True,
                    cwd=project_root,
                    capture_output=True,
                    timeout=60
                )
            except subprocess.TimeoutExpired:
                return None  # Skip if timeout
        
        if not modern_inp.exists():
            try:
                subprocess.run(
                    f"./build/find_pair_app data/pdb/{pdb_name}.pdb {pdb_name}_modern.inp",
                    shell=True,
                    cwd=project_root,
                    capture_output=True,
                    timeout=60
                )
            except subprocess.TimeoutExpired:
                return None  # Skip if timeout
    
    if not legacy_inp.exists() or not modern_inp.exists():
        return None
    
    legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True)
    modern_pairs = extract_base_pairs(modern_inp, is_legacy=False)
    
    legacy_norm = set(normalize_pair(p1, p2) for p1, p2 in legacy_pairs)
    modern_norm = set(normalize_pair(p1, p2) for p1, p2 in modern_pairs)
    
    common = legacy_norm & modern_norm
    only_legacy = legacy_norm - modern_norm
    only_modern = modern_norm - legacy_norm
    
    return {
        'pdb': pdb_name,
        'legacy_count': len(legacy_pairs),
        'modern_count': len(modern_pairs),
        'common_count': len(common),
        'only_legacy_count': len(only_legacy),
        'only_modern_count': len(only_modern),
        'overlap_ratio': len(common) / max(len(legacy_norm), len(modern_norm)) if max(len(legacy_norm), len(modern_norm)) > 0 else 0.0,
        'match': len(only_legacy) == 0 and len(only_modern) == 0 and len(legacy_norm) > 0
    }

def main():
    project_root = Path(__file__).parent.parent
    
    if len(sys.argv) > 1:
        num_pdbs = int(sys.argv[1])
    else:
        num_pdbs = 200  # Default to 200 PDBs
    
    print(f"Large-scale analysis: Testing {num_pdbs} PDBs...\n")
    
    pdbs = sorted([p.stem for p in (project_root / "data" / "pdb").glob("*.pdb")])[:num_pdbs]
    
    print(f"Processing {len(pdbs)} PDBs...")
    print("(This may take a while)\n")
    
    results = []
    categories = {
        'perfect_match': [],
        'modern_more': [],
        'legacy_more': [],
        'both_different': [],
        'legacy_zero': [],
        'modern_zero': [],
        'no_overlap': []
    }
    
    processed = 0
    for i, pdb in enumerate(pdbs):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(pdbs)}...")
        
        result = analyze_pdb(pdb, project_root, run_if_needed=True)
        
        if result is None:
            continue
        
        results.append(result)
        processed += 1
        
        # Categorize
        if result['match']:
            categories['perfect_match'].append(result)
        elif result['legacy_count'] == 0 and result['modern_count'] > 0:
            categories['legacy_zero'].append(result)
        elif result['modern_count'] == 0 and result['legacy_count'] > 0:
            categories['modern_zero'].append(result)
        elif result['only_modern_count'] > 0 and result['only_legacy_count'] == 0:
            categories['modern_more'].append(result)
        elif result['only_legacy_count'] > 0 and result['only_modern_count'] == 0:
            categories['legacy_more'].append(result)
        elif result['overlap_ratio'] == 0.0 and result['legacy_count'] > 0 and result['modern_count'] > 0:
            categories['no_overlap'].append(result)
            categories['both_different'].append(result)
        else:
            categories['both_different'].append(result)
    
    # Print summary
    print(f"\n{'='*80}")
    print("LARGE-SCALE ANALYSIS SUMMARY")
    print(f"{'='*80}\n")
    
    print(f"Total PDBs processed: {processed}")
    print(f"Successfully analyzed: {len(results)}\n")
    
    total = len(results)
    if total == 0:
        print("No results to analyze.")
        return
    
    print("Categories:")
    print(f"{'='*80}\n")
    
    for cat_name, cat_results in categories.items():
        if cat_results:
            count = len(cat_results)
            pct = (count / total) * 100
            print(f"{cat_name.replace('_', ' ').title()}: {count:4d} ({pct:5.1f}%)")
            
            if cat_name == 'perfect_match':
                avg_pairs = sum(r['legacy_count'] for r in cat_results) / count if count > 0 else 0
                print(f"  Average pairs: {avg_pairs:.1f}")
            elif cat_name == 'modern_more':
                avg_diff = sum(r['only_modern_count'] for r in cat_results) / count if count > 0 else 0
                print(f"  Average additional pairs: {avg_diff:.1f}")
            elif cat_name == 'legacy_more':
                avg_diff = sum(r['only_legacy_count'] for r in cat_results) / count if count > 0 else 0
                print(f"  Average legacy-only pairs: {avg_diff:.1f}")
            elif cat_name == 'no_overlap':
                print(f"  Cases with NO overlap (likely indexing differences)")
    
    # Statistics
    print(f"\n{'='*80}")
    print("STATISTICS")
    print(f"{'='*80}\n")
    
    total_legacy = sum(r['legacy_count'] for r in results)
    total_modern = sum(r['modern_count'] for r in results)
    total_common = sum(r['common_count'] for r in results)
    
    print(f"Total pairs found:")
    print(f"  Legacy: {total_legacy}")
    print(f"  Modern: {total_modern}")
    print(f"  Common: {total_common}")
    print(f"  Modern extra: {total_modern - total_common}")
    print(f"  Legacy extra: {total_legacy - total_common}")
    
    if total > 0:
        print(f"\nAverage pairs per PDB:")
        print(f"  Legacy: {total_legacy / total:.1f}")
        print(f"  Modern: {total_modern / total:.1f}")
    
    # Top differences
    print(f"\n{'='*80}")
    print("TOP 10 CASES: MODERN FINDS MORE")
    print(f"{'='*80}\n")
    
    modern_more_sorted = sorted(categories['modern_more'], 
                                key=lambda x: x['only_modern_count'], reverse=True)[:10]
    for r in modern_more_sorted:
        print(f"  {r['pdb']:8s}: Legacy={r['legacy_count']:3d}, Modern={r['modern_count']:3d}, "
              f"Extra={r['only_modern_count']:3d}")
    
    if categories['no_overlap']:
        print(f"\n{'='*80}")
        print(f"NO OVERLAP CASES (likely indexing differences): {len(categories['no_overlap'])}")
        print(f"{'='*80}\n")
        for r in categories['no_overlap'][:10]:
            print(f"  {r['pdb']:8s}: Legacy={r['legacy_count']:3d}, Modern={r['modern_count']:3d}")

if __name__ == "__main__":
    main()

