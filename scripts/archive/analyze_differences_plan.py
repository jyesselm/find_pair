#!/usr/bin/env python3
"""Plan and analyze differences between legacy and modern base pair detection."""

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

def analyze_pdb(pdb_name, project_root):
    """Analyze differences for a single PDB."""
    legacy_inp = project_root / f"{pdb_name}_legacy.inp"
    modern_inp = project_root / f"{pdb_name}_modern.inp"
    
    legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True) if legacy_inp.exists() else []
    modern_pairs = extract_base_pairs(modern_inp, is_legacy=False) if modern_inp.exists() else []
    
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
        'common_pairs': common,
        'only_legacy_pairs': only_legacy,
        'only_modern_pairs': only_modern,
        'match': len(only_legacy) == 0 and len(only_modern) == 0
    }

def main():
    if len(sys.argv) < 2:
        # Test on first 50 PDBs
        project_root = Path(__file__).parent.parent
        pdbs = sorted([p.stem for p in (project_root / "data" / "pdb").glob("*.pdb")])[:50]
        print(f"Testing {len(pdbs)} PDBs...")
    else:
        pdbs = sys.argv[1:]
    
    project_root = Path(__file__).parent.parent
    
    results = []
    categories = {
        'perfect_match': [],
        'modern_more': [],
        'legacy_more': [],
        'both_different': [],
        'legacy_zero': [],
        'modern_zero': []
    }
    
    for pdb in pdbs:
        # Run find_pair if needed
        legacy_inp = project_root / f"{pdb}_legacy.inp"
        modern_inp = project_root / f"{pdb}_modern.inp"
        
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
        
        result = analyze_pdb(pdb, project_root)
        results.append(result)
        
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
        else:
            categories['both_different'].append(result)
    
    # Print summary
    print(f"\n{'='*80}")
    print("DIFFERENCE ANALYSIS SUMMARY")
    print(f"{'='*80}\n")
    
    print(f"Total PDBs analyzed: {len(results)}\n")
    
    for cat_name, cat_results in categories.items():
        if cat_results:
            print(f"{cat_name.replace('_', ' ').title()}: {len(cat_results)}")
            
            # Show statistics
            if cat_name == 'perfect_match':
                avg_pairs = sum(r['legacy_count'] for r in cat_results) / len(cat_results)
                print(f"  Average pairs: {avg_pairs:.1f}")
            elif cat_name == 'modern_more':
                avg_diff = sum(r['only_modern_count'] for r in cat_results) / len(cat_results)
                print(f"  Average additional pairs: {avg_diff:.1f}")
                # Show top 5
                top5 = sorted(cat_results, key=lambda x: x['only_modern_count'], reverse=True)[:5]
                print(f"  Top differences: {[(r['pdb'], r['only_modern_count']) for r in top5]}")
            elif cat_name == 'legacy_more':
                avg_diff = sum(r['only_legacy_count'] for r in cat_results) / len(cat_results)
                print(f"  Average legacy-only pairs: {avg_diff:.1f}")
                # Show top 5
                top5 = sorted(cat_results, key=lambda x: x['only_legacy_count'], reverse=True)[:5]
                print(f"  Top differences: {[(r['pdb'], r['only_legacy_count']) for r in top5]}")
            elif cat_name == 'both_different':
                print(f"  Examples: {[r['pdb'] for r in cat_results[:5]]}")
    
    # Create investigation plan
    print(f"\n{'='*80}")
    print("INVESTIGATION PLAN")
    print(f"{'='*80}\n")
    
    print("1. PERFECT MATCHES (verify parameters match):")
    if categories['perfect_match']:
        print(f"   Test parameter comparison on: {[r['pdb'] for r in categories['perfect_match'][:5]]}")
    
    print("\n2. MODERN FINDS MORE (investigate why):")
    if categories['modern_more']:
        top_modern = sorted(categories['modern_more'], key=lambda x: x['only_modern_count'], reverse=True)[:5]
        print(f"   Priority cases: {[(r['pdb'], r['only_modern_count']) for r in top_modern]}")
        print("   Actions:")
        print("     - Check validation thresholds")
        print("     - Verify if additional pairs are valid")
        print("     - Compare hydrogen bond detection")
    
    print("\n3. LEGACY FINDS MORE (investigate why modern misses):")
    if categories['legacy_more']:
        top_legacy = sorted(categories['legacy_more'], key=lambda x: x['only_legacy_count'], reverse=True)[:5]
        print(f"   Priority cases: {[(r['pdb'], r['only_legacy_count']) for r in top_legacy]}")
        print("   Actions:")
        print("     - Check why modern validation rejects these pairs")
        print("     - Compare overlap calculations")
        print("     - Check frame calculation differences")
    
    print("\n4. LEGACY FINDS NONE BUT MODERN FINDS PAIRS:")
    if categories['legacy_zero']:
        print(f"   Cases: {[r['pdb'] for r in categories['legacy_zero'][:10]]}")
        print("   Actions:")
        print("     - Check PDB parsing differences")
        print("     - Verify if pairs are valid")
    
    print("\n5. NEXT STEPS:")
    print("   a) Run parameter comparison on perfect matches")
    print("   b) Deep dive into top 5 cases where modern finds more")
    print("   c) Deep dive into cases where legacy finds more")
    print("   d) Compare validation thresholds between legacy and modern")

if __name__ == "__main__":
    main()

