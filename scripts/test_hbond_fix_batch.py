#!/usr/bin/env python3
"""
Test the H-bond fix for modified nucleotides on a batch of PDBs
This script tests that the fix works correctly and doesn't cause regressions
Uses multiprocessing to run PDBs in parallel
"""

import json
import subprocess
import sys
import os
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial

def load_test_set(test_set_file):
    """Load PDB IDs from test set file"""
    with open(test_set_file) as f:
        data = json.load(f)
    return data.get('pdb_ids', [])

def process_single_pdb(args):
    """Process a single PDB - designed for parallel execution"""
    pdb_id, output_dir = args
    pdb_file = f"data/pdb/{pdb_id}.pdb"
    json_file = f"{output_dir}/{pdb_id}.json"
    
    if not os.path.exists(pdb_file):
        return {
            'pdb_id': pdb_id,
            'status': 'failed',
            'error': f"PDB file not found: {pdb_file}"
        }
    
    try:
        result = subprocess.run(
            ['build/generate_modern_json', pdb_file, json_file],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per PDB
        )
        if result.returncode != 0:
            return {
                'pdb_id': pdb_id,
                'status': 'failed',
                'error': result.stderr[:200]  # Truncate long errors
            }
        
        # Check for modified nucleotides
        modified = check_modified_nucleotides(json_file)
        has_modified = len(modified) > 0
        
        # Analyze H-bond types
        stats = analyze_hbond_types(json_file)
        if stats is None:
            return {
                'pdb_id': pdb_id,
                'status': 'no_hbonds'
            }
        
        if 'error' in stats:
            return {
                'pdb_id': pdb_id,
                'status': 'analysis_error',
                'error': stats['error']
            }
        
        status = 'success_modified' if has_modified else 'success'
        
        return {
            'pdb_id': pdb_id,
            'status': status,
            'stats': stats,
            'has_modified': has_modified,
            'modified_count': len(modified) if has_modified else 0,
            'modified': modified[:3] if has_modified else []
        }
    except subprocess.TimeoutExpired:
        return {
            'pdb_id': pdb_id,
            'status': 'failed',
            'error': 'Timeout (exceeded 5 minutes)'
        }
    except Exception as e:
        return {
            'pdb_id': pdb_id,
            'status': 'failed',
            'error': f"Exception: {str(e)}"
        }

def analyze_hbond_types(json_file):
    """Analyze H-bond types in the JSON file"""
    try:
        hbond_file = json_file.replace('.json', '_hbond_list.json')
        if not os.path.exists(hbond_file):
            return None
        
        with open(hbond_file) as f:
            data = json.load(f)
        
        stats = {
            'total_pairs': 0,
            'total_hbonds': 0,
            'standard_hbonds': 0,  # type='-'
            'non_standard_hbonds': 0,  # type='*'
            'invalid_hbonds': 0,  # type=' '
            'good_hbonds': 0,  # type='-' and dist in [2.5, 3.5]
            'modified_nucleotide_pairs': 0,
        }
        
        for record in data:
            if record.get('type') != 'hbond_list':
                continue
            
            stats['total_pairs'] += 1
            hbonds = record.get('hbonds', [])
            stats['total_hbonds'] += len(hbonds)
            
            for hb in hbonds:
                hb_type = hb.get('type', ' ')
                dist = hb.get('distance', 0)
                
                if hb_type == '-':
                    stats['standard_hbonds'] += 1
                    if 2.5 <= dist <= 3.5:
                        stats['good_hbonds'] += 1
                elif hb_type == '*':
                    stats['non_standard_hbonds'] += 1
                else:
                    stats['invalid_hbonds'] += 1
        
        return stats
    except Exception as e:
        return {'error': str(e)}

def check_modified_nucleotides(json_file):
    """Check if PDB has modified nucleotides"""
    try:
        base_frame_file = json_file.replace('.json', '_base_frame_calc.json')
        if not os.path.exists(base_frame_file):
            return []
        
        with open(base_frame_file) as f:
            data = json.load(f)
        
        modified = []
        for record in data:
            residue_name = record.get('residue_name', '').strip()
            base_type = record.get('base_type', '')
            
            # Check if it's a modified nucleotide (not standard A, C, G, T, U)
            standard_names = ['A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT', 'DU', 
                            'ADE', 'CYT', 'GUA', 'THY', 'URA']
            if residue_name not in standard_names and base_type in ['a', 'c', 'g', 't', 'u', 'A', 'C', 'G', 'T', 'U']:
                modified.append({
                    'residue_name': residue_name,
                    'base_type': base_type,
                    'residue_idx': record.get('residue_idx', -1)
                })
        
        return modified
    except Exception as e:
        return []

def main():
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)
    
    # Use test_set_50 for reasonable testing time
    test_set_file = project_root / "data" / "test_sets" / "test_set_50.json"
    if not test_set_file.exists():
        print(f"Test set file not found: {test_set_file}")
        sys.exit(1)
    
    pdb_ids = load_test_set(test_set_file)
    print(f"Testing H-bond fix on {len(pdb_ids)} PDBs from {test_set_file.name}")
    print("=" * 70)
    
    output_dir = "data/json_test_batch"
    os.makedirs(output_dir, exist_ok=True)
    
    # Use all available CPUs
    num_workers = cpu_count()
    print(f"Using {num_workers} parallel workers\n")
    
    # Prepare arguments for parallel processing
    args_list = [(pdb_id, output_dir) for pdb_id in pdb_ids]
    
    # Process PDBs in parallel
    print("Processing PDBs in parallel...")
    with Pool(processes=num_workers) as pool:
        results = pool.map(process_single_pdb, args_list)
    
    # Collect modified nucleotide PDBs
    modified_nucleotide_pdbs = []
    for result in results:
        if result.get('has_modified'):
            modified_nucleotide_pdbs.append({
                'pdb_id': result['pdb_id'],
                'modified_count': result.get('modified_count', 0),
                'modified': result.get('modified', [])
            })
    
    # Print results
    for i, result in enumerate(results, 1):
        pdb_id = result['pdb_id']
        status = result.get('status', 'unknown')
        
        if status == 'failed':
            error = result.get('error', 'Unknown error')
            print(f"[{i}/{len(results)}] {pdb_id}: ❌ FAILED - {error[:50]}")
        elif status == 'no_hbonds':
            print(f"[{i}/{len(results)}] {pdb_id}: ⚠️  No H-bond data")
        elif status in ['success', 'success_modified']:
            stats = result.get('stats', {})
            has_modified = result.get('has_modified', False)
            print(f"[{i}/{len(results)}] {pdb_id}: ✅ {stats.get('total_pairs', 0)} pairs, "
                  f"{stats.get('total_hbonds', 0)} H-bonds "
                  f"({stats.get('standard_hbonds', 0)} standard, {stats.get('good_hbonds', 0)} good)")
            if has_modified:
                print(f"   📝 {result.get('modified_count', 0)} modified nucleotides")
        else:
            print(f"[{i}/{len(results)}] {pdb_id}: ⚠️  {status}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    successful = [r for r in results if r.get('status') in ['success', 'success_modified']]
    failed = [r for r in results if r.get('status') == 'failed']
    no_hbonds = [r for r in results if r.get('status') == 'no_hbonds']
    
    print(f"\nTotal PDBs processed: {len(results)}")
    print(f"  ✅ Successful: {len(successful)}")
    print(f"  ❌ Failed: {len(failed)}")
    print(f"  ⚠️  No H-bonds: {len(no_hbonds)}")
    
    if modified_nucleotide_pdbs:
        print(f"\n📝 PDBs with modified nucleotides: {len(modified_nucleotide_pdbs)}")
        for pdb_info in modified_nucleotide_pdbs[:10]:  # Show first 10
            print(f"  - {pdb_info['pdb_id']}: {pdb_info['modified_count']} modified nucleotides")
            for mod in pdb_info['modified'][:2]:
                print(f"    {mod['residue_name']} -> base_type='{mod['base_type']}'")
    
    # Aggregate statistics
    if successful:
        total_pairs = sum(r['stats']['total_pairs'] for r in successful if 'stats' in r)
        total_hbonds = sum(r['stats']['total_hbonds'] for r in successful if 'stats' in r)
        total_standard = sum(r['stats']['standard_hbonds'] for r in successful if 'stats' in r)
        total_good = sum(r['stats']['good_hbonds'] for r in successful if 'stats' in r)
        
        print(f"\n📊 Aggregate Statistics:")
        print(f"  Total base pairs: {total_pairs}")
        print(f"  Total H-bonds: {total_hbonds}")
        print(f"  Standard H-bonds (type='-'): {total_standard} ({total_standard/total_hbonds*100:.1f}%)" if total_hbonds > 0 else "  Standard H-bonds: 0")
        print(f"  Good H-bonds (type='-' and dist in [2.5, 3.5]): {total_good}")
    
    # Save results
    results_file = "data/json_test_batch/test_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            'summary': {
                'total': len(results),
                'successful': len(successful),
                'failed': len(failed),
                'no_hbonds': len(no_hbonds),
                'modified_nucleotide_pdbs': len(modified_nucleotide_pdbs)
            },
            'results': results,
            'modified_nucleotide_pdbs': modified_nucleotide_pdbs
        }, f, indent=2)
    
    print(f"\n💾 Results saved to: {results_file}")
    print("=" * 70)
    
    return 0 if len(failed) == 0 else 1

if __name__ == '__main__':
    sys.exit(main())

