#!/usr/bin/env python3
"""
Validate frame calculations (base_frame_calc, ls_fitting, frame_calc) on all fast PDBs.

This script:
1. Loads legacy and modern JSON for frame calculations
2. Uses frame_comparison to compare (with automatic deduplication of legacy duplicates)
3. Reports success rate and saves results

Usage:
    python3 scripts/validate_frames_batch.py
"""

import json
import sys
from pathlib import Path
from typing import List, Dict

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from x3dna_json_compare.frame_comparison import compare_frames


def load_fast_pdbs() -> List[str]:
    """Load list of fast PDBs from JSON."""
    fast_pdbs_file = Path("data/valid_pdbs_fast.json")
    with open(fast_pdbs_file) as f:
        data = json.load(f)
    return data.get("valid_pdbs_with_atoms_and_frames", [])


def validate_single_pdb(pdb_id: str) -> Dict:
    """Validate frames for a single PDB."""
    result = {
        'pdb_id': pdb_id,
        'status': 'pending',
        'base_frame_calc': None,
        'ls_fitting': None,
        'frame_calc': None,
        'error': None
    }
    
    for frame_type in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
        legacy_file = Path(f"data/json_legacy/{frame_type}/{pdb_id}.json")
        modern_file = Path(f"data/json/{frame_type}/{pdb_id}.json")
        
        if not legacy_file.exists():
            result[frame_type] = 'no_legacy'
            continue
        
        if not modern_file.exists():
            result[frame_type] = 'no_modern'
            result['error'] = f"Missing modern {frame_type}"
            result['status'] = 'failed'
            continue
        
        try:
            with open(legacy_file) as f:
                legacy_data = json.load(f)
            with open(modern_file) as f:
                modern_data = json.load(f)
            
            # Compare (with automatic deduplication)
            comparison = compare_frames(legacy_data, modern_data, None, tolerance=1e-6)
            
            # Check if perfect match
            perfect = (comparison.matched == comparison.total_legacy and 
                      len(comparison.fp_only) == 0 and
                      len(comparison.mismatched) == 0)
            
            fp_only = len(comparison.fp_only) > 0 and len(comparison.mismatched) == 0
            
            if perfect:
                result[frame_type] = 'perfect'
            elif fp_only:
                result[frame_type] = 'fp_only'
            else:
                result[frame_type] = 'mismatch'
                result['error'] = f"{frame_type}: {comparison.matched}/{comparison.total_legacy} matched"
                result['status'] = 'failed'
                
        except Exception as e:
            result[frame_type] = 'error'
            result['error'] = f"{frame_type}: {str(e)}"
            result['status'] = 'failed'
    
    # Overall status
    if result['status'] == 'pending':
        # Check if all frame types passed
        all_passed = all(
            result[ft] in ['perfect', 'fp_only', 'no_legacy'] 
            for ft in ['base_frame_calc', 'ls_fitting', 'frame_calc']
        )
        result['status'] = 'passed' if all_passed else 'failed'
    
    return result


def main():
    print("Frame Validation on All Fast PDBs")
    print("="*60)
    
    # Load fast PDBs
    fast_pdbs = load_fast_pdbs()
    print(f"Total PDBs to validate: {len(fast_pdbs)}")
    
    results = []
    passed = 0
    failed = 0
    
    for i, pdb_id in enumerate(fast_pdbs, 1):
        result = validate_single_pdb(pdb_id)
        results.append(result)
        
        if result['status'] == 'passed':
            passed += 1
            if i <= 10 or passed % 100 == 0:
                print(f"✅ {pdb_id} ({i}/{len(fast_pdbs)})")
        else:
            failed += 1
            print(f"❌ {pdb_id}: {result.get('error', 'unknown error')} ({i}/{len(fast_pdbs)})")
        
        # Progress update every 500 PDBs
        if i % 500 == 0:
            success_rate = 100 * passed / i
            print(f"\n📊 Progress: {i}/{len(fast_pdbs)} ({success_rate:.2f}% passed)\n")
    
    # Final summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    print(f"Total PDBs: {len(results)}")
    print(f"Passed: {passed} ({100*passed/len(results):.2f}%)")
    print(f"Failed: {failed} ({100*failed/len(results):.2f}%)")
    
    # Save results
    output_file = Path("data/validation_results/frames_validation_full.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump({
            'total': len(results),
            'passed': passed,
            'failed': failed,
            'success_rate': passed / len(results),
            'results': results
        }, f, indent=2)
    
    print(f"\n✅ Results saved to {output_file}")
    
    # Show failed PDBs if any
    if failed > 0:
        print(f"\n❌ Failed PDBs ({failed}):")
        for r in results:
            if r['status'] == 'failed':
                print(f"   {r['pdb_id']}: {r.get('error', 'unknown')}")
    
    return 0 if passed == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())

