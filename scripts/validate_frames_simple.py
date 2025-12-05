#!/usr/bin/env python3
"""
Simple frame validation without complex dependencies.

Validates base_frame_calc and ls_fitting by:
1. Deduplicating legacy records by residue_idx
2. Comparing counts
3. Comparing values within tolerance
"""

import json
from pathlib import Path
from typing import List, Dict, Tuple


def deduplicate_legacy(records: List[Dict]) -> List[Dict]:
    """Deduplicate legacy records by residue_idx (keep first occurrence)."""
    seen = set()
    deduped = []
    for rec in records:
        residue_idx = rec.get('residue_idx')
        if residue_idx and residue_idx in seen:
            continue
        if residue_idx:
            seen.add(residue_idx)
        deduped.append(rec)
    return deduped


def compare_base_frame_calc(legacy: List[Dict], modern: List[Dict], tolerance: float = 1e-6) -> Tuple[bool, str]:
    """Compare base_frame_calc records."""
    legacy_dedup = deduplicate_legacy(legacy)
    
    if len(legacy_dedup) != len(modern):
        return False, f"Count mismatch: {len(legacy_dedup)} vs {len(modern)}"
    
    # Match by residue_idx
    legacy_by_idx = {r['residue_idx']: r for r in legacy_dedup}
    modern_by_idx = {r['residue_idx']: r for r in modern}
    
    for idx in legacy_by_idx:
        if idx not in modern_by_idx:
            return False, f"Missing residue_idx {idx} in modern"
        
        leg = legacy_by_idx[idx]
        mod = modern_by_idx[idx]
        
        # Compare RMS fit
        if abs(leg['rms_fit'] - mod['rms_fit']) > tolerance:
            return False, f"RMS mismatch at {idx}: {leg['rms_fit']} vs {mod['rms_fit']}"
        
        # Compare num_matched_atoms
        if leg['num_matched_atoms'] != mod['num_matched_atoms']:
            return False, f"Atom count mismatch at {idx}"
    
    return True, "Perfect match"


def compare_ls_fitting(legacy: List[Dict], modern: List[Dict], tolerance: float = 1e-6) -> Tuple[bool, str]:
    """Compare ls_fitting records."""
    legacy_dedup = deduplicate_legacy(legacy)
    
    if len(legacy_dedup) != len(modern):
        return False, f"Count mismatch: {len(legacy_dedup)} vs {len(modern)}"
    
    # Match by residue_idx
    legacy_by_idx = {r['residue_idx']: r for r in legacy_dedup}
    modern_by_idx = {r['residue_idx']: r for r in modern}
    
    for idx in legacy_by_idx:
        if idx not in modern_by_idx:
            return False, f"Missing residue_idx {idx} in modern"
        
        leg = legacy_by_idx[idx]
        mod = modern_by_idx[idx]
        
        # Compare RMS fit
        if abs(leg['rms_fit'] - mod['rms_fit']) > tolerance:
            return False, f"RMS mismatch at {idx}: {leg['rms_fit']} vs {mod['rms_fit']}"
        
        # Compare num_points
        if leg['num_points'] != mod['num_points']:
            return False, f"Point count mismatch at {idx}"
    
    return True, "Perfect match"


def validate_pdb(pdb_id: str) -> Dict:
    """Validate frames for a single PDB."""
    result = {
        'pdb_id': pdb_id,
        'status': 'passed',
        'base_frame_calc': 'not_tested',
        'ls_fitting': 'not_tested',
        'error': None
    }
    
    # Test base_frame_calc
    legacy_bf = Path(f"data/json_legacy/base_frame_calc/{pdb_id}.json")
    modern_bf = Path(f"data/json/base_frame_calc/{pdb_id}.json")
    
    if legacy_bf.exists() and modern_bf.exists():
        try:
            with open(legacy_bf) as f:
                legacy_data = json.load(f)
            with open(modern_bf) as f:
                modern_data = json.load(f)
            
            match, msg = compare_base_frame_calc(legacy_data, modern_data)
            if match:
                result['base_frame_calc'] = 'perfect'
            else:
                result['base_frame_calc'] = 'mismatch'
                result['error'] = f"base_frame_calc: {msg}"
                result['status'] = 'failed'
        except Exception as e:
            result['base_frame_calc'] = 'error'
            result['error'] = f"base_frame_calc: {str(e)}"
            result['status'] = 'failed'
    elif legacy_bf.exists():
        result['base_frame_calc'] = 'no_modern'
        result['error'] = 'Missing modern base_frame_calc'
        result['status'] = 'failed'
    
    # Test ls_fitting
    legacy_ls = Path(f"data/json_legacy/ls_fitting/{pdb_id}.json")
    modern_ls = Path(f"data/json/ls_fitting/{pdb_id}.json")
    
    if legacy_ls.exists() and modern_ls.exists():
        try:
            with open(legacy_ls) as f:
                legacy_data = json.load(f)
            with open(modern_ls) as f:
                modern_data = json.load(f)
            
            match, msg = compare_ls_fitting(legacy_data, modern_data)
            if match:
                result['ls_fitting'] = 'perfect'
            else:
                result['ls_fitting'] = 'mismatch'
                if not result['error']:
                    result['error'] = f"ls_fitting: {msg}"
                result['status'] = 'failed'
        except Exception as e:
            result['ls_fitting'] = 'error'
            if not result['error']:
                result['error'] = f"ls_fitting: {str(e)}"
            result['status'] = 'failed'
    elif legacy_ls.exists():
        result['ls_fitting'] = 'no_modern'
        if not result['error']:
            result['error'] = 'Missing modern ls_fitting'
        result['status'] = 'failed'
    
    return result


def main():
    print("Simple Frame Validation (All Fast PDBs)")
    print("="*60)
    
    # Load fast PDBs
    with open("data/valid_pdbs_fast.json") as f:
        data = json.load(f)
    fast_pdbs = data['valid_pdbs_with_atoms_and_frames']
    
    print(f"Total PDBs: {len(fast_pdbs)}")
    print("Validating...")
    
    results = []
    passed = 0
    failed = 0
    
    for i, pdb_id in enumerate(fast_pdbs, 1):
        result = validate_pdb(pdb_id)
        results.append(result)
        
        if result['status'] == 'passed':
            passed += 1
        else:
            failed += 1
            if failed <= 20:  # Show first 20 failures
                print(f"❌ {pdb_id}: {result['error']}")
        
        # Progress every 500
        if i % 500 == 0:
            print(f"Progress: {i}/{len(fast_pdbs)} ({100*passed/i:.1f}% passed)")
    
    # Summary
    print("\n" + "="*60)
    print(f"RESULTS:")
    print(f"  Total: {len(results)}")
    print(f"  Passed: {passed} ({100*passed/len(results):.2f}%)")
    print(f"  Failed: {failed} ({100*failed/len(results):.2f}%)")
    
    # Save results
    output_file = Path("data/validation_results/frames_validation_simple.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump({
            'total': len(results),
            'passed': passed,
            'failed': failed,
            'success_rate': passed / len(results) if results else 0,
            'results': results
        }, f, indent=2)
    
    print(f"\n✅ Results saved to {output_file}")
    
    if passed == len(results):
        print("\n🎉 100% MATCH ACHIEVED!")
    
    return 0 if passed == len(results) else 1


if __name__ == '__main__':
    import sys
    sys.exit(main())

