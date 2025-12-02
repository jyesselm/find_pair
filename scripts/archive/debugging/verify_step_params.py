#!/usr/bin/env python3
"""
Verify step parameters between legacy and modern code.
"""

import json
import os
import csv
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock


@dataclass
class StepResult:
    """Result for a single PDB step parameter comparison."""
    pdb_id: str
    status: str  # 'perfect', 'difference', 'error', 'missing'
    legacy_steps: int = 0
    modern_steps: int = 0
    matched_steps: int = 0
    max_shift_diff: float = 0.0
    max_slide_diff: float = 0.0
    max_rise_diff: float = 0.0
    max_tilt_diff: float = 0.0
    max_roll_diff: float = 0.0
    max_twist_diff: float = 0.0
    error_message: str = ""


def load_legacy_steps(filepath: str) -> Dict[Tuple[int, int], dict]:
    """Load legacy step parameters (only first occurrence of each bp_idx pair)."""
    if not os.path.exists(filepath):
        return {}
    
    with open(filepath) as f:
        data = json.load(f)
    
    steps = {}
    for rec in data:
        bp1 = rec.get('bp_idx1')
        bp2 = rec.get('bp_idx2')
        if bp1 is not None and bp2 is not None:
            key = (bp1, bp2)
            # Only keep first occurrence (legacy may have duplicates for multiple duplexes)
            if key not in steps:
                params = rec.get('params', {})
                steps[key] = {
                    'shift': params.get('Shift', 0),
                    'slide': params.get('Slide', 0),
                    'rise': params.get('Rise', 0),
                    'tilt': params.get('Tilt', 0),
                    'roll': params.get('Roll', 0),
                    'twist': params.get('Twist', 0),
                }
    return steps


def load_modern_steps(filepath: str) -> Dict[Tuple[int, int], dict]:
    """Load modern step parameters."""
    if not os.path.exists(filepath):
        return {}
    
    with open(filepath) as f:
        data = json.load(f)
    
    steps = {}
    for rec in data:
        bp1 = rec.get('bp_idx1')
        bp2 = rec.get('bp_idx2')
        if bp1 is not None and bp2 is not None:
            steps[(bp1, bp2)] = {
                'shift': rec.get('shift', 0),
                'slide': rec.get('slide', 0),
                'rise': rec.get('rise', 0),
                'tilt': rec.get('tilt', 0),
                'roll': rec.get('roll', 0),
                'twist': rec.get('twist', 0),
            }
    return steps


def compare_steps(legacy: Dict, modern: Dict) -> Tuple[int, Dict[str, float]]:
    """Compare step parameters and return (matched_count, max_diffs)."""
    common_keys = set(legacy.keys()) & set(modern.keys())
    
    max_diffs = {
        'shift': 0.0, 'slide': 0.0, 'rise': 0.0,
        'tilt': 0.0, 'roll': 0.0, 'twist': 0.0
    }
    
    for key in common_keys:
        l = legacy[key]
        m = modern[key]
        for param in max_diffs:
            diff = abs(l.get(param, 0) - m.get(param, 0))
            max_diffs[param] = max(max_diffs[param], diff)
    
    return len(common_keys), max_diffs


class StepVerifier:
    """Verifies step parameters between legacy and modern."""
    
    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        self.legacy_dir = os.path.join(base_dir, 'data/json_legacy/bpstep_params')
        self.modern_dir = os.path.join(base_dir, 'data/json/bpstep_params')
        self.lock = Lock()
        self.processed = 0
        self.total = 0
    
    def compare_pdb(self, pdb_id: str) -> StepResult:
        """Compare step parameters for a single PDB."""
        result = StepResult(pdb_id=pdb_id, status='error')
        
        try:
            legacy_path = os.path.join(self.legacy_dir, f'{pdb_id}.json')
            modern_path = os.path.join(self.modern_dir, f'{pdb_id}.json')
            
            if not os.path.exists(legacy_path):
                result.status = 'missing'
                result.error_message = 'Legacy file missing'
                return result
            
            if not os.path.exists(modern_path):
                result.status = 'missing'
                result.error_message = 'Modern file missing'
                return result
            
            legacy = load_legacy_steps(legacy_path)
            modern = load_modern_steps(modern_path)
            
            result.legacy_steps = len(legacy)
            result.modern_steps = len(modern)
            
            matched, max_diffs = compare_steps(legacy, modern)
            result.matched_steps = matched
            result.max_shift_diff = max_diffs['shift']
            result.max_slide_diff = max_diffs['slide']
            result.max_rise_diff = max_diffs['rise']
            result.max_tilt_diff = max_diffs['tilt']
            result.max_roll_diff = max_diffs['roll']
            result.max_twist_diff = max_diffs['twist']
            
            # Check if all params are within tolerance (0.001 for now)
            tol = 0.001
            all_match = all(v < tol for v in max_diffs.values())
            
            if all_match and result.legacy_steps == result.modern_steps:
                result.status = 'perfect'
            else:
                result.status = 'difference'
                
        except Exception as e:
            result.status = 'error'
            result.error_message = str(e)
        
        with self.lock:
            self.processed += 1
            if self.processed % 10 == 0 or self.processed == self.total:
                print(f"  Progress: {self.processed}/{self.total}")
        
        return result


def main():
    base_dir = '/Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2'
    
    # Find PDBs with both legacy and modern step params
    legacy_dir = os.path.join(base_dir, 'data/json_legacy/bpstep_params')
    modern_dir = os.path.join(base_dir, 'data/json/bpstep_params')
    
    legacy_pdbs = set(f.replace('.json', '') for f in os.listdir(legacy_dir) if f.endswith('.json'))
    modern_pdbs = set(f.replace('.json', '') for f in os.listdir(modern_dir) if f.endswith('.json'))
    common_pdbs = sorted(legacy_pdbs & modern_pdbs)
    
    print("=" * 80)
    print("STEP PARAMETER VERIFICATION")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Legacy step param files: {len(legacy_pdbs)}")
    print(f"Modern step param files: {len(modern_pdbs)}")
    print(f"PDBs with both: {len(common_pdbs)}")
    print("-" * 80)
    
    if not common_pdbs:
        print("No common PDBs found!")
        return
    
    verifier = StepVerifier(base_dir)
    verifier.total = len(common_pdbs)
    
    results = []
    print(f"\nComparing {len(common_pdbs)} PDBs...")
    
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(verifier.compare_pdb, pdb): pdb for pdb in common_pdbs}
        for future in as_completed(futures):
            results.append(future.result())
    
    # Categorize results
    perfect = [r for r in results if r.status == 'perfect']
    differences = [r for r in results if r.status == 'difference']
    errors = [r for r in results if r.status == 'error']
    missing = [r for r in results if r.status == 'missing']
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total compared: {len(results)}")
    print(f"  ✅ Perfect: {len(perfect)}")
    print(f"  ⚠️ Differences: {len(differences)}")
    print(f"  ❌ Errors: {len(errors)}")
    print(f"  ⏭️ Missing: {len(missing)}")
    
    # Calculate overall max differences
    all_results = perfect + differences
    if all_results:
        print("\n" + "-" * 80)
        print("MAX DIFFERENCES ACROSS ALL PDBs")
        print("-" * 80)
        print(f"  Shift: {max(r.max_shift_diff for r in all_results):.6f}")
        print(f"  Slide: {max(r.max_slide_diff for r in all_results):.6f}")
        print(f"  Rise:  {max(r.max_rise_diff for r in all_results):.6f}")
        print(f"  Tilt:  {max(r.max_tilt_diff for r in all_results):.6f}")
        print(f"  Roll:  {max(r.max_roll_diff for r in all_results):.6f}")
        print(f"  Twist: {max(r.max_twist_diff for r in all_results):.6f}")
    
    # Save results to CSV
    csv_path = os.path.join(base_dir, 'data/step_param_verification.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pdb_id', 'status', 'legacy_steps', 'modern_steps', 'matched_steps',
                        'max_shift', 'max_slide', 'max_rise', 'max_tilt', 'max_roll', 'max_twist', 'error'])
        for r in sorted(results, key=lambda x: x.pdb_id):
            writer.writerow([r.pdb_id, r.status, r.legacy_steps, r.modern_steps, r.matched_steps,
                            f'{r.max_shift_diff:.6f}', f'{r.max_slide_diff:.6f}', f'{r.max_rise_diff:.6f}',
                            f'{r.max_tilt_diff:.6f}', f'{r.max_roll_diff:.6f}', f'{r.max_twist_diff:.6f}',
                            r.error_message])
    
    print(f"\nResults saved to: {csv_path}")
    
    # Print differences
    if differences:
        print("\n" + "-" * 80)
        print("PDBs WITH DIFFERENCES")
        print("-" * 80)
        print(f"{'PDB':<10} {'L Steps':<10} {'M Steps':<10} {'Matched':<10} {'Max Diff':<15}")
        for r in differences:
            max_diff = max(r.max_shift_diff, r.max_slide_diff, r.max_rise_diff,
                          r.max_tilt_diff, r.max_roll_diff, r.max_twist_diff)
            print(f"{r.pdb_id:<10} {r.legacy_steps:<10} {r.modern_steps:<10} {r.matched_steps:<10} {max_diff:.6f}")


if __name__ == '__main__':
    main()

