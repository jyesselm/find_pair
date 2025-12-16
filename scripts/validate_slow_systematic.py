#!/usr/bin/env python3
"""Systematically validate slow PDBs through all stages."""

import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime

def main():
    # Load slow PDBs
    with open('data/slow_pdbs.json') as f:
        data = json.load(f)

    slow_pdbs = [p['pdb_id'] for p in data.get('slow_pdbs', [])]
    print(f"Total slow PDBs: {len(slow_pdbs)}")

    # Stage groups to test (in order)
    stage_groups = [
        ('atoms', '1', 'Atoms'),
        ('frames', 'frames', 'Frames (3-5)'),
        ('pairs', 'pairs', 'Pairs (6-10)'),
    ]

    results = {
        'start_time': datetime.now().isoformat(),
        'total': len(slow_pdbs),
        'stages': {}
    }

    for stage_key, stage_arg, stage_name in stage_groups:
        print(f"\n{'='*60}")
        print(f"VALIDATING STAGE: {stage_name}")
        print(f"{'='*60}")

        results['stages'][stage_key] = {
            'passed': [],
            'failed': [],
            'errors': []
        }

        for i, pdb_id in enumerate(slow_pdbs):
            # Check if modern JSON exists
            modern_json = Path(f'data/json/base_pair/{pdb_id}.json')
            legacy_json = Path(f'data/json_legacy/base_pair/{pdb_id}.json')

            if not modern_json.exists() or not legacy_json.exists():
                print(f"[{i+1}/{len(slow_pdbs)}] {pdb_id}: SKIP (missing JSON)")
                results['stages'][stage_key]['errors'].append(pdb_id)
                continue

            # Run validation
            try:
                result = subprocess.run(
                    ['fp2-validate', 'validate', stage_arg, '--pdb', pdb_id],
                    capture_output=True, text=True, timeout=60
                )

                if result.returncode == 0:
                    results['stages'][stage_key]['passed'].append(pdb_id)
                    status = "✅ PASSED"
                else:
                    results['stages'][stage_key]['failed'].append(pdb_id)
                    status = "❌ FAILED"

            except subprocess.TimeoutExpired:
                results['stages'][stage_key]['errors'].append(pdb_id)
                status = "⏱️ TIMEOUT"
            except Exception as e:
                results['stages'][stage_key]['errors'].append(pdb_id)
                status = f"⚠️ ERROR: {e}"

            passed = len(results['stages'][stage_key]['passed'])
            total = passed + len(results['stages'][stage_key]['failed'])
            rate = (passed / total * 100) if total > 0 else 0

            if (i+1) % 50 == 0 or i < 10:
                print(f"[{i+1}/{len(slow_pdbs)}] {pdb_id}: {status} ({rate:.1f}% pass rate)")

        # Stage summary
        stage_data = results['stages'][stage_key]
        passed = len(stage_data['passed'])
        failed = len(stage_data['failed'])
        errors = len(stage_data['errors'])
        total = passed + failed
        rate = (passed / total * 100) if total > 0 else 0

        print(f"\n{stage_name} Summary:")
        print(f"  Passed: {passed}/{total} ({rate:.1f}%)")
        print(f"  Failed: {failed}")
        print(f"  Errors/Skip: {errors}")

        # If too many failures at this stage, stop
        if rate < 80 and total > 20:
            print(f"\n⚠️ Pass rate below 80% at {stage_name} stage.")
            print("Stopping here - need to fix this stage before proceeding.")
            break

    # Save results
    results['end_time'] = datetime.now().isoformat()
    with open('data/slow_pdbs_systematic_validation.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("\nResults saved to data/slow_pdbs_systematic_validation.json")

if __name__ == '__main__':
    main()
