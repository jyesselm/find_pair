#!/usr/bin/env python3
"""Validate slow PDBs one at a time for pairs stage."""

import json
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

def main():
    # Load slow PDBs
    with open('data/slow_pdbs.json') as f:
        data = json.load(f)

    slow_pdbs = [p['pdb_id'] for p in data.get('slow_pdbs', [])]
    print(f"Total slow PDBs to validate: {len(slow_pdbs)}")

    results = {
        'start_time': datetime.now().isoformat(),
        'total': len(slow_pdbs),
        'passed': [],
        'failed': [],
        'errors': [],
        'current_idx': 0
    }

    results_file = Path('data/slow_pdbs_pairs_validation.json')

    # Check if we have existing progress
    if results_file.exists():
        with open(results_file) as f:
            existing = json.load(f)
        # Resume from where we left off
        done_pdbs = set(existing.get('passed', []) + existing.get('failed', []) + existing.get('errors', []))
        results['passed'] = existing.get('passed', [])
        results['failed'] = existing.get('failed', [])
        results['errors'] = existing.get('errors', [])
        print(f"Resuming from previous run: {len(done_pdbs)} already done")
    else:
        done_pdbs = set()

    env = os.environ.copy()
    env['X3DNA'] = '/Users/jyesselman2/local/installs/x3dna'

    for i, pdb_id in enumerate(slow_pdbs):
        if pdb_id in done_pdbs:
            continue

        results['current_idx'] = i
        print(f"\n[{i+1}/{len(slow_pdbs)}] Processing {pdb_id}...")

        pdb_file = Path(f'data/pdb/{pdb_id}.pdb')
        if not pdb_file.exists():
            print(f"  ⚠️  PDB file not found: {pdb_file}")
            results['errors'].append(pdb_id)
            continue

        # Check if JSON files exist, regenerate if not
        modern_json = Path(f'data/json/base_pair/{pdb_id}.json')
        legacy_json = Path(f'data/json_legacy/base_pair/{pdb_id}.json')

        try:
            # Regenerate modern JSON if missing
            if not modern_json.exists():
                print(f"  Generating modern JSON...")
                result = subprocess.run(
                    ['./build/generate_modern_json', str(pdb_file), 'data/json', '--stage=all'],
                    capture_output=True, text=True, timeout=300
                )
                if result.returncode != 0:
                    print(f"  ❌ Modern JSON generation failed")
                    results['errors'].append(pdb_id)
                    continue

            # Regenerate legacy JSON if missing
            if not legacy_json.exists():
                print(f"  Generating legacy JSON...")
                result = subprocess.run(
                    ['./org/build/bin/find_pair_analyze', str(pdb_file)],
                    capture_output=True, text=True, timeout=300, env=env,
                    cwd=str(Path.cwd())
                )
                if result.returncode != 0:
                    print(f"  ❌ Legacy JSON generation failed")
                    results['errors'].append(pdb_id)
                    continue

            # Validate pairs
            print(f"  Validating pairs...")
            result = subprocess.run(
                ['fp2-validate', 'validate', 'pairs', '--pdb', pdb_id],
                capture_output=True, text=True, timeout=120
            )

            if result.returncode == 0:
                print(f"  ✅ PASSED")
                results['passed'].append(pdb_id)
            else:
                print(f"  ❌ FAILED")
                results['failed'].append(pdb_id)

        except subprocess.TimeoutExpired:
            print(f"  ⏱️  Timeout")
            results['errors'].append(pdb_id)
        except Exception as e:
            print(f"  ❌ Error: {e}")
            results['errors'].append(pdb_id)

        # Save progress after each PDB
        results['last_update'] = datetime.now().isoformat()
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        # Print running stats
        total_done = len(results['passed']) + len(results['failed']) + len(results['errors'])
        pass_rate = len(results['passed']) / total_done * 100 if total_done > 0 else 0
        print(f"  Progress: {total_done}/{len(slow_pdbs)} ({pass_rate:.1f}% pass rate)")

    # Final summary
    results['end_time'] = datetime.now().isoformat()
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "="*60)
    print("FINAL RESULTS")
    print("="*60)
    print(f"Total: {len(slow_pdbs)}")
    print(f"Passed: {len(results['passed'])} ({len(results['passed'])/len(slow_pdbs)*100:.1f}%)")
    print(f"Failed: {len(results['failed'])} ({len(results['failed'])/len(slow_pdbs)*100:.1f}%)")
    print(f"Errors: {len(results['errors'])} ({len(results['errors'])/len(slow_pdbs)*100:.1f}%)")

    if results['failed']:
        print(f"\nFailed PDBs: {results['failed'][:20]}...")

    return 0 if not results['failed'] else 1

if __name__ == '__main__':
    sys.exit(main())
