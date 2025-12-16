#!/usr/bin/env python3
"""Validate slow PDBs that already have both modern and legacy JSON."""

import json
import subprocess
from pathlib import Path
from datetime import datetime

def main():
    # Load slow PDBs
    with open('data/slow_pdbs.json') as f:
        data = json.load(f)

    slow_pdbs = [p['pdb_id'] for p in data.get('slow_pdbs', [])]

    # Filter to those with both JSON files
    valid_pdbs = []
    for pdb in slow_pdbs:
        modern = Path(f'data/json/base_pair/{pdb}.json')
        legacy = Path(f'data/json_legacy/base_pair/{pdb}.json')
        if modern.exists() and legacy.exists():
            valid_pdbs.append(pdb)

    print(f"Slow PDBs with both JSON files: {len(valid_pdbs)}")

    # Validate atoms stage first
    print(f"\n{'='*60}")
    print("STAGE 1: ATOMS")
    print(f"{'='*60}")

    atoms_passed = []
    atoms_failed = []

    for i, pdb_id in enumerate(valid_pdbs):
        result = subprocess.run(
            ['fp2-validate', 'validate', '1', '--pdb', pdb_id],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0:
            atoms_passed.append(pdb_id)
            status = "✅"
        else:
            atoms_failed.append(pdb_id)
            status = "❌"

        if (i+1) % 20 == 0 or (i+1) == len(valid_pdbs):
            rate = len(atoms_passed) / (i+1) * 100
            print(f"Progress: {i+1}/{len(valid_pdbs)} ({rate:.1f}% pass rate)")

    print(f"\nAtoms: {len(atoms_passed)}/{len(valid_pdbs)} passed ({len(atoms_passed)/len(valid_pdbs)*100:.1f}%)")
    print(f"Failed: {atoms_failed[:10]}{'...' if len(atoms_failed) > 10 else ''}")

    # If atoms pass rate is high, continue to frames
    if len(atoms_passed) / len(valid_pdbs) >= 0.9:
        print(f"\n{'='*60}")
        print("STAGE 3-5: FRAMES")
        print(f"{'='*60}")

        frames_passed = []
        frames_failed = []

        for i, pdb_id in enumerate(atoms_passed):
            result = subprocess.run(
                ['fp2-validate', 'validate', 'frames', '--pdb', pdb_id],
                capture_output=True, text=True, timeout=60
            )
            if result.returncode == 0:
                frames_passed.append(pdb_id)
            else:
                frames_failed.append(pdb_id)

            if (i+1) % 20 == 0 or (i+1) == len(atoms_passed):
                rate = len(frames_passed) / (i+1) * 100
                print(f"Progress: {i+1}/{len(atoms_passed)} ({rate:.1f}% pass rate)")

        print(f"\nFrames: {len(frames_passed)}/{len(atoms_passed)} passed ({len(frames_passed)/len(atoms_passed)*100:.1f}%)")
        print(f"Failed: {frames_failed[:10]}{'...' if len(frames_failed) > 10 else ''}")

        # If frames pass rate is high, continue to pairs
        if len(frames_passed) / len(atoms_passed) >= 0.9:
            print(f"\n{'='*60}")
            print("STAGE 6-10: PAIRS")
            print(f"{'='*60}")

            pairs_passed = []
            pairs_failed = []

            for i, pdb_id in enumerate(frames_passed):
                result = subprocess.run(
                    ['fp2-validate', 'validate', 'pairs', '--pdb', pdb_id],
                    capture_output=True, text=True, timeout=60
                )
                if result.returncode == 0:
                    pairs_passed.append(pdb_id)
                else:
                    pairs_failed.append(pdb_id)

                if (i+1) % 20 == 0 or (i+1) == len(frames_passed):
                    rate = len(pairs_passed) / (i+1) * 100
                    print(f"Progress: {i+1}/{len(frames_passed)} ({rate:.1f}% pass rate)")

            print(f"\nPairs: {len(pairs_passed)}/{len(frames_passed)} passed ({len(pairs_passed)/len(frames_passed)*100:.1f}%)")
            print(f"Failed: {pairs_failed[:10]}{'...' if len(pairs_failed) > 10 else ''}")

    # Save results
    results = {
        'total_with_json': len(valid_pdbs),
        'atoms_passed': atoms_passed,
        'atoms_failed': atoms_failed,
    }
    if len(atoms_passed) / len(valid_pdbs) >= 0.9:
        results['frames_passed'] = frames_passed
        results['frames_failed'] = frames_failed
        if len(frames_passed) / len(atoms_passed) >= 0.9:
            results['pairs_passed'] = pairs_passed
            results['pairs_failed'] = pairs_failed

    with open('data/slow_pdbs_stage_validation.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("\nResults saved to data/slow_pdbs_stage_validation.json")

if __name__ == '__main__':
    main()
