#!/usr/bin/env python3
"""Regenerate modern JSON for slow PDBs that already have legacy JSON."""

import json
import subprocess
from pathlib import Path
from datetime import datetime
import sys

def main():
    # Load slow PDBs
    with open('data/slow_pdbs.json') as f:
        data = json.load(f)

    slow_pdbs = [p['pdb_id'] for p in data.get('slow_pdbs', [])]

    # Filter to those with existing modern JSON (these are faster)
    to_regenerate = []
    for pdb in slow_pdbs:
        modern = Path(f'data/json/base_pair/{pdb}.json')
        legacy = Path(f'data/json_legacy/base_pair/{pdb}.json')
        pdb_file = Path(f'data/pdb/{pdb}.pdb')
        if modern.exists() and legacy.exists() and pdb_file.exists():
            to_regenerate.append(pdb)

    print(f"Slow PDBs to regenerate: {len(to_regenerate)}")

    results = {
        'start_time': datetime.now().isoformat(),
        'total': len(to_regenerate),
        'success': [],
        'failed': [],
    }

    for i, pdb_id in enumerate(to_regenerate):
        print(f"[{i+1}/{len(to_regenerate)}] Regenerating {pdb_id}...", end=" ", flush=True)

        pdb_file = Path(f'data/pdb/{pdb_id}.pdb')

        try:
            result = subprocess.run(
                ['./build/generate_modern_json', str(pdb_file), 'data/json', '--stage=all'],
                capture_output=True, text=True, timeout=600
            )

            if result.returncode == 0:
                results['success'].append(pdb_id)
                print("✅")
            else:
                results['failed'].append(pdb_id)
                print("❌")
                if result.stderr:
                    print(f"  Error: {result.stderr[:200]}")

        except subprocess.TimeoutExpired:
            results['failed'].append(pdb_id)
            print("⏱️ TIMEOUT")
        except Exception as e:
            results['failed'].append(pdb_id)
            print(f"❌ Error: {e}")

    # Save results
    results['end_time'] = datetime.now().isoformat()
    with open('data/slow_pdbs_regeneration.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n{'='*60}")
    print("REGENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"Success: {len(results['success'])}/{len(to_regenerate)}")
    print(f"Failed: {len(results['failed'])}")

if __name__ == '__main__':
    main()
