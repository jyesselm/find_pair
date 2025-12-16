#!/usr/bin/env python3
"""Regenerate legacy JSON for slow PDBs with corrupted files."""

import json
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

def is_valid_json(filepath):
    """Check if JSON file is valid."""
    if not filepath.exists():
        return False
    try:
        with open(filepath) as f:
            json.load(f)
        return True
    except:
        return False

def main():
    # Load slow PDBs
    with open('data/slow_pdbs.json') as f:
        data = json.load(f)

    slow_pdbs = [p['pdb_id'] for p in data.get('slow_pdbs', [])]
    print(f"Total slow PDBs: {len(slow_pdbs)}")

    # Find PDBs with invalid legacy JSON
    to_regenerate = []
    for pdb in slow_pdbs:
        json_file = Path(f'data/json_legacy/base_pair/{pdb}.json')
        if not is_valid_json(json_file):
            to_regenerate.append(pdb)

    print(f"PDBs with invalid/missing legacy JSON: {len(to_regenerate)}")

    results = {
        'start_time': datetime.now().isoformat(),
        'total': len(to_regenerate),
        'success': [],
        'failed': [],
        'timeout': []
    }

    env = os.environ.copy()
    env['X3DNA'] = '/Users/jyesselman2/local/installs/x3dna'

    for i, pdb_id in enumerate(to_regenerate):
        print(f"\n[{i+1}/{len(to_regenerate)}] Regenerating {pdb_id}...")

        pdb_file = Path(f'data/pdb/{pdb_id}.pdb')
        if not pdb_file.exists():
            print(f"  ⚠️  PDB file not found")
            results['failed'].append(pdb_id)
            continue

        # Delete old corrupt files
        for stage in ['base_pair', 'pair_validation', 'distance_checks', 'frame_calc',
                      'base_frame_calc', 'hbond_list', 'bpstep_params', 'helical_params',
                      'find_bestpair_selection', 'pdb_atoms', 'residue_indices']:
            old_file = Path(f'data/json_legacy/{stage}/{pdb_id}.json')
            if old_file.exists():
                old_file.unlink()

        try:
            # Run legacy find_pair_analyze with 5 minute timeout
            result = subprocess.run(
                ['./org/build/bin/find_pair_analyze', str(pdb_file)],
                capture_output=True, text=True, timeout=300, env=env,
                cwd=str(Path.cwd())
            )

            # Check if output was generated and valid
            json_file = Path(f'data/json_legacy/base_pair/{pdb_id}.json')
            if is_valid_json(json_file):
                print(f"  ✅ Success")
                results['success'].append(pdb_id)
            else:
                print(f"  ❌ Generated but invalid JSON")
                results['failed'].append(pdb_id)

        except subprocess.TimeoutExpired:
            print(f"  ⏱️  Timeout (5 min)")
            results['timeout'].append(pdb_id)
        except Exception as e:
            print(f"  ❌ Error: {e}")
            results['failed'].append(pdb_id)

        # Progress
        done = len(results['success']) + len(results['failed']) + len(results['timeout'])
        print(f"  Progress: {done}/{len(to_regenerate)}")

    # Save results
    results['end_time'] = datetime.now().isoformat()
    with open('data/slow_pdbs_legacy_regen.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "="*60)
    print("REGENERATION COMPLETE")
    print("="*60)
    print(f"Success: {len(results['success'])}")
    print(f"Failed: {len(results['failed'])}")
    print(f"Timeout: {len(results['timeout'])}")

if __name__ == '__main__':
    main()
