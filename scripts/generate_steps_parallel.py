#!/usr/bin/env python3
"""Generate step params for all valid fast PDBs using parallel processing."""

import json
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

def process_pdb(pdb):
    try:
        result = subprocess.run(
            ['./build/generate_modern_json', f'data/pdb/{pdb}.pdb', 'data/json_test', '--stage=steps'],
            capture_output=True, timeout=60
        )
        return pdb, result.returncode == 0, None
    except Exception as e:
        return pdb, False, str(e)

def main():
    # Load valid PDBs
    with open('data/valid_pdbs_fast.json') as f:
        data = json.load(f)
    pdbs = data.get('valid_pdbs_with_atoms_and_frames', [])

    print(f"Total PDBs to process: {len(pdbs)}")

    # Check existing
    output_dir = Path('data/json_test/bpstep_params')
    existing = set(f.stem for f in output_dir.glob('*.json')) if output_dir.exists() else set()
    to_process = [p for p in pdbs if p not in existing]
    print(f"Already processed: {len(existing)}, To process: {len(to_process)}")

    if not to_process:
        print("Nothing to process")
        return

    success = 0
    with ProcessPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(process_pdb, pdb): pdb for pdb in to_process}
        for i, future in enumerate(as_completed(futures)):
            pdb, ok, err = future.result()
            if ok:
                success += 1
            if (i + 1) % 100 == 0:
                print(f"  Processed {i+1}/{len(to_process)}, success: {success}")
    print(f"\nDone! {success}/{len(to_process)} succeeded")

if __name__ == '__main__':
    main()
