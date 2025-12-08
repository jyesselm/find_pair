#!/usr/bin/env python3
"""
Parallel JSON generation using the C++ generate_modern_json tool.
"""

import subprocess
import sys
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple
import argparse


def process_pdb(args: Tuple[str, str, str, str]) -> Tuple[str, bool, str]:
    """Process a single PDB file."""
    pdb_id, pdb_dir, output_dir, tool_path = args
    pdb_file = Path(pdb_dir) / f"{pdb_id}.pdb"
    
    if not pdb_file.exists():
        return pdb_id, False, "PDB file not found"
    
    try:
        result = subprocess.run(
            [tool_path, str(pdb_file), output_dir, "--stage=all"],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per PDB
        )
        if result.returncode == 0:
            return pdb_id, True, ""
        else:
            return pdb_id, False, result.stderr[:200]
    except subprocess.TimeoutExpired:
        return pdb_id, False, "Timeout"
    except Exception as e:
        return pdb_id, False, str(e)[:200]


def main():
    parser = argparse.ArgumentParser(description="Generate modern JSON in parallel")
    parser.add_argument("--pdb-list", required=True, help="File with PDB IDs")
    parser.add_argument("--pdb-dir", default="data/pdb", help="PDB directory")
    parser.add_argument("--output-dir", default="data/json", help="Output directory")
    parser.add_argument("--workers", "-w", type=int, default=10, help="Number of workers")
    parser.add_argument("--tool", default="build/generate_modern_json", help="Path to C++ tool")
    parser.add_argument("--resume", action="store_true", help="Skip already processed")
    args = parser.parse_args()
    
    # Load PDB list
    with open(args.pdb_list) as f:
        pdb_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    # If resuming, filter out already processed
    if args.resume:
        existing = set()
        for subdir in ["base_pair", "find_bestpair_selection", "pair_validation"]:
            path = Path(args.output_dir) / subdir
            if path.exists():
                existing.update(f.stem for f in path.glob("*.json"))
        
        original_count = len(pdb_ids)
        pdb_ids = [p for p in pdb_ids if p not in existing]
        print(f"Resuming: {original_count - len(pdb_ids)} already done, {len(pdb_ids)} remaining")
    
    print(f"Processing {len(pdb_ids)} PDBs with {args.workers} workers")
    
    # Prepare arguments
    task_args = [(pdb_id, args.pdb_dir, args.output_dir, args.tool) for pdb_id in pdb_ids]
    
    succeeded = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_pdb, arg): arg[0] for arg in task_args}
        
        for i, future in enumerate(as_completed(futures), 1):
            pdb_id = futures[future]
            try:
                _, success, error = future.result()
                if success:
                    succeeded += 1
                else:
                    failed += 1
                    if error:
                        print(f"  {pdb_id}: {error}")
            except Exception as e:
                failed += 1
                print(f"  {pdb_id}: Exception: {e}")
            
            if i % 100 == 0 or i == len(pdb_ids):
                print(f"Progress: {i}/{len(pdb_ids)} ({100*i/len(pdb_ids):.1f}%) - {succeeded} ok, {failed} failed")
    
    print(f"\nComplete: {succeeded} succeeded, {failed} failed")


if __name__ == "__main__":
    main()

