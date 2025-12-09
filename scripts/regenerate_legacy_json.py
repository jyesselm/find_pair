#!/usr/bin/env python3
"""
Regenerate legacy JSON files using the updated code.

Uses the valid_pdbs_fast.json list and runs in parallel.
"""

import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple


PROJECT_ROOT = Path(__file__).parent.parent
FIND_PAIR_BINARY = PROJECT_ROOT / "org" / "build" / "bin" / "find_pair_analyze"
PDB_DIR = PROJECT_ROOT / "data" / "pdb"
VALID_PDBS_FILE = PROJECT_ROOT / "data" / "valid_pdbs_fast.json"


def get_valid_pdbs() -> List[str]:
    """Load the list of valid PDBs from JSON."""
    with open(VALID_PDBS_FILE) as f:
        data = json.load(f)
    return data.get("valid_pdbs_with_atoms_and_frames", [])


def run_find_pair(pdb_id: str) -> Tuple[str, bool, str]:
    """Run find_pair_analyze on a single PDB.
    
    Returns:
        Tuple of (pdb_id, success, error_message)
    """
    pdb_file = PDB_DIR / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        return (pdb_id, False, f"PDB file not found: {pdb_file}")
    
    try:
        result = subprocess.run(
            [str(FIND_PAIR_BINARY), str(pdb_file)],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per PDB
        )
        if result.returncode != 0:
            return (pdb_id, False, f"Exit code {result.returncode}")
        return (pdb_id, True, "")
    except subprocess.TimeoutExpired:
        return (pdb_id, False, "Timeout after 5 minutes")
    except Exception as e:
        return (pdb_id, False, str(e))


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Regenerate legacy JSON files")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers")
    parser.add_argument("--max-pdbs", type=int, help="Maximum number of PDBs to process")
    args = parser.parse_args()
    
    # Check binary exists
    if not FIND_PAIR_BINARY.exists():
        print(f"Error: Binary not found: {FIND_PAIR_BINARY}")
        print("Run: cd org/build && cmake --build . -j8")
        sys.exit(1)
    
    # Get PDB list
    pdbs = get_valid_pdbs()
    if args.max_pdbs:
        pdbs = pdbs[:args.max_pdbs]
    
    print(f"Regenerating JSON for {len(pdbs)} PDBs with {args.workers} workers...")
    print(f"Binary: {FIND_PAIR_BINARY}")
    print()
    
    # Process in parallel
    completed = 0
    failed = 0
    failed_pdbs = []
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_find_pair, pdb): pdb for pdb in pdbs}
        
        for future in as_completed(futures):
            pdb_id, success, error = future.result()
            completed += 1
            
            if success:
                status = "✓"
            else:
                status = "✗"
                failed += 1
                failed_pdbs.append((pdb_id, error))
            
            # Progress update every 100 or on failure
            if completed % 100 == 0 or not success:
                pct = completed / len(pdbs) * 100
                print(f"[{completed}/{len(pdbs)}] {pct:.1f}% - {pdb_id}: {status}")
                if not success:
                    print(f"  Error: {error}")
    
    # Summary
    print()
    print("=" * 60)
    print(f"COMPLETE: {completed - failed}/{completed} succeeded, {failed} failed")
    
    if failed_pdbs:
        print(f"\nFailed PDBs ({len(failed_pdbs)}):")
        for pdb_id, error in failed_pdbs[:10]:
            print(f"  {pdb_id}: {error}")
        if len(failed_pdbs) > 10:
            print(f"  ... and {len(failed_pdbs) - 10} more")


if __name__ == "__main__":
    main()

