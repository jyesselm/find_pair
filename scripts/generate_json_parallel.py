#!/usr/bin/env python3
"""Generate modern JSON for all PDBs using parallel processing."""

import subprocess
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import sys

PDB_DIR = Path("data/pdb")
OUTPUT_DIR = Path("data/json_test")
PROGRESS_FILE = Path("data/json_generation_helical_progress.json")
NUM_WORKERS = 20


def generate_single(pdb: str) -> dict:
    """Generate JSON for a single PDB."""
    pdb_file = PDB_DIR / f"{pdb}.pdb"
    if not pdb_file.exists():
        return {"pdb": pdb, "success": False, "error": "PDB file not found"}

    try:
        result = subprocess.run(
            ["./build/generate_modern_json", str(pdb_file), str(OUTPUT_DIR), "--stage=helical"],
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode == 0:
            return {"pdb": pdb, "success": True, "error": None}
        else:
            return {"pdb": pdb, "success": False, "error": result.stderr[:200]}
    except subprocess.TimeoutExpired:
        return {"pdb": pdb, "success": False, "error": "Timeout"}
    except Exception as e:
        return {"pdb": pdb, "success": False, "error": str(e)[:200]}


def main():
    # Get list of PDBs with legacy data
    legacy_dir = Path("data/json_legacy/helical_params")
    legacy_pdbs = {f.stem for f in legacy_dir.glob("*.json")}

    # Get already generated PDBs
    modern_dir = OUTPUT_DIR / "helical_params"
    modern_pdbs = {f.stem for f in modern_dir.glob("*.json")} if modern_dir.exists() else set()

    # Find missing PDBs
    missing = sorted(legacy_pdbs - modern_pdbs)
    print(f"Legacy: {len(legacy_pdbs)}, Modern: {len(modern_pdbs)}, Missing: {len(missing)}")

    if not missing:
        print("All PDBs already generated!")
        return

    # Load checkpoint if exists
    if PROGRESS_FILE.exists():
        with open(PROGRESS_FILE) as f:
            progress = json.load(f)
        completed = set(progress.get("completed", []))
        failed = progress.get("failed", [])
    else:
        completed = set()
        failed = []

    # Filter out already attempted
    to_process = [p for p in missing if p not in completed]
    print(f"To process: {len(to_process)} PDBs with {NUM_WORKERS} workers")
    print()

    if not to_process:
        print("All PDBs already processed!")
        return

    # Process in batches
    batch_size = 100
    start_time = time.time()
    success_count = 0
    fail_count = 0

    for batch_start in range(0, len(to_process), batch_size):
        batch = to_process[batch_start:batch_start + batch_size]
        batch_num = batch_start // batch_size + 1
        total_batches = (len(to_process) + batch_size - 1) // batch_size

        print(f"Processing batch {batch_num}/{total_batches} ({len(batch)} PDBs)...")

        with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
            futures = {executor.submit(generate_single, pdb): pdb for pdb in batch}

            for future in as_completed(futures):
                result = future.result()
                pdb = result["pdb"]
                completed.add(pdb)

                if result["success"]:
                    success_count += 1
                else:
                    fail_count += 1
                    failed.append({"pdb": pdb, "error": result["error"]})

        # Save checkpoint
        with open(PROGRESS_FILE, "w") as f:
            json.dump({"completed": list(completed), "failed": failed}, f)

        elapsed = time.time() - start_time
        done = batch_start + len(batch)
        rate = done / elapsed if elapsed > 0 else 0
        remaining = len(to_process) - done
        eta = remaining / rate if rate > 0 else 0

        print(f"  Done: {done}/{len(to_process)} | Success: {success_count} | Failed: {fail_count} | ETA: {eta/60:.1f}m")
        print()

    # Final summary
    print("=" * 60)
    print("GENERATION COMPLETE")
    print("=" * 60)
    print(f"Total processed: {len(to_process)}")
    print(f"Succeeded: {success_count}")
    print(f"Failed: {fail_count}")

    if failed:
        print(f"\nFirst 10 failures:")
        for f in failed[:10]:
            print(f"  {f['pdb']}: {f['error']}")


if __name__ == "__main__":
    main()
