#!/usr/bin/env python3
"""
Validate all PDBs for Stages 11 and 12 with checkpointing and parallel processing.
"""

import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import sys
import os

# Add tests_python to path
sys.path.insert(0, str(Path(__file__).parent.parent / "tests_python"))
from validation.runner import validate_pdb

JSON_DIR = Path("data/json_test")
CHECKPOINT_FILE = Path("data/validation_checkpoint_stage11_12.json")


def validate_single(pdb: str) -> dict:
    """Validate a single PDB for stages 11 and 12"""
    try:
        r11 = validate_pdb(pdb, 11, JSON_DIR, verbose=False)
        r12 = validate_pdb(pdb, 12, JSON_DIR, verbose=False)
        return {
            "pdb": pdb,
            "stage11_pass": r11.passed,
            "stage12_pass": r12.passed,
            "stage11_skip": r11.skipped,
            "stage12_skip": r12.skipped,
            "stage11_errors": len(r11.errors) if not r11.passed and not r11.skipped else 0,
            "stage12_errors": len(r12.errors) if not r12.passed and not r12.skipped else 0,
            "error": None,
        }
    except Exception as e:
        return {
            "pdb": pdb,
            "stage11_pass": False,
            "stage12_pass": False,
            "stage11_skip": False,
            "stage12_skip": False,
            "error": str(e),
        }


def main():
    # Load fast PDB list
    with open("data/valid_pdbs_fast.json") as f:
        data = json.load(f)
    pdbs = data.get("valid_pdbs_with_atoms_and_frames", [])

    print(f"Validating {len(pdbs)} PDBs for Stages 11 and 12 with 8 workers...")
    print()

    # Check if we have a checkpoint
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f:
            checkpoint = json.load(f)
        completed = set(checkpoint.get("completed", []))
        results = checkpoint.get("results", {})
        print(f"Resuming from checkpoint: {len(completed)} already completed")
    else:
        completed = set()
        results = {
            "stage11_pass": [],
            "stage11_fail": [],
            "stage11_skip": [],
            "stage12_pass": [],
            "stage12_fail": [],
            "stage12_skip": [],
        }

    # Filter to only uncompleted
    remaining = [p for p in pdbs if p not in completed]
    print(f"Remaining: {len(remaining)} PDBs")
    print()

    if not remaining:
        print("All PDBs already validated!")
    else:
        # Process in batches with checkpointing
        batch_size = 100
        start_time = time.time()

        for batch_start in range(0, len(remaining), batch_size):
            batch = remaining[batch_start : batch_start + batch_size]
            batch_num = batch_start // batch_size + 1
            total_batches = (len(remaining) + batch_size - 1) // batch_size

            print(f"Processing batch {batch_num}/{total_batches} ({len(batch)} PDBs)...")

            with ProcessPoolExecutor(max_workers=8) as executor:
                futures = {executor.submit(validate_single, pdb): pdb for pdb in batch}

                for future in as_completed(futures):
                    result = future.result()
                    pdb = result["pdb"]

                    if result["error"]:
                        print(f'  ERROR {pdb}: {result["error"]}')
                        # Save checkpoint before exiting
                        with open(CHECKPOINT_FILE, "w") as f:
                            json.dump({"completed": list(completed), "results": results}, f)
                        sys.exit(1)

                    completed.add(pdb)

                    if result.get("stage11_skip"):
                        results["stage11_skip"].append(pdb)
                    elif result["stage11_pass"]:
                        results["stage11_pass"].append(pdb)
                    else:
                        results["stage11_fail"].append(pdb)

                    if result.get("stage12_skip"):
                        results["stage12_skip"].append(pdb)
                    elif result["stage12_pass"]:
                        results["stage12_pass"].append(pdb)
                    else:
                        results["stage12_fail"].append(pdb)

            # Save checkpoint
            with open(CHECKPOINT_FILE, "w") as f:
                json.dump({"completed": list(completed), "results": results}, f)

            elapsed = time.time() - start_time
            done = len(completed)
            rate = done / elapsed if elapsed > 0 else 0
            remaining_count = len(pdbs) - done
            eta = remaining_count / rate if rate > 0 else 0

            s11_pass = len(results["stage11_pass"])
            s11_fail = len(results["stage11_fail"])
            s11_skip = len(results.get("stage11_skip", []))
            s12_pass = len(results["stage12_pass"])
            s12_fail = len(results["stage12_fail"])
            s12_skip = len(results.get("stage12_skip", []))

            print(
                f"  Completed: {done}/{len(pdbs)} | S11: {s11_pass}✓/{s11_fail}✗/{s11_skip}⊘ | S12: {s12_pass}✓/{s12_fail}✗/{s12_skip}⊘ | ETA: {eta/60:.1f}m"
            )
            print()

    # Final summary
    print("=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    total = len(pdbs)
    s11_pass = len(results["stage11_pass"])
    s11_fail = len(results["stage11_fail"])
    s11_skip = len(results.get("stage11_skip", []))
    s12_pass = len(results["stage12_pass"])
    s12_fail = len(results["stage12_fail"])
    s12_skip = len(results.get("stage12_skip", []))

    # Calculate pass rate excluding skips
    s11_tested = s11_pass + s11_fail
    s12_tested = s12_pass + s12_fail
    s11_pct = 100 * s11_pass / s11_tested if s11_tested > 0 else 0
    s12_pct = 100 * s12_pass / s12_tested if s12_tested > 0 else 0

    print(f"Stage 11 (Step Parameters):    {s11_pass}/{s11_tested} pass ({s11_pct:.1f}%), {s11_skip} skipped")
    print(f"Stage 12 (Helical Parameters): {s12_pass}/{s12_tested} pass ({s12_pct:.1f}%), {s12_skip} skipped")

    if results["stage11_fail"]:
        print(f'\nStage 11 failures ({len(results["stage11_fail"])}): first 10 = {results["stage11_fail"][:10]}')
    if results["stage12_fail"]:
        print(f'Stage 12 failures ({len(results["stage12_fail"])}): first 10 = {results["stage12_fail"][:10]}')


if __name__ == "__main__":
    main()

