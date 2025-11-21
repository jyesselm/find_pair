#!/usr/bin/env python3
"""
Compare pdb_atoms records between generated and legacy JSON files.
"""

import json
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


def compare_pdb_atoms(gen_file, leg_file):
    """Compare pdb_atoms records between generated and legacy JSON."""
    try:
        with open(gen_file) as f:
            gen_json = json.load(f)
        with open(leg_file) as f:
            leg_json = json.load(f)

        gen_atoms = [
            c for c in gen_json.get("calculations", []) if c.get("type") == "pdb_atoms"
        ]
        leg_atoms = [
            c for c in leg_json.get("calculations", []) if c.get("type") == "pdb_atoms"
        ]

        if not gen_atoms or not leg_atoms:
            return (gen_file.stem, False, "Missing pdb_atoms record")

        gen_rec = gen_atoms[0]
        leg_rec = leg_atoms[0]

        gen_atom_list = gen_rec.get("atoms", [])
        leg_atom_list = leg_rec.get("atoms", [])

        # Build maps using (chain_id, residue_seq, insertion, atom_name) as key
        gen_map = {}
        for atom in gen_atom_list:
            key = (
                atom.get("chain_id", ""),
                atom.get("residue_seq", 0),
                atom.get("insertion", ""),
                atom.get("atom_name", ""),
            )
            gen_map[key] = atom

        leg_map = {}
        for atom in leg_atom_list:
            key = (
                atom.get("chain_id", ""),
                atom.get("residue_seq", 0),
                atom.get("insertion", ""),
                atom.get("atom_name", ""),
            )
            leg_map[key] = atom

        gen_keys = set(gen_map.keys())
        leg_keys = set(leg_map.keys())

        missing = len(leg_keys - gen_keys)
        extra = len(gen_keys - leg_keys)

        # Check field mismatches
        common_keys = gen_keys & leg_keys
        mismatches = 0
        for key in common_keys:
            gen_atom = gen_map[key]
            leg_atom = leg_map[key]

            # Compare all fields
            for field in [
                "atom_name",
                "residue_name",
                "chain_id",
                "residue_seq",
                "record_type",
                "insertion",
            ]:
                if gen_atom.get(field) != leg_atom.get(field):
                    mismatches += 1
                    break

            # Check coordinates
            gen_xyz = gen_atom.get("xyz", [])
            leg_xyz = leg_atom.get("xyz", [])
            if len(gen_xyz) == 3 and len(leg_xyz) == 3:
                for j in range(3):
                    if abs(gen_xyz[j] - leg_xyz[j]) > 0.0001:
                        mismatches += 1
                        break

        if missing == 0 and extra == 0 and mismatches == 0:
            return (gen_file.stem, True, None)
        else:
            issues = []
            if missing > 0:
                issues.append(f"missing:{missing}")
            if extra > 0:
                issues.append(f"extra:{extra}")
            if mismatches > 0:
                issues.append(f"mismatches:{mismatches}")
            return (gen_file.stem, False, ", ".join(issues))

    except Exception as e:
        return (gen_file.stem, False, f"Error: {str(e)}")


def main():
    gen_dir = Path("data/json")
    leg_dir = Path("data/json_legacy")

    # Find all generated files
    gen_files = sorted(gen_dir.glob("*.json"))
    print(f"Found {len(gen_files)} generated JSON files")
    print(f"Comparing pdb_atoms records...\n")

    # Prepare comparison tasks
    tasks = []
    for gen_file in gen_files:
        leg_file = leg_dir / gen_file.name
        if leg_file.exists():
            tasks.append((gen_file, leg_file))

    print(f"Found {len(tasks)} file pairs to compare\n")

    # Compare in parallel using threads
    results = []
    num_workers = 24
    lock = threading.Lock()
    completed_count = [0]

    def update_progress():
        with lock:
            completed_count[0] += 1
            if completed_count[0] % 100 == 0:
                print(
                    f"Progress: {completed_count[0]}/{len(tasks)} files compared...",
                    flush=True,
                )

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(compare_pdb_atoms, gen, leg): (gen, leg)
            for gen, leg in tasks
        }

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            update_progress()

    # Analyze results
    passed = sum(1 for _, success, _ in results if success)
    failed = len(results) - passed

    print(f"\n{'='*70}")
    print(f"Comparison Results (pdb_atoms records only)")
    print(f"{'='*70}")
    print(f"Total files: {len(results)}")
    print(f"✓ Passed (exact match): {passed}")
    print(f"✗ Failed: {failed}")
    print(f"Success rate: {100*passed/len(results):.1f}%")

    if failed > 0:
        print(f"\nFailed files (showing first 30):")
        failed_results = [
            (name, reason) for name, success, reason in results if not success
        ]
        for name, reason in failed_results[:30]:
            print(f"  ✗ {name}: {reason}")
        if len(failed_results) > 30:
            print(f"  ... and {len(failed_results) - 30} more")

    print(f"{'='*70}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
