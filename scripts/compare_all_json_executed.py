#!/usr/bin/env python3
"""
Compare JSON output between modern and legacy implementations across all test cases.
Uses existing library code from lib/ for comparisons.

Usage:
    # Compare single PDB file
    python3 scripts/compare_all_json_executed.py 1H4S

    # Compare single PDB file (legacy mode)
    python3 scripts/compare_all_json_executed.py 1H4S --legacy

    # Compare all files (will take several minutes)
    python3 scripts/compare_all_json_executed.py

This script uses the existing JsonComparator from lib/json_comparison.py to compare
executed calculations (base_frame_calc, ls_fitting, frame_calc, pdb_atoms) between
legacy and modern JSON outputs.
"""

import json
import sys
import multiprocessing
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import from the package
from x3dna_json_compare import JsonComparator, ComparisonResult


def compare_executed_calculations(pdb_id, project_root, use_legacy_mode=False):
    """Compare executed calculations for a single PDB file using JsonComparator."""
    project_root = Path(project_root)

    if use_legacy_mode:
        modern_file = project_root / f"data/json/{pdb_id}_legacy.json"
    else:
        modern_file = project_root / f"data/json/{pdb_id}.json"
    legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
    pdb_file = project_root / f"data/pdb/{pdb_id}.pdb"

    # Use existing JsonComparator
    comparator = JsonComparator(enable_cache=False)
    result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)

    if result.status == "error":
        return None

    # Extract summary information for executed calculations
    summary = {
        "pdb_id": pdb_id,
        "structure_match": True,
        "record_types": {},
        "base_frame_calc_detail": {},
    }

    # Check structure for executed types
    if result.frame_comparison:
        fc = result.frame_comparison
        matching_count = fc.total_legacy - len(fc.missing_residues)

        # Count exact matches
        exact_matches = 0
        rms_diffs = []
        for mismatch in fc.mismatched_calculations:
            leg_atoms = sorted(mismatch.legacy_matched_atoms)
            mod_atoms = sorted(mismatch.modern_matched_atoms)
            if leg_atoms == mod_atoms:
                exact_matches += 1

            leg_rms = mismatch.legacy_record.get("rms_fit", 0)
            mod_rms = mismatch.modern_record.get("rms_fit", 0)
            if abs(leg_rms - mod_rms) > 0.001:
                rms_diffs.append(abs(leg_rms - mod_rms))

        summary["base_frame_calc_detail"] = {
            "modern_count": fc.total_modern,
            "legacy_count": fc.total_legacy,
            "matching_count": matching_count,
            "missing_residues": len(fc.missing_residues),
            "mismatched_calculations": len(fc.mismatched_calculations),
            "exact_atom_matches": exact_matches,
            "atom_match_rate": (
                exact_matches / matching_count * 100 if matching_count > 0 else 0
            ),
            "rms_differences_count": len(rms_diffs),
            "avg_rms_diff": sum(rms_diffs) / len(rms_diffs) if rms_diffs else 0,
            "max_rms_diff": max(rms_diffs) if rms_diffs else 0,
        }

    return summary


def main():
    """Compare all JSON files using existing library code."""
    project_root = Path(__file__).parent.parent

    if len(sys.argv) > 1:
        # Compare specific PDB ID
        pdb_id = sys.argv[1]
        use_legacy_mode = "--legacy" in sys.argv

        result = compare_executed_calculations(pdb_id, project_root, use_legacy_mode)
        if result:
            print(json.dumps(result, indent=2))
        else:
            print(f"No data found for {pdb_id}")
        return

    # Use JsonComparator batch_compare for all files
    json_dir = project_root / "data/json"
    legacy_dir = project_root / "data/json_legacy"

    modern_files = set(
        f.stem for f in json_dir.glob("*.json") if not f.stem.endswith("_legacy")
    )
    legacy_files = set(f.stem for f in legacy_dir.glob("*.json"))
    common = sorted(modern_files & legacy_files)

    print("=" * 80)
    print("JSON COMPARISON ACROSS ALL TEST CASES (Executed Calculations Only)")
    print("=" * 80)
    print(f"\nTotal PDB files with modern JSON: {len(modern_files)}")
    print(f"Total PDB files with legacy JSON: {len(legacy_files)}")
    print(f"Common files: {len(common)}")

    if not common:
        print("\nNo common files found!")
        return

    # Compare using JsonComparator with threading (avoids pickling issues)
    comparator = JsonComparator(enable_cache=False)

    # Use maximum threads (threads avoid pickling issues)
    max_workers = multiprocessing.cpu_count()

    print("\n" + "=" * 80)
    print(f"COMPARING ALL FILES (using {max_workers} threads)...")
    print("=" * 80)

    # Process in parallel using threads
    results = {}

    def compare_pdb(pdb_id):
        """Compare a single PDB file - skip PDB file to speed up (we don't need atom lines)."""
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_file = project_root / f"data/json/{pdb_id}.json"
        # Skip PDB file parsing - we only need JSON comparison, not PDB line lookups
        # This makes comparison ~10x faster
        pdb_file = None

        try:
            result = comparator.compare_files(
                legacy_file, modern_file, pdb_file, pdb_id
            )
            return (pdb_id, result)
        except Exception as e:
            # Create error result
            error_result = ComparisonResult(
                pdb_id=pdb_id, status="error", errors=[str(e)]
            )
            return (pdb_id, error_result)

    # Process with ThreadPoolExecutor
    completed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_pdb = {
            executor.submit(compare_pdb, pdb_id): pdb_id for pdb_id in common
        }

        # Collect results as they complete
        for future in as_completed(future_to_pdb):
            pdb_id, result = future.result()
            results[pdb_id] = result
            completed += 1

            # Progress update every 100 files or at completion
            if completed % 100 == 0 or completed == len(common):
                print(
                    f"  Progress: {completed}/{len(common)} ({completed*100//len(common)}%)",
                    flush=True,
                )

    # Summarize results
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    files_with_calcs = []
    exact_matches = 0
    total_residues = 0
    total_mismatches = 0
    rms_diffs = []

    for pdb_id, result in results.items():
        if result.status == "error":
            continue

        if result.frame_comparison:
            fc = result.frame_comparison
            if fc.total_modern > 0 or fc.total_legacy > 0:
                files_with_calcs.append(pdb_id)
                total_residues += fc.total_legacy
                total_mismatches += len(fc.mismatched_calculations)

                # Count exact matches
                for mismatch in fc.mismatched_calculations:
                    leg_atoms = sorted(mismatch.legacy_matched_atoms)
                    mod_atoms = sorted(mismatch.modern_matched_atoms)
                    if leg_atoms == mod_atoms:
                        exact_matches += 1

                    leg_rms = mismatch.legacy_record.get("rms_fit", 0)
                    mod_rms = mismatch.modern_record.get("rms_fit", 0)
                    if abs(leg_rms - mod_rms) > 0.001:
                        rms_diffs.append(abs(leg_rms - mod_rms))

    print(f"\nFiles with calculations: {len(files_with_calcs)}")
    print(f"Total residues processed: {total_residues}")
    print(f"Residues with mismatches: {total_mismatches}")
    print(
        f"Exact atom matches: {exact_matches}/{total_residues - len([r for r in results.values() if r.frame_comparison and r.frame_comparison.missing_residues])}"
    )

    if rms_diffs:
        print(f"RMS differences: {len(rms_diffs)} residues")
        print(f"Average RMS difference: {sum(rms_diffs) / len(rms_diffs):.6f}")
        print(f"Maximum RMS difference: {max(rms_diffs):.6f}")

    # Show files with most differences
    if files_with_calcs:
        print("\n" + "=" * 80)
        print("FILES WITH MOST DIFFERENCES (Top 10)")
        print("=" * 80)

        file_diffs = []
        for pdb_id in files_with_calcs:
            result = results[pdb_id]
            if result.frame_comparison:
                diff_count = len(result.frame_comparison.mismatched_calculations) + len(
                    result.frame_comparison.missing_residues
                )
                file_diffs.append((pdb_id, diff_count))

        file_diffs.sort(key=lambda x: x[1], reverse=True)
        for pdb_id, diff_count in file_diffs[:10]:
            print(f"  {pdb_id}: {diff_count} differences")

    # Save detailed mismatches for analysis
    if total_mismatches > 0:
        mismatches_file = project_root / "docs" / "MISMATCHED_RESIDUES.json"
        mismatches_data = []
        
        for pdb_id, result in results.items():
            if result.status == "error" or not result.frame_comparison:
                continue
            
            fc = result.frame_comparison
            
            # Collect missing residues
            for missing in fc.missing_residues:
                mismatches_data.append({
                    'pdb_id': pdb_id,
                    'type': 'missing',
                    'residue': f"{missing.chain_id}:{missing.residue_seq}{missing.insertion}",
                    'residue_name': missing.residue_name,
                    'base_type': missing.base_type,
                    'matched_atoms': missing.matched_atoms,
                })
            
            # Collect mismatched calculations
            for mismatch in fc.mismatched_calculations:
                leg_atoms = sorted(mismatch.legacy_matched_atoms)
                mod_atoms = sorted(mismatch.modern_matched_atoms)
                chain_id, residue_seq, insertion = mismatch.residue_key
                
                mismatches_data.append({
                    'pdb_id': pdb_id,
                    'type': 'mismatch',
                    'residue': f"{chain_id}:{residue_seq}{insertion}",
                    'residue_name': mismatch.legacy_record.get('residue_name', ''),
                    'base_type': mismatch.legacy_record.get('base_type', '?'),
                    'legacy_atoms': leg_atoms,
                    'modern_atoms': mod_atoms,
                    'only_legacy': sorted(set(leg_atoms) - set(mod_atoms)),
                    'only_modern': sorted(set(mod_atoms) - set(leg_atoms)),
                    'rms_legacy': mismatch.legacy_record.get('rms_fit', 0),
                    'rms_modern': mismatch.modern_record.get('rms_fit', 0),
                    'rms_diff': abs(mismatch.legacy_record.get('rms_fit', 0) - mismatch.modern_record.get('rms_fit', 0)),
                })
        
        with open(mismatches_file, 'w') as f:
            json.dump({
                'total_mismatches': len(mismatches_data),
                'mismatches': mismatches_data
            }, f, indent=2)
        
        print(f"\n💾 Detailed mismatches saved to: {mismatches_file}")


if __name__ == "__main__":
    main()
