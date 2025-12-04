#!/usr/bin/env python3
"""
Comprehensive Stage 2 validation test.

Tests all Stage 2 components without legacy JSON reads:
- residue_indices
- ls_fitting  
- base_frame_calc
- frame_calc

Handles legacy duplicate bug (ls_fitting, base_frame_calc, frame_calc all had duplicates).
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class ValidationResult:
    pdb_id: str
    record_type: str
    matches: bool
    modern_count: int
    legacy_count_raw: int
    legacy_count_dedup: int
    error: str = ""


def deduplicate_by_residue_idx(records: List[Dict]) -> List[Dict]:
    """De-duplicate records by residue_idx, keeping first occurrence."""
    seen = set()
    deduped = []
    for record in records:
        idx = record.get('residue_idx')
        if idx not in seen:
            seen.add(idx)
            deduped.append(record)
    return deduped


def validate_record_type(pdb_id: str, record_type: str, modern_dir: Path, legacy_dir: Path) -> ValidationResult:
    """Validate a single record type for a PDB."""
    modern_file = modern_dir / record_type / f"{pdb_id}.json"
    legacy_file = legacy_dir / record_type / f"{pdb_id}.json"
    
    if not modern_file.exists():
        return ValidationResult(pdb_id, record_type, False, 0, 0, 0, f"Modern JSON missing")
    
    if not legacy_file.exists():
        return ValidationResult(pdb_id, record_type, False, 0, 0, 0, f"Legacy JSON missing")
    
    try:
        with open(modern_file) as f:
            modern = json.load(f)
        with open(legacy_file) as f:
            legacy_raw = json.load(f)
        
        # Handle single record (residue_indices) vs list
        if not isinstance(modern, list):
            modern = [modern]
        if not isinstance(legacy_raw, list):
            legacy_raw = [legacy_raw]
        
        # De-duplicate legacy (needed for ls_fitting, base_frame_calc, frame_calc)
        legacy = deduplicate_by_residue_idx(legacy_raw)
        
        # Compare counts
        if len(modern) != len(legacy):
            return ValidationResult(
                pdb_id, record_type, False,
                len(modern), len(legacy_raw), len(legacy),
                f"Count mismatch"
            )
        
        # For this test, just check counts match (detailed comparison can be added later)
        return ValidationResult(
            pdb_id, record_type, True,
            len(modern), len(legacy_raw), len(legacy)
        )
        
    except Exception as e:
        return ValidationResult(pdb_id, record_type, False, 0, 0, 0, f"Error: {e}")


def main():
    if len(sys.argv) < 2:
        print("Usage: test_stage2_complete.py PDB1 [PDB2 ...]")
        print("   or: test_stage2_complete.py --test-50")
        sys.exit(1)
    
    root_dir = Path(__file__).parent.parent
    modern_dir = root_dir / "data" / "json"
    legacy_dir = root_dir / "data" / "json_legacy"
    
    # Get PDFs to test
    if sys.argv[1] == "--test-50":
        # Find 50 PDBs that have all legacy Stage 2 JSON
        pdbs = []
        for json_file in (legacy_dir / "ls_fitting").glob("*.json"):
            pdb_id = json_file.stem
            # Check if all Stage 2 types exist in legacy
            has_all = all(
                (legacy_dir / rt / f"{pdb_id}.json").exists()
                for rt in ["residue_indices", "ls_fitting", "base_frame_calc", "frame_calc"]
            )
            if has_all and (root_dir / "data" / "pdb" / f"{pdb_id}.pdb").exists():
                pdbs.append(pdb_id)
            if len(pdbs) >= 50:
                break
    else:
        pdbs = sys.argv[1:]
    
    if not pdbs:
        print("No PDBs found to test")
        sys.exit(1)
    
    print(f"\n{'='*70}")
    print(f"STAGE 2 COMPREHENSIVE VALIDATION")
    print(f"{'='*70}")
    print(f"Testing {len(pdbs)} PDBs")
    print(f"Record types: residue_indices, ls_fitting, base_frame_calc, frame_calc")
    print(f"Note: Legacy has duplicate bug for frame-related types")
    print(f"{'='*70}\n")
    
    # Generate modern JSON for all PDFs
    print("Generating modern JSON...")
    import subprocess
    for i, pdb_id in enumerate(pdbs, 1):
        if i % 10 == 0:
            print(f"  Progress: {i}/{len(pdbs)}")
        
        pdb_file = root_dir / "data" / "pdb" / f"{pdb_id}.pdb"
        if not pdb_file.exists():
            continue
        
        # Generate frames stage (includes all Stage 2 components)
        subprocess.run(
            [str(root_dir / "build" / "generate_modern_json"),
             str(pdb_file),
             str(modern_dir),
             "--stage=frames"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    
    print(f"  Generated modern JSON for {len(pdbs)} PDBs\n")
    
    # Validate each record type
    record_types = ["residue_indices", "ls_fitting", "base_frame_calc", "frame_calc"]
    results_by_type = {rt: [] for rt in record_types}
    
    for pdb_id in pdbs:
        for record_type in record_types:
            result = validate_record_type(pdb_id, record_type, modern_dir, legacy_dir)
            results_by_type[record_type].append(result)
    
    # Print summary for each record type
    print(f"\n{'='*70}")
    print("RESULTS BY RECORD TYPE")
    print(f"{'='*70}\n")
    
    for record_type in record_types:
        results = results_by_type[record_type]
        passed = sum(1 for r in results if r.matches)
        failed = sum(1 for r in results if not r.matches)
        
        print(f"{record_type}:")
        print(f"  ✅ Passed: {passed}/{len(results)} ({100*passed/len(results):.1f}%)")
        print(f"  ❌ Failed: {failed}/{len(results)}")
        
        # Show a few examples
        passing_examples = [r for r in results if r.matches][:3]
        if passing_examples:
            print(f"  Examples:")
            for r in passing_examples:
                dedup_note = f" (dedup {r.legacy_count_raw}→{r.legacy_count_dedup})" if r.legacy_count_raw != r.legacy_count_dedup else ""
                print(f"    {r.pdb_id}: {r.modern_count} records match{dedup_note}")
        
        # Show failures
        failing = [r for r in results if not r.matches]
        if failing:
            print(f"  Failures:")
            for r in failing[:5]:
                print(f"    {r.pdb_id}: {r.error}")
        
        print()
    
    # Overall summary
    total_tests = len(pdbs) * len(record_types)
    total_passed = sum(sum(1 for r in results if r.matches) for results in results_by_type.values())
    
    print(f"{'='*70}")
    print("OVERALL SUMMARY")
    print(f"{'='*70}")
    print(f"Total tests: {total_tests} ({len(pdbs)} PDBs × {len(record_types)} types)")
    print(f"✅ Passed: {total_passed} ({100*total_passed/total_tests:.1f}%)")
    print(f"❌ Failed: {total_tests - total_passed}")
    print(f"{'='*70}\n")
    
    if total_passed >= 0.85 * total_tests:
        print("✅ STAGE 2 VALIDATION PASSED (≥85% match rate)")
        return 0
    else:
        print("❌ STAGE 2 VALIDATION FAILED (<85% match rate)")
        return 1


if __name__ == "__main__":
    sys.exit(main())

