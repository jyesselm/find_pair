#!/usr/bin/env python3
"""
Test ls_fitting JSON generation and comparison.

Tests Stage 3: ls_fitting JSON generation after refactoring.
Uses x3dna_json_compare module for all comparison logic.
"""

import json
from pathlib import Path
from typing import Dict, Tuple, Optional
import sys

# Add project root to path to import comparison module
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
from x3dna_json_compare.frame_comparison import compare_frames
from scripts.test_utils import (
    find_executables,
    generate_legacy_json,
    generate_modern_json,
    cleanup_output_dir
)

def compare_ls_fitting_json(legacy_file: Path, modern_file: Path,
                           pdb_file: Path, project_root: Path) -> Tuple[bool, Dict]:
    """Compare legacy and modern ls_fitting JSON files using x3dna_json_compare."""
    try:
        # Load JSON files
        if not legacy_file.exists():
            return False, {"error": "Legacy file not found"}
        if not modern_file.exists():
            return False, {"error": "Modern file not found"}
        
        with open(legacy_file) as f:
            legacy_data = json.load(f)
        with open(modern_file) as f:
            modern_data = json.load(f)
        
        # Normalize to lists
        legacy_records = legacy_data if isinstance(legacy_data, list) else [legacy_data]
        modern_records = modern_data if isinstance(modern_data, list) else [modern_data]
        
        # Add type field to records for comparison function
        legacy_ls_records = []
        for r in legacy_records:
            if isinstance(r, dict):
                r_copy = r.copy()
                r_copy["type"] = "ls_fitting"
                legacy_ls_records.append(r_copy)
        
        modern_ls_records = []
        for r in modern_records:
            if isinstance(r, dict):
                r_copy = r.copy()
                r_copy["type"] = "ls_fitting"
                modern_ls_records.append(r_copy)
        
        # Load legacy atoms for atom_idx lookup
        legacy_atoms = []
        legacy_atoms_file = project_root / "data" / "json_legacy" / "pdb_atoms" / f"{pdb_file.stem}.json"
        if legacy_atoms_file.exists():
            try:
                with open(legacy_atoms_file) as f:
                    atoms_data = json.load(f)
                    if isinstance(atoms_data, list) and len(atoms_data) > 0:
                        if isinstance(atoms_data[0], dict) and "pdb_atoms" in atoms_data[0]:
                            legacy_atoms = atoms_data[0]["pdb_atoms"]
                        else:
                            legacy_atoms = atoms_data
                    elif isinstance(atoms_data, dict) and "pdb_atoms" in atoms_data:
                        legacy_atoms = atoms_data["pdb_atoms"]
            except Exception:
                pass  # Ignore errors loading atoms
        
        # Use x3dna_json_compare for comparison
        compare_result = compare_frames(
            legacy_ls_records, modern_ls_records, pdb_file, None, legacy_atoms
        )
        
        # Handle FrameComparison object - access attributes directly
        if hasattr(compare_result, "total_legacy"):
            total_legacy = compare_result.total_legacy
            total_modern = compare_result.total_modern
            missing_count = len(compare_result.missing_residues)
            mismatch_count = len(compare_result.mismatched_calculations)
            matched = total_legacy - missing_count - mismatch_count
            
            if total_legacy == 0 and total_modern == 0:
                return True, {"message": "No records to compare (both empty)"}
            
            if total_legacy != total_modern:
                return False, {
                    "error": f"Record count mismatch: legacy={total_legacy}, modern={total_modern}"
                }
            
            if missing_count > 0 or mismatch_count > 0:
                return False, {
                    "error": f"Not all records matched: {matched}/{total_legacy} (missing={missing_count}, mismatches={mismatch_count})",
                    "missing_residues": [str(r) for r in compare_result.missing_residues[:5]],
                    "mismatches": [str(m) for m in compare_result.mismatched_calculations[:5]]
                }
            
            return True, {
                "matched": matched,
                "total": total_legacy,
                "details": {
                    "total_legacy": total_legacy,
                    "total_modern": total_modern,
                    "matched": matched
                }
            }
        elif isinstance(compare_result, dict):
            # Handle dict result
            if compare_result.get("error"):
                return False, compare_result
            return True, compare_result
        else:
            return False, {
                "error": f"Unknown comparison result type: {type(compare_result)}"
            }
        
    except Exception as e:
        return False, {"error": f"Comparison error: {str(e)}"}

def test_single_pdb(pdb_id: str, pdb_file: Path, output_dir: Path,
                    legacy_exe: Optional[Path], modern_exe: Path,
                    project_root: Path) -> Dict:
    """Test ls_fitting generation for a single PDB."""
    result = {
        "pdb_id": pdb_id,
        "legacy_ok": False,
        "modern_ok": False,
        "compare_ok": False,
        "errors": []
    }
    
    # Generate legacy
    if legacy_exe:
        legacy_ok, legacy_msg = generate_legacy_json(
            pdb_id, pdb_file, output_dir, legacy_exe, project_root,
            record_types=["ls_fitting"]
        )
        result["legacy_ok"] = legacy_ok
        if not legacy_ok:
            result["errors"].append(f"Legacy: {legacy_msg}")
    else:
        result["errors"].append("Legacy executable not found")
    
    # Generate modern
    modern_ok, modern_msg = generate_modern_json(
        pdb_id, pdb_file, output_dir, modern_exe, stage="ls_fitting"
    )
    result["modern_ok"] = modern_ok
    if not modern_ok:
        result["errors"].append(f"Modern: {modern_msg}")
    
    # Compare if both generated successfully
    if legacy_ok and modern_ok:
        legacy_file = output_dir / "legacy" / "ls_fitting" / f"{pdb_id}.json"
        modern_file = output_dir / "modern" / "ls_fitting" / f"{pdb_id}.json"
        
        if legacy_file.exists() and modern_file.exists():
            compare_ok, compare_result = compare_ls_fitting_json(
                legacy_file, modern_file, pdb_file, project_root
            )
            result["compare_ok"] = compare_ok
            if compare_ok:
                result["matched"] = compare_result.get("matched", 0)
                result["total"] = compare_result.get("total", 0)
            else:
                if "error" in compare_result:
                    result["errors"].append(f"Compare: {compare_result['error']}")
                else:
                    result["errors"].append(f"Compare: {compare_result}")
        else:
            result["errors"].append("JSON files not found after generation")
    
    return result

def main():
    project_root = Path(".")
    
    # Find executables
    legacy_exe, modern_exe = find_executables(project_root)
    
    if not modern_exe:
        print("Error: Modern executable not found")
        print(f"Expected: {project_root / 'build' / 'generate_modern_json'}")
        sys.exit(1)
    
    # Test with a single PDB first
    import sys
    if len(sys.argv) > 1:
        pdb_id = sys.argv[1].upper()
    else:
        # Default to a small PDB for quick testing
        pdb_id = "1ABC"
    
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        print(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)
    
    output_dir = project_root / "tmp" / "test_ls_fitting"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Testing ls_fitting generation for {pdb_id}...")
    print(f"PDB file: {pdb_file}")
    print(f"Output dir: {output_dir}")
    print()
    
    try:
        result = test_single_pdb(
            pdb_id, pdb_file, output_dir, legacy_exe, modern_exe, project_root
        )
        
        # Print results
        print(f"Results for {pdb_id}:")
        print(f"  Legacy generation: {'✓' if result['legacy_ok'] else '✗'}")
        print(f"  Modern generation: {'✓' if result['modern_ok'] else '✗'}")
        print(f"  Comparison: {'✓' if result['compare_ok'] else '✗'}")
        
        if result.get("matched") is not None:
            print(f"  Matched: {result['matched']}/{result.get('total', 0)} records")
        
        if result["errors"]:
            print("\nErrors:")
            for error in result["errors"]:
                print(f"  - {error}")
        
        success = result["legacy_ok"] and result["modern_ok"] and result["compare_ok"]
        
        if success:
            print("\n✓ All tests passed!")
            return_code = 0
        else:
            print("\n✗ Some tests failed")
            return_code = 1
    finally:
        # Clean up output directory
        cleanup_output_dir(output_dir)
        print(f"\nCleaned up output directory: {output_dir}")
    
    return return_code

if __name__ == "__main__":
    sys.exit(main())

