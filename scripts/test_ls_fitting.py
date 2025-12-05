#!/usr/bin/env python3
"""
Test ls_fitting JSON generation and comparison.

Tests Stage 3: ls_fitting JSON generation after refactoring.
"""

import subprocess
import json
from pathlib import Path
from typing import Dict, Tuple, Optional
import sys

# Add parent directory to path to import comparison module
sys.path.insert(0, str(Path(__file__).parent.parent))
from x3dna_json_compare.frame_comparison import compare_frames

def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    modern_exe = project_root / "build" / "generate_modern_json"
    
    if not legacy_exe.exists():
        legacy_exe = None
    if not modern_exe.exists():
        modern_exe = None
    
    return legacy_exe, modern_exe

def generate_legacy_ls_fitting(pdb_id: str, pdb_file: Path, output_dir: Path,
                                legacy_exe: Path, project_root: Path) -> Tuple[bool, str]:
    """Generate legacy ls_fitting JSON - check existing files first."""
    try:
        legacy_output = output_dir / "legacy"
        legacy_output.mkdir(parents=True, exist_ok=True)
        
        # Check if legacy ls_fitting JSON already exists
        ls_fitting_file = project_root / "data" / "json_legacy" / "ls_fitting" / f"{pdb_id}.json"
        
        if ls_fitting_file.exists():
            # Copy existing file
            (legacy_output / "ls_fitting").mkdir(exist_ok=True)
            import shutil
            shutil.copy2(ls_fitting_file, legacy_output / "ls_fitting" / f"{pdb_id}.json")
            return True, "OK (existing file)"
        
        # File doesn't exist - need to generate it
        # Legacy generates all frames together, so we need base_frame_calc and frame_calc too
        # But we'll just check if ls_fitting exists after running legacy
        if pdb_file.is_absolute():
            try:
                pdb_file_rel = pdb_file.relative_to(project_root)
            except ValueError:
                pdb_file_rel = pdb_file
        else:
            pdb_file_rel = pdb_file
        
        # Build path for org/ directory
        if str(pdb_file_rel).startswith("data/"):
            pdb_path_for_org = "../" + str(pdb_file_rel)
        else:
            pdb_path_for_org = str(pdb_file_rel)
        
        cmd = [str(legacy_exe.resolve()), pdb_path_for_org]
        
        # Run with timeout
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=300,  # 5 minute timeout
                cwd=str(project_root / "org")
            )
            
            if result.returncode != 0:
                return False, f"Legacy generation failed (return code {result.returncode})"
        except subprocess.TimeoutExpired:
            return False, "Legacy generation timeout (exceeded 5 minutes)"
        
        # Check if file was created
        if ls_fitting_file.exists():
            # Copy to output directory
            (legacy_output / "ls_fitting").mkdir(exist_ok=True)
            import shutil
            shutil.copy2(ls_fitting_file, legacy_output / "ls_fitting" / f"{pdb_id}.json")
            return True, "OK"
        else:
            return False, "Legacy ls_fitting JSON file not created"
        
    except Exception as e:
        return False, f"Legacy error: {str(e)}"

def generate_modern_ls_fitting(pdb_id: str, pdb_file: Path, output_dir: Path,
                               modern_exe: Path) -> Tuple[bool, str]:
    """Generate modern ls_fitting JSON."""
    try:
        modern_output = output_dir / "modern"
        modern_output.mkdir(parents=True, exist_ok=True)
        
        cmd = [str(modern_exe), str(pdb_file), str(modern_output), "--stage=ls_fitting"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout
        )
        
        if result.returncode != 0:
            return False, f"Modern generation failed: {result.stderr[:200]}"
        
        # Check if modern JSON file was created
        ls_fitting_file = modern_output / "ls_fitting" / f"{pdb_id}.json"
        
        if not ls_fitting_file.exists():
            return False, "Modern ls_fitting JSON file not created"
        
        return True, "OK"
        
    except subprocess.TimeoutExpired:
        return False, "Modern generation timeout"
    except Exception as e:
        return False, f"Modern error: {str(e)}"

def compare_ls_fitting_json(legacy_file: Path, modern_file: Path,
                           pdb_file: Path, project_root: Path) -> Tuple[bool, Dict]:
    """Compare legacy and modern ls_fitting JSON files."""
    try:
        # Load legacy JSON
        if not legacy_file.exists():
            return False, {"error": "Legacy file not found"}
        
        with open(legacy_file) as f:
            legacy_data = json.load(f)
        
        # Load modern JSON
        if not modern_file.exists():
            return False, {"error": "Modern file not found"}
        
        with open(modern_file) as f:
            modern_data = json.load(f)
        
        # Normalize to lists
        legacy_records = legacy_data if isinstance(legacy_data, list) else [legacy_data]
        modern_records = modern_data if isinstance(modern_data, list) else [modern_data]
        
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
        
        # Add type field to records for comparison function
        # Legacy records don't have type, modern records might not either
        legacy_ls_records = []
        for r in legacy_records:
            if isinstance(r, dict):
                r_copy = r.copy()
                r_copy["type"] = "ls_fitting"  # Add type field
                legacy_ls_records.append(r_copy)
        
        modern_ls_records = []
        for r in modern_records:
            if isinstance(r, dict):
                r_copy = r.copy()
                r_copy["type"] = "ls_fitting"  # Add type field
                modern_ls_records.append(r_copy)
        
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
        legacy_ok, legacy_msg = generate_legacy_ls_fitting(
            pdb_id, pdb_file, output_dir, legacy_exe, project_root
        )
        result["legacy_ok"] = legacy_ok
        if not legacy_ok:
            result["errors"].append(f"Legacy: {legacy_msg}")
    else:
        result["errors"].append("Legacy executable not found")
    
    # Generate modern
    modern_ok, modern_msg = generate_modern_ls_fitting(
        pdb_id, pdb_file, output_dir, modern_exe
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
        import shutil
        if output_dir.exists():
            shutil.rmtree(output_dir)
            print(f"\nCleaned up output directory: {output_dir}")
    
    return return_code

if __name__ == "__main__":
    sys.exit(main())

