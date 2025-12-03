#!/usr/bin/env python3
"""
Test frames JSON generation and comparison in batches with parallel processing.

Generates both legacy and modern frames JSON (base_frame_calc, frame_calc, ls_fitting),
compares them, and cleans up after each batch.
"""

import subprocess
import json
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional
import time
import sys

# Add parent directory to path to import comparison module
sys.path.insert(0, str(Path(__file__).parent.parent))
from x3dna_json_compare.frame_comparison import compare_frames
from x3dna_json_compare.pdb_utils import PdbFileReader

def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    modern_exe = project_root / "build" / "generate_modern_json"
    
    if not legacy_exe.exists():
        legacy_exe = None
    if not modern_exe.exists():
        modern_exe = None
    
    return legacy_exe, modern_exe

def generate_legacy_frames(pdb_id: str, pdb_file: Path, output_dir: Path, 
                           legacy_exe: Path, project_root: Path) -> Tuple[bool, str]:
    """Generate legacy frames JSON - check existing files first, regenerate if needed."""
    try:
        legacy_output = output_dir / "legacy"
        legacy_output.mkdir(parents=True, exist_ok=True)
        
        # Check if legacy frames JSON already exist
        base_frame_file = project_root / "data" / "json_legacy" / "base_frame_calc" / f"{pdb_id}.json"
        frame_calc_file = project_root / "data" / "json_legacy" / "frame_calc" / f"{pdb_id}.json"
        ls_fitting_file = project_root / "data" / "json_legacy" / "ls_fitting" / f"{pdb_id}.json"
        
        # Check for ref_frame as alternative to frame_calc
        if not frame_calc_file.exists():
            ref_frame_file = project_root / "data" / "json_legacy" / "ref_frame" / f"{pdb_id}.json"
            if ref_frame_file.exists():
                frame_calc_file = ref_frame_file
        
        all_exist = base_frame_file.exists() and frame_calc_file.exists() and ls_fitting_file.exists()
        
        if all_exist:
            # Copy existing files
            (legacy_output / "base_frame_calc").mkdir(exist_ok=True)
            (legacy_output / "frame_calc").mkdir(exist_ok=True)
            (legacy_output / "ls_fitting").mkdir(exist_ok=True)
            shutil.copy2(base_frame_file, legacy_output / "base_frame_calc" / f"{pdb_id}.json")
            shutil.copy2(frame_calc_file, legacy_output / "frame_calc" / f"{pdb_id}.json")
            shutil.copy2(ls_fitting_file, legacy_output / "ls_fitting" / f"{pdb_id}.json")
            return True, "OK (existing files)"
        
        # Files don't exist - regenerate using legacy executable
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
        
        # Check if files were created
        if base_frame_file.exists() and frame_calc_file.exists() and ls_fitting_file.exists():
            # Copy to output directory
            (legacy_output / "base_frame_calc").mkdir(exist_ok=True)
            (legacy_output / "frame_calc").mkdir(exist_ok=True)
            (legacy_output / "ls_fitting").mkdir(exist_ok=True)
            shutil.copy2(base_frame_file, legacy_output / "base_frame_calc" / f"{pdb_id}.json")
            shutil.copy2(frame_calc_file, legacy_output / "frame_calc" / f"{pdb_id}.json")
            shutil.copy2(ls_fitting_file, legacy_output / "ls_fitting" / f"{pdb_id}.json")
            return True, "OK (regenerated)"
        else:
            missing = []
            if not base_frame_file.exists():
                missing.append("base_frame_calc")
            if not frame_calc_file.exists():
                missing.append("frame_calc")
            if not ls_fitting_file.exists():
                missing.append("ls_fitting")
            return False, f"Legacy generation completed but files missing: {', '.join(missing)}"
        
    except Exception as e:
        return False, f"Legacy error: {str(e)}"

def generate_modern_frames(pdb_id: str, pdb_file: Path, output_dir: Path,
                          modern_exe: Path) -> Tuple[bool, str]:
    """Generate modern frames JSON."""
    try:
        modern_output = output_dir / "modern"
        modern_output.mkdir(parents=True, exist_ok=True)
        
        cmd = [str(modern_exe), str(pdb_file), str(modern_output), "--stage=frames"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout for frames
        )
        
        if result.returncode != 0:
            return False, f"Modern generation failed: {result.stderr[:200]}"
        
        # Check if modern JSON files were created
        base_frame_file = modern_output / "base_frame_calc" / f"{pdb_id}.json"
        frame_calc_file = modern_output / "frame_calc" / f"{pdb_id}.json"
        ls_fitting_file = modern_output / "ls_fitting" / f"{pdb_id}.json"
        
        if not base_frame_file.exists():
            return False, "Modern base_frame_calc JSON file not created"
        if not frame_calc_file.exists():
            return False, "Modern frame_calc JSON file not created"
        if not ls_fitting_file.exists():
            return False, "Modern ls_fitting JSON file not created"
        
        return True, "OK"
        
    except subprocess.TimeoutExpired:
        return False, "Modern generation timeout"
    except Exception as e:
        return False, f"Modern error: {str(e)}"

def compare_frames_json(legacy_base_frame: Path, legacy_frame_calc: Path, legacy_ls_fitting: Path,
                       modern_base_frame: Path, modern_frame_calc: Path, modern_ls_fitting: Path,
                       pdb_file: Path, project_root: Path) -> Tuple[bool, Dict]:
    """Compare legacy and modern frames JSON files."""
    try:
        # Load legacy JSON files
        legacy_records = []
        for file_path in [legacy_base_frame, legacy_frame_calc, legacy_ls_fitting]:
            if not file_path.exists():
                continue
            with open(file_path) as f:
                data = json.load(f)
                if isinstance(data, list):
                    legacy_records.extend(data)
                elif isinstance(data, dict):
                    if "type" in data:
                        legacy_records.append(data)
                    elif "calculations" in data:
                        legacy_records.extend(data.get("calculations", []))
                    else:
                        legacy_records.append(data)
        
        # Load modern JSON files
        modern_records = []
        for file_path in [modern_base_frame, modern_frame_calc, modern_ls_fitting]:
            if not file_path.exists():
                continue
            with open(file_path) as f:
                data = json.load(f)
                if isinstance(data, list):
                    modern_records.extend(data)
                elif isinstance(data, dict):
                    if "type" in data:
                        modern_records.append(data)
                    elif "calculations" in data:
                        modern_records.extend(data.get("calculations", []))
                    else:
                        modern_records.append(data)
        
        # Load legacy atoms for atom_idx lookup (optional but helpful)
        legacy_atoms = []
        legacy_atoms_file = project_root / "data" / "json_legacy" / "pdb_atoms" / f"{pdb_file.stem}.json"
        if legacy_atoms_file.exists():
            try:
                with open(legacy_atoms_file) as f:
                    atoms_data = json.load(f)
                    if isinstance(atoms_data, list) and len(atoms_data) > 0:
                        # Legacy atoms JSON is a list with a single pdb_atoms record
                        if isinstance(atoms_data[0], dict) and "atoms" in atoms_data[0]:
                            legacy_atoms = atoms_data[0]["atoms"]
                    elif isinstance(atoms_data, dict) and "atoms" in atoms_data:
                        legacy_atoms = atoms_data["atoms"]
            except Exception:
                pass  # If we can't load atoms, comparison will still work
        
        # Use the comparison module
        pdb_reader = PdbFileReader(pdb_file) if pdb_file.exists() else None
        comparison_result = compare_frames(legacy_records, modern_records, pdb_file, pdb_reader, legacy_atoms)
        
        if (comparison_result.missing_residues or 
            comparison_result.mismatched_calculations):
            errors = []
            if comparison_result.missing_residues:
                errors.append(f"Missing {len(comparison_result.missing_residues)} residues in modern")
            if comparison_result.mismatched_calculations:
                errors.append(f"{len(comparison_result.mismatched_calculations)} calculation mismatches")
            
            return False, {
                "errors": errors,
                "missing_residues": list(comparison_result.missing_residues)[:10],
                "mismatched_calculations": comparison_result.mismatched_calculations[:5]
            }
        
        return True, {
            "matched": comparison_result.total_legacy,
            "total_legacy": comparison_result.total_legacy,
            "total_modern": comparison_result.total_modern
        }
        
    except Exception as e:
        return False, {"error": f"Comparison error: {str(e)}"}

def process_pdb(pdb_id: str, pdb_file: str, output_dir: str,
                legacy_exe: Optional[str], modern_exe: str,
                project_root: str, verbose: bool = False) -> Dict:
    """Process a single PDB: generate both JSONs and compare.
    
    Note: Arguments are strings (not Path objects) for pickling compatibility with ProcessPoolExecutor.
    """
    # Convert string arguments back to Path objects for internal use
    pdb_file = Path(pdb_file)
    output_dir = Path(output_dir)
    legacy_exe = Path(legacy_exe) if legacy_exe else None
    modern_exe = Path(modern_exe)
    project_root = Path(project_root)
    
    result = {
        "pdb_id": pdb_id,
        "legacy_ok": False,
        "modern_ok": False,
        "compare_ok": False,
        "errors": []
    }
    
    # Generate legacy JSON
    if legacy_exe:
        legacy_ok, legacy_msg = generate_legacy_frames(
            pdb_id, pdb_file, output_dir, legacy_exe, project_root
        )
        result["legacy_ok"] = legacy_ok
        if not legacy_ok:
            result["errors"].append(f"Legacy: {legacy_msg}")
    else:
        result["errors"].append("Legacy executable not found")
        return result
    
    # Generate modern JSON
    modern_ok, modern_msg = generate_modern_frames(
        pdb_id, pdb_file, output_dir, modern_exe
    )
    result["modern_ok"] = modern_ok
    if not modern_ok:
        result["errors"].append(f"Modern: {modern_msg}")
        return result
    
    # Compare if both generated successfully
    if legacy_ok and modern_ok:
        legacy_base_frame = output_dir / "legacy" / "base_frame_calc" / f"{pdb_id}.json"
        legacy_frame_calc = output_dir / "legacy" / "frame_calc" / f"{pdb_id}.json"
        legacy_ls_fitting = output_dir / "legacy" / "ls_fitting" / f"{pdb_id}.json"
        
        modern_base_frame = output_dir / "modern" / "base_frame_calc" / f"{pdb_id}.json"
        modern_frame_calc = output_dir / "modern" / "frame_calc" / f"{pdb_id}.json"
        modern_ls_fitting = output_dir / "modern" / "ls_fitting" / f"{pdb_id}.json"
        
        if (legacy_base_frame.exists() and legacy_frame_calc.exists() and legacy_ls_fitting.exists() and
            modern_base_frame.exists() and modern_frame_calc.exists() and modern_ls_fitting.exists()):
            compare_ok, compare_result = compare_frames_json(
                legacy_base_frame, legacy_frame_calc, legacy_ls_fitting,
                modern_base_frame, modern_frame_calc, modern_ls_fitting,
                pdb_file, project_root
            )
            result["compare_ok"] = compare_ok
            if compare_ok:
                result["frames_matched"] = compare_result.get("total_legacy", 0)
            else:
                # Format compare errors nicely
                if "error" in compare_result:
                    result["errors"].append(f"Compare: {compare_result['error']}")
                elif "errors" in compare_result:
                    result["errors"].extend([f"Compare: {e}" for e in compare_result["errors"]])
                else:
                    result["errors"].append(f"Compare: {compare_result}")
        else:
            result["errors"].append("JSON files not found after generation")
    
    return result

def main():
    project_root = Path(".")
    pdb_dir = project_root / "data" / "pdb"
    output_base = Path("/tmp/test_frames_batch")
    results_file = project_root / "data" / "frames_test_results.json"
    
    # Find executables
    legacy_exe, modern_exe = find_executables(project_root)
    
    if not modern_exe:
        print("❌ Modern executable not found!")
        return 1
    
    if not legacy_exe:
        print("⚠️  Legacy executable not found, will only test modern generation")
    
    # Load valid PDBs list (prefer fast PDBs if available)
    fast_pdbs_file = project_root / "data" / "valid_pdbs_fast.json"
    valid_pdbs_file = project_root / "data" / "valid_pdbs.json"
    
    if fast_pdbs_file.exists():
        with open(fast_pdbs_file) as f:
            pdbs_data = json.load(f)
            pdb_list = pdbs_data.get("valid_pdbs_with_atoms_and_frames", [])
        print(f"📋 Using {fast_pdbs_file.name} ({len(pdb_list)} PDBs)")
    elif valid_pdbs_file.exists():
        with open(valid_pdbs_file) as f:
            pdbs_data = json.load(f)
            pdb_list = pdbs_data.get("valid_pdbs_with_atoms_and_frames", [])
        print(f"📋 Using {valid_pdbs_file.name} ({len(pdb_list)} PDBs)")
    else:
        print("❌ No valid PDBs file found!")
        return 1
    
    # Load existing results if any
    existing_results = {}
    if results_file.exists():
        try:
            with open(results_file) as f:
                existing_data = json.load(f)
                existing_results = {r["pdb_id"]: r for r in existing_data.get("results", {}).get("all", [])}
            print(f"📂 Loaded {len(existing_results)} existing results")
        except Exception as e:
            print(f"⚠️  Could not load existing results: {e}")
    
    # Filter out already completed PDBs
    remaining_pdbs = [pdb_id for pdb_id in pdb_list 
                     if pdb_id not in existing_results or 
                     not existing_results[pdb_id].get("compare_ok", False)]
    
    if not remaining_pdbs:
        print("✅ All PDBs already tested!")
        return 0
    
    print(f"🔄 Testing {len(remaining_pdbs)} PDBs (skipping {len(pdb_list) - len(remaining_pdbs)} already completed)")
    
    # Batch processing settings
    batch_size = 50
    num_threads = 10
    
    # Use 10 processes for frames testing
    effective_workers = min(num_threads, 10) if legacy_exe else num_threads
    
    print(f"⚙️  Batch size: {batch_size}, Workers: {effective_workers}")
    print()
    
    all_results = list(existing_results.values())
    batch_num = 0
    
    for i in range(0, len(remaining_pdbs), batch_size):
        batch_num += 1
        batch = remaining_pdbs[i:i+batch_size]
        batch_output = output_base / f"batch_{batch_num}"
        batch_output.mkdir(parents=True, exist_ok=True)
        
        print(f"📦 Batch {batch_num}: Processing {len(batch)} PDBs...")
        start_time = time.time()
        
        batch_results = []
        
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            futures = {}
            for pdb_id in batch:
                pdb_file = pdb_dir / f"{pdb_id}.pdb"
                if not pdb_file.exists():
                    batch_results.append({
                        "pdb_id": pdb_id,
                        "legacy_ok": False,
                        "modern_ok": False,
                        "compare_ok": False,
                        "errors": [f"PDB file not found: {pdb_file}"]
                    })
                    continue
                
                future = executor.submit(
                    process_pdb,
                    pdb_id,
                    str(pdb_file),
                    str(batch_output),
                    str(legacy_exe) if legacy_exe else None,
                    str(modern_exe),
                    str(project_root)
                )
                futures[future] = pdb_id
            
            for future in as_completed(futures):
                pdb_id = futures[future]
                try:
                    result = future.result()
                    batch_results.append(result)
                    
                    # Status: ✅ if modern OK and (legacy OK or not needed) and compare OK
                    if result["modern_ok"]:
                        if result.get("legacy_ok") and result.get("compare_ok"):
                            status = "✅"
                        elif not result.get("legacy_ok") and not legacy_exe:
                            status = "✅"  # No legacy exe, modern OK is success
                        elif not result.get("legacy_ok"):
                            status = "⚠️ "  # Modern OK but legacy missing
                        else:
                            status = "❌"
                    else:
                        status = "❌"
                    print(f"  {status} {pdb_id}")
                    if result["errors"]:
                        for err in result["errors"][:2]:  # Show first 2 errors
                            print(f"     {err}")
                except Exception as e:
                    error_msg = f"Exception: {str(e)}"
                    batch_results.append({
                        "pdb_id": pdb_id,
                        "legacy_ok": False,
                        "modern_ok": False,
                        "compare_ok": False,
                        "errors": [error_msg]
                    })
                    print(f"  ❌ {pdb_id}: {error_msg}")
        
        elapsed = time.time() - start_time
        print(f"  ⏱️  Batch {batch_num} completed in {elapsed:.1f}s ({elapsed/len(batch):.2f}s/PDB)")
        
        # Categorize results
        passed = [r for r in batch_results if r.get("modern_ok") and r.get("legacy_ok") and r.get("compare_ok")]
        partial = [r for r in batch_results if r.get("modern_ok") and (not r.get("legacy_ok") or not r.get("compare_ok"))]
        failed = [r for r in batch_results if not r.get("modern_ok")]
        
        print(f"  📊 Passed: {len(passed)}, Partial: {len(partial)}, Failed: {len(failed)}")
        print()
        
        # Add to all results
        all_results.extend(batch_results)
        
        # Save results after each batch
        results_data = {
            "total_pdbs": len(all_results),
            "passed": len([r for r in all_results if r.get("modern_ok") and r.get("legacy_ok") and r.get("compare_ok")]),
            "partial": len([r for r in all_results if r.get("modern_ok") and (not r.get("legacy_ok") or not r.get("compare_ok"))]),
            "failed": len([r for r in all_results if not r.get("modern_ok")]),
            "results": {
                "all": all_results,
                "passed": [r for r in all_results if r.get("modern_ok") and r.get("legacy_ok") and r.get("compare_ok")],
                "partial": [r for r in all_results if r.get("modern_ok") and (not r.get("legacy_ok") or not r.get("compare_ok"))],
                "failed": [r for r in all_results if not r.get("modern_ok")]
            }
        }
        
        with open(results_file, 'w') as f:
            json.dump(results_data, f, indent=2)
        
        # Clean up batch output directory
        if batch_output.exists():
            shutil.rmtree(batch_output)
        
        # Check for critical errors that should stop processing
        critical_errors = [r for r in batch_results 
                          if "executable not found" in str(r.get("errors", [])) or
                          "timeout" in str(r.get("errors", [])).lower()]
        
        if critical_errors and len(critical_errors) == len(batch_results):
            print("⚠️  Critical error detected in entire batch - stopping for inspection")
            break
    
    # Final summary
    print()
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Total tested: {len(all_results)}")
    print(f"✅ Passed: {len([r for r in all_results if r.get('modern_ok') and r.get('legacy_ok') and r.get('compare_ok')])}")
    print(f"⚠️  Partial: {len([r for r in all_results if r.get('modern_ok') and (not r.get('legacy_ok') or not r.get('compare_ok'))])}")
    print(f"❌ Failed: {len([r for r in all_results if not r.get('modern_ok')])}")
    print()
    print(f"Results saved to: {results_file}")
    print("=" * 70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

