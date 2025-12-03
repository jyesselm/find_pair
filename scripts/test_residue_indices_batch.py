#!/usr/bin/env python3
"""
Test residue_indices JSON generation and comparison in batches with parallel processing.

Generates both legacy and modern residue_indices JSON, compares them, and cleans up after each batch.
Residue indices map residues to atom ranges (seidx array).
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
from x3dna_json_compare.residue_indices_comparison import (
    compare_residue_indices,
    parse_seidx_array
)

def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    modern_exe = project_root / "build" / "generate_modern_json"
    
    if not legacy_exe.exists():
        legacy_exe = None
    if not modern_exe.exists():
        modern_exe = None
    
    return legacy_exe, modern_exe

def generate_legacy_residue_indices(pdb_id: str, pdb_file: Path, output_dir: Path, 
                                   legacy_exe: Path, project_root: Path) -> Tuple[bool, str]:
    """Generate legacy residue_indices JSON - check existing file first, regenerate if needed."""
    try:
        legacy_output = output_dir / "legacy"
        legacy_output.mkdir(parents=True, exist_ok=True)
        
        # Check if legacy residue_indices JSON already exists
        legacy_json = project_root / "data" / "json_legacy" / "residue_indices" / f"{pdb_id}.json"
        if legacy_json.exists():
            # Copy existing file
            shutil.copy2(legacy_json, legacy_output / f"{pdb_id}.json")
            return True, "OK (existing file)"
        
        # File doesn't exist - regenerate it using legacy executable
        # Legacy executable runs from org/ directory and expects PDB file path relative to org/
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
        
        # Check if the file was created
        if legacy_json.exists():
            # Copy to output directory
            shutil.copy2(legacy_json, legacy_output / f"{pdb_id}.json")
            return True, "OK (regenerated)"
        else:
            return False, "Legacy generation completed but residue_indices JSON not found"
        
    except Exception as e:
        return False, f"Legacy error: {str(e)}"

def generate_modern_residue_indices(pdb_id: str, pdb_file: Path, output_dir: Path,
                                   modern_exe: Path) -> Tuple[bool, str]:
    """Generate modern residue_indices JSON."""
    try:
        modern_output = output_dir / "modern"
        modern_output.mkdir(parents=True, exist_ok=True)
        
        cmd = [str(modern_exe), str(pdb_file), str(modern_output), "--stage=residue_indices"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode != 0:
            return False, f"Modern generation failed: {result.stderr[:200]}"
        
        # Check if modern JSON was created
        modern_json = modern_output / "residue_indices" / f"{pdb_id}.json"
        if not modern_json.exists():
            return False, "Modern JSON file not created"
        
        return True, "OK"
        
    except subprocess.TimeoutExpired:
        return False, "Modern generation timeout"
    except Exception as e:
        return False, f"Modern error: {str(e)}"

def compare_residue_indices_json(legacy_file: Path, modern_file: Path) -> Tuple[bool, Dict]:
    """Compare legacy and modern residue_indices JSON files."""
    try:
        # Load legacy JSON
        with open(legacy_file) as f:
            legacy_data = json.load(f)
        
        # Handle legacy JSON format (may be list or dict)
        if isinstance(legacy_data, list):
            if len(legacy_data) > 0 and isinstance(legacy_data[0], dict):
                legacy_records = legacy_data
            else:
                return False, {"error": "Legacy JSON is array but first element is not a dict"}
        elif isinstance(legacy_data, dict):
            if "type" in legacy_data and legacy_data.get("type") == "residue_indices":
                legacy_records = [legacy_data]
            elif "calculations" in legacy_data:
                # Find residue_indices in calculations array
                legacy_records = []
                for record in legacy_data.get("calculations", []):
                    if isinstance(record, dict) and record.get("type") == "residue_indices":
                        legacy_records.append(record)
                        break
            else:
                legacy_records = [legacy_data]
        else:
            return False, {"error": "Legacy JSON format not recognized"}
        
        # Load modern JSON
        with open(modern_file) as f:
            modern_data = json.load(f)
        
        # Handle modern JSON format (may be list or dict)
        if isinstance(modern_data, list):
            if len(modern_data) > 0 and isinstance(modern_data[0], dict):
                modern_records = modern_data
            else:
                return False, {"error": "Modern JSON is array but first element is not a dict"}
        elif isinstance(modern_data, dict):
            if "type" in modern_data and modern_data.get("type") == "residue_indices":
                modern_records = [modern_data]
            elif "calculations" in modern_data:
                # Find residue_indices in calculations array
                modern_records = []
                for record in modern_data.get("calculations", []):
                    if isinstance(record, dict) and record.get("type") == "residue_indices":
                        modern_records.append(record)
                        break
            else:
                modern_records = [modern_data]
        else:
            return False, {"error": "Modern JSON format not recognized"}
        
        # Use the comparison module
        comparison_result = compare_residue_indices(legacy_records, modern_records)
        
        if (comparison_result.missing_in_modern or 
            comparison_result.extra_in_modern or 
            not comparison_result.num_residue_match or
            comparison_result.mismatched_entries):
            errors = []
            if comparison_result.missing_in_modern:
                errors.append("Missing in modern")
            if comparison_result.extra_in_modern:
                errors.append("Extra in modern")
            if not comparison_result.num_residue_match:
                errors.append(f"num_residue mismatch: legacy={comparison_result.legacy_num_residue}, modern={comparison_result.modern_num_residue}")
            if comparison_result.mismatched_entries:
                errors.append(f"{len(comparison_result.mismatched_entries)} seidx entry mismatches")
            
            return False, {
                "errors": errors,
                "mismatched_entries": comparison_result.mismatched_entries[:10],
                "total_mismatches": len(comparison_result.mismatched_entries)
            }
        
        return True, {
            "matched": comparison_result.legacy_num_residue,
            "num_residue": comparison_result.legacy_num_residue
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
        legacy_ok, legacy_msg = generate_legacy_residue_indices(
            pdb_id, pdb_file, output_dir, legacy_exe, project_root
        )
        result["legacy_ok"] = legacy_ok
        if not legacy_ok:
            result["errors"].append(f"Legacy: {legacy_msg}")
    else:
        result["errors"].append("Legacy executable not found")
        return result
    
    # Generate modern JSON
    modern_ok, modern_msg = generate_modern_residue_indices(
        pdb_id, pdb_file, output_dir, modern_exe
    )
    result["modern_ok"] = modern_ok
    if not modern_ok:
        result["errors"].append(f"Modern: {modern_msg}")
        return result
    
    # Compare if both generated successfully
    if legacy_ok and modern_ok:
        legacy_file = output_dir / "legacy" / f"{pdb_id}.json"
        modern_file = output_dir / "modern" / "residue_indices" / f"{pdb_id}.json"
        
        if legacy_file.exists() and modern_file.exists():
            compare_ok, compare_result = compare_residue_indices_json(legacy_file, modern_file)
            result["compare_ok"] = compare_ok
            if compare_ok:
                result["residues_matched"] = compare_result.get("num_residue", 0)
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
    output_base = Path("/tmp/test_residue_indices_batch")
    results_file = project_root / "data" / "residue_indices_test_results.json"
    
    # Find executables
    legacy_exe, modern_exe = find_executables(project_root)
    
    if not modern_exe:
        print("❌ Modern executable not found!")
        return 1
    
    if not legacy_exe:
        print("⚠️  Legacy executable not found, will only test modern generation")
    
    # Load valid PDBs list - MUST use valid_pdbs_fast.json
    fast_pdbs_file = project_root / "data" / "valid_pdbs_fast.json"
    
    if not fast_pdbs_file.exists():
        print(f"❌ {fast_pdbs_file} not found!")
        print("   This script requires valid_pdbs_fast.json for testing.")
        return 1
    
    print(f"  Using {fast_pdbs_file.name} (excludes slow PDBs)")
    
    valid_pdb_ids = set()
    try:
        with open(fast_pdbs_file) as f:
                valid_data = json.load(f)
                # Handle different possible formats
                if isinstance(valid_data, list):
                    valid_pdb_ids = set(valid_data)
                elif isinstance(valid_data, dict):
                    # Check for common keys
                    if "valid_pdbs_with_atoms_and_frames" in valid_data:
                        valid_pdb_ids = set(valid_data["valid_pdbs_with_atoms_and_frames"])
                    elif "valid_pdbs" in valid_data:
                        valid_pdb_ids = set(valid_data["valid_pdbs"])
                    else:
                        # Assume keys are PDB IDs
                        valid_pdb_ids = set(valid_data.keys())
        except Exception as e:
            print(f"⚠️  Could not load {fast_pdbs_file.name}: {e}")
            print("   Will test all PDB files")
    
    # Get all PDB files, filtered by valid_pdbs if available
    all_pdb_files = sorted(list(pdb_dir.glob("*.pdb")))
    if valid_pdb_ids:
        pdb_files = [f for f in all_pdb_files if f.stem in valid_pdb_ids]
        print(f"  Filtered to {len(pdb_files)} valid PDBs (from {len(all_pdb_files)} total)")
    else:
        pdb_files = all_pdb_files
    
    if not pdb_files:
        print("❌ No PDB files found!")
        return 1
    
    # Batch processing settings
    batch_size = 10
    num_threads = 10
    total_batches = (len(pdb_files) + batch_size - 1) // batch_size
    
    print("=" * 60)
    print("Testing residue_indices JSON generation and comparison")
    print(f"  PDBs: {len(pdb_files)}")
    if valid_pdb_ids:
        print(f"  (Filtered from {len(all_pdb_files)} total using {fast_pdbs_file.name})")
    print(f"  Batches: {total_batches} (batch size: {batch_size})")
    print(f"  Threads: {num_threads}")
    print(f"  Legacy exe: {'✅' if legacy_exe else '❌'} (using existing files)")
    print(f"  Modern exe: {'✅' if modern_exe else '❌'}")
    print()
    
    # Load previous results if they exist
    previous_results = {}
    if results_file.exists():
        try:
            with open(results_file) as f:
                previous_results = json.load(f)
        except:
            pass
    
    # Track critical errors that should stop processing
    critical_errors = []
    stop_on_error = True  # Stop if critical error detected
    
    start_time = time.time()
    all_results = []
    
    for batch_num in range(total_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(pdb_files))
        batch_files = pdb_files[start_idx:end_idx]
        
        print(f"Batch {batch_num + 1}/{total_batches} (PDBs {start_idx + 1}-{end_idx})...")
        
        # Create batch output directory
        batch_output = output_base / f"batch_{batch_num + 1}"
        batch_output.mkdir(parents=True, exist_ok=True)
        
        # Process batch in parallel using processes (not threads)
        batch_results = []
        
        # Use 10 processes for residue_indices testing (faster than atoms)
        effective_workers = min(num_threads, 10) if legacy_exe else num_threads
        
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            futures = {
                executor.submit(
                    process_pdb,
                    pdb_file.stem,
                    str(pdb_file),
                    str(batch_output / pdb_file.stem),
                    str(legacy_exe) if legacy_exe else None,
                    str(modern_exe),
                    str(project_root)
                ): pdb_file.stem
                for pdb_file in batch_files
            }
            
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
                    print(f"  ❌ {pdb_id}: {error_msg}")
                    batch_results.append({
                        "pdb_id": pdb_id,
                        "legacy_ok": False,
                        "modern_ok": False,
                        "compare_ok": False,
                        "errors": [error_msg]
                    })
                    # Check if this is a critical error
                    if any(keyword in str(e).lower() for keyword in ["not found", "permission", "cannot", "no such file", "executable"]):
                        critical_errors.append(f"{pdb_id}: {error_msg}")
        
        all_results.extend(batch_results)
        
        # Check if we're stopping due to mismatch - if so, preserve output
        has_mismatch = any(r.get("modern_ok") and r.get("legacy_ok") and not r.get("compare_ok") 
                          for r in batch_results)
        # Check for critical errors (excluding legacy generation failures - those are expected)
        has_critical_error = any(
            any(keyword in err.lower() for keyword in ["permission", "cannot", "executable", "timeout", "no such file or directory"])
            and not ("legacy" in err.lower() and (
                "not found" in err.lower() or 
                "missing" in err.lower() or
                "generation failed" in err.lower() or
                "not available" in err.lower()
            ))
            for r in batch_results for err in r.get("errors", [])
        )
        
        # Only clean up if we're continuing (no mismatches or critical errors)
        if not (has_mismatch or has_critical_error):
            try:
                if batch_output.exists():
                    shutil.rmtree(batch_output)
            except Exception as e:
                print(f"  ⚠️  Warning: Could not clean up batch output: {e}")
        
        # Report batch summary
        batch_success = sum(1 for r in batch_results 
                           if r["modern_ok"] and 
                           (not legacy_exe or r["legacy_ok"]) and
                           r["compare_ok"])
        batch_failed = len(batch_results) - batch_success
        
        print(f"  Batch summary: {batch_success}/{len(batch_results)} succeeded")
        
        # Check for critical errors and stop if needed
        if stop_on_error and critical_errors:
            print(f"\n❌ CRITICAL ERROR DETECTED - Stopping processing")
            print(f"   Critical errors:")
            for err in critical_errors[:5]:
                print(f"     - {err}")
            break
        
        # Check for comparison mismatches - stop to inspect
        comparison_mismatches = [r for r in batch_results 
                                if r.get("modern_ok") and r.get("legacy_ok") and not r.get("compare_ok")]
        if comparison_mismatches and stop_on_error:
            print(f"\n❌ COMPARISON MISMATCH DETECTED - Stopping to allow inspection")
            print(f"   Found {len(comparison_mismatches)} PDB(s) with comparison mismatches in this batch:")
            for r in comparison_mismatches:
                print(f"     - {r['pdb_id']}: {', '.join(r.get('errors', [])[:1])}")
            print(f"\n   Output files preserved for inspection:")
            for r in comparison_mismatches:
                legacy_file = batch_output / r["pdb_id"] / "legacy" / f"{r['pdb_id']}.json"
                modern_file = batch_output / r["pdb_id"] / "modern" / "residue_indices" / f"{r['pdb_id']}.json"
                if legacy_file.exists():
                    print(f"     Legacy: {legacy_file}")
                if modern_file.exists():
                    print(f"     Modern: {modern_file}")
            print(f"\n   Batch output directory: {batch_output}")
            print(f"   (Not cleaned up - inspect files above)")
            break
        
        # Check if all in batch failed with critical errors
        if batch_failed == len(batch_results) and len(batch_results) > 0:
            # Check if it's a consistent critical error pattern
            critical_error_types = {}
            for r in batch_results:
                for err in r.get("errors", []):
                    err_lower = err.lower()
                    # Exclude legacy-related failures
                    if "legacy" in err_lower and (
                        "not found" in err_lower or 
                        "missing" in err_lower or
                        "generation failed" in err_lower or
                        "not available" in err_lower
                    ):
                        continue
                    # Only consider actual critical errors
                    if any(keyword in err_lower for keyword in ["permission", "cannot", "executable", "timeout", "no such file or directory"]):
                        error_key = err.split(":")[0] if ":" in err else err[:50]
                        critical_error_types[error_key] = critical_error_types.get(error_key, 0) + 1
            
            # If all have the same critical error, might be a systemic issue
            if len(critical_error_types) == 1 and list(critical_error_types.values())[0] == len(batch_results):
                print(f"  ⚠️  All PDBs in batch failed with same critical error - might indicate systemic issue")
                print(f"     Error: {list(critical_error_types.keys())[0]}")
                if stop_on_error:
                    print(f"  Stopping processing due to consistent critical error")
                    print(f"  Batch output directory: {batch_output} (preserved for inspection)")
                    break
        
        # Save results after each batch (for progress tracking and crash recovery)
        _save_results(results_file, all_results, pdb_files, start_time, legacy_exe, 
                      previous_results, critical_errors)
        
        print()  # Empty line between batches
    
    # Final summary (results already saved after each batch)
    elapsed = time.time() - start_time
    
    # Re-categorize for final summary display
    passed = []
    failed = []
    partial = []
    
    for r in all_results:
        modern_ok = r["modern_ok"]
        legacy_ok = r.get("legacy_ok", False)
        compare_ok = r.get("compare_ok", False)
        
        is_full_success = (modern_ok and 
                          (not legacy_exe or legacy_ok) and
                          compare_ok)
        is_partial = (modern_ok and not legacy_ok and legacy_exe)
        
        if is_full_success:
            passed.append(r)
        elif is_partial:
            partial.append(r)
        else:
            failed.append(r)
    
    total_success = len(passed)
    modern_only_success = len(partial)
    
    print("=" * 60)
    print("Final Results:")
    print(f"  Total processed: {len(all_results)} PDBs")
    if len(all_results) < len(pdb_files):
        print(f"  ⚠️  Stopped early (processed {len(all_results)}/{len(pdb_files)})")
    print(f"  ✅ Full success (modern + legacy + compare): {total_success}")
    print(f"  ⚠️  Partial (modern OK, legacy missing): {modern_only_success}")
    print(f"  ❌ Failed: {len(all_results) - total_success - modern_only_success}")
    print(f"  Time: {elapsed:.1f} seconds ({elapsed/len(all_results):.2f} sec/PDB if all_results else 0)")
    print(f"  Results saved to: {results_file}")
    
    # Check if we stopped due to mismatch or critical error
    stopped_due_to_mismatch = any(r.get("modern_ok") and r.get("legacy_ok") and not r.get("compare_ok") 
                                  for r in all_results)
    stopped_due_to_critical = len(critical_errors) > 0
    
    # Final cleanup - only if no mismatches or critical errors
    if not (stopped_due_to_mismatch or stopped_due_to_critical):
        try:
            if output_base.exists():
                shutil.rmtree(output_base)
                print(f"  ✅ Cleaned up temporary output directory")
        except Exception as e:
            print(f"  ⚠️  Warning: Could not clean up temporary output: {e}")
    else:
        print(f"  ⚠️  Output preserved for inspection: {output_base}")
    
    return 0

def _save_results(results_file: Path, all_results: List[Dict], pdb_files: List[Path],
                 start_time: float, legacy_exe: Optional[Path], previous_results: Dict,
                 critical_errors: List[str]) -> None:
    """Save results to JSON file. Called after each batch to preserve progress."""
    elapsed = time.time() - start_time
    
    # Reload existing results file to get results from previous batches in this run
    existing_results = {}
    if results_file.exists():
        try:
            with open(results_file) as f:
                existing_results = json.load(f)
        except:
            pass
    
    # Create lookup of all current results (from this run so far)
    current_results_lookup = {r["pdb_id"]: r for r in all_results}
    
    # Merge: start with existing results, then update with current results
    merged_results = {
        "passed": [],
        "partial": [],
        "failed": [],
        "issues": []
    }
    
    # Add existing results (excluding ones we're updating)
    if existing_results and "results" in existing_results:
        for category in ["passed", "partial", "failed", "issues"]:
            for entry in existing_results["results"].get(category, []):
                pdb_id = entry.get("pdb_id")
                if pdb_id not in current_results_lookup:
                    merged_results[category].append(entry)
    
    # Add current results
    passed = []
    failed = []
    partial = []
    issues = []
    
    for r in all_results:
        modern_ok = r["modern_ok"]
        legacy_ok = r.get("legacy_ok", False)
        compare_ok = r.get("compare_ok", False)
        
        is_full_success = (modern_ok and 
                          (not legacy_exe or legacy_ok) and
                          compare_ok)
        is_partial = (modern_ok and not legacy_ok and legacy_exe)
        
        result_entry = {
            "pdb_id": r["pdb_id"],
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "modern_ok": modern_ok,
            "legacy_ok": legacy_ok,
            "compare_ok": compare_ok,
            "residues_matched": r.get("residues_matched", 0),
            "errors": r.get("errors", [])
        }
        
        if is_full_success:
            passed.append(result_entry)
        elif is_partial:
            partial.append(result_entry)
        else:
            failed.append(result_entry)
            if r.get("errors"):
                issues.append({
                    "pdb_id": r["pdb_id"],
                    "errors": r["errors"]
                })
    
    # Add current batch results to merged results
    merged_results["passed"].extend(passed)
    merged_results["partial"].extend(partial)
    merged_results["failed"].extend(failed)
    merged_results["issues"].extend(issues)
    
    # Calculate totals from merged results
    total_merged = (len(merged_results["passed"]) + 
                   len(merged_results["partial"]) + 
                   len(merged_results["failed"]))
    
    # Save merged results to JSON file
    results_data = {
        "test_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_pdbs": len(all_results),
        "total_available": len(pdb_files),
        "stopped_early": len(all_results) < len(pdb_files),
        "passed": len(merged_results["passed"]),
        "partial": len(merged_results["partial"]),
        "failed": len(merged_results["failed"]),
        "elapsed_seconds": round(elapsed, 2),
        "critical_errors": critical_errors,
        "summary": {
            "total": total_merged,
            "passed": len(merged_results["passed"]),
            "partial": len(merged_results["partial"]),
            "failed": len(merged_results["failed"]),
            "modern_ok": sum(1 for r in all_results if r["modern_ok"]),
            "legacy_ok": sum(1 for r in all_results if r.get("legacy_ok", False)),
            "compare_ok": sum(1 for r in all_results if r.get("compare_ok", False)),
            "success_rate": round(len(merged_results["passed"]) / total_merged * 100, 2) if total_merged else 0
        },
        "results": {
            "passed": sorted(merged_results["passed"], key=lambda x: x["pdb_id"]),
            "partial": sorted(merged_results["partial"], key=lambda x: x["pdb_id"]),
            "failed": sorted(merged_results["failed"], key=lambda x: x["pdb_id"]),
            "issues": sorted(merged_results["issues"], key=lambda x: x["pdb_id"])
        }
    }
    
    # Write results file
    results_file.parent.mkdir(parents=True, exist_ok=True)
    with open(results_file, 'w') as f:
        json.dump(results_data, f, indent=2)

if __name__ == "__main__":
    exit(main())

