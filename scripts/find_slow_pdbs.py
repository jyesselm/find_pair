#!/usr/bin/env python3
"""
Identify PDBs that take more than 30 seconds to generate legacy output.

These PDBs can be excluded from batch testing to avoid timeouts and slowdowns.
"""

import subprocess
import json
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional

def _save_progress(output_file: Path, slow_pdbs: List[Dict], fast_pdbs: List[Dict],
                  failed_pdbs: List[Dict], total: int, start_time: float, timeout: int):
    """Save progress incrementally."""
    elapsed_total = time.time() - start_time
    results = {
        "test_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "timeout_threshold_seconds": timeout,
        "total_tested": len(slow_pdbs) + len(fast_pdbs) + len(failed_pdbs),
        "total_available": total,
        "slow_pdbs_count": len(slow_pdbs),
        "fast_pdbs_count": len(fast_pdbs),
        "failed_count": len(failed_pdbs),
        "total_time_seconds": round(elapsed_total, 2),
        "in_progress": True,
        "slow_pdbs": sorted(slow_pdbs, key=lambda x: x["elapsed_seconds"], reverse=True),
        "fast_pdbs_sample": sorted(fast_pdbs, key=lambda x: x["elapsed_seconds"], reverse=True)[:10]
    }
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

def test_legacy_generation_time(pdb_id: str, pdb_file: str, legacy_exe: str, 
                                project_root: str, timeout_seconds: int = 30) -> Tuple[str, bool, float, str]:
    """Test how long legacy generation takes for a single PDB.
    
    Returns:
        (pdb_id, is_slow, elapsed_time, status_message)
    """
    pdb_file_path = Path(pdb_file)
    legacy_exe_path = Path(legacy_exe)
    project_root_path = Path(project_root)
    
    # Convert PDB path to relative path for org/ directory
    if pdb_file_path.is_absolute():
        try:
            pdb_file_rel = pdb_file_path.relative_to(project_root_path)
        except ValueError:
            pdb_file_rel = pdb_file_path
    else:
        pdb_file_rel = pdb_file_path
    
    # Build path for org/ directory
    if str(pdb_file_rel).startswith("data/"):
        pdb_path_for_org = "../" + str(pdb_file_rel)
    else:
        pdb_path_for_org = str(pdb_file_rel)
    
    cmd = [str(legacy_exe_path.resolve()), pdb_path_for_org]
    
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout_seconds,
            cwd=str(project_root_path / "org")
        )
        elapsed = time.time() - start_time
        
        if result.returncode == 0:
            if elapsed > timeout_seconds:
                return (pdb_id, True, elapsed, f"Completed but took {elapsed:.1f}s (> {timeout_seconds}s)")
            else:
                return (pdb_id, False, elapsed, f"OK ({elapsed:.1f}s)")
        else:
            return (pdb_id, True, elapsed, f"Failed (return code {result.returncode})")
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        return (pdb_id, True, elapsed, f"Timeout after {elapsed:.1f}s")
    except Exception as e:
        elapsed = time.time() - start_time
        return (pdb_id, True, elapsed, f"Error: {str(e)}")

def main():
    project_root = Path(".")
    valid_pdbs_file = project_root / "data" / "valid_pdbs.json"
    output_file = project_root / "data" / "slow_pdbs.json"
    
    # Find legacy executable
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    if not legacy_exe.exists():
        print(f"❌ Legacy executable not found: {legacy_exe}")
        return 1
    
    # Load valid PDBs
    if not valid_pdbs_file.exists():
        print(f"❌ Valid PDBs file not found: {valid_pdbs_file}")
        return 1
    
    with open(valid_pdbs_file) as f:
        valid_data = json.load(f)
    
    # Handle different possible structures of valid_pdbs.json
    if isinstance(valid_data, list):
        valid_pdb_ids = valid_data
    elif isinstance(valid_data, dict):
        # Try different possible keys
        valid_pdb_ids = (valid_data.get("valid_pdbs", []) or 
                        valid_data.get("valid_pdbs_with_atoms_and_frames", []) or
                        valid_data.get("valid_pdbs_atoms_only", []))
    else:
        print("❌ Invalid format in valid_pdbs.json")
        return 1
    
    if not valid_pdb_ids:
        print("❌ No valid PDBs found in file")
        return 1
    
    print(f"Testing {len(valid_pdb_ids)} PDBs for generation time...")
    print(f"Timeout threshold: 30 seconds")
    print(f"Legacy executable: {legacy_exe}")
    print()
    
    # Get PDB files
    pdb_dir = project_root / "data" / "pdb"
    pdb_files = []
    missing = []
    
    for pdb_id in valid_pdb_ids:
        pdb_file = pdb_dir / f"{pdb_id}.pdb"
        if pdb_file.exists():
            pdb_files.append((pdb_id, str(pdb_file)))
        else:
            missing.append(pdb_id)
    
    if missing:
        print(f"⚠️  {len(missing)} PDBs not found in {pdb_dir}")
    
    print(f"Testing {len(pdb_files)} PDBs...")
    print()
    
    # Test in parallel (use fewer workers to avoid overwhelming system)
    num_workers = min(5, len(pdb_files))
    slow_pdbs = []
    fast_pdbs = []
    failed_pdbs = []
    
    start_time = time.time()
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                test_legacy_generation_time,
                pdb_id,
                pdb_file,
                str(legacy_exe),
                str(project_root),
                15  # 30 second timeout
            ): pdb_id
            for pdb_id, pdb_file in pdb_files
        }
        
        completed = 0
        for future in as_completed(futures):
            pdb_id = futures[future]
            try:
                result_pdb_id, is_slow, elapsed, status = future.result()
                
                if is_slow:
                    slow_pdbs.append({
                        "pdb_id": result_pdb_id,
                        "elapsed_seconds": round(elapsed, 2),
                        "status": status
                    })
                    print(f"  ⏱️  {result_pdb_id}: {status}")
                else:
                    fast_pdbs.append({
                        "pdb_id": result_pdb_id,
                        "elapsed_seconds": round(elapsed, 2)
                    })
                    if completed % 50 == 0:  # Print every 50th fast one
                        print(f"  ✅ {result_pdb_id}: {status}")
                
                completed += 1
                
                # Save progress every 50 PDBs to avoid losing work if interrupted
                if completed % 50 == 0:
                    _save_progress(output_file, slow_pdbs, fast_pdbs, failed_pdbs, 
                                  len(pdb_files), start_time, 30)
                
                if completed % 10 == 0:  # Print progress more frequently
                    elapsed_so_far = time.time() - start_time
                    rate = completed / elapsed_so_far if elapsed_so_far > 0 else 0
                    remaining = (len(pdb_files) - completed) / rate if rate > 0 else 0
                    print(f"  Progress: {completed}/{len(pdb_files)} ({completed/len(pdb_files)*100:.1f}%) | "
                          f"Slow: {len(slow_pdbs)} | "
                          f"ETA: {remaining/60:.1f} min")
            except Exception as e:
                failed_pdbs.append({
                    "pdb_id": pdb_id,
                    "error": str(e)
                })
                print(f"  ❌ {pdb_id}: Exception - {e}")
                completed += 1
    
    elapsed_total = time.time() - start_time
    
    # Save final results
    results = {
        "test_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "timeout_threshold_seconds": 30,
        "total_tested": len(pdb_files),
        "total_available": len(pdb_files),
        "slow_pdbs_count": len(slow_pdbs),
        "fast_pdbs_count": len(fast_pdbs),
        "failed_count": len(failed_pdbs),
        "total_time_seconds": round(elapsed_total, 2),
        "in_progress": False,
        "slow_pdbs": sorted(slow_pdbs, key=lambda x: x["elapsed_seconds"], reverse=True),
        "fast_pdbs_sample": sorted(fast_pdbs, key=lambda x: x["elapsed_seconds"], reverse=True)[:10]  # Top 10 fastest
    }
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print()
    print("=" * 60)
    print("Results:")
    print(f"  Total tested: {len(pdb_files)}")
    print(f"  ⏱️  Slow (>30s): {len(slow_pdbs)}")
    print(f"  ✅ Fast (≤30s): {len(fast_pdbs)}")
    print(f"  ❌ Failed: {len(failed_pdbs)}")
    print(f"  Total time: {elapsed_total:.1f} seconds")
    print(f"  Results saved to: {output_file}")
    print()
    
    if slow_pdbs:
        print(f"Slow PDBs (will be excluded from batch testing):")
        for pdb in slow_pdbs[:20]:  # Show first 20
            print(f"  - {pdb['pdb_id']}: {pdb['status']}")
        if len(slow_pdbs) > 20:
            print(f"  ... and {len(slow_pdbs) - 20} more")
    
    return 0

if __name__ == "__main__":
    exit(main())

