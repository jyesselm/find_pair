#!/usr/bin/env python3
"""
Comprehensive index validation for all PDBs.

Tests ResidueTracker on all PDBs in batches, tracking results in CSV.
Stops on first failure for investigation.
Uses threading for parallel processing.

Features:
- Auto-resumes: skips PDBs that already passed, retests failures
- Cleans up: removes generated JSON files after successful validation (default)
- Parallel processing: uses multiple threads for speed

Usage:
    # Validate all untested PDBs (auto-resume)
    python3 scripts/validate_all_indices.py
    
    # Validate with custom settings
    python3 scripts/validate_all_indices.py --threads 8 --batch-size 100
    
    # Resume from a specific PDB
    python3 scripts/validate_all_indices.py --start-from 6V9Q
    
    # Revalidate all PDBs (including those that passed)
    python3 scripts/validate_all_indices.py --revalidate
    
    # Keep generated files (don't clean up)
    python3 scripts/validate_all_indices.py --no-clean
    
Options:
    --threads N         Number of parallel threads (default: 4)
    --batch-size N      Process N PDBs at a time (default: 50)
    --start-from PDB    Resume from this PDB ID
    --revalidate        Revalidate all PDBs, including those that already passed
    --clean             Remove generated files after validation (default: True)
    --no-clean          Keep generated files after validation
"""

import sys
import subprocess
import csv
import json
from pathlib import Path
from typing import Optional
import argparse
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Thread-safe lock for CSV writing
csv_lock = threading.Lock()

def get_valid_pdb_ids(project_root: Path) -> list[str]:
    """Get valid PDB IDs from valid_pdbs.json."""
    valid_pdbs_file = project_root / "data" / "valid_pdbs.json"
    
    if not valid_pdbs_file.exists():
        print("Warning: valid_pdbs.json not found, using all PDBs from data/pdb/")
        pdb_dir = project_root / "data" / "pdb"
        if not pdb_dir.exists():
            return []
        
        pdb_ids = []
        for pdb_file in pdb_dir.glob("*.pdb"):
            pdb_ids.append(pdb_file.stem)
        
        return sorted(pdb_ids)
    
    try:
        with open(valid_pdbs_file) as f:
            data = json.load(f)
            # Try different possible keys
            for key in ['valid_pdbs_with_atoms_and_frames', 'valid_pdbs', 'pdb_ids']:
                if key in data:
                    return data[key]
            return []
    except Exception as e:
        print(f"Error loading valid_pdbs.json: {e}")
        return []

def check_index_validation(pdb_id: str, project_root: Path, clean: bool = False) -> dict:
    """Run index validation on a single PDB and return results."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_dir = project_root / "data" / "json"
    legacy_json = project_root / "data" / "json_legacy" / "base_frame_calc" / f"{pdb_id}.json"
    
    result = {
        "pdb_id": pdb_id,
        "status": "UNKNOWN",
        "num_modern": 0,
        "num_legacy": 0,
        "num_matched": 0,
        "num_filtered": 0,
        "error": None,
        "timestamp": datetime.now().isoformat()
    }
    
    # Check if legacy JSON exists
    if not legacy_json.exists():
        result["status"] = "SKIP"
        result["error"] = "No legacy JSON found"
        return result
    
    # Check if PDB file exists
    if not pdb_file.exists():
        result["status"] = "SKIP"
        result["error"] = "PDB file not found"
        return result
    
    try:
        # Run generate_modern_json (this does validation internally)
        cmd = [
            str(project_root / "build" / "generate_modern_json"),
            str(pdb_file),
            str(json_dir)
        ]
        
        output = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout per PDB
        )
        
        # Check validation output
        if "Status: ❌ FAIL" in output.stdout or "Status: ✅ PASS" in output.stdout:
            if "Status: ❌ FAIL" in output.stdout:
                result["status"] = "FAIL"
                # Extract error details
                if "Count mismatch" in output.stdout:
                    result["error"] = "Count mismatch"
                elif "has no legacy index" in output.stdout:
                    result["error"] = "Missing legacy index"
                else:
                    result["error"] = "Validation failed (check output)"
                
                # Extract numbers from validation output
                for line in output.stdout.split('\n'):
                    if "Modern indices:" in line:
                        try:
                            result["num_modern"] = int(line.split(':')[1].strip())
                        except:
                            pass
                    elif "Legacy indices:" in line:
                        try:
                            result["num_legacy"] = int(line.split(':')[1].strip())
                        except:
                            pass
                    elif "Matched:" in line:
                        try:
                            result["num_matched"] = int(line.split(':')[1].strip())
                        except:
                            pass
                    elif "Filtered out:" in line:
                        try:
                            result["num_filtered"] = int(line.split(':')[1].strip())
                        except:
                            pass
                
                return result
            
            elif "Status: ✅ PASS" in output.stdout:
                result["status"] = "PASS"
                # Extract numbers
                for line in output.stdout.split('\n'):
                    if "Modern indices:" in line:
                        try:
                            result["num_modern"] = int(line.split(':')[1].strip())
                        except:
                            pass
                    elif "Legacy indices:" in line:
                        try:
                            result["num_legacy"] = int(line.split(':')[1].strip())
                        except:
                            pass
                    elif "Matched:" in line:
                        try:
                            result["num_matched"] = int(line.split(':')[1].strip())
                        except:
                            pass
                    elif "Filtered out:" in line:
                        try:
                            result["num_filtered"] = int(line.split(':')[1].strip())
                        except:
                            pass
                
                # AGGRESSIVE CLEANUP: Remove ALL generated files for this PDB if clean flag is set
                if clean:
                    json_dir = project_root / "data" / "json"
                    for subdir in json_dir.glob("*"):
                        if subdir.is_dir():
                            pdb_file = subdir / f"{pdb_id}.json"
                            if pdb_file.exists():
                                pdb_file.unlink()
                    
                    # Remove mapping file (only keep for failures)
                    mapping_file = project_root / "data" / "index_mapping" / f"{pdb_id}.json"
                    if mapping_file.exists():
                        mapping_file.unlink()
                
                return result
            else:
                result["status"] = "UNKNOWN"
                result["error"] = "Could not parse validation status"
                return result
        else:
            # Check exit code
            if output.returncode != 0:
                result["status"] = "ERROR"
                result["error"] = f"Process failed with exit code {output.returncode}"
                if output.stderr:
                    result["error"] += f": {output.stderr[:200]}"
            else:
                result["status"] = "SKIP"
                result["error"] = "No validation output found"
            
            return result
    
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
        result["error"] = "Process timed out (>120s)"
        return result
    except Exception as e:
        result["status"] = "ERROR"
        result["error"] = str(e)
        return result

def load_status_csv(csv_path: Path) -> dict[str, dict]:
    """Load existing status CSV."""
    if not csv_path.exists():
        return {}
    
    status = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert from CSV format to internal format
            pdb_id = row['pdb_id']
            status[pdb_id] = {
                'pdb_id': pdb_id,
                'status': row.get('match_status', 'UNKNOWN'),
                'num_modern': int(row.get('num_modern', 0)),
                'num_legacy': int(row.get('num_legacy', 0)),
                'num_matched': int(row.get('num_modern', 0)),  # Same as num_modern if passed
                'num_filtered': int(row.get('num_filtered', 0)),
                'error': row.get('notes', ''),
                'timestamp': ''
            }
    
    return status

def save_status_csv(csv_path: Path, status: dict[str, dict]):
    """Save status to CSV."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Match existing CSV format: pdb_id, num_residues_read, num_legacy, num_modern, num_filtered, match_status, notes
    fieldnames = ['pdb_id', 'num_residues_read', 'num_legacy', 'num_modern', 'num_filtered', 'match_status', 'notes']
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for pdb_id in sorted(status.keys()):
            row = status[pdb_id].copy()
            # Convert to existing CSV format
            output_row = {
                'pdb_id': row.get('pdb_id', pdb_id),
                'num_residues_read': row.get('num_modern', 0),
                'num_legacy': row.get('num_legacy', 0),
                'num_modern': row.get('num_modern', 0),
                'num_filtered': row.get('num_filtered', 0),
                'match_status': row.get('status', 'UNKNOWN'),
                'notes': row.get('error', '') if row.get('status') != 'PASS' else 'All indices match perfectly'
            }
            writer.writerow(output_row)

def process_batch_parallel(batch: list[str], project_root: Path, status: dict, 
                          status_csv: Path, clean: bool, threads: int) -> Optional[str]:
    """Process a batch of PDBs in parallel. Returns PDB ID of first failure or None."""
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all tasks
        future_to_pdb = {
            executor.submit(check_index_validation, pdb_id, project_root, clean): pdb_id
            for pdb_id in batch
        }
        
        # Process as they complete
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                result = future.result()
                
                # Thread-safe status update
                with csv_lock:
                    status[pdb_id] = result
                    save_status_csv(status_csv, status)
                
                # Print result
                if result["status"] == "PASS":
                    print(f"  ✅ {pdb_id}: PASS ({result['num_modern']} nucleotides)")
                elif result["status"] == "FAIL":
                    print(f"  ❌ {pdb_id}: FAIL - {result['error']}")
                    return pdb_id  # Return first failure
                elif result["status"] == "SKIP":
                    print(f"  ⏭️  {pdb_id}: SKIP - {result['error']}")
                else:
                    print(f"  ⚠️  {pdb_id}: {result['status']} - {result['error']}")
                    
            except Exception as e:
                print(f"  ❌ {pdb_id}: ERROR - {e}")
                with csv_lock:
                    status[pdb_id] = {
                        "pdb_id": pdb_id,
                        "status": "ERROR",
                        "error": str(e),
                        "num_modern": 0,
                        "num_legacy": 0,
                        "num_matched": 0,
                        "num_filtered": 0,
                        "timestamp": datetime.now().isoformat()
                    }
                    save_status_csv(status_csv, status)
    
    return None  # No failures

def main():
    parser = argparse.ArgumentParser(description='Validate indices for all PDBs')
    parser.add_argument('--threads', type=int, default=4,
                       help='Number of parallel threads (default: 4)')
    parser.add_argument('--batch-size', type=int, default=50, 
                       help='Process N PDBs at a time (default: 50)')
    parser.add_argument('--start-from', type=str, default=None,
                       help='Resume from this PDB ID')
    parser.add_argument('--clean', action='store_true', default=True,
                       help='Remove generated files after successful validation (default: True)')
    parser.add_argument('--no-clean', dest='clean', action='store_false',
                       help='Keep generated files after successful validation')
    parser.add_argument('--revalidate', action='store_true',
                       help='Revalidate all PDBs, including those that already passed')
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    status_csv = project_root / "data" / "index_validation_status.csv"
    
    # Load existing status
    status = load_status_csv(status_csv)
    
    # Get valid PDB IDs from valid_pdbs.json
    all_pdb_ids = get_valid_pdb_ids(project_root)
    
    if not all_pdb_ids:
        print("Error: No valid PDB IDs found")
        return 1
    
    print(f"Found {len(all_pdb_ids)} valid PDB IDs")
    
    # Filter to untested or start from specified PDB
    if args.start_from:
        start_idx = all_pdb_ids.index(args.start_from) if args.start_from in all_pdb_ids else 0
        pdb_ids_to_test = all_pdb_ids[start_idx:]
        print(f"Starting from {args.start_from} ({len(pdb_ids_to_test)} PDBs remaining)")
    elif args.revalidate:
        pdb_ids_to_test = all_pdb_ids
        print(f"Revalidating all {len(pdb_ids_to_test)} PDBs")
    else:
        # Auto-resume: skip only those that PASSED, retest FAIL/SKIP/ERROR/TIMEOUT
        pdb_ids_to_test = [pdb for pdb in all_pdb_ids 
                          if pdb not in status or status[pdb].get('status') != 'PASS']
        already_passed = len(all_pdb_ids) - len(pdb_ids_to_test)
        print(f"Auto-resuming: {len(pdb_ids_to_test)} PDBs to test ({already_passed} already passed)")
        if already_passed > 0:
            print(f"  (Use --revalidate to retest all PDBs)")
    
    if not pdb_ids_to_test:
        print("All PDBs already tested!")
        return 0
    
    # Process in batches
    total = len(pdb_ids_to_test)
    batch_size = args.batch_size
    
    cleanup_msg = "cleaning up" if args.clean else "keeping"
    print(f"\nProcessing in batches of {batch_size} with {args.threads} threads ({cleanup_msg} generated files)...\n")
    
    for batch_start in range(0, total, batch_size):
        batch_end = min(batch_start + batch_size, total)
        batch = pdb_ids_to_test[batch_start:batch_end]
        
        print(f"Batch {batch_start//batch_size + 1}: Testing PDBs {batch_start+1}-{batch_end} of {total}")
        
        # Process batch in parallel
        failed_pdb = process_batch_parallel(batch, project_root, status, status_csv, 
                                           args.clean, args.threads)
        
        if failed_pdb:
            print(f"\n⚠️  STOPPING: Found failure at {failed_pdb}")
            print(f"   Modern: {status[failed_pdb]['num_modern']}, Legacy: {status[failed_pdb]['num_legacy']}")
            print(f"   Error: {status[failed_pdb]['error']}")
            print(f"   See {status_csv} for details")
            print(f"   To continue after fixing: --start-from {failed_pdb}")
            return 1
        
        # Print batch summary
        batch_statuses = [status[pdb]['status'] for pdb in batch if pdb in status]
        pass_count = batch_statuses.count('PASS')
        fail_count = batch_statuses.count('FAIL')
        skip_count = batch_statuses.count('SKIP')
        
        print(f"\n  Batch summary: {pass_count} PASS, {fail_count} FAIL, {skip_count} SKIP\n")
    
    # Final summary
    all_statuses = [s['status'] for s in status.values()]
    total_pass = all_statuses.count('PASS')
    total_fail = all_statuses.count('FAIL')
    total_skip = all_statuses.count('SKIP')
    total_error = all_statuses.count('ERROR')
    total_timeout = all_statuses.count('TIMEOUT')
    
    print(f"\n{'='*60}")
    print(f"VALIDATION COMPLETE")
    print(f"{'='*60}")
    print(f"Total tested: {len(status)}")
    print(f"✅ PASS: {total_pass}")
    print(f"❌ FAIL: {total_fail}")
    print(f"⏭️  SKIP: {total_skip}")
    if total_error > 0:
        print(f"⚠️  ERROR: {total_error}")
    if total_timeout > 0:
        print(f"⏱️  TIMEOUT: {total_timeout}")
    print(f"\nResults saved to: {status_csv}")
    
    return 0 if total_fail == 0 else 1

if __name__ == '__main__':
    sys.exit(main())

