#!/usr/bin/env python3
"""
Generate missing legacy JSON files for PDBs that were skipped.

This script:
1. Finds PDBs that were skipped due to missing legacy JSON
2. Runs the legacy find_pair_original code to generate JSON
3. The legacy code automatically writes JSON to data/json_legacy/

Usage:
    python3 scripts/generate_missing_legacy_json.py [--threads N] [--limit N]
    
Example:
    # Generate legacy JSON for all skipped PDBs
    python3 scripts/generate_missing_legacy_json.py --threads 8
    
    # Test on first 10 PDBs
    python3 scripts/generate_missing_legacy_json.py --limit 10
    
    # Generate for specific PDB
    python3 scripts/generate_missing_legacy_json.py --pdb 1FFK
"""

import sys
import subprocess
import csv
import json
from pathlib import Path
from typing import List
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from datetime import datetime

csv_lock = threading.Lock()


def get_skipped_pdbs(status_csv: Path) -> List[str]:
    """Get list of PDBs that were skipped due to missing legacy JSON."""
    skipped = []
    
    if not status_csv.exists():
        return skipped
    
    with open(status_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['match_status'] == 'SKIP' and 'No legacy JSON' in row.get('notes', ''):
                skipped.append(row['pdb_id'])
    
    return sorted(skipped)


def generate_legacy_json(pdb_id: str, project_root: Path) -> dict:
    """Run legacy find_pair_original to generate JSON for one PDB."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    legacy_bin = project_root / "org" / "build" / "bin" / "find_pair_original"
    legacy_json = project_root / "data" / "json_legacy" / "base_frame_calc" / f"{pdb_id}.json"
    
    result = {
        "pdb_id": pdb_id,
        "status": "UNKNOWN",
        "error": None,
        "timestamp": datetime.now().isoformat()
    }
    
    if not pdb_file.exists():
        result["status"] = "ERROR"
        result["error"] = "PDB file not found"
        return result
    
    if not legacy_bin.exists():
        result["status"] = "ERROR"
        result["error"] = "Legacy binary not found (run: cd org/build && make)"
        return result
    
    try:
        # Run legacy find_pair_original
        # The legacy code automatically writes JSON to data/json_legacy/
        # We need to redirect stdout to avoid cluttering output
        outfile = project_root / "temp" / f"{pdb_id}.inp"
        outfile.parent.mkdir(exist_ok=True)
        
        cmd = [str(legacy_bin), str(pdb_file), str(outfile)]
        
        output = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minutes for large structures
            cwd=project_root
        )
        
        # Clean up temp file
        if outfile.exists():
            outfile.unlink()
        
        # Check if JSON was generated
        if legacy_json.exists():
            result["status"] = "SUCCESS"
            result["error"] = None
        elif output.returncode != 0:
            result["status"] = "ERROR"
            # Extract meaningful error from stderr
            error_msg = output.stderr.strip()
            if "No atoms found" in error_msg:
                result["error"] = "No atoms found in PDB"
            elif "No nucleotides" in error_msg:
                result["error"] = "No nucleotides found"
            else:
                result["error"] = f"Failed: {error_msg[:200]}"
        else:
            result["status"] = "ERROR"
            result["error"] = "No JSON generated (unknown reason)"
        
    except subprocess.TimeoutExpired:
        result["status"] = "TIMEOUT"
        result["error"] = "Process timed out (>120s)"
    except Exception as e:
        result["status"] = "ERROR"
        result["error"] = str(e)
    
    return result


def process_batch(batch: List[str], project_root: Path, threads: int) -> dict:
    """Process a batch of PDBs in parallel."""
    results = {"SUCCESS": 0, "ERROR": 0, "TIMEOUT": 0}
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_pdb = {
            executor.submit(generate_legacy_json, pdb_id, project_root): pdb_id
            for pdb_id in batch
        }
        
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                result = future.result()
                status = result["status"]
                results[status] = results.get(status, 0) + 1
                
                if status == "SUCCESS":
                    print(f"  ✅ {pdb_id}: Generated legacy JSON")
                elif status == "ERROR":
                    print(f"  ❌ {pdb_id}: {result['error']}")
                elif status == "TIMEOUT":
                    print(f"  ⏱️  {pdb_id}: Timeout")
            except Exception as e:
                print(f"  ❌ {pdb_id}: Exception - {e}")
                results["ERROR"] = results.get("ERROR", 0) + 1
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Generate missing legacy JSON files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Generate for all skipped PDBs with 8 threads
  python3 scripts/generate_missing_legacy_json.py --threads 8
  
  # Test on first 10 PDBs
  python3 scripts/generate_missing_legacy_json.py --limit 10
  
  # Generate for specific PDB
  python3 scripts/generate_missing_legacy_json.py --pdb 1FFK
        '''
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of parallel threads (default: 4)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='Only process first N PDBs (for testing)'
    )
    parser.add_argument(
        '--pdb',
        type=str,
        help='Process only this specific PDB'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=50,
        help='Process N PDBs per batch (default: 50)'
    )
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    status_csv = project_root / "data" / "index_validation_status.csv"
    
    # Check if legacy binary exists
    legacy_bin = project_root / "org" / "build" / "bin" / "find_pair_original"
    if not legacy_bin.exists():
        print(f"❌ Error: Legacy binary not found: {legacy_bin}")
        print("Build it with: cd org/build && cmake .. && make")
        return 1
    
    print(f"✅ Legacy binary found: {legacy_bin}")
    
    # Get list of skipped PDBs
    if args.pdb:
        pdb_ids = [args.pdb]
        print(f"\nProcessing single PDB: {args.pdb}")
    else:
        pdb_ids = get_skipped_pdbs(status_csv)
        print(f"\nFound {len(pdb_ids)} PDBs missing legacy JSON")
        
        if args.limit:
            pdb_ids = pdb_ids[:args.limit]
            print(f"Limited to first {args.limit} PDBs")
    
    if not pdb_ids:
        print("No PDBs to process")
        return 0
    
    print(f"\nGenerating legacy JSON with {args.threads} threads...\n")
    
    # Process in batches
    total_results = {"SUCCESS": 0, "ERROR": 0, "TIMEOUT": 0}
    batch_size = args.batch_size
    total = len(pdb_ids)
    
    for batch_start in range(0, total, batch_size):
        batch_end = min(batch_start + batch_size, total)
        batch = pdb_ids[batch_start:batch_end]
        
        print(f"Batch {batch_start//batch_size + 1}: Processing PDBs {batch_start+1}-{batch_end} of {total}")
        
        batch_results = process_batch(batch, project_root, args.threads)
        
        # Aggregate results
        for status, count in batch_results.items():
            total_results[status] = total_results.get(status, 0) + count
        
        print(f"  Batch: {batch_results.get('SUCCESS', 0)} SUCCESS, "
              f"{batch_results.get('ERROR', 0)} ERROR, "
              f"{batch_results.get('TIMEOUT', 0)} TIMEOUT\n")
    
    # Final summary
    print(f"\n{'='*60}")
    print("LEGACY JSON GENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"✅ SUCCESS: {total_results.get('SUCCESS', 0)}")
    print(f"❌ ERROR:   {total_results.get('ERROR', 0)}")
    print(f"⏱️  TIMEOUT: {total_results.get('TIMEOUT', 0)}")
    print(f"📊 Total:   {total}")
    print(f"\nGenerated JSON files are in: data/json_legacy/base_frame_calc/")
    
    return 0 if total_results.get('ERROR', 0) == 0 else 1


if __name__ == '__main__':
    sys.exit(main())

