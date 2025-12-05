#!/usr/bin/env python3
"""
Download PDB files from RCSB PDB and place them in data/pdb/.

This script:
1. Downloads PDB files from RCSB PDB (https://files.rcsb.org/download/)
2. Saves them to data/pdb/{PDB_ID}.pdb
3. Validates downloaded files
4. Skips existing files
5. Supports batch downloads from test sets or PDB ID lists

Usage:
    # Download PDBs from test set
    python3 scripts/download_pdbs.py --test-set 100

    # Download specific PDB IDs
    python3 scripts/download_pdbs.py 1H4S 2BNA 3DNA

    # Download from file containing PDB IDs
    python3 scripts/download_pdbs.py --from-file pdb_list.txt
"""

import argparse
import sys
import urllib.request
import urllib.error
from pathlib import Path
from typing import List, Optional, Set
import json
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from scripts.compare_json import load_test_set


# RCSB PDB download URL format
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Alternative download URL (if primary fails)
RCSB_ALTERNATIVE_URL = "https://files.rcsb.org/view/{pdb_id}.pdb"


def download_pdb(pdb_id: str, output_path: Path, retry: int = 3) -> bool:
    """
    Download a PDB file from RCSB PDB.
    
    Args:
        pdb_id: PDB ID (e.g., "1H4S")
        output_path: Where to save the file
        retry: Number of retry attempts
        
    Returns:
        True if successful, False otherwise
    """
    pdb_id_upper = pdb_id.upper()
    
    # Try primary URL first
    urls = [
        RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id_upper),
        RCSB_ALTERNATIVE_URL.format(pdb_id=pdb_id_upper),
    ]
    
    for attempt in range(retry):
        for url in urls:
            try:
                # Create request with user agent (some servers require this)
                req = urllib.request.Request(url)
                req.add_header('User-Agent', 'Mozilla/5.0 (Python PDB Downloader)')
                
                # Download file
                with urllib.request.urlopen(req, timeout=30) as response:
                    # Check if we got valid content
                    content = response.read()
                    
                    # Check if it's an error page (often starts with <!DOCTYPE)
                    if content.startswith(b'<!DOCTYPE'):
                        continue  # Try next URL
                    
                    # Check if it's actually a PDB file (should start with HEADER, REMARK, ATOM, etc.)
                    content_str = content.decode('utf-8', errors='ignore')
                    if not any(content_str.strip().startswith(prefix) 
                              for prefix in ['HEADER', 'REMARK', 'ATOM', 'HETATM', 'CRYST1']):
                        continue  # Try next URL
                    
                    # Write to file
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(output_path, 'wb') as f:
                        f.write(content)
                    
                    # Validate file size (PDB files should be > 0 bytes)
                    if output_path.stat().st_size == 0:
                        output_path.unlink()
                        continue
                    
                    return True
                    
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    # File not found - try next URL or return False
                    continue
                elif e.code == 429:
                    # Rate limited - wait and retry
                    if attempt < retry - 1:
                        time.sleep(2 ** attempt)  # Exponential backoff
                        continue
                else:
                    # Other HTTP error - try next URL
                    continue
            except urllib.error.URLError as e:
                # Network error - wait and retry
                if attempt < retry - 1:
                    time.sleep(2 ** attempt)
                    continue
            except Exception as e:
                # Other error - try next URL
                continue
        
        # If all URLs failed, wait before retrying
        if attempt < retry - 1:
            time.sleep(2 ** attempt)
    
    return False


def validate_pdb_file(pdb_path: Path) -> bool:
    """
    Validate that a downloaded file is a valid PDB file.
    
    Args:
        pdb_path: Path to PDB file
        
    Returns:
        True if valid, False otherwise
    """
    if not pdb_path.exists():
        return False
    
    # Check file size
    if pdb_path.stat().st_size == 0:
        return False
    
    # Check file starts with valid PDB record types
    try:
        with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            first_lines = [f.readline().strip() for _ in range(10)]
            
        # Check for valid PDB record type prefixes
        valid_prefixes = ['HEADER', 'REMARK', 'ATOM', 'HETATM', 'CRYST1', 'TITLE', 'COMPND', 'SOURCE']
        for line in first_lines:
            if line and any(line.startswith(prefix) for prefix in valid_prefixes):
                return True
        
        return False
    except Exception:
        return False


def load_pdb_list(
    project_root: Path,
    test_set: Optional[int],
    from_file: Optional[Path],
    pdb_ids: List[str]
) -> List[str]:
    """Load list of PDB IDs to download."""
    all_ids = set()
    
    # Load from test set
    if test_set:
        test_pdb_ids = load_test_set(project_root, test_set)
        if test_pdb_ids is None:
            print(f"Error: Test set of size {test_set} not found.", file=sys.stderr)
            sys.exit(1)
        all_ids.update(test_pdb_ids)
        print(f"Loaded {len(test_pdb_ids)} PDB IDs from test set {test_set}")
    
    # Load from file
    if from_file:
        if not from_file.exists():
            print(f"Error: File not found: {from_file}", file=sys.stderr)
            sys.exit(1)
        with open(from_file) as f:
            file_ids = [line.strip().upper() for line in f if line.strip()]
        all_ids.update(file_ids)
        print(f"Loaded {len(file_ids)} PDB IDs from file")
    
    # Add command-line PDB IDs
    if pdb_ids:
        all_ids.update(pdb_id.upper() for pdb_id in pdb_ids)
    
    if not all_ids:
        print("Error: No PDB IDs specified. Use --test-set, --from-file, or provide PDB IDs as arguments.", file=sys.stderr)
        sys.exit(1)
    
    return sorted(all_ids)


def download_single_pdb(
    pdb_id: str,
    output_dir: Path,
    args,
    counter_lock: threading.Lock,
    counter: list  # [current, total]
) -> tuple:
    """
    Download a single PDB file (thread-safe wrapper).
    
    Returns:
        (pdb_id, status, message) where status is 'success', 'skipped', or 'failed'
    """
    output_path = output_dir / f"{pdb_id}.pdb"
    
    # Check if file exists
    if args.skip_existing and output_path.exists():
        if args.validate and not validate_pdb_file(output_path):
            # File exists but is invalid, re-download
            pass
        else:
            with counter_lock:
                counter[0] += 1
                current = counter[0]
                total = counter[1]
            return (pdb_id, 'skipped', f"[{current}/{total}] ⏭️  {pdb_id}: Already exists, skipping")
    
    # Download file
    with counter_lock:
        counter[0] += 1
        current = counter[0]
        total = counter[1]
    
    success = download_pdb(pdb_id, output_path, retry=args.retry)
    
    if success:
        # Validate if requested
        if args.validate:
            if validate_pdb_file(output_path):
                file_size = output_path.stat().st_size
                return (pdb_id, 'success', f"[{current}/{total}] ✅ {pdb_id}: ({file_size:,} bytes)")
            else:
                output_path.unlink()
                return (pdb_id, 'failed', f"[{current}/{total}] ❌ {pdb_id}: Downloaded file is invalid")
        else:
            file_size = output_path.stat().st_size
            return (pdb_id, 'success', f"[{current}/{total}] ✅ {pdb_id}: ({file_size:,} bytes)")
    else:
        return (pdb_id, 'failed', f"[{current}/{total}] ❌ {pdb_id}: Failed")


def main():
    parser = argparse.ArgumentParser(
        description="Download PDB files from RCSB PDB",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download PDBs from test set
  python3 scripts/download_pdbs.py --test-set 100

  # Download specific PDB IDs
  python3 scripts/download_pdbs.py 1H4S 2BNA 3DNA

  # Download from file containing PDB IDs (one per line)
  python3 scripts/download_pdbs.py --from-file pdb_list.txt

  # Download with threading (faster)
  python3 scripts/download_pdbs.py --test-set 100 --threads 10

  # Download from test set, but skip existing files
  python3 scripts/download_pdbs.py --test-set 100 --skip-existing

  # Download with progress and resume failed downloads
  python3 scripts/download_pdbs.py --test-set 500 --continue-on-error
"""
    )
    
    parser.add_argument(
        "pdb_ids",
        nargs="*",
        help="PDB IDs to download (optional if using --test-set or --from-file)"
    )
    parser.add_argument(
        "--test-set",
        type=int,
        choices=[10, 50, 100, 500, 1000],
        help="Download PDBs from test set"
    )
    parser.add_argument(
        "--from-file",
        type=Path,
        help="Load PDB IDs from file (one per line)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: data/pdb/)"
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip files that already exist"
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue downloading even if some fail"
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.5,
        help="Delay between downloads in seconds (default: 0.5)"
    )
    parser.add_argument(
        "--retry",
        type=int,
        default=3,
        help="Number of retry attempts (default: 3)"
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        default=True,
        help="Validate downloaded files (default: True)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of parallel download threads (default: 5, set to 0 for single-threaded)"
    )
    
    args = parser.parse_args()
    
    # Set default threads if not specified
    if args.threads is None:
        args.threads = 5  # Reasonable default for RCSB PDB
    elif args.threads == 0:
        args.threads = 1  # Single-threaded mode
    
    # Determine output directory
    output_dir = args.output_dir or (project_root / "data" / "pdb")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load PDB list
    pdb_ids = load_pdb_list(project_root, args.test_set, args.from_file, args.pdb_ids)
    
    print(f"\nDownloading {len(pdb_ids)} PDB files to {output_dir}")
    print(f"Threads: {args.threads}")
    if args.threads > 1:
        print(f"Delay between downloads: {args.delay}s (per thread)")
    else:
        print(f"Delay between downloads: {args.delay}s")
    if args.skip_existing:
        print("Will skip existing files")
    print()
    
    # Track results
    successful = []
    skipped = []
    failed = []
    
    # Thread-safe counter for progress
    counter_lock = threading.Lock()
    counter = [0, len(pdb_ids)]  # [current, total]
    
    start_time = datetime.now()
    
    if args.threads == 1:
        # Single-threaded mode (original behavior)
        for i, pdb_id in enumerate(pdb_ids, 1):
            output_path = output_dir / f"{pdb_id}.pdb"
            
            # Check if file exists
            if args.skip_existing and output_path.exists():
                if args.validate and not validate_pdb_file(output_path):
                    print(f"[{i}/{len(pdb_ids)}] ❌ {pdb_id}: Existing file is invalid, re-downloading...")
                else:
                    print(f"[{i}/{len(pdb_ids)}] ⏭️  {pdb_id}: Already exists, skipping")
                    skipped.append(pdb_id)
                    continue
            
            # Download file
            print(f"[{i}/{len(pdb_ids)}] ⬇️  {pdb_id}: Downloading...", end=" ", flush=True)
            
            success = download_pdb(pdb_id, output_path, retry=args.retry)
            
            if success:
                # Validate if requested
                if args.validate:
                    if validate_pdb_file(output_path):
                        file_size = output_path.stat().st_size
                        print(f"✅ ({file_size:,} bytes)")
                        successful.append(pdb_id)
                    else:
                        print(f"❌ Downloaded file is invalid")
                        output_path.unlink()
                        failed.append(pdb_id)
                        if not args.continue_on_error:
                            print(f"\nError: Failed to download {pdb_id}. Use --continue-on-error to continue.")
                            sys.exit(1)
                else:
                    file_size = output_path.stat().st_size
                    print(f"✅ ({file_size:,} bytes)")
                    successful.append(pdb_id)
            else:
                print(f"❌ Failed")
                failed.append(pdb_id)
                if not args.continue_on_error:
                    print(f"\nError: Failed to download {pdb_id}. Use --continue-on-error to continue.")
                    sys.exit(1)
            
            # Delay between downloads (be nice to RCSB servers)
            if i < len(pdb_ids):
                time.sleep(args.delay)
    else:
        # Multi-threaded mode
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            # Submit all download tasks
            future_to_pdb = {
                executor.submit(download_single_pdb, pdb_id, output_dir, args, counter_lock, counter): pdb_id
                for pdb_id in pdb_ids
            }
            
            # Process completed downloads
            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result_pdb_id, status, message = future.result()
                    
                    # Print result
                    print(message)
                    
                    # Track results
                    if status == 'success':
                        successful.append(result_pdb_id)
                    elif status == 'skipped':
                        skipped.append(result_pdb_id)
                    else:  # failed
                        failed.append(result_pdb_id)
                        if not args.continue_on_error:
                            print(f"\nError: Failed to download {result_pdb_id}. Use --continue-on-error to continue.")
                            # Cancel remaining futures
                            for f in future_to_pdb:
                                f.cancel()
                            sys.exit(1)
                    
                    # Small delay to avoid overwhelming RCSB servers
                    # (downloads happen in parallel, so use minimal delay)
                    if args.delay > 0 and args.threads > 1:
                        time.sleep(args.delay * 0.1)  # Reduced delay for parallel mode
                    
                except Exception as e:
                    failed.append(pdb_id)
                    with counter_lock:
                        counter[0] += 1
                        current = counter[0]
                        total = counter[1]
                    print(f"[{current}/{total}] ❌ {pdb_id}: Exception: {e}")
                    if not args.continue_on_error:
                        # Cancel remaining futures
                        for f in future_to_pdb:
                            f.cancel()
                        sys.exit(1)
    
    # Print summary
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print("\n" + "=" * 80)
    print("DOWNLOAD SUMMARY")
    print("=" * 80)
    print(f"Total: {len(pdb_ids)}")
    print(f"Successful: {len(successful)}")
    print(f"Skipped: {len(skipped)}")
    print(f"Failed: {len(failed)}")
    print(f"Duration: {duration:.1f} seconds ({duration/60:.1f} minutes)")
    print()
    
    if failed:
        print("Failed PDB IDs:")
        for pdb_id in failed:
            print(f"  - {pdb_id}")
        print()
        print("You can retry failed downloads by running:")
        print(f"  python3 scripts/download_pdbs.py {' '.join(failed)} --continue-on-error")
    else:
        print("✅ All downloads completed successfully!")
    
    # Save download log
    log_file = output_dir / ".download_log.json"
    log_data = {
        "timestamp": datetime.now().isoformat(),
        "successful": successful,
        "skipped": skipped,
        "failed": failed,
        "total": len(pdb_ids),
        "duration_seconds": duration,
    }
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=2)
    print(f"\n📋 Download log saved to: {log_file}")


if __name__ == "__main__":
    main()

