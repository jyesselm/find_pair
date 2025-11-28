#!/usr/bin/env python3
"""
Cleanup script to remove optional/debug JSON record types to save storage space.

This script removes JSON record types that are only needed for debugging,
keeping only the essential record types needed for core functionality validation.

Essential record types (kept):
- find_bestpair_selection (CRITICAL - final output)
- base_pair (geometric properties)
- hbond_list (H-bond details)
- base_frame_calc (frame metadata)
- frame_calc (frame matrices)

Optional record types (removed):
- pair_validation (debug only - LARGE)
- distance_checks (debug only)
- pdb_atoms (skip if already verified - LARGE)
- ls_fitting (redundant with frame_calc)

Usage:
    # Dry run (show what would be deleted)
    python3 scripts/cleanup_optional_json.py --dry-run

    # Actually delete files
    python3 scripts/cleanup_optional_json.py

    # Delete specific record types only
    python3 scripts/cleanup_optional_json.py --types pair_validation distance_checks

    # Keep pdb_atoms (don't delete them)
    python3 scripts/cleanup_optional_json.py --keep-atoms
"""

import sys
import argparse
from pathlib import Path
from typing import List, Set

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


# Optional record types that can be deleted
OPTIONAL_RECORD_TYPES = {
    "pair_validation",      # Debug only - LARGE (~2 MB per PDB)
    "distance_checks",      # Debug only (~500 KB per PDB)
    "pdb_atoms",            # Skip if verified - LARGE (~5 MB per PDB)
    "ls_fitting",           # Redundant with frame_calc (~100 KB per PDB)
}

# Essential record types (never delete)
ESSENTIAL_RECORD_TYPES = {
    "find_bestpair_selection",  # CRITICAL - final output
    "base_pair",                 # Geometric properties
    "hbond_list",                # H-bond details
    "base_frame_calc",           # Frame metadata
    "frame_calc",                # Frame matrices
    "ref_frame",                 # Alternative frame format
}


def find_json_files(base_dir: Path, record_type: str) -> List[Path]:
    """Find all JSON files for a record type."""
    files = []
    
    # Check new structure: <record_type>/<PDB_ID>.json
    record_dir = base_dir / record_type
    if record_dir.exists():
        files.extend(record_dir.glob("*.json"))
    
    # Check old structure: <PDB_ID>_<record_type>.json
    pattern = f"*_{record_type}.json"
    files.extend(base_dir.glob(pattern))
    
    return files


def get_file_size_mb(file_path: Path) -> float:
    """Get file size in MB."""
    try:
        return file_path.stat().st_size / (1024 * 1024)
    except:
        return 0.0


def cleanup_record_type(
    base_dir: Path,
    record_type: str,
    dry_run: bool = False,
    verbose: bool = False
) -> tuple[int, float]:
    """Clean up a specific record type. Returns (file_count, total_size_mb)."""
    files = find_json_files(base_dir, record_type)
    
    if not files:
        return 0, 0.0
    
    total_size = sum(get_file_size_mb(f) for f in files)
    
    if dry_run:
        print(f"  [DRY RUN] Would delete {len(files)} files ({total_size:.2f} MB)")
        if verbose:
            for f in sorted(files)[:10]:  # Show first 10
                size_mb = get_file_size_mb(f)
                print(f"    - {f} ({size_mb:.2f} MB)")
            if len(files) > 10:
                print(f"    ... and {len(files) - 10} more files")
    else:
        deleted = 0
        for f in files:
            try:
                f.unlink()
                deleted += 1
            except Exception as e:
                print(f"  Error deleting {f}: {e}")
        
        if deleted > 0:
            print(f"  ✓ Deleted {deleted}/{len(files)} files ({total_size:.2f} MB)")
    
    return len(files), total_size


def main():
    parser = argparse.ArgumentParser(
        description="Clean up optional JSON record types to save storage space"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without actually deleting"
    )
    parser.add_argument(
        "--types",
        nargs="+",
        choices=list(OPTIONAL_RECORD_TYPES),
        default=list(OPTIONAL_RECORD_TYPES),
        help="Specific record types to delete (default: all optional types)"
    )
    parser.add_argument(
        "--keep-atoms",
        action="store_true",
        help="Keep pdb_atoms files (don't delete them)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output"
    )
    parser.add_argument(
        "--modern-only",
        action="store_true",
        help="Only clean modern JSON files (keep legacy)"
    )
    parser.add_argument(
        "--legacy-only",
        action="store_true",
        help="Only clean legacy JSON files (keep modern)"
    )
    
    args = parser.parse_args()
    
    # Determine which record types to delete
    record_types_to_delete = set(args.types)
    if args.keep_atoms:
        record_types_to_delete.discard("pdb_atoms")
    
    # Validate record types
    invalid_types = record_types_to_delete - OPTIONAL_RECORD_TYPES
    if invalid_types:
        print(f"Error: Invalid record types: {invalid_types}")
        print(f"Valid optional types: {OPTIONAL_RECORD_TYPES}")
        return 1
    
    # Determine which directories to clean
    dirs_to_clean = []
    if not args.legacy_only:
        modern_dir = project_root / "data" / "json"
        if modern_dir.exists():
            dirs_to_clean.append(("Modern", modern_dir))
    
    if not args.modern_only:
        legacy_dir = project_root / "data" / "json_legacy"
        if legacy_dir.exists():
            dirs_to_clean.append(("Legacy", legacy_dir))
    
    if not dirs_to_clean:
        print("Error: No JSON directories found")
        return 1
    
    # Show what will be deleted
    mode = "DRY RUN" if args.dry_run else "CLEANUP"
    print(f"=== {mode} ===")
    print(f"Record types to delete: {', '.join(sorted(record_types_to_delete))}")
    print(f"Directories: {', '.join([d[0] for d in dirs_to_clean])}")
    print()
    
    # Clean up each directory
    total_files = 0
    total_size = 0.0
    
    for dir_name, base_dir in dirs_to_clean:
        print(f"{dir_name} JSON ({base_dir}):")
        
        dir_files = 0
        dir_size = 0.0
        
        for record_type in sorted(record_types_to_delete):
            files, size = cleanup_record_type(
                base_dir, record_type, args.dry_run, args.verbose
            )
            dir_files += files
            dir_size += size
        
        total_files += dir_files
        total_size += dir_size
        
        if dir_files > 0:
            print(f"  Total: {dir_files} files, {dir_size:.2f} MB")
        else:
            print(f"  No files found to delete")
        print()
    
    # Summary
    print("=" * 60)
    if args.dry_run:
        print(f"DRY RUN SUMMARY:")
        print(f"  Would delete: {total_files} files")
        print(f"  Would free: {total_size:.2f} MB")
        print()
        print("Run without --dry-run to actually delete files")
    else:
        print(f"CLEANUP SUMMARY:")
        print(f"  Deleted: {total_files} files")
        print(f"  Freed: {total_size:.2f} MB")
        print()
        print("✓ Cleanup complete!")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

