#!/usr/bin/env python3
"""
Script to remove corrupted or empty JSON files from json_legacy directory.

This script:
1. Scans all JSON files in data/json_legacy/
2. Identifies files that are empty or contain invalid JSON
3. Removes those files (with optional dry-run mode)
"""

import sys
from pathlib import Path

# Import reusable validation from lib
from x3dna_json_compare import JsonValidator


def clean_json_legacy(directory, dry_run=True):
    """
    Remove corrupted or empty JSON files from the specified directory.
    
    Args:
        directory: Path to the json_legacy directory
        dry_run: If True, only report what would be deleted without actually deleting
    """
    directory_path = Path(directory)
    
    if not directory_path.exists():
        print(f"Error: Directory '{directory}' does not exist")
        return
    
    if not directory_path.is_dir():
        print(f"Error: '{directory}' is not a directory")
        return
    
    # Find all JSON files
    json_files = list(directory_path.glob("*.json"))
    
    if not json_files:
        print(f"No JSON files found in '{directory}'")
        return
    
    print(f"Found {len(json_files)} JSON files in '{directory}'")
    print(f"Mode: {'DRY RUN (no files will be deleted)' if dry_run else 'LIVE (files will be deleted)'}")
    print("-" * 80)
    
    corrupted_files = []
    empty_files = []
    valid_files = 0
    
    for json_file in json_files:
        is_valid, error_msg = JsonValidator.validate_file(json_file)
        
        if not is_valid:
            if "empty" in error_msg.lower():
                empty_files.append((json_file, error_msg))
            else:
                corrupted_files.append((json_file, error_msg))
        else:
            valid_files += 1
    
    # Report results
    print(f"\nSummary:")
    print(f"  Valid files: {valid_files}")
    print(f"  Empty files: {len(empty_files)}")
    print(f"  Corrupted files: {len(corrupted_files)}")
    print(f"  Total to remove: {len(empty_files) + len(corrupted_files)}")
    
    if empty_files:
        print(f"\nEmpty files ({len(empty_files)}):")
        for file_path, error_msg in empty_files[:10]:  # Show first 10
            print(f"  - {file_path.name}: {error_msg}")
        if len(empty_files) > 10:
            print(f"  ... and {len(empty_files) - 10} more")
    
    if corrupted_files:
        print(f"\nCorrupted files ({len(corrupted_files)}):")
        for file_path, error_msg in corrupted_files[:10]:  # Show first 10
            print(f"  - {file_path.name}: {error_msg[:60]}...")
        if len(corrupted_files) > 10:
            print(f"  ... and {len(corrupted_files) - 10} more")
    
    # Remove files if not in dry-run mode
    if not dry_run and (empty_files or corrupted_files):
        print(f"\nRemoving {len(empty_files) + len(corrupted_files)} files...")
        removed_count = 0
        
        for file_path, _ in empty_files + corrupted_files:
            try:
                file_path.unlink()
                removed_count += 1
            except Exception as e:
                print(f"  Error removing {file_path.name}: {e}")
        
        print(f"Successfully removed {removed_count} files")
    elif dry_run and (empty_files or corrupted_files):
        print(f"\n[DRY RUN] Would remove {len(empty_files) + len(corrupted_files)} files")
        print("Run with --execute to actually remove files")
    else:
        print("\nNo files to remove - all JSON files are valid!")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Remove corrupted or empty JSON files from json_legacy directory"
    )
    parser.add_argument(
        "--directory",
        "-d",
        default="data/json_legacy",
        help="Path to json_legacy directory (default: data/json_legacy)"
    )
    parser.add_argument(
        "--execute",
        "-e",
        action="store_true",
        help="Actually remove files (default is dry-run mode)"
    )
    
    args = parser.parse_args()
    
    # Get the script's directory and resolve relative paths
    script_dir = Path(__file__).parent.parent
    directory = script_dir / args.directory
    
    clean_json_legacy(directory, dry_run=not args.execute)


if __name__ == "__main__":
    main()

