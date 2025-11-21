#!/usr/bin/env python3
"""
Script to remove empty JSON files from a directory.

This script:
1. Scans all JSON files in the specified directory (recursively by default)
2. Identifies files that are empty (0 bytes)
3. Removes those files (with optional dry-run mode)
"""

import os
import sys
from pathlib import Path


def find_empty_json_files(directory, recursive=True):
    """
    Find all empty JSON files in the specified directory.
    
    Args:
        directory: Path to the directory to search
        recursive: If True, search recursively in subdirectories
        
    Returns:
        list: List of Path objects for empty JSON files
    """
    directory_path = Path(directory)
    
    if not directory_path.exists():
        print(f"Error: Directory '{directory}' does not exist")
        return []
    
    if not directory_path.is_dir():
        print(f"Error: '{directory}' is not a directory")
        return []
    
    # Find all JSON files
    if recursive:
        json_files = list(directory_path.rglob("*.json"))
    else:
        json_files = list(directory_path.glob("*.json"))
    
    empty_files = []
    
    for json_file in json_files:
        try:
            if os.path.getsize(json_file) == 0:
                empty_files.append(json_file)
        except OSError as e:
            print(f"Warning: Could not check size of {json_file}: {e}")
    
    return empty_files


def remove_empty_json_files(directory, recursive=True, dry_run=True):
    """
    Remove empty JSON files from the specified directory.
    
    Args:
        directory: Path to the directory to search
        recursive: If True, search recursively in subdirectories
        dry_run: If True, only report what would be deleted without actually deleting
    """
    empty_files = find_empty_json_files(directory, recursive)
    
    if not empty_files:
        print(f"No empty JSON files found in '{directory}'")
        return
    
    print(f"Found {len(empty_files)} empty JSON file(s)")
    print(f"Mode: {'DRY RUN (no files will be deleted)' if dry_run else 'LIVE (files will be deleted)'}")
    print("-" * 80)
    
    # Show files that will be removed
    print("\nEmpty JSON files:")
    for file_path in empty_files[:20]:  # Show first 20
        print(f"  - {file_path}")
    if len(empty_files) > 20:
        print(f"  ... and {len(empty_files) - 20} more")
    
    # Remove files if not in dry-run mode
    if not dry_run:
        print(f"\nRemoving {len(empty_files)} file(s)...")
        removed_count = 0
        
        for file_path in empty_files:
            try:
                file_path.unlink()
                removed_count += 1
            except Exception as e:
                print(f"  Error removing {file_path}: {e}")
        
        print(f"Successfully removed {removed_count} file(s)")
    else:
        print(f"\n[DRY RUN] Would remove {len(empty_files)} file(s)")
        print("Run with --execute to actually remove files")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Remove empty JSON files from a directory"
    )
    parser.add_argument(
        "--directory",
        "-d",
        default="data/json_legacy",
        help="Path to directory to search (default: data/json_legacy)"
    )
    parser.add_argument(
        "--execute",
        "-e",
        action="store_true",
        help="Actually remove files (default is dry-run mode)"
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Do not search recursively in subdirectories"
    )
    
    args = parser.parse_args()
    
    # Get the script's directory and resolve relative paths
    script_dir = Path(__file__).parent.parent
    directory = script_dir / args.directory
    
    remove_empty_json_files(
        directory,
        recursive=not args.no_recursive,
        dry_run=not args.execute
    )


if __name__ == "__main__":
    main()

