#!/usr/bin/env python3
"""
Clean up JSON data directory structure.

This script removes all existing JSON files to start fresh with the new structure:
- data/json/<record_type>/ - modern JSON files organized by record type
  - pdb_atoms/, base_frame_calc/, frame_calc/, base_pair/, etc.
- data/json_legacy/<record_type>/ - legacy JSON files organized by record type
  - pdb_atoms/, base_frame_calc/, frame_calc/, base_pair/, etc.

Usage:
    # Dry run (show what would be deleted)
    python3 scripts/cleanup_json_data.py --dry-run

    # Actually remove files
    python3 scripts/cleanup_json_data.py
"""

import sys
from pathlib import Path
import click


def cleanup_json_data(project_root: Path, dry_run: bool = False):
    """Clean up JSON data directories."""
    
    # Directories to clean
    dirs_to_clean = [
        project_root / "data" / "json",
        project_root / "data" / "json_legacy",
    ]
    
    # Directories to remove completely
    dirs_to_remove = [
        project_root / "data" / "json_final",
        project_root / "data" / "json_test",
        project_root / "data" / "json_test_batch",
        project_root / "data" / "residue_ordering",
        project_root / "data" / "residue_ordering_legacy",
    ]
    
    total_files = 0
    total_dirs = 0
    
    # Clean JSON files from main directories and all subdirectories
    for json_dir in dirs_to_clean:
        if json_dir.exists():
            # Remove JSON files from root directory
            json_files = list(json_dir.glob("*.json"))
            total_files += len(json_files)
            
            if dry_run:
                if json_files:
                    print(f"Would remove {len(json_files)} files from {json_dir}")
            else:
                for json_file in json_files:
                    json_file.unlink()
                if json_files:
                    print(f"Removed {len(json_files)} files from {json_dir}")
            
            # Remove JSON files from all subdirectories
            for subdir in json_dir.iterdir():
                if subdir.is_dir():
                    subdir_json_files = list(subdir.glob("*.json"))
                    total_files += len(subdir_json_files)
                    
                    if dry_run:
                        if subdir_json_files:
                            print(f"Would remove {len(subdir_json_files)} files from {subdir}")
                    else:
                        for json_file in subdir_json_files:
                            json_file.unlink()
                        if subdir_json_files:
                            print(f"Removed {len(subdir_json_files)} files from {subdir}")
    
    # Remove entire directories
    for dir_to_remove in dirs_to_remove:
        if dir_to_remove.exists():
            # Count files recursively
            files_in_dir = list(dir_to_remove.rglob("*"))
            files_count = sum(1 for f in files_in_dir if f.is_file())
            total_files += files_count
            total_dirs += 1
            
            if dry_run:
                print(f"Would remove directory {dir_to_remove} ({files_count} files)")
            else:
                import shutil
                shutil.rmtree(dir_to_remove)
                print(f"Removed directory {dir_to_remove} ({files_count} files)")
    
    # Create new directory structure (one directory per record type)
    record_types = [
        "pdb_atoms",
        "base_frame_calc",
        "frame_calc",
        "base_pair",
        "pair_validation",
        "distance_checks",
        "hbond_list",
        "find_bestpair_selection",
        "atoms",
    ]
    
    new_dirs = []
    for record_type in record_types:
        new_dirs.append(project_root / "data" / "json" / record_type)
        new_dirs.append(project_root / "data" / "json_legacy" / record_type)
    
    for new_dir in new_dirs:
        if dry_run:
            print(f"Would create directory {new_dir}")
        else:
            new_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created directory {new_dir}")
    
    if dry_run:
        print(f"\nDry run: Would remove {total_files} files and {total_dirs} directories")
        print("Run without --dry-run to actually perform cleanup")
    else:
        print(f"\nCleanup complete: Removed {total_files} files and {total_dirs} directories")
        print("New directory structure created")
        print("\nTo regenerate JSON files:")
        print("  python3 scripts/rebuild_json.py regenerate --modern-only")
        print("  python3 scripts/rebuild_json.py regenerate --legacy-only")


@click.command()
@click.option("--dry-run", is_flag=True, help="Show what would be deleted without actually deleting")
def main(dry_run: bool):
    """Clean up JSON data directory structure."""
    project_root = Path(__file__).parent.parent
    
    if not dry_run:
        response = input("This will delete all JSON files. Are you sure? (yes/no): ")
        if response.lower() != "yes":
            print("Aborted")
            sys.exit(1)
    
    cleanup_json_data(project_root, dry_run=dry_run)


if __name__ == "__main__":
    main()

