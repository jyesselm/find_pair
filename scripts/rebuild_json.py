#!/usr/bin/env python3
"""
Unified JSON Rebuilding Script

This script provides comprehensive functionality for rebuilding, regenerating,
validating, and cleaning JSON files (both legacy and modern).

Usage:
    # Regenerate all JSON files (legacy and modern)
    python3 scripts/rebuild_json.py regenerate

    # Regenerate only legacy JSON files
    python3 scripts/rebuild_json.py regenerate --legacy-only

    # Regenerate only modern JSON files
    python3 scripts/rebuild_json.py regenerate --modern-only

    # Regenerate specific PDB file(s)
    python3 scripts/rebuild_json.py regenerate 1H4S
    python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA

    # Validate existing JSON files
    python3 scripts/rebuild_json.py validate

    # Clean invalid/empty JSON files
    python3 scripts/rebuild_json.py clean --execute

    # Reorganize JSON files (array to grouped format)
    python3 scripts/rebuild_json.py reorganize <json_file>
"""

import json
import os
import sys
import subprocess
import click
import random
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple
import time
from collections import defaultdict
from datetime import datetime

from x3dna_json_compare import JsonValidator


def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    if not legacy_exe.exists():
        legacy_exe = None

    modern_exe = project_root / "build" / "generate_modern_json"
    if not modern_exe.exists():
        modern_exe = None

    return legacy_exe, modern_exe


def get_all_pdb_ids(project_root: Path) -> List[str]:
    """Get all PDB IDs from the data/pdb directory."""
    pdb_dir = project_root / "data" / "pdb"
    if not pdb_dir.exists():
        return []

    pdb_ids = []
    for pdb_file in pdb_dir.glob("*.pdb"):
        pdb_ids.append(pdb_file.stem)

    return sorted(pdb_ids)


def get_test_set_path(project_root: Path, size: int) -> Path:
    """Get the path to a test set file."""
    test_sets_dir = project_root / "data" / "test_sets"
    test_sets_dir.mkdir(parents=True, exist_ok=True)
    return test_sets_dir / f"test_set_{size}.json"


def load_test_set(project_root: Path, size: int) -> Optional[List[str]]:
    """Load a test set of the specified size."""
    test_set_file = get_test_set_path(project_root, size)
    if not test_set_file.exists():
        return None

    try:
        with open(test_set_file, "r") as f:
            data = json.load(f)
            return data.get("pdb_ids", [])
    except Exception:
        return None


def save_test_set(project_root: Path, size: int, pdb_ids: List[str]) -> bool:
    """Save a test set to disk."""
    test_set_file = get_test_set_path(project_root, size)
    try:
        data = {
            "size": size,
            "pdb_ids": pdb_ids,
            "generated": datetime.now().isoformat(),
        }
        with open(test_set_file, "w") as f:
            json.dump(data, f, indent=2)
        return True
    except Exception:
        return False


def get_valid_pdbs_with_atoms(project_root: Path) -> List[str]:
    """Get list of PDBs that have both atom data AND frame calculation data in legacy JSON."""
    valid_file = project_root / "data" / "valid_pdbs.json"
    if valid_file.exists():
        try:
            with open(valid_file) as f:
                data = json.load(f)
                # Prefer PDBs with both atoms and frames
                both = data.get("valid_pdbs_with_atoms_and_frames", [])
                if both:
                    return both
                # Fallback to atoms only if no both list exists
                return data.get("valid_pdbs_with_atoms", [])
        except Exception:
            pass

    # Fallback: check existing JSON files
    legacy_json_dir = project_root / "data" / "json_legacy"
    valid_pdbs = []

    for json_file in legacy_json_dir.glob("*.json"):
        if "_" in json_file.stem:
            continue

        pdb_id = json_file.stem
        has_atoms = False

        # Check split file
        atoms_file = legacy_json_dir / f"{pdb_id}_pdb_atoms.json"
        if atoms_file.exists():
            try:
                with open(atoms_file) as af:
                    atoms_data = json.load(af)
                    if isinstance(atoms_data, list) and len(atoms_data) > 0:
                        atoms = atoms_data[0].get("atoms", [])
                        if atoms:
                            has_atoms = True
            except Exception:
                pass

        # Check main file
        if not has_atoms:
            try:
                with open(json_file) as f:
                    data = json.load(f)
                    calc = data.get("calculations", [])
                    if isinstance(calc, list):
                        for c in calc:
                            if c.get("type") == "pdb_atoms" and "atoms" in c:
                                atoms = c.get("atoms", [])
                                if atoms:
                                    has_atoms = True
                                    break
            except Exception:
                pass

        # Also check for frame calculations (base_frame_calc, frame_calc, ls_fitting)
        has_frames = False
        base_frame_file = legacy_json_dir / f"{pdb_id}_base_frame_calc.json"
        frame_calc_file = legacy_json_dir / f"{pdb_id}_frame_calc.json"
        ls_fitting_file = legacy_json_dir / f"{pdb_id}_ls_fitting.json"

        for frame_file in [base_frame_file, frame_calc_file, ls_fitting_file]:
            if frame_file.exists():
                try:
                    with open(frame_file) as f:
                        frame_data = json.load(f)
                        if isinstance(frame_data, list) and len(frame_data) > 0:
                            has_frames = True
                            break
                except Exception:
                    pass

        # Only include if has both atoms AND frames
        if has_atoms and has_frames:
            valid_pdbs.append(pdb_id)

    return sorted(valid_pdbs)


def generate_test_set(
    project_root: Path, size: int, seed: Optional[int] = None
) -> List[str]:
    """Generate a test set by randomly selecting PDB files that have atom data."""
    valid_pdbs = get_valid_pdbs_with_atoms(project_root)

    if not valid_pdbs:
        click.echo(
            "Warning: No valid PDBs with atom data found. Falling back to all PDB files.",
            err=True,
        )
        valid_pdbs = get_all_pdb_ids(project_root)

    if not valid_pdbs:
        return []

    if size > len(valid_pdbs):
        click.echo(
            f"Warning: Requested size {size} is larger than available valid PDBs ({len(valid_pdbs)}). Using all available.",
            err=True,
        )
        size = len(valid_pdbs)

    if seed is not None:
        random.seed(seed)

    selected = random.sample(valid_pdbs, size)
    return sorted(selected)


def regenerate_legacy_json(
    pdb_id: str, executable_path: str, project_root: str
) -> Dict:
    """Regenerate legacy JSON for a single PDB with validation."""
    if isinstance(project_root, str):
        project_root = Path(project_root)
    if isinstance(executable_path, str):
        executable_path = Path(executable_path)

    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_output_dir = project_root / "data" / "json_legacy"  # Output directory

    if not pdb_file.exists():
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": f"PDB file not found: {pdb_file}",
        }

    try:
        # Legacy code outputs split files to json_legacy directory
        # It will create files like: <PDB_ID>_pdb_atoms.json, <PDB_ID>_base_pair.json, etc.
        cmd = [str(executable_path), str(pdb_file)]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(project_root / "org"),
        )

        # Check if any JSON files were created (legacy outputs split files)
        # Check for at least one record type file
        record_types = ["pdb_atoms", "base_frame_calc", "base_pair"]
        files_created = False
        for record_type in record_types:
            # Legacy outputs to root: <PDB_ID>_<record_type>.json
            json_file = json_output_dir / f"{pdb_id}_{record_type}.json"
            if json_file.exists():
                files_created = True
                # Validate the file
                is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(
                    json_file, remove_invalid=True
                )
                if not is_valid:
                    return {
                        "pdb_id": pdb_id,
                        "status": "error",
                        "message": f"Invalid JSON removed: {error_msg}",
                        "removed": was_removed,
                    }

        if not files_created:
            return {
                "pdb_id": pdb_id,
                "status": "error",
                "message": f"No JSON files created in {json_output_dir}",
                "stderr": result.stderr[-500:] if result.stderr else "",
            }

        return {
            "pdb_id": pdb_id,
            "status": "success",
            "message": "JSON regenerated and validated successfully",
        }

    except subprocess.TimeoutExpired:
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": "Timeout after 5 minutes",
            "removed": False,
        }
    except Exception as e:
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": str(e),
            "removed": False,
        }


def regenerate_modern_json(
    pdb_id: str,
    executable_path: Path,
    project_root: Path,
    use_legacy_mode: bool = False,
) -> Dict:
    """Regenerate modern JSON for a single PDB."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_output_dir = project_root / "data" / "json"  # Output directory, not file

    if not pdb_file.exists():
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": f"PDB file not found: {pdb_file}",
        }

    try:
        # Create output directory if needed
        json_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Pass directory to executable (not a file path)
        cmd = [str(executable_path), str(pdb_file), str(json_output_dir)]
        if use_legacy_mode:
            cmd.append("--legacy")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        # Check if any JSON files were created in the output directory
        # Look for at least one record type file
        record_types = ["pdb_atoms", "base_frame_calc", "base_pair"]
        files_created = False
        for record_type in record_types:
            json_file = json_output_dir / record_type / f"{pdb_id}.json"
            if json_file.exists():
                files_created = True
                # Validate the file
                is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(
                    json_file, remove_invalid=True
                )
                if not is_valid:
                    return {
                        "pdb_id": pdb_id,
                        "status": "error",
                        "message": f"Invalid JSON removed: {error_msg}",
                        "removed": was_removed,
                    }

        if files_created:
            return {
                "pdb_id": pdb_id,
                "status": "success",
                "message": "JSON regenerated and validated successfully",
            }
        else:
            return {
                "pdb_id": pdb_id,
                "status": "error",
                "message": f"No JSON files created in {json_output_dir}",
            }

    except subprocess.TimeoutExpired:
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": "Timeout during generation",
            "removed": False,
        }
    except Exception as e:
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": str(e),
            "removed": False,
        }


def validate_existing_json(json_file: Path) -> Dict:
    """Validate an existing JSON file and remove if invalid."""
    pdb_id = json_file.stem

    is_valid, error_msg, was_removed = JsonValidator.validate_and_clean(
        json_file, remove_invalid=True
    )

    if not is_valid:
        return {
            "pdb_id": pdb_id,
            "status": "removed" if was_removed else "error",
            "message": f"Invalid JSON removed: {error_msg}",
            "removed": was_removed,
        }

    return {
        "pdb_id": pdb_id,
        "status": "valid",
        "message": "JSON is valid",
    }


def reorganize_json_file(json_file: Path) -> bool:
    """Reorganize JSON file to group calculations by type."""
    try:
        with open(json_file, "r") as f:
            data = json.load(f)

        if isinstance(data.get("calculations"), dict):
            return True

        calculations = data.get("calculations", [])
        if not isinstance(calculations, list):
            return True

        grouped = defaultdict(list)
        for calc in calculations:
            calc_type = calc.get("type")
            if calc_type:
                calc_copy = {k: v for k, v in calc.items() if k != "type"}
                grouped[calc_type].append(calc_copy)

        data["calculations"] = dict(grouped)

        with open(json_file, "w") as f:
            json.dump(data, f, indent=2)

        return True
    except Exception as e:
        click.echo(f"Error reorganizing JSON: {e}", err=True)
        return False


@click.group()
def cli():
    """Rebuild, regenerate, validate, and clean JSON files."""
    pass


@cli.command()
@click.option(
    "--legacy-only",
    is_flag=True,
    help="Regenerate only legacy JSON files",
)
@click.option(
    "--modern-only",
    is_flag=True,
    help="Regenerate only modern JSON files",
)
@click.option(
    "--legacy-mode",
    is_flag=True,
    help="Use legacy mode for modern JSON generation",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=None,
    help="Number of threads (default: CPU count)",
)
@click.option(
    "--test-set",
    type=click.Choice(["10", "50", "100", "500", "1000"], case_sensitive=False),
    help="Use a saved test set of the specified size (10, 50, 100, 500, or 1000 PDBs)",
)
@click.argument("pdb_ids", nargs=-1)
def regenerate(legacy_only, modern_only, legacy_mode, threads, test_set, pdb_ids):
    """Regenerate JSON files for specified PDB(s) or all available."""
    project_root = Path(__file__).parent.parent
    legacy_exe, modern_exe = find_executables(project_root)

    if test_set:
        test_pdb_ids = load_test_set(project_root, int(test_set))
        if test_pdb_ids is None:
            click.echo(
                f"Error: Test set of size {test_set} not found. Generate it first with 'generate-test-sets' command.",
                err=True,
            )
            return
        pdb_ids = test_pdb_ids
        click.echo(f"Using test set of size {test_set}: {len(pdb_ids)} PDB files")
    elif not pdb_ids:
        pdb_ids = get_all_pdb_ids(project_root)
        if not pdb_ids:
            click.echo("No PDB files found!", err=True)
            return
        click.echo(f"Found {len(pdb_ids)} PDB files to process")

    max_workers = threads or os.cpu_count() or 4

    # Regenerate legacy JSON
    if not modern_only and legacy_exe:
        click.echo("=" * 80)
        click.echo("REGENERATING LEGACY JSON FILES")
        click.echo("=" * 80)
        click.echo(f"Processing {len(pdb_ids)} PDBs...")
        click.echo()

        start_time = time.time()
        results = []

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {
                executor.submit(
                    regenerate_legacy_json, pdb_id, str(legacy_exe), str(project_root)
                ): pdb_id
                for pdb_id in pdb_ids
            }

            completed = 0
            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result = future.result()
                    results.append(result)
                    completed += 1

                    status = "✓" if result["status"] == "success" else "✗"
                    click.echo(
                        f"{status} [{completed}/{len(pdb_ids)}] {pdb_id}: {result['message']}"
                    )

                except Exception as e:
                    results.append(
                        {
                            "pdb_id": pdb_id,
                            "status": "error",
                            "message": f"Exception: {str(e)}",
                        }
                    )
                    completed += 1
                    click.echo(
                        f"✗ [{completed}/{len(pdb_ids)}] {pdb_id}: Exception: {str(e)}"
                    )

        elapsed = time.time() - start_time
        success_count = sum(1 for r in results if r["status"] == "success")
        error_count = sum(1 for r in results if r["status"] == "error")

        click.echo()
        click.echo("=" * 80)
        click.echo("LEGACY JSON REGENERATION SUMMARY")
        click.echo("=" * 80)
        click.echo(f"Total processed: {len(results)}")
        click.echo(f"Success: {success_count}")
        click.echo(f"Errors: {error_count}")
        click.echo(f"Time elapsed: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        click.echo()
    elif not modern_only:
        click.echo(
            "Warning: Legacy executable not found, skipping legacy regeneration",
            err=True,
        )

    # Regenerate modern JSON
    if not legacy_only and modern_exe:
        click.echo("=" * 80)
        click.echo("REGENERATING MODERN JSON FILES")
        click.echo("=" * 80)
        click.echo(f"Processing {len(pdb_ids)} PDBs...")
        click.echo()

        start_time = time.time()
        results = []

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {
                executor.submit(
                    regenerate_modern_json,
                    pdb_id,
                    modern_exe,
                    project_root,
                    legacy_mode,
                ): pdb_id
                for pdb_id in pdb_ids
            }

            completed = 0
            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result = future.result()
                    results.append(result)
                    completed += 1

                    status = "✓" if result["status"] == "success" else "✗"
                    if completed % 100 == 0 or completed == len(pdb_ids):
                        click.echo(
                            f"{status} [{completed}/{len(pdb_ids)}] {pdb_id}: {result['message']}"
                        )

                except Exception as e:
                    results.append(
                        {
                            "pdb_id": pdb_id,
                            "status": "error",
                            "message": f"Exception: {str(e)}",
                        }
                    )
                    completed += 1
                    click.echo(
                        f"✗ [{completed}/{len(pdb_ids)}] {pdb_id}: Exception: {str(e)}"
                    )

        elapsed = time.time() - start_time
        success_count = sum(1 for r in results if r["status"] == "success")
        error_count = sum(1 for r in results if r["status"] == "error")

        click.echo()
        click.echo("=" * 80)
        click.echo("MODERN JSON REGENERATION SUMMARY")
        click.echo("=" * 80)
        click.echo(f"Total processed: {len(results)}")
        click.echo(f"Success: {success_count}")
        click.echo(f"Errors: {error_count}")
        click.echo(f"Time elapsed: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        click.echo()
    elif not legacy_only:
        click.echo(
            "Warning: Modern executable not found, skipping modern regeneration",
            err=True,
        )


@cli.command()
@click.option(
    "--legacy-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory containing legacy JSON files",
)
@click.option(
    "--modern-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory containing modern JSON files",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=None,
    help="Number of threads (default: CPU count)",
)
def validate(legacy_dir, modern_dir, threads):
    """Validate existing JSON files and remove invalid ones."""
    project_root = Path(__file__).parent.parent

    if legacy_dir is None:
        legacy_dir = project_root / "data" / "json_legacy"
    if modern_dir is None:
        modern_dir = project_root / "data" / "json"

    max_workers = threads or os.cpu_count() or 4

    # Validate legacy JSON files
    if legacy_dir.exists():
        click.echo("Validating legacy JSON files...")
        legacy_json_files = list(legacy_dir.glob("*.json"))
        click.echo(f"Found {len(legacy_json_files)} legacy JSON files")

        removed_count = 0
        valid_count = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {
                executor.submit(validate_existing_json, json_file): json_file
                for json_file in legacy_json_files
            }

            for future in as_completed(future_to_file):
                result = future.result()
                if result["status"] == "removed":
                    removed_count += 1
                    click.echo(f"  ✗ Removed: {result['pdb_id']} - {result['message']}")
                elif result["status"] == "valid":
                    valid_count += 1

        click.echo(f"\nLegacy JSON: {valid_count} valid, {removed_count} removed")

    # Validate modern JSON files
    if modern_dir.exists():
        click.echo("\nValidating modern JSON files...")
        modern_json_files = list(modern_dir.glob("*.json"))
        click.echo(f"Found {len(modern_json_files)} modern JSON files")

        removed_count = 0
        valid_count = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {
                executor.submit(validate_existing_json, json_file): json_file
                for json_file in modern_json_files
            }

            for future in as_completed(future_to_file):
                result = future.result()
                if result["status"] == "removed":
                    removed_count += 1
                    click.echo(f"  ✗ Removed: {result['pdb_id']} - {result['message']}")
                elif result["status"] == "valid":
                    valid_count += 1

        click.echo(f"\nModern JSON: {valid_count} valid, {removed_count} removed")


@cli.command()
@click.option(
    "--directory",
    "-d",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory to clean (default: both legacy and modern)",
)
@click.option(
    "--legacy-only",
    is_flag=True,
    help="Clean only legacy JSON files",
)
@click.option(
    "--modern-only",
    is_flag=True,
    help="Clean only modern JSON files",
)
@click.option(
    "--execute",
    "-e",
    is_flag=True,
    help="Actually remove files (default is dry-run mode)",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=None,
    help="Number of threads (default: CPU count)",
)
def clean(directory, legacy_only, modern_only, execute, threads):
    """Remove corrupted or empty JSON files."""
    project_root = Path(__file__).parent.parent

    if directory:
        directories = [directory]
    else:
        directories = []
        if not modern_only:
            legacy_dir = project_root / "data" / "json_legacy"
            if legacy_dir.exists():
                directories.append(legacy_dir)
        if not legacy_only:
            modern_dir = project_root / "data" / "json"
            if modern_dir.exists():
                directories.append(modern_dir)

    if not directories:
        click.echo("No directories to clean!", err=True)
        return

    max_workers = threads or os.cpu_count() or 4

    for dir_path in directories:
        click.echo(f"\nCleaning: {dir_path}")
        click.echo(
            f"Mode: {'DRY RUN (no files will be deleted)' if not execute else 'LIVE (files will be deleted)'}"
        )

        json_files = list(dir_path.glob("*.json"))
        if not json_files:
            click.echo("No JSON files found")
            continue

        click.echo(f"Found {len(json_files)} JSON files")

        corrupted_files = []
        empty_files = []
        valid_files = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {
                executor.submit(validate_existing_json, json_file): json_file
                for json_file in json_files
            }

            for future in as_completed(future_to_file):
                result = future.result()
                if result["status"] == "removed":
                    if "empty" in result["message"].lower():
                        empty_files.append(result["pdb_id"])
                    else:
                        corrupted_files.append(result["pdb_id"])
                elif result["status"] == "valid":
                    valid_files += 1

        click.echo(f"\nSummary:")
        click.echo(f"  Valid files: {valid_files}")
        click.echo(f"  Empty files: {len(empty_files)}")
        click.echo(f"  Corrupted files: {len(corrupted_files)}")
        click.echo(f"  Total to remove: {len(empty_files) + len(corrupted_files)}")

        if execute and (empty_files or corrupted_files):
            removed = 0
            for pdb_id in empty_files + corrupted_files:
                json_file = dir_path / f"{pdb_id}.json"
                if json_file.exists():
                    try:
                        json_file.unlink()
                        removed += 1
                    except Exception as e:
                        click.echo(f"  Error removing {pdb_id}: {e}", err=True)
            click.echo(f"\nSuccessfully removed {removed} files")
        elif not execute and (empty_files or corrupted_files):
            click.echo(
                f"\n[DRY RUN] Would remove {len(empty_files) + len(corrupted_files)} files"
            )
            click.echo("Run with --execute to actually remove files")


@cli.command()
@click.argument("json_file", type=click.Path(path_type=Path, exists=True))
def reorganize(json_file):
    """Reorganize JSON file to group calculations by type."""
    if reorganize_json_file(json_file):
        click.echo(f"✓ Successfully reorganized: {json_file}")
    else:
        click.echo(f"✗ Failed to reorganize: {json_file}", err=True)
        sys.exit(1)


@cli.command("generate-test-sets")
@click.option(
    "--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)"
)
@click.option(
    "--force", is_flag=True, help="Regenerate test sets even if they already exist"
)
def generate_test_sets(seed, force):
    """Generate test sets of 10, 50, 100, 500, and 1000 PDBs from available PDB files."""
    project_root = Path(__file__).parent.parent

    all_pdbs = get_all_pdb_ids(project_root)
    if not all_pdbs:
        click.echo("Error: No PDB files found in data/pdb/", err=True)
        return

    click.echo(f"Found {len(all_pdbs)} total PDB files")
    click.echo(f"Using random seed: {seed}")
    click.echo()

    test_sizes = [10, 50, 100, 500, 1000]

    for size in test_sizes:
        test_set_file = get_test_set_path(project_root, size)

        if test_set_file.exists() and not force:
            click.echo(
                f"Test set of size {size} already exists. Use --force to regenerate."
            )
            continue

        if size > len(all_pdbs):
            click.echo(
                f"Skipping size {size}: Only {len(all_pdbs)} PDB files available"
            )
            continue

        click.echo(f"Generating test set of size {size}...")
        pdb_ids = generate_test_set(project_root, size, seed)

        if save_test_set(project_root, size, pdb_ids):
            click.echo(f"✓ Saved test set of size {size} to {test_set_file}")
        else:
            click.echo(f"✗ Failed to save test set of size {size}", err=True)
        click.echo()

    click.echo("Test set generation complete!")
    click.echo(f"Test sets saved in: {project_root / 'data' / 'test_sets'}")


def main():
    """Main entry point."""
    cli()


if __name__ == "__main__":
    main()
