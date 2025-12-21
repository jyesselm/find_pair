#!/usr/bin/env python3
"""
JSON Rebuilding and Regeneration Module

Provides functionality for rebuilding, regenerating, validating, and cleaning
JSON files (both legacy and modern).

This module is used by:
- The validation runner (auto-regenerate before comparison)
- The CLI (manual regeneration commands)
- The legacy scripts/rebuild_json.py (backwards compatibility)
"""

import json
import os
import subprocess
import random
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple
import time
from collections import defaultdict
from datetime import datetime

import click

from . import JsonValidator

# All record types that may be generated
ALL_RECORD_TYPES = [
    "pdb_atoms",
    "base_frame_calc",
    "frame_calc",
    "ls_fitting",
    "base_pair",
    "hbond_list",
    "step_params",
    "helical_params",
    "helix_organization",
]

MAX_RETRY_ATTEMPTS = 3


def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    if not legacy_exe.exists():
        legacy_exe = None

    modern_exe = project_root / "build" / "generate_modern_json"
    if not modern_exe.exists():
        modern_exe = None

    return legacy_exe, modern_exe


def validate_json_file(json_file: Path) -> Tuple[bool, Optional[str]]:
    """
    Validate a single JSON file for corruption.

    Checks for:
    - File exists and is readable
    - File is not empty or truncated (minimum 2 bytes for "[]" or "{}")
    - Valid JSON syntax
    - Content is array or object (not primitive)

    Returns:
        Tuple of (is_valid, error_message or None)
    """
    try:
        # Check file exists
        if not json_file.exists():
            return False, "File does not exist"

        # Check file size (minimum valid JSON is "[]" or "{}" = 2 bytes)
        file_size = json_file.stat().st_size
        if file_size < 2:
            return False, f"File too small ({file_size} bytes), likely truncated"

        # Read and parse JSON
        with open(json_file, 'r') as f:
            content = f.read()

        # Check for obvious truncation (doesn't end with ] or })
        content_stripped = content.rstrip()
        if content_stripped and content_stripped[-1] not in ']})':
            return False, f"File appears truncated (ends with '{content_stripped[-20:]}')"

        # Parse JSON
        data = json.loads(content)

        # Accept both arrays and objects as valid JSON
        if isinstance(data, (list, dict)):
            # Additional check: if it's a list, make sure it's not empty for split files
            # (split files should have at least one record)
            return True, None
        else:
            return False, f"Expected array or object, got {type(data).__name__}"

    except json.JSONDecodeError as e:
        return False, f"JSON parse error: {str(e)}"
    except Exception as e:
        return False, f"Error reading file: {str(e)}"


def validate_all_json_files(
    json_dir: Path, pdb_id: str, record_types: List[str] = None
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate all JSON files for a PDB.

    Returns:
        Tuple of (all_valid, valid_files, invalid_files)
    """
    from .json_file_finder import find_json_file

    if record_types is None:
        record_types = ALL_RECORD_TYPES

    valid_files = []
    invalid_files = []

    for record_type in record_types:
        json_file = find_json_file(json_dir, pdb_id, record_type)
        if json_file and json_file.exists():
            is_valid, error_msg = validate_json_file(json_file)
            if is_valid:
                valid_files.append(str(json_file))
            else:
                invalid_files.append(f"{json_file}: {error_msg}")

    return len(invalid_files) == 0, valid_files, invalid_files


def delete_json_files_for_pdb(json_dir: Path, pdb_id: str) -> int:
    """Delete all JSON files for a PDB before regeneration."""
    from .json_file_finder import find_json_file

    deleted = 0
    for record_type in ALL_RECORD_TYPES:
        json_file = find_json_file(json_dir, pdb_id, record_type)
        if json_file and json_file.exists():
            try:
                json_file.unlink()
                deleted += 1
            except Exception:
                pass
    return deleted


def regenerate_modern_json(
    pdb_id: str,
    project_root: Path,
    use_legacy_mode: bool = False,
) -> Dict:
    """Regenerate modern JSON for a single PDB with validation and retry.

    Args:
        pdb_id: The PDB ID to regenerate
        project_root: Path to project root
        use_legacy_mode: Use legacy mode for generation

    Returns:
        Dict with keys: pdb_id, status ('success' or 'error'), message, retries
    """
    if isinstance(project_root, str):
        project_root = Path(project_root)

    _, modern_exe = find_executables(project_root)
    if not modern_exe:
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": "Modern executable not found at build/generate_modern_json",
        }

    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_output_dir = project_root / "data" / "json"

    if not pdb_file.exists():
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": f"PDB file not found: {pdb_file}",
        }

    last_error = None
    retry_count = 0

    for attempt in range(MAX_RETRY_ATTEMPTS):
        try:
            # Delete existing files before regeneration (clean slate)
            if attempt > 0:
                delete_json_files_for_pdb(json_output_dir, pdb_id)

            # Create output directory if needed
            json_output_dir.mkdir(parents=True, exist_ok=True)

            cmd = [str(modern_exe), str(pdb_file), str(json_output_dir)]
            if use_legacy_mode:
                cmd.append("--legacy")

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            # Validate ALL generated JSON files
            all_valid, valid_files, invalid_files = validate_all_json_files(
                json_output_dir, pdb_id
            )

            if not valid_files:
                last_error = f"No JSON files created in {json_output_dir}"
                if result.stderr:
                    last_error += f" (stderr: {result.stderr[-200:]})"
                retry_count += 1
                continue

            if not all_valid:
                last_error = f"Invalid JSON files: {', '.join(invalid_files)}"
                retry_count += 1
                continue

            # All files validated successfully
            message = f"JSON regenerated and validated ({len(valid_files)} files)"
            if retry_count > 0:
                message += f" after {retry_count} retry(s)"

            return {
                "pdb_id": pdb_id,
                "status": "success",
                "message": message,
                "retries": retry_count,
            }

        except subprocess.TimeoutExpired:
            last_error = "Timeout after 5 minutes"
            retry_count += 1
        except Exception as e:
            last_error = str(e)
            retry_count += 1

    # All retries exhausted
    return {
        "pdb_id": pdb_id,
        "status": "error",
        "message": f"Failed after {MAX_RETRY_ATTEMPTS} attempts: {last_error}",
        "retries": retry_count,
    }


def regenerate_legacy_json(
    pdb_id: str, executable_path: str, project_root: str
) -> Dict:
    """Regenerate legacy JSON for a single PDB with validation and retry."""
    if isinstance(project_root, str):
        project_root = Path(project_root)
    if isinstance(executable_path, str):
        executable_path = Path(executable_path)

    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_output_dir = project_root / "data" / "json_legacy"

    if not pdb_file.exists():
        return {
            "pdb_id": pdb_id,
            "status": "error",
            "message": f"PDB file not found: {pdb_file}",
        }

    last_error = None
    retry_count = 0

    for attempt in range(MAX_RETRY_ATTEMPTS):
        try:
            # Delete existing files before regeneration (clean slate)
            if attempt > 0:
                delete_json_files_for_pdb(json_output_dir, pdb_id)

            cmd = [str(executable_path), str(pdb_file)]
            # Legacy code requires X3DNA environment variable
            env = os.environ.copy()
            if 'X3DNA' not in env:
                # Try to find X3DNA installation
                possible_paths = [
                    Path.home() / "local" / "installs" / "x3dna",
                    Path("/usr/local/x3dna"),
                    Path("/opt/x3dna"),
                ]
                for p in possible_paths:
                    if p.exists():
                        env['X3DNA'] = str(p)
                        break
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=str(project_root / "org"),
                env=env,
            )

            # Validate ALL generated JSON files
            all_valid, valid_files, invalid_files = validate_all_json_files(
                json_output_dir, pdb_id
            )

            if not valid_files:
                last_error = f"No JSON files created in {json_output_dir}"
                if result.stderr:
                    last_error += f" (stderr: {result.stderr[-200:]})"
                retry_count += 1
                continue

            if not all_valid:
                last_error = f"Invalid JSON files: {', '.join(invalid_files)}"
                retry_count += 1
                continue

            # All files validated successfully
            message = f"JSON regenerated and validated ({len(valid_files)} files)"
            if retry_count > 0:
                message += f" after {retry_count} retry(s)"

            return {
                "pdb_id": pdb_id,
                "status": "success",
                "message": message,
                "retries": retry_count,
            }

        except subprocess.TimeoutExpired:
            last_error = "Timeout after 5 minutes"
            retry_count += 1
        except Exception as e:
            last_error = str(e)
            retry_count += 1

    # All retries exhausted
    return {
        "pdb_id": pdb_id,
        "status": "error",
        "message": f"Failed after {MAX_RETRY_ATTEMPTS} attempts: {last_error}",
        "retries": retry_count,
    }


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
    resources_dir = project_root / "resources" / "test_sets"
    data_dir = project_root / "data" / "test_sets"

    test_sets_dir = resources_dir if resources_dir.exists() else data_dir
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
                both = data.get("valid_pdbs_with_atoms_and_frames", [])
                if both:
                    return both
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

        # Also check for frame calculations
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


# =============================================================================
# CLI Commands
# =============================================================================

@click.group()
def rebuild_cli():
    """Rebuild, regenerate, validate, and clean JSON files."""
    pass


@rebuild_cli.command('regenerate')
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
    type=click.Choice(["10", "50", "100", "500", "1000", "fast"], case_sensitive=False),
    help="Use a test set: 10/50/100/500/1000 or 'fast' (3602 validated PDBs)",
)
@click.option(
    "--project-root",
    type=click.Path(path_type=Path, exists=True),
    help="Project root directory",
)
@click.argument("pdb_ids", nargs=-1)
def regenerate_cmd(legacy_only, modern_only, legacy_mode, threads, test_set, project_root, pdb_ids):
    """Regenerate JSON files for specified PDB(s) or all available."""
    if project_root is None:
        project_root = Path.cwd()

    legacy_exe, modern_exe = find_executables(project_root)

    if test_set:
        if test_set.lower() == "fast":
            fast_file = project_root / "data" / "valid_pdbs_fast.json"
            if not fast_file.exists():
                click.echo(f"Error: {fast_file} not found!", err=True)
                return
            with open(fast_file) as f:
                data = json.load(f)
                pdb_ids = data.get("valid_pdbs_with_atoms_and_frames", [])
            click.echo(f"Using fast PDB list: {len(pdb_ids)} PDB files")
        else:
            test_pdb_ids = load_test_set(project_root, int(test_set))
            if test_pdb_ids is None:
                click.echo(
                    f"Error: Test set of size {test_set} not found.",
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

                    status = "+" if result["status"] == "success" else "x"
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
                        f"x [{completed}/{len(pdb_ids)}] {pdb_id}: Exception: {str(e)}"
                    )

        elapsed = time.time() - start_time
        success_count = sum(1 for r in results if r["status"] == "success")
        error_count = sum(1 for r in results if r["status"] == "error")
        retry_count = sum(r.get("retries", 0) for r in results)
        retried_count = sum(1 for r in results if r.get("retries", 0) > 0)

        click.echo()
        click.echo("=" * 80)
        click.echo("LEGACY JSON REGENERATION SUMMARY")
        click.echo("=" * 80)
        click.echo(f"Total processed: {len(results)}")
        click.echo(f"Success: {success_count}")
        click.echo(f"Errors: {error_count}")
        if retried_count > 0:
            click.echo(f"Files requiring retry: {retried_count} ({retry_count} total retries)")
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

                    status = "+" if result["status"] == "success" else "x"
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
                        f"x [{completed}/{len(pdb_ids)}] {pdb_id}: Exception: {str(e)}"
                    )

        elapsed = time.time() - start_time
        success_count = sum(1 for r in results if r["status"] == "success")
        error_count = sum(1 for r in results if r["status"] == "error")
        retry_count = sum(r.get("retries", 0) for r in results)
        retried_count = sum(1 for r in results if r.get("retries", 0) > 0)

        click.echo()
        click.echo("=" * 80)
        click.echo("MODERN JSON REGENERATION SUMMARY")
        click.echo("=" * 80)
        click.echo(f"Total processed: {len(results)}")
        click.echo(f"Success: {success_count}")
        click.echo(f"Errors: {error_count}")
        if retried_count > 0:
            click.echo(f"Files requiring retry: {retried_count} ({retry_count} total retries)")
        click.echo(f"Time elapsed: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        click.echo()
    elif not legacy_only:
        click.echo(
            "Warning: Modern executable not found, skipping modern regeneration",
            err=True,
        )


@rebuild_cli.command('validate-json')
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
@click.option(
    "--project-root",
    type=click.Path(path_type=Path, exists=True),
    help="Project root directory",
)
def validate_json_cmd(legacy_dir, modern_dir, threads, project_root):
    """Validate existing JSON files and remove invalid ones."""
    if project_root is None:
        project_root = Path.cwd()

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
                    click.echo(f"  x Removed: {result['pdb_id']} - {result['message']}")
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
                    click.echo(f"  x Removed: {result['pdb_id']} - {result['message']}")
                elif result["status"] == "valid":
                    valid_count += 1

        click.echo(f"\nModern JSON: {valid_count} valid, {removed_count} removed")


@rebuild_cli.command('clean')
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
@click.option(
    "--project-root",
    type=click.Path(path_type=Path, exists=True),
    help="Project root directory",
)
def clean_cmd(directory, legacy_only, modern_only, execute, threads, project_root):
    """Remove corrupted or empty JSON files."""
    if project_root is None:
        project_root = Path.cwd()

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


@rebuild_cli.command('reorganize')
@click.argument("json_file", type=click.Path(path_type=Path, exists=True))
def reorganize_cmd(json_file):
    """Reorganize JSON file to group calculations by type."""
    if reorganize_json_file(json_file):
        click.echo(f"+ Successfully reorganized: {json_file}")
    else:
        click.echo(f"x Failed to reorganize: {json_file}", err=True)
        raise SystemExit(1)


@rebuild_cli.command("generate-test-sets")
@click.option(
    "--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)"
)
@click.option(
    "--force", is_flag=True, help="Regenerate test sets even if they already exist"
)
@click.option(
    "--project-root",
    type=click.Path(path_type=Path, exists=True),
    help="Project root directory",
)
def generate_test_sets_cmd(seed, force, project_root):
    """Generate test sets of 10, 50, 100, 500, and 1000 PDBs from available PDB files."""
    if project_root is None:
        project_root = Path.cwd()

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
            click.echo(f"+ Saved test set of size {size} to {test_set_file}")
        else:
            click.echo(f"x Failed to save test set of size {size}", err=True)
        click.echo()

    click.echo("Test set generation complete!")
    test_sets_dir = project_root / "resources" / "test_sets"
    click.echo(f"Test sets saved in: {test_sets_dir}")
