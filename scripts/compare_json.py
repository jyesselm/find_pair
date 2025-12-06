#!/usr/bin/env python3
"""
Unified JSON Comparison Script

⚠️  DEPRECATION NOTICE ⚠️
This script is being replaced by the unified fp2-validate CLI tool.
Please use 'fp2-validate compare' instead for all new work.

This script provides a comprehensive, multithreaded tool for comparing JSON files
between legacy and modern implementations. It generates detailed human-readable reports.

Usage:
    # Compare all available files (atoms and frames)
    python3 scripts/compare_json.py compare

    # Compare specific PDB file(s)
    python3 scripts/compare_json.py compare 1H4S
    python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

    # Compare only atoms
    python3 scripts/compare_json.py atoms 1H4S

    # Compare only frames
    python3 scripts/compare_json.py frames 1H4S

    # List available PDB files
    python3 scripts/compare_json.py list

    # Use legacy mode JSON files
    python3 scripts/compare_json.py compare --legacy-mode

    # Save report to file
    python3 scripts/compare_json.py compare --output report.md

    # Verbose output
    python3 scripts/compare_json.py compare --verbose

    # Show only files with differences
    python3 scripts/compare_json.py compare --diff-only
"""

import json
import sys
import click
import subprocess
import random
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import multiprocessing
from datetime import datetime

from x3dna_json_compare import JsonComparator, ComparisonResult, JsonValidator
from x3dna_json_compare.config import load_config, get_comparison_flags
from x3dna_json_compare.verbose_reporter import (
    VerboseReporter, 
    RecordComparison,
    create_record_comparison_from_dicts
)


def get_all_pdb_files(project_root: Path) -> List[str]:
    """Get all PDB IDs from the data/pdb directory."""
    pdb_dir = project_root / "data" / "pdb"
    if not pdb_dir.exists():
        return []

    pdb_ids = []
    for pdb_file in pdb_dir.glob("*.pdb"):
        pdb_ids.append(pdb_file.stem)

    return sorted(pdb_ids)


def find_available_pdbs(project_root: Path, use_legacy_mode: bool = False) -> List[str]:
    """Find all PDB IDs that have both legacy and modern JSON files."""
    from x3dna_json_compare.json_file_finder import find_json_file
    
    legacy_dir = project_root / "data" / "json_legacy"
    modern_dir = project_root / "data" / "json"

    # Find PDBs with files in new structure: <record_type>/<PDB_ID>.json
    # Check pdb_atoms directory as indicator
    legacy_pdb_atoms_dir = legacy_dir / "pdb_atoms"
    modern_pdb_atoms_dir = modern_dir / "pdb_atoms"
    
    legacy_files = set()
    if legacy_pdb_atoms_dir.exists():
        legacy_files = {f.stem for f in legacy_pdb_atoms_dir.glob("*.json")}
    else:
        # Fall back to old structure: <PDB_ID>_*.json in root
        legacy_files = {f.stem.split("_")[0] for f in legacy_dir.glob("*_pdb_atoms.json")}
        # Also check for main files
        legacy_files.update({f.stem for f in legacy_dir.glob("*.json") if not "_" in f.stem})

    modern_files = set()
    if modern_pdb_atoms_dir.exists():
        modern_files = {f.stem for f in modern_pdb_atoms_dir.glob("*.json")}
    else:
        # Fall back to old structure
        modern_files = {f.stem.split("_")[0] for f in modern_dir.glob("*_pdb_atoms.json")}
        # Also check for main files
        for f in modern_dir.glob("*.json"):
            stem = f.stem
            if use_legacy_mode:
                if stem.endswith("_legacy"):
                    modern_files.add(stem[:-7])
            else:
                if not stem.endswith("_legacy"):
                    modern_files.add(stem)

    common = sorted(legacy_files & modern_files)
    return common


def get_test_set_path(project_root: Path, size: int) -> Path:
    """Get the path to a test set file."""
    # Check resources/test_sets first (new location), fall back to data/test_sets (old location)
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
        with open(test_set_file, 'r') as f:
            data = json.load(f)
            return data.get('pdb_ids', [])
    except Exception:
        return None


def save_test_set(project_root: Path, size: int, pdb_ids: List[str]) -> bool:
    """Save a test set to disk."""
    test_set_file = get_test_set_path(project_root, size)
    try:
        data = {
            'size': size,
            'pdb_ids': pdb_ids,
            'generated': datetime.now().isoformat()
        }
        with open(test_set_file, 'w') as f:
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


def generate_test_set(project_root: Path, size: int, seed: Optional[int] = None) -> List[str]:
    """Generate a test set by randomly selecting PDB files that have atom data."""
    valid_pdbs = get_valid_pdbs_with_atoms(project_root)
    
    if not valid_pdbs:
        click.echo("Warning: No valid PDBs with atom data found. Falling back to all PDB files.", err=True)
        valid_pdbs = get_all_pdb_files(project_root)
    
    if not valid_pdbs:
        return []
    
    if size > len(valid_pdbs):
        click.echo(f"Warning: Requested size {size} is larger than available valid PDBs ({len(valid_pdbs)}). Using all available.", err=True)
        size = len(valid_pdbs)
    
    if seed is not None:
        random.seed(seed)
    
    selected = random.sample(valid_pdbs, size)
    return sorted(selected)


def get_modern_file_path(
    project_root: Path, pdb_id: str, use_legacy_mode: bool = False
) -> Optional[Path]:
    """Get a dummy path to modern JSON file (actual files are in record-type directories)."""
    modern_dir = project_root / "data" / "json"

    if use_legacy_mode:
        # Legacy mode files have _legacy suffix
        return modern_dir / f"{pdb_id}_legacy.json"
    else:
        # Return dummy path - JsonComparator will find split files in subdirectories
        return modern_dir / f"{pdb_id}.json"


def find_executables(project_root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find legacy and modern executables."""
    # Legacy executable (org code)
    legacy_exe = project_root / "org" / "build" / "bin" / "find_pair_analyze"
    if not legacy_exe.exists():
        legacy_exe = None

    # Modern executable (new code)
    modern_exe = project_root / "build" / "generate_modern_json"
    if not modern_exe.exists():
        modern_exe = None

    return legacy_exe, modern_exe


def regenerate_legacy_json(
    pdb_id: str, project_root: Path, legacy_exe: Optional[Path] = None
) -> bool:
    """Regenerate legacy JSON file using org code."""
    if legacy_exe is None:
        legacy_exe, _ = find_executables(project_root)
        if legacy_exe is None:
            return False

    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"

    if not pdb_file.exists():
        return False

    try:
        # Ensure output directory exists
        json_file.parent.mkdir(parents=True, exist_ok=True)

        # Run find_pair_analyze from org directory
        cmd = [str(legacy_exe), str(pdb_file)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(project_root / "org"),
        )

        if json_file.exists():
            # Validate the generated JSON
            is_valid, _, _ = JsonValidator.validate_and_clean(
                json_file, remove_invalid=True
            )
            return is_valid
        return False
    except Exception:
        return False


def regenerate_modern_json(
    pdb_id: str,
    project_root: Path,
    modern_exe: Optional[Path] = None,
    use_legacy_mode: bool = False,
) -> bool:
    """Regenerate modern JSON file using new code."""
    if modern_exe is None:
        _, modern_exe = find_executables(project_root)
        if modern_exe is None:
            return False

    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_output_dir = project_root / "data" / "json"  # Output directory, not file

    if not pdb_file.exists():
        return False

    try:
        # Ensure output directory exists
        json_output_dir.mkdir(parents=True, exist_ok=True)

        # Run generate_modern_json (expects directory, not file)
        cmd = [str(modern_exe), str(pdb_file), str(json_output_dir)]
        if use_legacy_mode:
            cmd.append("--legacy")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        # Check if segmented files were created (new structure)
        from x3dna_json_compare.json_file_finder import find_json_file
        json_file_check = find_json_file(json_output_dir, pdb_id, "pdb_atoms")
        if json_file_check and json_file_check.exists():
            return True  # Segmented files created successfully
        
        # Also check for hbond_list as a fallback indicator
        json_file_check = find_json_file(json_output_dir, pdb_id, "hbond_list")
        if json_file_check and json_file_check.exists():
            return True  # Segmented files created successfully
        
        return False
    except Exception:
        return False


def create_comparator(config_path: Optional[Path] = None, project_root: Optional[Path] = None) -> JsonComparator:
    """Create a JsonComparator from YAML config file."""
    if project_root is None:
        project_root = Path.cwd()
    
    # Load config
    if config_path is None:
        # Try to find config in project root
        config_path = project_root / "comparison_config.yaml"
    
    config = load_config(config_path)
    flags = get_comparison_flags(config)
    
    # Create comparator with config values (cache removed for reliability)
    return JsonComparator(
        tolerance=config.get("tolerance", 1e-6),
        **flags
    )


def compare_single_pdb(
    pdb_id: str,
    project_root: Path,
    use_legacy_mode: bool = False,
    comparator: Optional[JsonComparator] = None,
    regenerate: bool = False,
    config_path: Optional[Path] = None,
) -> ComparisonResult:
    """Compare a single PDB file."""
    if comparator is None:
        comparator = create_comparator(config_path, project_root)

    # Use a dummy file path for comparison - the JsonComparator will find split files
    # The actual files are in record-type-specific directories
    legacy_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
    modern_file = project_root / "data" / "json" / f"{pdb_id}.json"
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"

    # Only regenerate if --regenerate flag is set
    if regenerate:
        if not legacy_file.exists():
            click.echo(
                f"  Regenerating legacy JSON for {pdb_id} (from org code)...", err=True
            )
            if regenerate_legacy_json(pdb_id, project_root):
                click.echo(f"  ✓ Legacy JSON regenerated", err=True)
            else:
                return ComparisonResult(
                    pdb_id=pdb_id,
                    status="error",
                    errors=[
                        f"Legacy file not found and could not be regenerated: {legacy_file}"
                    ],
                )
        else:
            # Force regeneration if --regenerate flag is set
            click.echo(
                f"  Regenerating legacy JSON for {pdb_id} (--regenerate flag)...", err=True
            )
            if regenerate_legacy_json(pdb_id, project_root):
                click.echo(f"  ✓ Legacy JSON regenerated", err=True)

        if not modern_file or not modern_file.exists():
            click.echo(
                f"  Regenerating modern JSON for {pdb_id} (from new code)...", err=True
            )
            if regenerate_modern_json(
                pdb_id, project_root, use_legacy_mode=use_legacy_mode
            ):
                click.echo(f"  ✓ Modern JSON regenerated", err=True)
                # Re-get the file path after regeneration
                modern_file = get_modern_file_path(project_root, pdb_id, use_legacy_mode)
            else:
                return ComparisonResult(
                    pdb_id=pdb_id,
                    status="error",
                    errors=[
                        f"Modern file not found and could not be regenerated for {pdb_id}"
                    ],
                )
        else:
            # Force regeneration if --regenerate flag is set
            click.echo(
                f"  Regenerating modern JSON for {pdb_id} (--regenerate flag)...", err=True
            )
            if regenerate_modern_json(
                pdb_id, project_root, use_legacy_mode=use_legacy_mode
            ):
                click.echo(f"  ✓ Modern JSON regenerated", err=True)
                # Re-get the file path after regeneration
                modern_file = get_modern_file_path(project_root, pdb_id, use_legacy_mode)
    
    # Check if files exist in new directory structure
    # Check for at least one record type file in new structure, or fall back to old structure
    from x3dna_json_compare.json_file_finder import find_json_file
    
    legacy_dir = legacy_file.parent
    modern_dir = modern_file.parent
    
    # Check if legacy data exists (new or old structure)
    has_legacy_data = legacy_file.exists()
    if not has_legacy_data:
        # Check new structure: <record_type>/<PDB_ID>.json
        # Try essential record types (since pdb_atoms may have been deleted)
        # Also include step parameter types for step comparisons
        record_types_to_check = ["find_bestpair_selection", "base_pair", "hbond_list", "base_frame_calc", "pdb_atoms", "bpstep_params", "helical_params"]
        for record_type in record_types_to_check:
            legacy_file_check = find_json_file(legacy_dir, pdb_id, record_type)
            if legacy_file_check and legacy_file_check.exists():
                has_legacy_data = True
                break
    
    if not has_legacy_data:
        return ComparisonResult(
            pdb_id=pdb_id,
            status="error",
            errors=[
                f"Legacy JSON files not found for {pdb_id}. Use --regenerate to create them."
            ],
        )
    
    # Check if modern data exists (new or old structure)
    # Skip if modern_file is a directory (segmented JSON structure)
    has_modern_data = modern_file.exists() and modern_file.is_file()
    if not has_modern_data:
        # Check new structure: <record_type>/<PDB_ID>.json
        # Try essential record types (since pdb_atoms may have been deleted)
        # Also include step parameter types for step comparisons
        record_types_to_check = ["find_bestpair_selection", "base_pair", "hbond_list", "base_frame_calc", "pdb_atoms", "bpstep_params", "helical_params"]
        for record_type in record_types_to_check:
            modern_file_check = find_json_file(modern_dir, pdb_id, record_type)
            if modern_file_check and modern_file_check.exists():
                has_modern_data = True
                break
    
    if not has_modern_data:
        return ComparisonResult(
            pdb_id=pdb_id,
            status="error",
            errors=[
                f"Modern JSON files not found for {pdb_id}. Use --regenerate to create them."
            ],
        )

    result = comparator.compare_files(legacy_file, modern_file, pdb_file, pdb_id)
    return result


def analyze_results(results: Dict[str, ComparisonResult]) -> Dict:
    """Analyze comparison results and generate comprehensive statistics."""
    stats = {
        "total_pdbs": len(results),
        "error_count": 0,
        "match_count": 0,
        "diff_count": 0,
        "total_residues": 0,
        "missing_residues": 0,
        "mismatched_calculations": 0,
        "exact_atom_matches": 0,
        "atom_set_matches": 0,
        "real_atom_differences": 0,
        "rms_differences": [],
        "base_type_stats": defaultdict(
            lambda: {
                "total": 0,
                "exact_matches": 0,
                "set_matches": 0,
                "real_diffs": 0,
                "missing": 0,
            }
        ),
        "files_with_differences": [],
        # Atom comparison stats
        "atom_total_legacy": 0,
        "atom_total_modern": 0,
        "atom_common_count": 0,
        "atom_missing_in_modern": 0,
        "atom_extra_in_modern": 0,
        "atom_mismatched_fields": 0,
        "atom_count_difference": False,
        # Step parameter comparison stats
        "step_total_legacy": 0,
        "step_total_modern": 0,
        "step_missing_steps": 0,
        "step_mismatched_steps": 0,
        "helical_total_legacy": 0,
        "helical_total_modern": 0,
        "helical_missing_steps": 0,
        "helical_mismatched_steps": 0,
        "find_bestpair_total_legacy": 0,
        "find_bestpair_total_modern": 0,
        "find_bestpair_missing_in_modern": 0,
        "find_bestpair_extra_in_modern": 0,
        "base_pair_total_legacy": 0,
        "base_pair_total_modern": 0,
        "base_pair_missing_in_modern": 0,
        "base_pair_extra_in_modern": 0,
        "base_pair_mismatched": 0,
        "pair_validation_total_legacy": 0,
        "pair_validation_total_modern": 0,
        "pair_validation_missing_in_modern": 0,
        "pair_validation_extra_in_modern": 0,
        "pair_validation_mismatched": 0,
        "distance_checks_total_legacy": 0,
        "distance_checks_total_modern": 0,
        "distance_checks_missing_in_modern": 0,
        "distance_checks_extra_in_modern": 0,
        "distance_checks_mismatched": 0,
        # H-bond list comparison stats
        "hbond_list_total_legacy": 0,
        "hbond_list_total_modern": 0,
        "hbond_list_missing_in_modern": 0,
        "hbond_list_extra_in_modern": 0,
        "hbond_list_mismatched_pairs": 0,
        # Residue indices comparison stats
        "residue_indices_missing_in_modern": 0,
        "residue_indices_extra_in_modern": 0,
        "residue_indices_num_residue_match": True,
        "residue_indices_mismatched_entries": 0,
    }

    for pdb_id, result in results.items():
        if result.status == "error":
            stats["error_count"] += 1
            continue

        if not result.has_differences():
            stats["match_count"] += 1
        else:
            stats["diff_count"] += 1
            stats["files_with_differences"].append(pdb_id)

        # Process atom comparison
        if result.atom_comparison:
            ac = result.atom_comparison
            stats["atom_total_legacy"] += ac.total_legacy
            stats["atom_total_modern"] += ac.total_modern
            stats["atom_common_count"] += ac.common_count
            stats["atom_missing_in_modern"] += len(ac.missing_in_modern)
            stats["atom_extra_in_modern"] += len(ac.extra_in_modern)
            stats["atom_mismatched_fields"] += len(ac.mismatched_fields)
            if ac.count_difference:
                stats["atom_count_difference"] = True

        if result.frame_comparison:
            fc = result.frame_comparison
            stats["total_residues"] += fc.total_legacy
            stats["missing_residues"] += len(fc.missing_residues)
            stats["mismatched_calculations"] += len(fc.mismatched_calculations)

            matching_residues = fc.total_legacy - len(fc.missing_residues)

            for mismatch in fc.mismatched_calculations:
                leg_atoms = set(mismatch.legacy_matched_atoms)
                mod_atoms = set(mismatch.modern_matched_atoms)
                leg_atoms_list = sorted(mismatch.legacy_matched_atoms)
                mod_atoms_list = sorted(mismatch.modern_matched_atoms)

                if leg_atoms_list == mod_atoms_list:
                    stats["exact_atom_matches"] += 1
                    stats["atom_set_matches"] += 1
                elif leg_atoms == mod_atoms:
                    stats["atom_set_matches"] += 1
                else:
                    stats["real_atom_differences"] += 1

                # Handle nested mismatches structure (base_frame_calc, ls_fitting, or frame_calc)
                mismatches_dict = mismatch.mismatches
                if 'base_frame_calc' in mismatches_dict:
                    mismatches_dict = mismatches_dict['base_frame_calc']
                elif 'ls_fitting' in mismatches_dict:
                    mismatches_dict = mismatches_dict['ls_fitting']
                elif 'frame_calc' in mismatches_dict:
                    mismatches_dict = mismatches_dict['frame_calc']
                
                # Check RMS from mismatches dict if available, otherwise from records
                if 'rms' in mismatches_dict:
                    rms_diff = mismatches_dict['rms'].get('diff', 0.0)
                    if rms_diff > 0.001:
                        stats["rms_differences"].append(rms_diff)
                else:
                    leg_rms = mismatch.legacy_record.get("rms_fit", 0.0)
                    mod_rms = mismatch.modern_record.get("rms_fit", 0.0)
                    if abs(leg_rms - mod_rms) > 0.001:
                        stats["rms_differences"].append(abs(leg_rms - mod_rms))

                base_type = mismatch.legacy_record.get("base_type", "?")
                base_stats = stats["base_type_stats"][base_type]
                base_stats["total"] += 1

                if leg_atoms_list == mod_atoms_list:
                    base_stats["exact_matches"] += 1
                    base_stats["set_matches"] += 1
                elif leg_atoms == mod_atoms:
                    base_stats["set_matches"] += 1
                else:
                    base_stats["real_diffs"] += 1

            for missing in fc.missing_residues:
                base_type = missing.base_type
                stats["base_type_stats"][base_type]["missing"] += 1
        
        # Process step comparison (bpstep_params)
        if result.step_comparison:
            sc = result.step_comparison
            stats["step_total_legacy"] += sc.total_legacy
            stats["step_total_modern"] += sc.total_modern
            stats["step_missing_steps"] += len(sc.missing_steps)
            stats["step_mismatched_steps"] += len(sc.mismatched_steps)
        
        # Process helical comparison (helical_params)
        if result.helical_comparison:
            hc = result.helical_comparison
            stats["helical_total_legacy"] += hc.total_legacy
            stats["helical_total_modern"] += hc.total_modern
            stats["helical_missing_steps"] += len(hc.missing_steps)
            stats["helical_mismatched_steps"] += len(hc.mismatched_steps)
        
        # Find bestpair selection comparison statistics (PRIMARY - actual selected pairs)
        if hasattr(result, 'find_bestpair_comparison') and result.find_bestpair_comparison:
            fbc = result.find_bestpair_comparison
            stats["find_bestpair_total_legacy"] += fbc.total_legacy
            stats["find_bestpair_total_modern"] += fbc.total_modern
            stats["find_bestpair_missing_in_modern"] += len(fbc.missing_in_modern)
            stats["find_bestpair_extra_in_modern"] += len(fbc.extra_in_modern)
        
        # Base pair comparison statistics (pairs that pass all checks from calculate_more_bppars)
        if hasattr(result, 'base_pair_comparison') and result.base_pair_comparison:
            bpc = result.base_pair_comparison
            stats["base_pair_total_legacy"] += bpc.total_legacy
            stats["base_pair_total_modern"] += bpc.total_modern
            stats["base_pair_missing_in_modern"] += len(bpc.missing_in_modern)
            stats["base_pair_extra_in_modern"] += len(bpc.extra_in_modern)
            stats["base_pair_mismatched"] += len(bpc.mismatched_pairs)
        
        # Pair validation comparison statistics (SECONDARY - for debugging)
        if hasattr(result, 'pair_validation_comparison') and result.pair_validation_comparison:
            pvc = result.pair_validation_comparison
            stats["pair_validation_total_legacy"] += pvc.total_legacy
            stats["pair_validation_total_modern"] += pvc.total_modern
            stats["pair_validation_missing_in_modern"] += len(pvc.missing_in_modern)
            stats["pair_validation_extra_in_modern"] += len(pvc.extra_in_modern)
            stats["pair_validation_mismatched"] += len(pvc.mismatched_validations)
        
        # Distance checks comparison statistics
        if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
            dcc = result.distance_checks_comparison
            stats["distance_checks_total_legacy"] += dcc.total_legacy
            stats["distance_checks_total_modern"] += dcc.total_modern
            stats["distance_checks_missing_in_modern"] += len(dcc.missing_in_modern)
            stats["distance_checks_extra_in_modern"] += len(dcc.extra_in_modern)
            stats["distance_checks_mismatched"] += len(dcc.mismatched_checks)
        
        # H-bond list comparison statistics
        if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
            hlc = result.hbond_list_comparison
            stats["hbond_list_total_legacy"] += hlc.total_legacy
            stats["hbond_list_total_modern"] += hlc.total_modern
            stats["hbond_list_missing_in_modern"] += len(hlc.missing_in_modern)
            stats["hbond_list_extra_in_modern"] += len(hlc.extra_in_modern)
            stats["hbond_list_mismatched_pairs"] += len(hlc.mismatched_pairs)
        
        # Residue indices comparison statistics
        if hasattr(result, 'residue_indices_comparison') and result.residue_indices_comparison:
            ric = result.residue_indices_comparison
            if ric.missing_in_modern:
                stats["residue_indices_missing_in_modern"] += 1
            if ric.extra_in_modern:
                stats["residue_indices_extra_in_modern"] += 1
            if not ric.num_residue_match:
                stats["residue_indices_num_residue_match"] = False
            stats["residue_indices_mismatched_entries"] += len(ric.mismatched_entries)

    return stats


def generate_verbose_report(
    pdb_id: str,
    result: ComparisonResult,
    project_root: Path,
    tolerance: float = 1e-6
) -> str:
    """Generate a detailed verbose report for a single PDB."""
    reporter = VerboseReporter(
        tolerance=tolerance,
        show_provenance=False,
        show_related=True,
        diff_only=False
    )
    
    # Add header
    stages = []
    if result.atom_comparison:
        stages.append("atoms")
    if result.frame_comparison:
        stages.append("frames")
    if result.step_comparison:
        stages.append("steps")
    if result.helical_comparison:
        stages.append("helical")
    if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
        stages.append("distance_checks")
    if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
        stages.append("hbond_list")
    if hasattr(result, 'pair_validation_comparison') and result.pair_validation_comparison:
        stages.append("pair_validation")
    if hasattr(result, 'find_bestpair_comparison') and result.find_bestpair_comparison:
        stages.append("find_bestpair_selection")
    if hasattr(result, 'base_pair_comparison') and result.base_pair_comparison:
        stages.append("base_pair")
    
    reporter.add_header(pdb_id, stages)
    
    # Add distance checks comparison (if available)
    if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
        dcc = result.distance_checks_comparison
        reporter.add_stage_header("distance_checks", 3)
        reporter.add_stage_summary(
            dcc.total_legacy,
            dcc.total_modern,
            dcc.total_legacy - len(dcc.missing_in_modern),
            len(dcc.missing_in_modern),
            len(dcc.extra_in_modern),
            len(dcc.mismatched_checks)
        )
        
        # Show mismatched checks
        for mismatch in dcc.mismatched_checks[:reporter.max_mismatches_per_stage]:
            pair_key = (mismatch.get('base_i'), mismatch.get('base_j'))
            legacy_rec = mismatch.get('legacy_record', {})
            modern_rec = mismatch.get('modern_record', {})
            
            fields = ['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area']
            
            rec_comp = create_record_comparison_from_dicts(
                record_key=pair_key,
                record_type="distance_checks",
                legacy_record=legacy_rec,
                modern_record=modern_rec,
                fields_to_compare=fields,
                tolerance=tolerance,
                legacy_source=f"data/json_legacy/distance_checks/{pdb_id}.json",
                modern_source=f"data/json/distance_checks/{pdb_id}.json"
            )
            reporter.add_record_comparison(rec_comp)
    
    # Add h-bond list comparison (if available)
    if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
        hlc = result.hbond_list_comparison
        reporter.add_stage_header("hbond_list", 4)
        reporter.add_stage_summary(
            hlc.total_legacy,
            hlc.total_modern,
            hlc.total_legacy - len(hlc.missing_in_modern),
            len(hlc.missing_in_modern),
            len(hlc.extra_in_modern),
            len(hlc.mismatched_pairs)
        )
        
        # Show mismatched pairs
        for mismatch in hlc.mismatched_pairs[:reporter.max_mismatches_per_stage]:
            pair_key = (mismatch.get('base_i'), mismatch.get('base_j'))
            legacy_rec = mismatch.get('legacy_record', {})
            modern_rec = mismatch.get('modern_record', {})
            
            fields = ['num_hbonds']
            
            rec_comp = create_record_comparison_from_dicts(
                record_key=pair_key,
                record_type="hbond_list",
                legacy_record=legacy_rec,
                modern_record=modern_rec,
                fields_to_compare=fields,
                tolerance=tolerance,
                legacy_source=f"data/json_legacy/hbond_list/{pdb_id}.json",
                modern_source=f"data/json/hbond_list/{pdb_id}.json"
            )
            reporter.add_record_comparison(rec_comp)
    
    # Add frame comparison (if available)
    if result.frame_comparison:
        fc = result.frame_comparison
        reporter.add_stage_header("frames", 2)
        reporter.add_stage_summary(
            fc.total_legacy,
            fc.total_modern,
            fc.total_legacy - len(fc.missing_residues),
            len(fc.missing_residues),
            0,  # No extra in modern for frames
            len(fc.mismatched_calculations)
        )
        
        # Show mismatched calculations
        for mismatch in fc.mismatched_calculations[:reporter.max_mismatches_per_stage]:
            residue_key = mismatch.residue_key
            legacy_rec = mismatch.legacy_record
            modern_rec = mismatch.modern_record
            
            # Determine which type of frame record this is
            record_type = "unknown"
            mismatches_dict = mismatch.mismatches
            if 'base_frame_calc' in mismatches_dict:
                record_type = "base_frame_calc"
                fields = ['rms_fit', 'num_matched_atoms', 'base_type']
            elif 'ls_fitting' in mismatches_dict:
                record_type = "ls_fitting"
                fields = ['rms_fit', 'num_points']
            elif 'frame_calc' in mismatches_dict:
                record_type = "frame_calc"
                fields = ['rms_fit', 'num_matched_atoms']
            else:
                fields = ['rms_fit']
            
            rec_comp = create_record_comparison_from_dicts(
                record_key=residue_key,
                record_type=record_type,
                legacy_record=legacy_rec,
                modern_record=modern_rec,
                fields_to_compare=fields,
                tolerance=tolerance,
                legacy_source=f"data/json_legacy/{record_type}/{pdb_id}.json",
                modern_source=f"data/json/{record_type}/{pdb_id}.json"
            )
            reporter.add_record_comparison(rec_comp)
    
    # Add summary
    stages_compared = len(stages)
    stages_with_diffs = 0
    diff_details = []
    
    if result.frame_comparison and len(result.frame_comparison.mismatched_calculations) > 0:
        stages_with_diffs += 1
        diff_details.append(f"frames: {len(result.frame_comparison.mismatched_calculations)} mismatches")
    
    if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
        if len(result.distance_checks_comparison.mismatched_checks) > 0:
            stages_with_diffs += 1
            diff_details.append(f"distance_checks: {len(result.distance_checks_comparison.mismatched_checks)} mismatches")
    
    if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
        if len(result.hbond_list_comparison.mismatched_pairs) > 0:
            stages_with_diffs += 1
            diff_details.append(f"hbond_list: {len(result.hbond_list_comparison.mismatched_pairs)} mismatches")
    
    reporter.add_summary(
        stages_compared=stages_compared,
        perfect_matches=stages_compared - stages_with_diffs,
        stages_with_diffs=stages_with_diffs,
        diff_details=diff_details
    )
    
    return reporter.generate_report()


def generate_report(
    stats: Dict,
    results: Dict[str, ComparisonResult],
    verbose: bool = False,
    diff_only: bool = False,
) -> str:
    """Generate a detailed human-readable report."""
    lines = []

    lines.append("=" * 80)
    lines.append("JSON COMPARISON REPORT")
    lines.append("=" * 80)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    lines.append("SUMMARY")
    lines.append("-" * 80)
    lines.append(f"Total PDB files compared: {stats['total_pdbs']}")
    lines.append(f"  ✅ Perfect matches: {stats['match_count']}")
    lines.append(f"  ⚠️  Files with differences: {stats['diff_count']}")
    lines.append(f"  ❌ Errors: {stats['error_count']}")
    lines.append("")
    
    # Display errors if any
    if stats['error_count'] > 0:
        lines.append("=" * 80)
        lines.append("ERRORS")
        lines.append("=" * 80)
        for pdb_id, result in results.items():
            if result.status == "error" and result.errors:
                lines.append(f"\n{pdb_id}:")
                for error in result.errors:
                    lines.append(f"  ❌ {error}")
        lines.append("")

    # Atom comparison statistics
    if stats["atom_total_legacy"] > 0 or stats["atom_total_modern"] > 0:
        lines.append("ATOM COMPARISON STATISTICS")
        lines.append("-" * 80)
        lines.append(f"Total legacy atoms: {stats['atom_total_legacy']}")
        lines.append(f"Total modern atoms: {stats['atom_total_modern']}")
        lines.append(f"Common atoms: {stats['atom_common_count']}")
        if stats["atom_count_difference"]:
            lines.append(f"⚠️  Atom count difference detected")
        if stats["atom_missing_in_modern"] > 0:
            lines.append(f"Missing in modern: {stats['atom_missing_in_modern']}")
        if stats["atom_extra_in_modern"] > 0:
            lines.append(f"Extra in modern: {stats['atom_extra_in_modern']}")
        if stats["atom_mismatched_fields"] > 0:
            lines.append(f"Mismatched fields: {stats['atom_mismatched_fields']}")
        lines.append("")

    if stats["total_residues"] > 0:
        lines.append("FRAME CALCULATION STATISTICS")
        lines.append("-" * 80)
        lines.append(f"Total residues: {stats['total_residues']}")
        lines.append(f"Missing residues: {stats['missing_residues']}")
        lines.append(f"Mismatched calculations: {stats['mismatched_calculations']}")
        lines.append("")
    
    # Step parameter comparison statistics (bpstep_params)
    if stats["step_total_legacy"] > 0 or stats["step_total_modern"] > 0:
        lines.append("STEP PARAMETER STATISTICS (bpstep_params)")
        lines.append("-" * 80)
        lines.append(f"Total legacy steps: {stats['step_total_legacy']}")
        lines.append(f"Total modern steps: {stats['step_total_modern']}")
        lines.append(f"Missing steps: {stats['step_missing_steps']}")
        lines.append(f"Mismatched steps: {stats['step_mismatched_steps']}")
        lines.append("")
    
    # Helical parameter comparison statistics (helical_params)
    if stats["helical_total_legacy"] > 0 or stats["helical_total_modern"] > 0:
        lines.append("HELICAL PARAMETER STATISTICS (helical_params)")
        lines.append("-" * 80)
        lines.append(f"Total legacy helical: {stats['helical_total_legacy']}")
        lines.append(f"Total modern helical: {stats['helical_total_modern']}")
        lines.append(f"Missing helical: {stats['helical_missing_steps']}")
        lines.append(f"Mismatched helical: {stats['helical_mismatched_steps']}")
        lines.append("")
    
    # Find bestpair selection statistics (PRIMARY - actual selected pairs)
    if stats["find_bestpair_total_legacy"] > 0 or stats["find_bestpair_total_modern"] > 0:
        lines.append("FIND_BESTPAIR SELECTION STATISTICS (Actual Selected Pairs)")
        lines.append("-" * 80)
        lines.append(f"Total legacy selected pairs: {stats['find_bestpair_total_legacy']}")
        lines.append(f"Total modern selected pairs: {stats['find_bestpair_total_modern']}")
        lines.append(f"Common pairs: {stats['find_bestpair_total_legacy'] - stats['find_bestpair_missing_in_modern']}")
        lines.append(f"Missing in modern: {stats['find_bestpair_missing_in_modern']}")
        lines.append(f"Extra in modern: {stats['find_bestpair_extra_in_modern']}")
        lines.append("")
    
    # Base pair comparison statistics (pairs that pass all checks from calculate_more_bppars)
    if stats["base_pair_total_legacy"] > 0 or stats["base_pair_total_modern"] > 0:
        lines.append("BASE PAIR STATISTICS (Pairs Passing All Checks)")
        lines.append("-" * 80)
        lines.append(f"Total legacy base pairs: {stats['base_pair_total_legacy']}")
        lines.append(f"Total modern base pairs: {stats['base_pair_total_modern']}")
        lines.append(f"Common pairs: {stats['base_pair_total_legacy'] - stats['base_pair_missing_in_modern']}")
        lines.append(f"Missing in modern: {stats['base_pair_missing_in_modern']}")
        lines.append(f"Extra in modern: {stats['base_pair_extra_in_modern']}")
        lines.append(f"Mismatched pairs: {stats['base_pair_mismatched']}")
        lines.append("")
    
    # Pair validation comparison statistics (SECONDARY - for debugging)
    if stats["pair_validation_total_legacy"] > 0 or stats["pair_validation_total_modern"] > 0:
        lines.append("PAIR VALIDATION STATISTICS (All Validated Pairs)")
        lines.append("-" * 80)
        lines.append(f"Total legacy validations: {stats['pair_validation_total_legacy']}")
        lines.append(f"Total modern validations: {stats['pair_validation_total_modern']}")
        lines.append(f"Missing in modern: {stats['pair_validation_missing_in_modern']}")
        lines.append(f"Extra in modern: {stats['pair_validation_extra_in_modern']}")
        lines.append(f"Mismatched validations: {stats['pair_validation_mismatched']}")
        lines.append("")
    
    # Distance checks comparison statistics
    if stats["distance_checks_total_legacy"] > 0 or stats["distance_checks_total_modern"] > 0:
        lines.append("DISTANCE CHECKS STATISTICS")
        lines.append("-" * 80)
        lines.append(f"Total legacy checks: {stats['distance_checks_total_legacy']}")
        lines.append(f"Total modern checks: {stats['distance_checks_total_modern']}")
        lines.append(f"Missing in modern: {stats['distance_checks_missing_in_modern']}")
        lines.append(f"Extra in modern: {stats['distance_checks_extra_in_modern']}")
        lines.append(f"Mismatched checks: {stats['distance_checks_mismatched']}")
        lines.append("")
    
    # H-bond list comparison statistics
    if stats["hbond_list_total_legacy"] > 0 or stats["hbond_list_total_modern"] > 0:
        lines.append("H-BOND LIST STATISTICS")
        lines.append("-" * 80)
        lines.append(f"Total legacy H-bond lists: {stats['hbond_list_total_legacy']}")
        lines.append(f"Total modern H-bond lists: {stats['hbond_list_total_modern']}")
        lines.append(f"Common pairs: {stats['hbond_list_total_legacy'] - stats['hbond_list_missing_in_modern']}")
        lines.append(f"Missing in modern: {stats['hbond_list_missing_in_modern']}")
        lines.append(f"Extra in modern: {stats['hbond_list_extra_in_modern']}")
        lines.append(f"Mismatched pairs: {stats['hbond_list_mismatched_pairs']}")
        lines.append("")
    
    # Residue indices comparison statistics
    # Show statistics if residue indices comparisons were performed
    # Check if any results have residue_indices_comparison set (not None means comparison was performed)
    has_residue_indices_comparison = False
    total_residue_indices_pdbs = 0
    for result in results.values():
        if result.residue_indices_comparison is not None:
            has_residue_indices_comparison = True
            total_residue_indices_pdbs += 1
    
    if has_residue_indices_comparison:
        lines.append("RESIDUE INDICES STATISTICS")
        lines.append("-" * 80)
        
        # Show summary even if everything matches
        if (stats["residue_indices_missing_in_modern"] == 0 and 
            stats["residue_indices_extra_in_modern"] == 0 and
            stats["residue_indices_num_residue_match"] and
            stats["residue_indices_mismatched_entries"] == 0):
            lines.append(f"✅ All residue indices match perfectly ({total_residue_indices_pdbs} PDBs compared)")
        else:
            if stats["residue_indices_missing_in_modern"] > 0:
                lines.append(f"Missing in modern: {stats['residue_indices_missing_in_modern']}")
            if stats["residue_indices_extra_in_modern"] > 0:
                lines.append(f"Extra in modern: {stats['residue_indices_extra_in_modern']}")
            if not stats["residue_indices_num_residue_match"]:
                lines.append("⚠️  Residue count mismatch detected")
            if stats["residue_indices_mismatched_entries"] > 0:
                lines.append(f"Mismatched entries: {stats['residue_indices_mismatched_entries']}")
        lines.append("")
    
    matching_residues = stats["total_residues"] - stats["missing_residues"]
    if matching_residues > 0:
            exact_rate = stats["exact_atom_matches"] / matching_residues * 100
            set_rate = stats["atom_set_matches"] / matching_residues * 100
            real_diff_rate = stats["real_atom_differences"] / matching_residues * 100

            lines.append("ATOM MATCHING STATISTICS")
            lines.append("-" * 80)
            lines.append(
                f"Exact atom matches: {stats['exact_atom_matches']}/{matching_residues} "
                f"({exact_rate:.1f}%)"
            )
            lines.append(
                f"Atom set matches: {stats['atom_set_matches']}/{matching_residues} "
                f"({set_rate:.1f}%)"
            )
            lines.append(
                f"Real atom differences: {stats['real_atom_differences']}/{matching_residues} "
                f"({real_diff_rate:.1f}%)"
            )
            lines.append("")
    else:
        # Check if atom comparison was skipped
        has_atom_comparison = any(
            result.atom_comparison is not None
            for result in results.values()
            if hasattr(result, "atom_comparison")
        )
        if not has_atom_comparison:
            # Atom comparison was skipped, don't show atom matching stats
            pass

        if stats["rms_differences"]:
            avg_rms = sum(stats["rms_differences"]) / len(stats["rms_differences"])
            max_rms = max(stats["rms_differences"])
            lines.append("RMS FIT DIFFERENCES")
            lines.append("-" * 80)
            lines.append(
                f"Residues with RMS differences: {len(stats['rms_differences'])}"
            )
            lines.append(f"Average RMS difference: {avg_rms:.6f} Å")
            lines.append(f"Maximum RMS difference: {max_rms:.6f} Å")
            lines.append("")

    if stats["base_type_stats"]:
        lines.append("BASE TYPE STATISTICS")
        lines.append("-" * 80)
        for base_type in sorted(stats["base_type_stats"].keys()):
            base_stats = stats["base_type_stats"][base_type]
            if base_stats["total"] > 0 or base_stats["missing"] > 0:
                lines.append(f"{base_type}: {base_stats['total']} residues")
                if base_stats["total"] > 0:
                    exact_rate = base_stats["exact_matches"] / base_stats["total"] * 100
                    set_rate = base_stats["set_matches"] / base_stats["total"] * 100
                    diff_rate = base_stats["real_diffs"] / base_stats["total"] * 100
                    lines.append(
                        f"  Exact matches: {base_stats['exact_matches']} ({exact_rate:.1f}%)"
                    )
                    lines.append(
                        f"  Set matches: {base_stats['set_matches']} ({set_rate:.1f}%)"
                    )
                    lines.append(
                        f"  Real differences: {base_stats['real_diffs']} ({diff_rate:.1f}%)"
                    )
                if base_stats["missing"] > 0:
                    lines.append(f"  Missing: {base_stats['missing']}")
                lines.append("")

    if stats["diff_count"] > 0:
        lines.append("=" * 80)
        lines.append("FILES WITH MOST DIFFERENCES (Top 20)")
        lines.append("=" * 80)
        lines.append("")

        file_diffs = []
        for pdb_id in stats["files_with_differences"]:
            result = results[pdb_id]
            if result.status == "error" or not result.has_differences():
                continue

            diff_count = 0
            if result.frame_comparison:
                diff_count = len(result.frame_comparison.missing_residues) + len(
                    result.frame_comparison.mismatched_calculations
                )

            if diff_count > 0:
                file_diffs.append((pdb_id, diff_count))

        file_diffs.sort(key=lambda x: x[1], reverse=True)
        for pdb_id, diff_count in file_diffs[:20]:
            lines.append(f"  {pdb_id}: {diff_count} differences")
        lines.append("")

    if verbose and diff_only and len(results) == 1:
        pdb_id = list(results.keys())[0]
        result = results[pdb_id]
        if result.frame_comparison:
            fc = result.frame_comparison

            if fc.missing_residues:
                lines.append("=" * 80)
                lines.append("MISSING RESIDUES")
                lines.append("=" * 80)
                for missing in fc.missing_residues[:20]:
                    idx_str = f" [legacy_residue_idx={missing.legacy_residue_idx}]" if missing.legacy_residue_idx else ""
                    lines.append(
                        f"  {missing.base_type} {missing.chain_id}:{missing.residue_seq}"
                        f"{missing.insertion}{idx_str}"
                    )
                if len(fc.missing_residues) > 20:
                    lines.append(f"  ... and {len(fc.missing_residues) - 20} more")
                lines.append("")

            if fc.mismatched_calculations:
                lines.append("=" * 80)
                lines.append("MISMATCHED CALCULATIONS (First 20)")
                lines.append("=" * 80)
                for mismatch in fc.mismatched_calculations[:20]:
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    base_type = mismatch.legacy_record.get("base_type", "?")
                    leg_atoms = mismatch.legacy_matched_atoms
                    mod_atoms = mismatch.modern_matched_atoms
                    # Show legacy_residue_idx if available
                    legacy_residue_idx = mismatch.modern_record.get('legacy_residue_idx')
                    idx_str = f" [legacy_residue_idx={legacy_residue_idx}]" if legacy_residue_idx else ""

                    # Handle nested structure for record type (base_frame_calc, ls_fitting, or frame_calc)
                    mismatches_dict = mismatch.mismatches
                    record_type = None
                    if 'base_frame_calc' in mismatches_dict:
                        record_type = 'base_frame_calc'
                        mismatches_dict = mismatches_dict['base_frame_calc']
                    elif 'ls_fitting' in mismatches_dict:
                        record_type = 'ls_fitting'
                        mismatches_dict = mismatches_dict['ls_fitting']
                    elif 'frame_calc' in mismatches_dict:
                        record_type = 'frame_calc'
                        mismatches_dict = mismatches_dict['frame_calc']
                    
                    type_str = f" [{record_type}]" if record_type else ""
                    lines.append(
                        f"\n  {base_type} {chain_id}:{residue_seq}{insertion}:{idx_str}{type_str}"
                    )
                    if "rms" in mismatches_dict:
                        lines.append(
                            f"    RMS: legacy={mismatches_dict['rms']['legacy']:.6f}, "
                            f"modern={mismatches_dict['rms']['modern']:.6f}"
                        )
                    if "matched_atoms" in mismatches_dict:
                        lines.append(
                            f"    Legacy only: {mismatches_dict['matched_atoms']['only_legacy']}"
                        )
                        lines.append(
                            f"    Modern only: {mismatches_dict['matched_atoms']['only_modern']}"
                        )
                    if "residue_idx" in mismatches_dict:
                        lines.append(
                            f"    residue_idx: legacy={mismatches_dict['residue_idx']['legacy']}, "
                            f"modern={mismatches_dict['residue_idx']['modern']}"
                        )
                    if "legacy_residue_idx_mismatch" in mismatches_dict:
                        mismatch_info = mismatches_dict['legacy_residue_idx_mismatch']
                        lines.append(
                            f"    legacy_residue_idx mismatch: {mismatch_info.get('note', '')} "
                            f"(legacy={mismatch_info.get('legacy_residue_idx')}, "
                            f"modern={mismatch_info.get('modern_legacy_residue_idx')})"
                        )
                    if "base_type" in mismatches_dict:
                        lines.append(
                            f"    base_type: legacy='{mismatches_dict['base_type']['legacy']}', "
                            f"modern='{mismatches_dict['base_type']['modern']}'"
                        )
                    if "residue_name" in mismatches_dict:
                        lines.append(
                            f"    residue_name: legacy='{mismatches_dict['residue_name']['legacy']}', "
                            f"modern='{mismatches_dict['residue_name']['modern']}'"
                        )
                    if "num_points" in mismatches_dict:
                        lines.append(
                            f"    num_points: legacy={mismatches_dict['num_points']['legacy']}, "
                            f"modern={mismatches_dict['num_points']['modern']}"
                        )
                    if "num_matched_atoms" in mismatches_dict:
                        lines.append(
                            f"    num_matched_atoms: legacy={mismatches_dict['num_matched_atoms']['legacy']}, "
                            f"modern={mismatches_dict['num_matched_atoms']['modern']}"
                        )
                    if "rotation_matrix" in mismatches_dict:
                        rot_info = mismatches_dict['rotation_matrix']
                        if 'max_diff' in rot_info:
                            lines.append(
                                f"    rotation_matrix: max_diff={rot_info['max_diff']:.6e}"
                            )
                        else:
                            lines.append("    rotation_matrix: dimension mismatch")
                    if "translation" in mismatches_dict:
                        trans_info = mismatches_dict['translation']
                        if 'max_diff' in trans_info:
                            lines.append(
                                f"    translation: max_diff={trans_info['max_diff']:.6e}"
                            )
                        else:
                            lines.append("    translation: dimension mismatch")
                    if "template_file" in mismatches_dict:
                        lines.append(
                            f"    template_file: legacy={mismatches_dict['template_file']['legacy']}, "
                            f"modern={mismatches_dict['template_file']['modern']}"
                        )
                    if "matched_coordinates" in mismatches_dict:
                        coord_info = mismatches_dict['matched_coordinates']
                        if 'note' in coord_info:
                            lines.append(
                                f"    matched_coordinates: {coord_info['note']} "
                                f"(legacy={coord_info['legacy_count']}, modern={coord_info['modern_count']})"
                            )
                        elif 'mismatched_pairs' in coord_info:
                            lines.append(
                                f"    matched_coordinates: {coord_info['mismatched_pairs']}/{coord_info['total_pairs']} pairs differ"
                            )
                            if coord_info['differences']:
                                for diff in coord_info['differences'][:3]:  # Show first 3
                                    lines.append(
                                        f"      Pair {diff['index']} (atom_idx={diff['atom_idx']}): {list(diff['mismatches'].keys())}"
                                    )
                                if len(coord_info['differences']) > 3:
                                    lines.append(f"      ... and {len(coord_info['differences']) - 3} more")
                    elif set(leg_atoms) != set(mod_atoms):
                        leg_only = sorted(set(leg_atoms) - set(mod_atoms))
                        mod_only = sorted(set(mod_atoms) - set(leg_atoms))
                        if leg_only:
                            lines.append(f"    Legacy only: {leg_only}")
                        if mod_only:
                            lines.append(f"    Modern only: {mod_only}")

                    # Display atom information with legacy indices
                    if mismatch.atom_pdb_lines:
                        lines.append("    Atoms:")
                        # Create lookup for legacy atom indices
                        atom_idx_map = {
                            info.atom_name: info.legacy_atom_idx
                            for info in mismatch.atom_pdb_lines
                            if info.legacy_atom_idx is not None
                        }

                        # Show legacy atoms with indices
                        if leg_atoms:
                            lines.append("      Legacy matched atoms:")
                            for atom_name in sorted(leg_atoms):
                                atom_idx = atom_idx_map.get(atom_name, "?")
                                lines.append(f"        {atom_name} (idx={atom_idx})")

                        # Show modern atoms with indices
                        if mod_atoms:
                            lines.append("      Modern matched atoms:")
                            for atom_name in sorted(mod_atoms):
                                atom_idx = atom_idx_map.get(atom_name, "?")
                                lines.append(f"        {atom_name} (idx={atom_idx})")
                lines.append("")

    # Atom comparison details (when verbose mode is enabled)
    if verbose and len(results) == 1:
        pdb_id = list(results.keys())[0]
        result = results[pdb_id]
        if result.atom_comparison:
            ac = result.atom_comparison
            if ac.missing_in_modern or ac.extra_in_modern or ac.mismatched_fields:
                lines.append("=" * 80)
                lines.append("ATOM COMPARISON DETAILS")
                lines.append("=" * 80)

                if ac.missing_in_modern:
                    lines.append(
                        f"\nMISSING IN MODERN ({len(ac.missing_in_modern)} atoms):"
                    )
                    for atom_info in ac.missing_in_modern[:50]:  # Limit to 50
                        idx_str = (
                            f"idx={atom_info.legacy_atom_idx}"
                            if atom_info.legacy_atom_idx
                            else ""
                        )
                        lines.append(
                            f"  {atom_info.atom_name} {atom_info.chain_id}:{atom_info.residue_seq}"
                            f"{atom_info.insertion} {atom_info.residue_name} {idx_str}"
                        )
                        if atom_info.pdb_line:
                            line_num = (
                                f"Line {atom_info.line_number}"
                                if atom_info.line_number
                                else ""
                            )
                            lines.append(
                                f"    {line_num}: {atom_info.pdb_line.strip()}"
                            )
                    if len(ac.missing_in_modern) > 50:
                        lines.append(f"  ... and {len(ac.missing_in_modern) - 50} more")
                    lines.append("")

                if ac.extra_in_modern:
                    lines.append(
                        f"\nEXTRA IN MODERN ({len(ac.extra_in_modern)} atoms):"
                    )
                    for atom_info in ac.extra_in_modern[:50]:  # Limit to 50
                        idx_str = (
                            f"idx={atom_info.legacy_atom_idx}"
                            if atom_info.legacy_atom_idx
                            else ""
                        )
                        lines.append(
                            f"  {atom_info.atom_name} {atom_info.chain_id}:{atom_info.residue_seq}"
                            f"{atom_info.insertion} {atom_info.residue_name} {idx_str}"
                        )
                        if atom_info.pdb_line:
                            line_num = (
                                f"Line {atom_info.line_number}"
                                if atom_info.line_number
                                else ""
                            )
                            lines.append(
                                f"    {line_num}: {atom_info.pdb_line.strip()}"
                            )
                    if len(ac.extra_in_modern) > 50:
                        lines.append(f"  ... and {len(ac.extra_in_modern) - 50} more")
                    lines.append("")

                if ac.mismatched_fields:
                    lines.append(
                        f"\nMISMATCHED FIELDS ({len(ac.mismatched_fields)} atoms):"
                    )
                    for mismatch in ac.mismatched_fields[:50]:  # Limit to 50
                        chain_id, residue_seq, insertion, atom_name = mismatch.atom_key
                        lines.append(
                            f"  {atom_name} {chain_id}:{residue_seq}{insertion}"
                        )
                        lines.append(f"    Field: {mismatch.field_name}")
                        if isinstance(mismatch.legacy_value, dict):
                            for field, values in mismatch.legacy_value.items():
                                lines.append(
                                    f"      {field}: legacy={values.get('legacy')}, modern={values.get('modern')}"
                                )
                        else:
                            lines.append(f"      Legacy: {mismatch.legacy_value}")
                            lines.append(f"      Modern: {mismatch.modern_value}")
                        if mismatch.legacy_pdb_line:
                            line_num = (
                                f"Line {mismatch.legacy_line_number}"
                                if mismatch.legacy_line_number
                                else ""
                            )
                            lines.append(
                                f"    Legacy PDB {line_num}: {mismatch.legacy_pdb_line.strip()}"
                            )
                        if mismatch.modern_pdb_line:
                            line_num = (
                                f"Line {mismatch.modern_line_number}"
                                if mismatch.modern_line_number
                                else ""
                            )
                            lines.append(
                                f"    Modern PDB {line_num}: {mismatch.modern_pdb_line.strip()}"
                            )
                        lines.append("")
                    if len(ac.mismatched_fields) > 50:
                        lines.append(f"  ... and {len(ac.mismatched_fields) - 50} more")
                    lines.append("")

    # Step parameter comparison details (bpstep_params) - when verbose mode is enabled
    if verbose and len(results) == 1:
        pdb_id = list(results.keys())[0]
        result = results[pdb_id]
        if result.step_comparison:
            sc = result.step_comparison
            if sc.missing_steps or sc.mismatched_steps:
                lines.append("=" * 80)
                lines.append("STEP PARAMETER COMPARISON DETAILS (bpstep_params)")
                lines.append("=" * 80)
                
                if sc.missing_steps:
                    lines.append(f"\nMISSING IN MODERN ({len(sc.missing_steps)} steps):")
                    for missing in sc.missing_steps[:20]:
                        bp_idx1 = missing.get('bp_idx1', '?')
                        bp_idx2 = missing.get('bp_idx2', '?')
                        lines.append(f"  Step bp_idx1={bp_idx1}, bp_idx2={bp_idx2}")
                    if len(sc.missing_steps) > 20:
                        lines.append(f"  ... and {len(sc.missing_steps) - 20} more")
                    lines.append("")
                
                if sc.mismatched_steps:
                    lines.append(f"\nMISMATCHED STEPS ({len(sc.mismatched_steps)} steps):")
                    for mismatch in sc.mismatched_steps[:20]:
                        bp_idx1 = mismatch.get('bp_idx1', '?')
                        bp_idx2 = mismatch.get('bp_idx2', '?')
                        mismatches = mismatch.get('mismatches', {})
                        lines.append(f"\n  Step bp_idx1={bp_idx1}, bp_idx2={bp_idx2}:")
                        for field, diff_info in list(mismatches.items())[:10]:
                            leg_val = diff_info.get('legacy', '?')
                            mod_val = diff_info.get('modern', '?')
                            diff = diff_info.get('diff', 0.0)
                            lines.append(f"    {field}: legacy={leg_val}, modern={mod_val}, diff={diff:.6e}")
                        if len(mismatches) > 10:
                            lines.append(f"    ... and {len(mismatches) - 10} more parameter differences")
                    if len(sc.mismatched_steps) > 20:
                        lines.append(f"  ... and {len(sc.mismatched_steps) - 20} more mismatched steps")
                    lines.append("")
        
        # Helical parameter comparison details (helical_params) - when verbose mode is enabled
        if result.helical_comparison:
            hc = result.helical_comparison
            if hc.missing_steps or hc.mismatched_steps:
                lines.append("=" * 80)
                lines.append("HELICAL PARAMETER COMPARISON DETAILS (helical_params)")
                lines.append("=" * 80)
                
                if hc.missing_steps:
                    lines.append(f"\nMISSING IN MODERN ({len(hc.missing_steps)} helical steps):")
                    for missing in hc.missing_steps[:20]:
                        bp_idx1 = missing.get('bp_idx1', '?')
                        bp_idx2 = missing.get('bp_idx2', '?')
                        lines.append(f"  Helical step bp_idx1={bp_idx1}, bp_idx2={bp_idx2}")
                    if len(hc.missing_steps) > 20:
                        lines.append(f"  ... and {len(hc.missing_steps) - 20} more")
                    lines.append("")
                
                if hc.mismatched_steps:
                    lines.append(f"\nMISMATCHED HELICAL STEPS ({len(hc.mismatched_steps)} steps):")
                    for mismatch in hc.mismatched_steps[:20]:
                        bp_idx1 = mismatch.get('bp_idx1', '?')
                        bp_idx2 = mismatch.get('bp_idx2', '?')
                        mismatches = mismatch.get('mismatches', {})
                        lines.append(f"\n  Helical step bp_idx1={bp_idx1}, bp_idx2={bp_idx2}:")
                        for field, diff_info in list(mismatches.items())[:10]:
                            leg_val = diff_info.get('legacy', '?')
                            mod_val = diff_info.get('modern', '?')
                            diff = diff_info.get('diff', 0.0)
                            lines.append(f"    {field}: legacy={leg_val}, modern={mod_val}, diff={diff:.6e}")
                        if len(mismatches) > 10:
                            lines.append(f"    ... and {len(mismatches) - 10} more parameter differences")
                    if len(hc.mismatched_steps) > 20:
                        lines.append(f"  ... and {len(hc.mismatched_steps) - 20} more mismatched helical steps")
                    lines.append("")

    return "\n".join(lines)


# Shared options
def common_options(f):
    """Common options for comparison commands."""
    f = click.option(
        "--legacy-mode", is_flag=True, help="Use modern JSON files with _legacy suffix"
    )(f)
    f = click.option(
        "--output", "-o", type=click.Path(path_type=Path), help="Save report to file"
    )(f)
    f = click.option(
        "--threads",
        "-t",
        type=int,
        default=None,
        help="Number of threads (default: CPU count)",
    )(f)
    f = click.option(
        "--regenerate",
        "-r",
        is_flag=True,
        help="Regenerate JSON files if missing (legacy from org code, modern from new code)",
    )(f)
    f = click.option(
        "--test-set",
        type=click.Choice(['10', '50', '100', '500', '1000'], case_sensitive=False),
        help="Use a saved test set of the specified size (10, 50, 100, 500, or 1000 PDBs)",
    )(f)
    f = click.option(
        "--config",
        type=click.Path(path_type=Path, exists=True),
        default=None,
        help="Path to YAML configuration file (default: comparison_config.yaml in project root)",
    )(f)
    return f


def _compare_pdb_worker(args):
    """Worker function for parallel processing (needs to be picklable)."""
    pdb_id, project_root_str, use_legacy_mode, regenerate = args
    project_root = Path(project_root_str)
    
    # Create new comparator for this process (loads config automatically)
    comparator = create_comparator(None, project_root)
    
    return compare_single_pdb(
        pdb_id, project_root, use_legacy_mode, comparator, regenerate
    )


def run_comparison(
    pdb_ids: List[str],
    project_root: Path,
    use_legacy_mode: bool,
    comparator: JsonComparator,
    max_workers: int,
    regenerate: bool = False,
) -> Dict[str, ComparisonResult]:
    """Run comparison for given PDB IDs."""
    results = {}

    if len(pdb_ids) == 1:
        result = compare_single_pdb(
            pdb_ids[0], project_root, use_legacy_mode, comparator, regenerate
        )
        results[pdb_ids[0]] = result
    else:
        completed = 0
        # Use ProcessPoolExecutor for CPU-bound work (better than threads due to GIL)
        # Processes allow true parallelism for CPU-bound JSON parsing/comparison
        # Prepare arguments for worker function (must be picklable)
        worker_args = [
            (
                pdb_id,
                str(project_root),  # Convert Path to string for pickling
                use_legacy_mode,
                regenerate,
            )
            for pdb_id in pdb_ids
        ]
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {
                executor.submit(_compare_pdb_worker, args): args[0]
                for args in worker_args
            }

            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    result = future.result()
                    results[pdb_id] = result
                except Exception as e:
                    results[pdb_id] = ComparisonResult(
                        pdb_id=pdb_id,
                        status="error",
                        errors=[f"Exception: {str(e)}"],
                    )

                completed += 1
                if completed % 100 == 0 or completed == len(pdb_ids):
                    click.echo(
                        f"  Progress: {completed}/{len(pdb_ids)} "
                        f"({completed*100//len(pdb_ids)}%)",
                        err=True,
                    )

    return results


@click.group()
def cli():
    """Compare JSON files between legacy and modern X3DNA implementations.

    Use subcommands to specify what to compare:
    - compare: Compare atoms, frames, and step parameters (default)
    - atoms: Compare only atom records
    - frames: Compare only frame calculations
    - steps: Compare only step parameters (bpstep_params, helical_params)
    - list: List available PDB files
    """
    pass


@cli.command()
@click.option("--verbose", "-v", is_flag=True, help="Show detailed output")
@click.option("--diff-only", is_flag=True, help="Show only files with differences (default: True)")
@click.option("--show-all", is_flag=True, help="Show all files, including perfect matches")
@common_options
@click.argument("pdb_ids", nargs=-1)
def compare(verbose, diff_only, show_all, legacy_mode, output, threads, regenerate, test_set, config, pdb_ids):
    """Compare atoms and frames for specified PDB(s) or all available."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()
    
    # Default to diff_only=True unless --show-all is set
    if not show_all:
        diff_only = True

    if test_set:
        test_pdb_ids = load_test_set(project_root, int(test_set))
        if test_pdb_ids is None:
            click.echo(f"Error: Test set of size {test_set} not found. Generate it first with 'generate-test-sets' command.", err=True)
            return
        pdb_ids = test_pdb_ids
        click.echo(f"Using test set of size {test_set}: {len(pdb_ids)} PDB files")
    elif not pdb_ids:
        click.echo("Finding available PDB files...")
        pdb_ids = find_available_pdbs(project_root, legacy_mode)
        if not pdb_ids:
            click.echo("No common PDB files found!", err=True)
            return
        click.echo(f"Found {len(pdb_ids)} PDB files to compare")

    click.echo(f"Comparing {len(pdb_ids)} file(s) using {max_workers} thread(s)...")
    if regenerate:
        click.echo(
            "  (Will regenerate JSON files if missing or with --regenerate flag)"
        )
    
    # Load config from YAML file
    config_path = config or (project_root / "comparison_config.yaml")
    comparator = create_comparator(config_path, project_root)
    
    # Show which comparisons are enabled
    enabled = []
    if comparator.compare_atoms:
        enabled.append("atoms")
    if comparator.compare_frames:
        enabled.append("frames")
    if comparator.compare_steps:
        enabled.append("steps")
    if comparator.compare_pairs:
        enabled.append("pairs")
    if comparator.compare_hbond_list:
        enabled.append("hbond_list")
    if comparator.compare_residue_indices:
        enabled.append("residue_indices")
    if enabled:
        click.echo(f"  Enabled comparisons: {', '.join(enabled)}")
    results = run_comparison(
        pdb_ids, project_root, legacy_mode, comparator, max_workers, regenerate
    )

    # Generate report based on mode
    if verbose and len(pdb_ids) == 1:
        # Verbose mode for single PDB
        pdb_id = pdb_ids[0]
        result = results.get(pdb_id)
        if result:
            report = generate_verbose_report(pdb_id, result, project_root, comparator.tolerance)
        else:
            report = f"Error: No comparison result for {pdb_id}"
    else:
        # Standard summary report
        if verbose and len(pdb_ids) > 1:
            click.echo("Warning: Verbose mode only supported for single PDB. Using standard report.", err=True)
        stats = analyze_results(results)
        report = generate_report(stats, results, verbose, diff_only)

    if output:
        output.write_text(report)
        click.echo(f"Report saved to: {output}")
    else:
        click.echo(report)

    if diff_only:
        diff_files = [
            pdb_id for pdb_id, result in results.items() if result.has_differences()
        ]
        if diff_files:
            click.echo(f"\nFiles with differences: {', '.join(sorted(diff_files))}")
        else:
            click.echo("\n✅ No files with differences found!")


@cli.command()
@click.option("--verbose", "-v", is_flag=True, help="Show detailed output")
@click.option("--diff-only", is_flag=True, help="Show only files with differences (default: True)")
@click.option("--show-all", is_flag=True, help="Show all files, including perfect matches")
@common_options
@click.argument("pdb_ids", nargs=-1)
def atoms(verbose, diff_only, show_all, legacy_mode, output, threads, regenerate, test_set, config, pdb_ids):
    """Compare only atom records (skip frame calculations)."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()
    
    # Default to diff_only=True unless --show-all is set
    if not show_all:
        diff_only = True

    if test_set:
        test_pdb_ids = load_test_set(project_root, int(test_set))
        if test_pdb_ids is None:
            click.echo(f"Error: Test set of size {test_set} not found. Generate it first with 'generate-test-sets' command.", err=True)
            return
        pdb_ids = test_pdb_ids
        click.echo(f"Using test set of size {test_set}: {len(pdb_ids)} PDB files")
    elif not pdb_ids:
        click.echo("Finding available PDB files...")
        pdb_ids = find_available_pdbs(project_root, legacy_mode)
        if not pdb_ids:
            click.echo("No common PDB files found!", err=True)
            return
        click.echo(f"Found {len(pdb_ids)} PDB files to compare")

    if regenerate:
        click.echo("Regenerating JSON files if missing...")

    click.echo(
        f"Comparing atoms for {len(pdb_ids)} file(s) using {max_workers} thread(s)..."
    )

    comparator = JsonComparator(
        compare_atoms=True, compare_frames=False, compare_steps=False
    )
    results = run_comparison(
        pdb_ids, project_root, legacy_mode, comparator, max_workers, regenerate
    )

    stats = analyze_results(results)
    report = generate_report(stats, results, verbose, diff_only)

    if output:
        output.write_text(report)
        click.echo(f"Report saved to: {output}")
    else:
        click.echo(report)


@cli.command()
@click.option("--verbose", "-v", is_flag=True, help="Show detailed output")
@click.option("--diff-only", is_flag=True, help="Show only files with differences (default: True)")
@click.option("--show-all", is_flag=True, help="Show all files, including perfect matches")
@common_options
@click.argument("pdb_ids", nargs=-1)
def frames(verbose, diff_only, show_all, legacy_mode, output, threads, regenerate, test_set, config, pdb_ids):
    """Compare only frame calculations (skip atom records)."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()
    
    # Default to diff_only=True unless --show-all is set
    if not show_all:
        diff_only = True

    if test_set:
        test_pdb_ids = load_test_set(project_root, int(test_set))
        if test_pdb_ids is None:
            click.echo(f"Error: Test set of size {test_set} not found. Generate it first with 'generate-test-sets' command.", err=True)
            return
        pdb_ids = test_pdb_ids
        click.echo(f"Using test set of size {test_set}: {len(pdb_ids)} PDB files")
    elif not pdb_ids:
        click.echo("Finding available PDB files...")
        pdb_ids = find_available_pdbs(project_root, legacy_mode)
        if not pdb_ids:
            click.echo("No common PDB files found!", err=True)
            return
        click.echo(f"Found {len(pdb_ids)} PDB files to compare")

    if regenerate:
        click.echo("Regenerating JSON files if missing...")

    click.echo(
        f"Comparing frames for {len(pdb_ids)} file(s) using {max_workers} thread(s)..."
    )

    comparator = JsonComparator(
        compare_atoms=False, compare_frames=True, compare_steps=False
    )
    results = run_comparison(
        pdb_ids, project_root, legacy_mode, comparator, max_workers, regenerate
    )

    stats = analyze_results(results)
    report = generate_report(stats, results, verbose, diff_only)

    if output:
        output.write_text(report)
        click.echo(f"Report saved to: {output}")
    else:
        click.echo(report)


@cli.command()
@click.option("--verbose", "-v", is_flag=True, help="Show detailed output")
@click.option("--diff-only", is_flag=True, help="Show only files with differences (default: True)")
@click.option("--show-all", is_flag=True, help="Show all files, including perfect matches")
@common_options
@click.argument("pdb_ids", nargs=-1)
def steps(verbose, diff_only, show_all, legacy_mode, output, threads, regenerate, test_set, config, pdb_ids):
    """Compare step parameters (bpstep_params, helical_params). 
    
    IMPORTANT: Frames are verified first since step parameters depend on frames.
    If frames don't match, step parameter comparison may be unreliable.
    """
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()
    
    # Default to diff_only=True unless --show-all is set
    if not show_all:
        diff_only = True

    if test_set:
        test_pdb_ids = load_test_set(project_root, int(test_set))
        if test_pdb_ids is None:
            click.echo(f"Error: Test set of size {test_set} not found. Generate it first with 'generate-test-sets' command.", err=True)
            return
        pdb_ids = test_pdb_ids
        click.echo(f"Using test set of size {test_set}: {len(pdb_ids)} PDB files")
    elif not pdb_ids:
        click.echo("Finding available PDB files...")
        pdb_ids = find_available_pdbs(project_root, legacy_mode)
        if not pdb_ids:
            click.echo("No common PDB files found!", err=True)
            return
        click.echo(f"Found {len(pdb_ids)} PDB files to compare")

    if regenerate:
        click.echo("Regenerating JSON files if missing...")

    click.echo(
        f"Comparing step parameters for {len(pdb_ids)} file(s) using {max_workers} thread(s)..."
    )
    click.echo("  NOTE: Frames will be verified first (step parameters depend on frames)")

    # Compare frames FIRST, then steps (frames are required for step parameter comparison)
    comparator = JsonComparator(
        compare_atoms=False, compare_frames=True, compare_steps=True
    )
    results = run_comparison(
        pdb_ids, project_root, legacy_mode, comparator, max_workers, regenerate
    )

    stats = analyze_results(results)
    report = generate_report(stats, results, verbose, diff_only)

    if output:
        output.write_text(report)
        click.echo(f"Report saved to: {output}")
    else:
        click.echo(report)


@cli.command()
@click.option("--legacy-mode", is_flag=True, help="Include _legacy suffix files")
def list(legacy_mode):
    """List all available PDB files for comparison."""
    project_root = Path(__file__).parent.parent
    pdb_ids = find_available_pdbs(project_root, legacy_mode)

    if not pdb_ids:
        click.echo("No common PDB files found!")
        return

    click.echo(f"Found {len(pdb_ids)} PDB files:")
    for pdb_id in pdb_ids:
        click.echo(f"  {pdb_id}")


@cli.command("ring-atoms")
@click.option("--verbose", "-v", is_flag=True, help="Show detailed output")
@click.option("--diff-only", is_flag=True, help="Show only residues with differences (default: True)")
@click.option("--show-all", is_flag=True, help="Show all residues, including perfect matches")
@click.option("--show-pdb-lines", is_flag=True, help="Show PDB lines for atoms that differ")
@click.option("--output-dir", "-o", type=click.Path(path_type=Path), help="Directory to save output files")
@common_options
@click.argument("pdb_ids", nargs=-1)
def ring_atoms(verbose, diff_only, show_all, show_pdb_lines, output_dir, legacy_mode, output, threads, regenerate, test_set, config, pdb_ids):
    """Compare ring atom matching between legacy and modern code."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()
    
    # Default to diff_only=True unless --show-all is set
    if not show_all:
        diff_only = True
    
    # Import ring atom comparison functionality
    import sys
    sys.path.insert(0, str(project_root))
    from scripts.compare_ring_atoms import compare_ring_atoms
    
    if test_set:
        test_pdb_ids = load_test_set(project_root, int(test_set))
        if test_pdb_ids is None:
            click.echo(f"Error: Test set of size {test_set} not found. Generate it first with 'generate-test-sets' command.", err=True)
            return
        pdb_ids = test_pdb_ids
        click.echo(f"Using test set of size {test_set}: {len(pdb_ids)} PDB files")
    elif not pdb_ids:
        click.echo("Finding available PDB files...")
        pdb_ids = find_available_pdbs(project_root, legacy_mode)
        if not pdb_ids:
            click.echo("No common PDB files found!", err=True)
            return
        click.echo(f"Found {len(pdb_ids)} PDB files to compare")
    
    if len(pdb_ids) == 1:
        # Single PDB comparison
        pdb_id = pdb_ids[0]
        legacy_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
        modern_file = get_modern_file_path(project_root, pdb_id, legacy_mode)
        pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
        
        if not legacy_file.exists() or not modern_file or not modern_file.exists():
            click.echo(f"Error: Missing JSON files for {pdb_id}", err=True)
            return
        
        results = compare_ring_atoms(legacy_file, modern_file, pdb_file if show_pdb_lines else None)
        
        # Generate report
        lines = []
        lines.append("=" * 80)
        lines.append("RING ATOM MATCHING COMPARISON")
        lines.append("=" * 80)
        lines.append(f"PDB ID: {pdb_id}")
        lines.append("")
        
        matches = sum(1 for r in results if r.is_match)
        mismatches = len(results) - matches
        
        lines.append("SUMMARY")
        lines.append(f"  Total residues: {len(results)}")
        lines.append(f"  Matches: {matches}")
        lines.append(f"  Mismatches: {mismatches}")
        lines.append("")
        
        if mismatches > 0 or not diff_only:
            lines.append("=" * 80)
            lines.append("DETAILED RESULTS")
            lines.append("=" * 80)
            lines.append("")
            
            for result in results:
                if diff_only and result.is_match:
                    continue
                
                chain_id, residue_seq, insertion = result.residue_key
                status = "✓" if result.is_match else "✗"
                
                lines.append(f"{status} {result.residue_name} {chain_id}:{residue_seq}{insertion} ({result.base_type})")
                lines.append(f"  Legacy: {result.num_matched_legacy} atoms - {result.legacy_matched_atoms}")
                lines.append(f"  Modern: {result.num_matched_modern} atoms - {result.modern_matched_atoms}")
                
                if result.legacy_atom_indices or result.modern_atom_indices:
                    lines.append(f"  Legacy indices: {result.legacy_atom_indices}")
                    lines.append(f"  Modern indices: {result.modern_atom_indices}")
                
                if result.differences:
                    lines.append(f"  Differences:")
                    for diff in result.differences:
                        lines.append(f"    - {diff}")
                
                if show_pdb_lines and (result.legacy_atom_details or result.modern_atom_details):
                    leg_set = set(result.legacy_atom_indices)
                    mod_set = set(result.modern_atom_indices)
                    only_legacy = leg_set - mod_set
                    only_modern = mod_set - leg_set
                    
                    if only_legacy or only_modern:
                        lines.append(f"  PDB Lines for differing atoms:")
                        if only_legacy:
                            lines.append(f"    Only in Legacy:")
                            for idx in sorted(only_legacy):
                                if idx in result.legacy_atom_details:
                                    details = result.legacy_atom_details[idx]
                                    atom_name = details.get('atom_name', '?')
                                    line_num = details.get('line_number', 0)
                                    pdb_line = details.get('pdb_line', '')
                                    if pdb_line:
                                        lines.append(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: {pdb_line[:80]}")
                                    else:
                                        lines.append(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: (PDB line not available)")
                        if only_modern:
                            lines.append(f"    Only in Modern:")
                            for idx in sorted(only_modern):
                                if idx in result.modern_atom_details:
                                    details = result.modern_atom_details[idx]
                                    atom_name = details.get('atom_name', '?')
                                    line_num = details.get('line_number', 0)
                                    pdb_line = details.get('pdb_line', '')
                                    if pdb_line:
                                        lines.append(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: {pdb_line[:80]}")
                                    else:
                                        lines.append(f"      legacy_atom_idx {idx} ({atom_name}) line {line_num}: (PDB line not available)")
                lines.append("")
        
        report = "\n".join(lines)
        
        if output:
            output.write_text(report)
            click.echo(f"Report saved to: {output}")
        else:
            click.echo(report)
    else:
        # Batch comparison
        if output_dir is None:
            output_dir = project_root / "ring_atom_comparisons"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # If --show-all is set, override diff_only
        if show_all:
            diff_only = False
        
        click.echo(f"Comparing ring atoms for {len(pdb_ids)} file(s) using {max_workers} thread(s)...")
        click.echo(f"Output directory: {output_dir}")
        
        def compare_single_ring(pdb_id: str):
            legacy_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
            modern_file = get_modern_file_path(project_root, pdb_id, legacy_mode)
            pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
            
            if not legacy_file.exists() or not modern_file or not modern_file.exists():
                return (pdb_id, None, "Missing JSON files")
            
            try:
                results = compare_ring_atoms(legacy_file, modern_file, None)
                return (pdb_id, results, None)
            except Exception as e:
                return (pdb_id, None, str(e))
        
        completed = 0
        successful = 0
        failed = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {
                executor.submit(compare_single_ring, pdb_id): pdb_id
                for pdb_id in pdb_ids
            }
            
            for future in as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    pdb_id_result, results, error = future.result()
                    completed += 1
                    
                    if error:
                        failed += 1
                        click.echo(f"✗ {pdb_id_result}: {error}", err=True)
                    else:
                        successful += 1
                        # Write output to file
                        output_file = output_dir / f"{pdb_id_result}_ring_atoms.txt"
                        matches = sum(1 for r in results if r.is_match) if results else 0
                        mismatches = len(results) - matches if results else 0
                        
                        with open(output_file, 'w') as f:
                            f.write(f"RING ATOM COMPARISON: {pdb_id_result}\n")
                            f.write(f"Total residues: {len(results)}\n")
                            f.write(f"Matches: {matches}\n")
                            f.write(f"Mismatches: {mismatches}\n\n")
                            
                            for result in results:
                                if diff_only and result.is_match:
                                    continue
                                chain_id, residue_seq, insertion = result.residue_key
                                status = "✓" if result.is_match else "✗"
                                f.write(f"{status} {result.residue_name} {chain_id}:{residue_seq}{insertion} ({result.base_type})\n")
                                f.write(f"  Legacy: {result.legacy_matched_atoms}\n")
                                f.write(f"  Modern: {result.modern_matched_atoms}\n")
                                if result.differences:
                                    for diff in result.differences:
                                        f.write(f"  - {diff}\n")
                                f.write("\n")
                        
                        if completed % 100 == 0 or completed == len(pdb_ids):
                            click.echo(f"  Progress: {completed}/{len(pdb_ids)} ({completed*100//len(pdb_ids)}%)", err=True)
                except Exception as e:
                    failed += 1
                    completed += 1
                    click.echo(f"✗ {pdb_id}: Exception: {e}", err=True)
        
        click.echo("\n" + "=" * 80)
        click.echo("SUMMARY")
        click.echo("=" * 80)
        click.echo(f"Total PDBs: {len(pdb_ids)}")
        click.echo(f"Successful: {successful}")
        click.echo(f"Failed: {failed}")
        click.echo(f"Output directory: {output_dir}")


@cli.command("generate-test-sets")
@click.option("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)")
@click.option("--force", is_flag=True, help="Regenerate test sets even if they already exist")
def generate_test_sets(seed, force):
    """Generate test sets of 10, 50, 100, 500, and 1000 PDBs from available PDB files."""
    project_root = Path(__file__).parent.parent
    
    all_pdbs = get_all_pdb_files(project_root)
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
            click.echo(f"Test set of size {size} already exists. Use --force to regenerate.")
            continue
        
        if size > len(all_pdbs):
            click.echo(f"Skipping size {size}: Only {len(all_pdbs)} PDB files available")
            continue
        
        click.echo(f"Generating test set of size {size}...")
        pdb_ids = generate_test_set(project_root, size, seed)
        
        if save_test_set(project_root, size, pdb_ids):
            click.echo(f"✓ Saved test set of size {size} to {test_set_file}")
        else:
            click.echo(f"✗ Failed to save test set of size {size}", err=True)
        click.echo()
    
    click.echo("Test set generation complete!")
    test_sets_dir = project_root / "resources" / "test_sets"
    click.echo(f"Test sets saved in: {test_sets_dir}")


def main():
    """Main entry point."""
    cli()


if __name__ == "__main__":
    main()

