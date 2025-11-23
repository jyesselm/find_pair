#!/usr/bin/env python3
"""
Unified JSON Comparison Script

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
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
from datetime import datetime

from x3dna_json_compare import JsonComparator, ComparisonResult, JsonValidator


def find_available_pdbs(project_root: Path, use_legacy_mode: bool = False) -> List[str]:
    """Find all PDB IDs that have both legacy and modern JSON files."""
    legacy_dir = project_root / "data" / "json_legacy"
    modern_dir = project_root / "data" / "json"

    legacy_files = {f.stem for f in legacy_dir.glob("*.json")}

    modern_files = set()
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


def get_modern_file_path(
    project_root: Path, pdb_id: str, use_legacy_mode: bool = False
) -> Optional[Path]:
    """Get the path to modern JSON file."""
    modern_dir = project_root / "data" / "json"

    if use_legacy_mode:
        legacy_mode_file = modern_dir / f"{pdb_id}_legacy.json"
        if legacy_mode_file.exists():
            return legacy_mode_file

    regular_file = modern_dir / f"{pdb_id}.json"
    if regular_file.exists():
        return regular_file

    return None


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
    json_file = project_root / "data" / "json" / f"{pdb_id}.json"

    if not pdb_file.exists():
        return False

    try:
        # Ensure output directory exists
        json_file.parent.mkdir(parents=True, exist_ok=True)

        # Run generate_modern_json
        cmd = [str(modern_exe), str(pdb_file), str(json_file)]
        if use_legacy_mode:
            cmd.append("--legacy")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if json_file.exists():
            # Validate the generated JSON
            is_valid, _, _ = JsonValidator.validate_and_clean(
                json_file, remove_invalid=True
            )
            return is_valid
        return False
    except Exception:
        return False


def compare_single_pdb(
    pdb_id: str,
    project_root: Path,
    use_legacy_mode: bool = False,
    comparator: Optional[JsonComparator] = None,
    regenerate: bool = False,
) -> ComparisonResult:
    """Compare a single PDB file."""
    if comparator is None:
        comparator = JsonComparator(enable_cache=False)

    legacy_file = project_root / "data" / "json_legacy" / f"{pdb_id}.json"
    modern_file = get_modern_file_path(project_root, pdb_id, use_legacy_mode)
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"

    # Regenerate files if missing (always try to regenerate missing files)
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
    elif regenerate:
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
    elif regenerate:
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

    return stats


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
    elif any(
        result.frame_comparison is None
        for result in results.values()
        if hasattr(result, "frame_comparison")
    ):
        # Frame comparison was skipped
        pass

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

                    lines.append(
                        f"\n  {base_type} {chain_id}:{residue_seq}{insertion}:{idx_str}"
                    )
                    if "rms" in mismatch.mismatches:
                        lines.append(
                            f"    RMS: legacy={mismatch.mismatches['rms']['legacy']:.6f}, "
                            f"modern={mismatch.mismatches['rms']['modern']:.6f}"
                        )
                    if "matched_atoms" in mismatch.mismatches:
                        lines.append(
                            f"    Legacy only: {mismatch.mismatches['matched_atoms']['only_legacy']}"
                        )
                        lines.append(
                            f"    Modern only: {mismatch.mismatches['matched_atoms']['only_modern']}"
                        )
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
    return f


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
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdb = {
                executor.submit(
                    compare_single_pdb,
                    pdb_id,
                    project_root,
                    use_legacy_mode,
                    comparator,
                    regenerate,
                ): pdb_id
                for pdb_id in pdb_ids
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
    - compare: Compare both atoms and frames (default)
    - atoms: Compare only atom records
    - frames: Compare only frame calculations
    - list: List available PDB files
    """
    pass


@cli.command()
@click.option("--verbose", "-v", is_flag=True, help="Show detailed output")
@click.option("--diff-only", is_flag=True, help="Show only files with differences")
@common_options
@click.argument("pdb_ids", nargs=-1)
def compare(verbose, diff_only, legacy_mode, output, threads, regenerate, pdb_ids):
    """Compare atoms and frames for specified PDB(s) or all available."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()

    if not pdb_ids:
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

    comparator = JsonComparator(
        enable_cache=False, compare_atoms=True, compare_frames=True
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
@click.option("--diff-only", is_flag=True, help="Show only files with differences")
@common_options
@click.argument("pdb_ids", nargs=-1)
def atoms(verbose, diff_only, legacy_mode, output, threads, regenerate, pdb_ids):
    """Compare only atom records (skip frame calculations)."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()

    if not pdb_ids:
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
        enable_cache=False, compare_atoms=True, compare_frames=False
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
@click.option("--diff-only", is_flag=True, help="Show only files with differences")
@common_options
@click.argument("pdb_ids", nargs=-1)
def frames(verbose, diff_only, legacy_mode, output, threads, regenerate, pdb_ids):
    """Compare only frame calculations (skip atom records)."""
    project_root = Path(__file__).parent.parent
    max_workers = threads or multiprocessing.cpu_count()

    if not pdb_ids:
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
        enable_cache=False, compare_atoms=False, compare_frames=True
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


def main():
    """Main entry point."""
    cli()


if __name__ == "__main__":
    main()
