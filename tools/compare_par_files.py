#!/usr/bin/env python3
"""
Compare modern .par files against legacy X3DNA output.

Compares step parameters (bp_step.par) and helical parameters (bp_helical.par)
by matching unique step labels (e.g., "GC/CG"). Order doesn't matter.

Usage:
    python tools/compare_par_files.py <pdb_file> [--tolerance 0.01]
    python tools/compare_par_files.py <pdb_file> --use-json  # Compare against legacy JSON

Example:
    python tools/compare_par_files.py data/pdb/1EHZ.pdb
    python tools/compare_par_files.py data/pdb/1EHZ.pdb --use-json
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class StepParams:
    """Step parameters for a single step."""
    label: str
    values: List[float]
    line_num: int  # For debugging
    bp_idx1: int = 0  # For JSON matching
    bp_idx2: int = 0


def parse_par_file(filepath: Path, param_type: str = "step", is_legacy: bool = False) -> List[StepParams]:
    """
    Parse a .par file and return parameters as a list in order.

    Args:
        filepath: Path to .par file
        param_type: "step" or "helical"
        is_legacy: If True, parse legacy X3DNA format (12 columns: bp params + step params)

    Returns:
        List of StepParams in order of appearance
    """
    params = []

    if not filepath.exists():
        return params

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line_num, line in enumerate(lines, 1):
        line = line.strip()

        # Skip header and empty lines
        if not line or line.startswith('#'):
            continue

        # Skip metadata lines (e.g., "30 # base-pairs")
        if '# base-pairs' in line or '# ***' in line:
            continue

        parts = line.split()

        # Legacy format has 13 columns (label + 6 bp params + 6 step params)
        # Modern format has 7 columns (label + 6 step params)
        if is_legacy:
            if len(parts) < 13:  # label + 12 parameters
                continue
            label = parts[0]
            try:
                # Step parameters are in columns 7-12 (indices 7-12)
                values = [float(x) for x in parts[7:13]]
            except ValueError:
                continue
            # Skip rows with all-zero step params (first bp has no step)
            if all(abs(v) < 0.001 for v in values):
                continue
        else:
            if len(parts) < 7:  # label + 6 parameters
                continue
            label = parts[0]
            try:
                values = [float(x) for x in parts[1:7]]
            except ValueError:
                continue

        step = StepParams(label=label, values=values, line_num=line_num)
        params.append(step)

    return params


def parse_par_file_by_label(filepath: Path, param_type: str = "step") -> Dict[str, List[StepParams]]:
    """
    Parse a .par file and return parameters grouped by step label.

    Args:
        filepath: Path to .par file
        param_type: "step" or "helical"

    Returns:
        Dict mapping step label -> list of StepParams (multiple entries possible)
    """
    params_list = parse_par_file(filepath, param_type)
    params = {}
    for step in params_list:
        if step.label not in params:
            params[step.label] = []
        params[step.label].append(step)
    return params


def load_legacy_json(json_dir: Path, pdb_name: str, param_type: str) -> Dict[str, List[StepParams]]:
    """
    Load legacy parameters from JSON files.

    Args:
        json_dir: Directory containing legacy JSON files
        pdb_name: PDB name (e.g., "1EHZ")
        param_type: "bpstep_params" or "helical_params"

    Returns:
        Dict mapping (bp_idx1, bp_idx2) tuple -> list of StepParams
    """
    params = {}

    if param_type == "bpstep_params":
        json_file = json_dir / "bpstep_params" / f"{pdb_name}.json"
    else:
        json_file = json_dir / "helical_params" / f"{pdb_name}.json"

    if not json_file.exists():
        return params

    with open(json_file, 'r') as f:
        data = json.load(f)

    for i, record in enumerate(data):
        if record.get('type') != param_type:
            continue

        bp_idx1 = record.get('bp_idx1', 0)
        bp_idx2 = record.get('bp_idx2', 0)

        if param_type == "bpstep_params":
            p = record.get('params', {})
            values = [
                p.get('Shift', 0),
                p.get('Slide', 0),
                p.get('Rise', 0),
                p.get('Tilt', 0),
                p.get('Roll', 0),
                p.get('Twist', 0)
            ]
        else:  # helical_params
            p = record.get('params', [0, 0, 0, 0, 0, 0])
            values = list(p[:6]) if len(p) >= 6 else [0, 0, 0, 0, 0, 0]

        # Use bp_idx pair as key for matching
        key = (bp_idx1, bp_idx2)
        step = StepParams(
            label=f"bp{bp_idx1}-{bp_idx2}",
            values=values,
            line_num=i,
            bp_idx1=bp_idx1,
            bp_idx2=bp_idx2
        )

        if key not in params:
            params[key] = []
        params[key].append(step)

    return params


def load_modern_json(json_dir: Path, pdb_name: str, param_type: str) -> Dict[str, List[StepParams]]:
    """
    Load modern parameters from JSON files.

    Args:
        json_dir: Directory containing modern JSON files
        pdb_name: PDB name (e.g., "1EHZ")
        param_type: "bpstep_params" or "helical_params"

    Returns:
        Dict mapping (bp_idx1, bp_idx2) tuple -> list of StepParams
    """
    params = {}

    if param_type == "bpstep_params":
        json_file = json_dir / "bpstep_params" / f"{pdb_name}.json"
    else:
        json_file = json_dir / "helical_params" / f"{pdb_name}.json"

    if not json_file.exists():
        return params

    with open(json_file, 'r') as f:
        data = json.load(f)

    for i, record in enumerate(data):
        if record.get('type') != param_type:
            continue

        bp_idx1 = record.get('bp_idx1', 0)
        bp_idx2 = record.get('bp_idx2', 0)

        if param_type == "bpstep_params":
            values = [
                record.get('shift', 0),
                record.get('slide', 0),
                record.get('rise', 0),
                record.get('tilt', 0),
                record.get('roll', 0),
                record.get('twist', 0)
            ]
        else:  # helical_params
            values = [
                record.get('x_displacement', 0),
                record.get('y_displacement', 0),
                record.get('rise', 0),
                record.get('inclination', 0),
                record.get('tip', 0),
                record.get('twist', 0)
            ]

        key = (bp_idx1, bp_idx2)
        step = StepParams(
            label=f"bp{bp_idx1}-{bp_idx2}",
            values=values,
            line_num=i,
            bp_idx1=bp_idx1,
            bp_idx2=bp_idx2
        )

        if key not in params:
            params[key] = []
        params[key].append(step)

    return params


def run_legacy_x3dna(pdb_file: Path, work_dir: Path, project_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Run legacy X3DNA find_pair and analyze to generate .par files.
    Uses the local org/build/bin binaries.

    Returns:
        Tuple of (step_par_path, helical_par_path) or (None, None) on error
    """
    # Use local legacy build
    org_bin_dir = project_dir / 'org' / 'build' / 'bin'
    find_pair = org_bin_dir / 'find_pair_original'
    analyze = org_bin_dir / 'analyze_original'

    if not find_pair.exists():
        print(f"Error: find_pair_original not found at {find_pair}")
        print("  Run 'make org-release' to build the legacy binaries")
        return None, None

    if not analyze.exists():
        print(f"Error: analyze_original not found at {analyze}")
        print("  Run 'make org-release' to build the legacy binaries")
        return None, None

    # X3DNA env var is needed for legacy code to find config files
    x3dna_dir = os.environ.get('X3DNA')
    if not x3dna_dir:
        print("Warning: X3DNA environment variable not set, legacy code may fail")
        print("  Set it with: export X3DNA=/path/to/x3dna")

    # Set up environment for subprocess - must include X3DNA
    env = os.environ.copy()
    if x3dna_dir:
        env['X3DNA'] = x3dna_dir

    # Create symlink to data directory for templates
    data_link = work_dir / 'data'
    if not data_link.exists():
        data_link.symlink_to(project_dir / 'data')

    # Run find_pair_original
    inp_file = work_dir / 'legacy.inp'
    try:
        result = subprocess.run(
            [str(find_pair), str(pdb_file.absolute()), str(inp_file)],
            cwd=work_dir,
            capture_output=True,
            timeout=60,
            env=env
        )
        if result.returncode != 0:
            print(f"Warning: find_pair_original returned {result.returncode}")
            if result.stderr:
                print(f"  stderr: {result.stderr.decode()[:200]}")
    except Exception as e:
        print(f"Error running find_pair_original: {e}")
        return None, None

    if not inp_file.exists():
        print("Error: find_pair_original did not create .inp file")
        return None, None

    # Run analyze_original
    try:
        result = subprocess.run(
            [str(analyze), str(inp_file)],
            cwd=work_dir,
            capture_output=True,
            timeout=60,
            env=env
        )
        if result.returncode != 0:
            print(f"Warning: analyze_original returned {result.returncode}")
            if result.stderr:
                stderr_text = result.stderr.decode()
                print(f"  stderr: {stderr_text[:200]}")
                if 'version' in stderr_text.lower() or 'config' in stderr_text.lower():
                    print(f"  Note: X3DNA={x3dna_dir}")
                    print("  Ensure $X3DNA/config/version exists")
    except Exception as e:
        print(f"Error running analyze_original: {e}")
        return None, None

    step_par = work_dir / 'bp_step.par'
    helical_par = work_dir / 'bp_helical.par'

    return (
        step_par if step_par.exists() else None,
        helical_par if helical_par.exists() else None
    )


def run_modern_find_pair(pdb_file: Path, work_dir: Path, build_dir: Path, project_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Run modern find_pair_app to generate .par files.

    Returns:
        Tuple of (step_par_path, helical_par_path) or (None, None) on error
    """
    find_pair_app = build_dir / 'find_pair_app'

    if not find_pair_app.exists():
        print(f"Error: find_pair_app not found at {find_pair_app}")
        return None, None

    # Create symlink to data directory for templates
    data_link = work_dir / 'data'
    if not data_link.exists():
        data_link.symlink_to(project_dir / 'data')

    try:
        result = subprocess.run(
            [str(find_pair_app), str(pdb_file.absolute())],
            cwd=work_dir,
            capture_output=True,
            timeout=60
        )
        if result.returncode != 0:
            print(f"Warning: find_pair_app returned {result.returncode}")
            if result.stderr:
                print(result.stderr.decode())
    except Exception as e:
        print(f"Error running find_pair_app: {e}")
        return None, None

    step_par = work_dir / 'bp_step.par'
    helical_par = work_dir / 'bp_helical.par'

    return (
        step_par if step_par.exists() else None,
        helical_par if helical_par.exists() else None
    )


def compare_params_by_label(
    legacy_params: Dict[str, List[StepParams]],
    modern_params: Dict[str, List[StepParams]],
    param_names: List[str],
    tolerance: float = 0.01
) -> Tuple[int, int, List[str]]:
    """
    Compare parameters between legacy and modern by label.

    Returns:
        Tuple of (matches, mismatches, error_messages)
    """
    matches = 0
    mismatches = 0
    errors = []

    # Get all unique labels
    all_labels = set(legacy_params.keys()) | set(modern_params.keys())

    for label in sorted(all_labels):
        legacy_list = legacy_params.get(label, [])
        modern_list = modern_params.get(label, [])

        if not legacy_list:
            errors.append(f"  {label}: only in modern ({len(modern_list)} occurrences)")
            mismatches += len(modern_list)
            continue

        if not modern_list:
            errors.append(f"  {label}: only in legacy ({len(legacy_list)} occurrences)")
            mismatches += len(legacy_list)
            continue

        # Match by finding closest values
        # For each modern entry, find the best matching legacy entry
        used_legacy = set()

        for modern in modern_list:
            best_legacy_idx = None
            best_diff = float('inf')

            for i, legacy in enumerate(legacy_list):
                if i in used_legacy:
                    continue

                # Calculate total absolute difference
                total_diff = sum(abs(l - m) for l, m in zip(legacy.values, modern.values))
                if total_diff < best_diff:
                    best_diff = total_diff
                    best_legacy_idx = i

            if best_legacy_idx is not None:
                used_legacy.add(best_legacy_idx)
                legacy = legacy_list[best_legacy_idx]

                # Check each parameter
                param_errors = []
                for i, (lv, mv) in enumerate(zip(legacy.values, modern.values)):
                    diff = abs(lv - mv)
                    if diff > tolerance:
                        param_errors.append(f"{param_names[i]}: {lv:.2f} vs {mv:.2f} (diff={diff:.4f})")

                if param_errors:
                    errors.append(f"  {label}: " + ", ".join(param_errors))
                    mismatches += 1
                else:
                    matches += 1
            else:
                errors.append(f"  {label}: no matching legacy entry")
                mismatches += 1

        # Check for unmatched legacy entries
        for i, legacy in enumerate(legacy_list):
            if i not in used_legacy:
                errors.append(f"  {label}: unmatched legacy entry")
                mismatches += 1

    return matches, mismatches, errors


def compare_params_by_order(
    legacy_params: List[StepParams],
    modern_params: List[StepParams],
    param_names: List[str],
    tolerance: float = 0.01
) -> Tuple[int, int, List[str]]:
    """
    Compare parameters between legacy and modern by position/order.
    This works when label formats differ but step order is the same.

    Returns:
        Tuple of (matches, mismatches, error_messages)
    """
    matches = 0
    mismatches = 0
    errors = []

    # Compare by position
    min_len = min(len(legacy_params), len(modern_params))

    for i in range(min_len):
        legacy = legacy_params[i]
        modern = modern_params[i]

        # Check each parameter
        param_errors = []
        for j, (lv, mv) in enumerate(zip(legacy.values, modern.values)):
            diff = abs(lv - mv)
            if diff > tolerance:
                param_errors.append(f"{param_names[j]}: {lv:.2f} vs {mv:.2f} (diff={diff:.4f})")

        if param_errors:
            errors.append(f"  Step {i+1} [{legacy.label} vs {modern.label}]: " + ", ".join(param_errors))
            mismatches += 1
        else:
            matches += 1

    # Report count differences
    if len(legacy_params) > min_len:
        extra = len(legacy_params) - min_len
        errors.append(f"  {extra} extra steps in legacy")
        mismatches += extra

    if len(modern_params) > min_len:
        extra = len(modern_params) - min_len
        errors.append(f"  {extra} extra steps in modern")
        mismatches += extra

    return matches, mismatches, errors


def main():
    parser = argparse.ArgumentParser(
        description='Compare modern .par files against legacy X3DNA output'
    )
    parser.add_argument('pdb_file', type=Path, help='PDB file to analyze')
    parser.add_argument('--tolerance', type=float, default=0.01,
                        help='Tolerance for parameter comparison (default: 0.01)')
    parser.add_argument('--build-dir', type=Path, default=None,
                        help='Path to build directory (default: auto-detect)')
    parser.add_argument('--use-json', action='store_true',
                        help='Compare against legacy JSON files instead of running legacy binary')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show detailed comparison')

    args = parser.parse_args()

    # Verify PDB file exists
    if not args.pdb_file.exists():
        print(f"Error: PDB file not found: {args.pdb_file}")
        sys.exit(1)

    # Find build directory and project directory
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
    build_dir = args.build_dir or project_dir / 'build'

    if not build_dir.exists():
        print(f"Error: Build directory not found: {build_dir}")
        sys.exit(1)

    pdb_name = args.pdb_file.stem

    print(f"Comparing parameters for {pdb_name}")
    print(f"Tolerance: {args.tolerance}")
    if args.use_json:
        print("Mode: JSON comparison")
    print()

    if args.use_json:
        # Compare using existing JSON files
        legacy_json_dir = project_dir / 'data' / 'json_legacy'
        modern_json_dir = project_dir / 'data' / 'json'

        # Compare step parameters
        print("=" * 60)
        print("STEP PARAMETERS (from JSON)")
        print("=" * 60)

        legacy_params = load_legacy_json(legacy_json_dir, pdb_name, "bpstep_params")
        modern_params = load_modern_json(modern_json_dir, pdb_name, "bpstep_params")

        if legacy_params and modern_params:
            step_names = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']
            matches, mismatches, errors = compare_params_by_label(
                legacy_params, modern_params, step_names, args.tolerance
            )

            total = matches + mismatches
            if total > 0:
                pct = 100 * matches / total
                print(f"Matches: {matches}/{total} ({pct:.1f}%)")

                if mismatches > 0 and args.verbose:
                    print("\nMismatches:")
                    for err in errors:
                        print(err)
            else:
                print("No steps to compare")
        else:
            print("Could not compare - missing JSON files")
            if not legacy_params:
                print(f"  - Legacy JSON not found in {legacy_json_dir}")
            if not modern_params:
                print(f"  - Modern JSON not found in {modern_json_dir}")

        print()

        # Compare helical parameters
        print("=" * 60)
        print("HELICAL PARAMETERS (from JSON)")
        print("=" * 60)

        legacy_params = load_legacy_json(legacy_json_dir, pdb_name, "helical_params")
        modern_params = load_modern_json(modern_json_dir, pdb_name, "helical_params")

        if legacy_params and modern_params:
            helical_names = ['X-disp', 'Y-disp', 'Rise', 'Incl', 'Tip', 'Twist']
            matches, mismatches, errors = compare_params_by_label(
                legacy_params, modern_params, helical_names, args.tolerance
            )

            total = matches + mismatches
            if total > 0:
                pct = 100 * matches / total
                print(f"Matches: {matches}/{total} ({pct:.1f}%)")

                if mismatches > 0 and args.verbose:
                    print("\nMismatches:")
                    for err in errors:
                        print(err)
            else:
                print("No helical params to compare")
        else:
            print("Could not compare - missing JSON files")

        print()
        return

    # Create temp directories for legacy and modern output
    with tempfile.TemporaryDirectory() as legacy_dir, \
         tempfile.TemporaryDirectory() as modern_dir:

        legacy_path = Path(legacy_dir)
        modern_path = Path(modern_dir)

        # Run legacy X3DNA (using org/build/bin binaries)
        print("Running legacy X3DNA (org/build/bin)...")
        legacy_step, legacy_helical = run_legacy_x3dna(args.pdb_file, legacy_path, project_dir)

        # Run modern find_pair_app
        print("Running modern find_pair_app...")
        modern_step, modern_helical = run_modern_find_pair(args.pdb_file, modern_path, build_dir, project_dir)

        print()

        # Compare step parameters
        print("=" * 60)
        print("STEP PARAMETERS (bp_step.par)")
        print("=" * 60)

        if legacy_step and modern_step:
            # Legacy format has 12 columns (bp params + step params), step params are columns 7-12
            legacy_params = parse_par_file(legacy_step, "step", is_legacy=True)
            modern_params = parse_par_file(modern_step, "step", is_legacy=False)

            print(f"Legacy: {len(legacy_params)} steps, Modern: {len(modern_params)} steps")

            step_names = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']
            matches, mismatches, errors = compare_params_by_order(
                legacy_params, modern_params, step_names, args.tolerance
            )

            total = matches + mismatches
            if total > 0:
                pct = 100 * matches / total
                print(f"Matches: {matches}/{total} ({pct:.1f}%)")

                if mismatches > 0 and args.verbose:
                    print("\nMismatches:")
                    for err in errors:
                        print(err)
            else:
                print("No steps to compare")
        else:
            print("Could not compare - missing files")
            if not legacy_step:
                print("  - Legacy bp_step.par not found")
            if not modern_step:
                print("  - Modern bp_step.par not found")

        print()

        # Compare helical parameters
        print("=" * 60)
        print("HELICAL PARAMETERS (bp_helical.par)")
        print("=" * 60)

        if legacy_helical and modern_helical:
            # Legacy format has 12 columns (bp params + helical params), helical params are columns 7-12
            legacy_params = parse_par_file(legacy_helical, "helical", is_legacy=True)
            modern_params = parse_par_file(modern_helical, "helical", is_legacy=False)

            print(f"Legacy: {len(legacy_params)} entries, Modern: {len(modern_params)} entries")

            helical_names = ['X-disp', 'Y-disp', 'Rise', 'Incl', 'Tip', 'Twist']
            matches, mismatches, errors = compare_params_by_order(
                legacy_params, modern_params, helical_names, args.tolerance
            )

            total = matches + mismatches
            if total > 0:
                pct = 100 * matches / total
                print(f"Matches: {matches}/{total} ({pct:.1f}%)")

                if mismatches > 0 and args.verbose:
                    print("\nMismatches:")
                    for err in errors:
                        print(err)
            else:
                print("No helical params to compare")
        else:
            print("Could not compare - missing files")
            if not legacy_helical:
                print("  - Legacy bp_helical.par not found")
            if not modern_helical:
                print("  - Modern bp_helical.par not found")

        print()


if __name__ == '__main__':
    main()
