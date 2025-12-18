#!/usr/bin/env python3
"""
Generate C++ header with constexpr values from parameters.json.

This script reads the single source of truth (parameters.json) and generates
a C++ header with compile-time constants.

Usage:
    python scripts/generate_parameters.py

Output:
    include/x3dna/config/parameters_generated.hpp

CMake integration:
    This script is run automatically during build when parameters.json changes.
"""

import json
import sys
from pathlib import Path
from datetime import datetime


# Manual mapping for cleaner constant names
NAME_MAPPINGS = {
    # Validation
    'validation_distance_min_dorg': 'MIN_DORG',
    'validation_distance_max_dorg': 'MAX_DORG',
    'validation_distance_min_dv': 'MIN_DV',
    'validation_distance_max_dv': 'MAX_DV',
    'validation_distance_min_dNN': 'MIN_DNN',
    'validation_distance_max_dNN': 'MAX_DNN',
    'validation_angle_min_plane_angle': 'MIN_PLANE_ANGLE',
    'validation_angle_max_plane_angle': 'MAX_PLANE_ANGLE',
    'validation_overlap_threshold': 'OVERLAP_THRESHOLD',

    # H-bond
    'hydrogen_bond_detection_hb_lower': 'HB_LOWER',
    'hydrogen_bond_detection_hb_dist1': 'HB_DIST1',
    'hydrogen_bond_detection_hb_dist2': 'HB_DIST2',
    'hydrogen_bond_detection_hb_atoms': 'HB_ATOMS',
    'hydrogen_bond_thresholds_good_min': 'HB_GOOD_MIN',
    'hydrogen_bond_thresholds_good_max': 'HB_GOOD_MAX',
    'hydrogen_bond_thresholds_filter_max': 'HB_FILTER_MAX',
    'hydrogen_bond_thresholds_nonstandard_min': 'HB_NONSTANDARD_MIN',
    'hydrogen_bond_thresholds_nonstandard_max': 'HB_NONSTANDARD_MAX',
    'hydrogen_bond_thresholds_default_dist2': 'HB_DEFAULT_DIST2',
    'hydrogen_bond_min_base_hb': 'MIN_BASE_HB',
    'hydrogen_bond_linkage_conflict': 'HB_LINKAGE_CONFLICT',

    # Quality score
    'quality_score_d_v_weight': 'D_V_WEIGHT',
    'quality_score_plane_angle_divisor': 'PLANE_ANGLE_DIVISOR',
    'quality_score_wc_bonus': 'WC_QUALITY_BONUS',

    # Nucleotide
    'nucleotide_rmsd_cutoff': 'NT_RMSD_CUTOFF',
    'nucleotide_dnn_fallback': 'DNN_FALLBACK',
    'nucleotide_bond_distance': 'BOND_DISTANCE',
    'nucleotide_min_atom_distance': 'MIN_ATOM_DISTANCE',

    # Helix
    'helix_helix_break': 'HELIX_BREAK',
    'helix_end_stack_xang': 'END_STACK_XANG',
    'helix_std_curved': 'STD_CURVED',

    # Misc
    'misc_alt_list': 'ALT_LIST',
    'misc_o3p_dist': 'O3P_DIST',
    'misc_xbig': 'XBIG',
    'misc_gamut': 'GAMUT',
}


def flatten_dict(d: dict, parent_key: str = '', sep: str = '_') -> dict:
    """Flatten nested dict into single-level dict with compound keys."""
    items = []
    for k, v in d.items():
        if k.startswith('_'):  # Skip description fields
            continue
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def cpp_value(val) -> str:
    """Convert Python value to C++ literal."""
    if isinstance(val, bool):
        return "true" if val else "false"
    elif isinstance(val, int):
        return str(val)
    elif isinstance(val, float):
        if val >= 1e10:
            return f"{val:.0e}".replace("e+", "e")
        elif val == int(val):
            return f"{val:.1f}"
        else:
            return str(val)
    elif isinstance(val, str):
        return f'"{val}"'
    else:
        return str(val)


def cpp_type(val) -> str:
    """Get C++ type for value."""
    if isinstance(val, bool):
        return "bool"
    elif isinstance(val, int):
        return "int"
    elif isinstance(val, float):
        return "double"
    elif isinstance(val, str):
        return "const char*"
    else:
        return "auto"


def generate_header(params: dict, timestamp: str = None) -> str:
    """Generate C++ header content."""

    if timestamp is None:
        timestamp = datetime.now().isoformat()

    flat = flatten_dict(params)

    # Map to clean names
    constants = {}
    for key, val in flat.items():
        if key in NAME_MAPPINGS:
            constants[NAME_MAPPINGS[key]] = val
        else:
            # Fallback: uppercase the key
            constants[key.upper()] = val

    # Group by prefix for organization
    groups = {
        'Validation thresholds': ['MIN_DORG', 'MAX_DORG', 'MIN_DV', 'MAX_DV', 'MIN_DNN', 'MAX_DNN',
                                   'MIN_PLANE_ANGLE', 'MAX_PLANE_ANGLE', 'OVERLAP_THRESHOLD'],
        'Hydrogen bond': ['HB_LOWER', 'HB_DIST1', 'HB_DIST2', 'HB_ATOMS', 'HB_GOOD_MIN', 'HB_GOOD_MAX',
                          'HB_FILTER_MAX', 'HB_NONSTANDARD_MIN', 'HB_NONSTANDARD_MAX', 'HB_DEFAULT_DIST2',
                          'MIN_BASE_HB', 'HB_LINKAGE_CONFLICT'],
        'Quality score': ['D_V_WEIGHT', 'PLANE_ANGLE_DIVISOR', 'WC_QUALITY_BONUS'],
        'Nucleotide identification': ['NT_RMSD_CUTOFF', 'DNN_FALLBACK', 'BOND_DISTANCE', 'MIN_ATOM_DISTANCE'],
        'Helix organization': ['HELIX_BREAK', 'END_STACK_XANG', 'STD_CURVED'],
        'Miscellaneous': ['ALT_LIST', 'O3P_DIST', 'XBIG', 'GAMUT'],
    }

    lines = [
        "/**",
        " * @file parameters_generated.hpp",
        " * @brief Auto-generated compile-time constants from parameters.json",
        " *",
        f" * Generated: {timestamp}",
        " * Source: resources/config/parameters.json",
        " *",
        " * DO NOT EDIT THIS FILE DIRECTLY!",
        " * Edit parameters.json and rebuild (CMake runs generate_parameters.py automatically)",
        " */",
        "",
        "#pragma once",
        "",
        "namespace x3dna {",
        "namespace config {",
        "namespace params {",
        "",
    ]

    used = set()
    for group_name, names in groups.items():
        group_constants = [(n, constants[n]) for n in names if n in constants]
        if not group_constants:
            continue

        lines.append(f"// {group_name}")
        for name, val in group_constants:
            ctype = cpp_type(val)
            cval = cpp_value(val)
            lines.append(f'constexpr {ctype} {name} = {cval};')
            used.add(name)
        lines.append("")

    # Add any unmapped constants
    unmapped = [(k, v) for k, v in constants.items() if k not in used]
    if unmapped:
        lines.append("// Other")
        for name, val in sorted(unmapped):
            ctype = cpp_type(val)
            cval = cpp_value(val)
            lines.append(f'constexpr {ctype} {name} = {cval};')
        lines.append("")

    lines.extend([
        "} // namespace params",
        "} // namespace config",
        "} // namespace x3dna",
        "",
    ])

    return '\n'.join(lines)


def main():
    # Support being called with explicit paths (for CMake)
    if len(sys.argv) == 3:
        json_path = Path(sys.argv[1])
        output_path = Path(sys.argv[2])
    else:
        script_dir = Path(__file__).parent
        project_root = script_dir.parent
        json_path = project_root / "resources" / "config" / "parameters.json"
        output_path = project_root / "include" / "x3dna" / "config" / "parameters_generated.hpp"

    # Load JSON
    with open(json_path) as f:
        params = json.load(f)

    # Generate header
    header = generate_header(params)

    # Only write if content changed (avoid unnecessary rebuilds)
    if output_path.exists():
        existing = output_path.read_text()
        # Compare ignoring timestamp line
        new_lines = [l for l in header.split('\n') if not l.startswith(' * Generated:')]
        old_lines = [l for l in existing.split('\n') if not l.startswith(' * Generated:')]
        if new_lines == old_lines:
            print(f"No changes to {output_path}")
            return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(header)

    print(f"Generated: {output_path}")


if __name__ == "__main__":
    main()
