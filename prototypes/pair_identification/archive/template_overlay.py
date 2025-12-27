#!/usr/bin/env python3
"""Overlay base pair templates with target pairs for PyMOL visualization.

This tool aligns idealized templates to target base pairs and generates
PyMOL scripts or session files for visual comparison.

Usage:
    # Overlay template on a specific pair from a PDB
    python template_overlay.py 1EHZ A.G1 A.C72 --lw cWW

    # All cWW pairs from DSSR
    python template_overlay.py 1EHZ --from-dssr --lw cWW

    # Non-standard pairs only (exclude GC, CG, AU, UA cWW)
    python template_overlay.py 1EHZ --from-dssr --lw cWW --exclude-standard

    # Limit to 100 pairs max from one PDB
    python template_overlay.py 1EHZ --from-dssr --lw cWW --max-pairs 100

    # Scan ALL PDBs for first 100 AG-cWW pairs
    python template_overlay.py --scan-all --lw cWW --seq AG --max-pairs 100

    # All tSH pairs and save as PyMOL session
    python template_overlay.py 1EHZ --from-dssr --lw tSH --save-pse

    # List available LW types in a PDB
    python template_overlay.py 1EHZ --list-types
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import Counter

import numpy as np

from core.residue import Residue, Atom, normalize_base_type
from core.pdb_parser import parse_template_pdb, write_pair_pdb
from core.alignment import kabsch_align, compute_rmsd, align_atom_dicts


# Standard Watson-Crick pairs to exclude with --exclude-standard
STANDARD_WC_PAIRS = {"GC", "CG", "AU", "UA", "AT", "TA", "GU", "UG", "GT", "TG"}


def parse_pdb_residues(pdb_path: Path) -> Dict[str, Residue]:
    """Parse PDB file and return residues keyed by DSSR-style ID.

    Uses core.pdb_parser but converts res_id format from canonical
    "chain-base-seq" to DSSR-style "chain.basename+seq" for compatibility.
    """
    residues: Dict[str, Residue] = {}

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21].strip() or "A"
            res_seq = line[22:26].strip()
            ins_code = line[26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) > 76 else atom_name[0]

            # Build DSSR-style key: "chain.residue_name+seq[ins]"
            seq_with_ins = res_seq.strip() + ins_code if ins_code else res_seq.strip()
            dssr_key = f"{chain}.{res_name}{seq_with_ins}"

            # Create or update residue
            if dssr_key not in residues:
                base_type = normalize_base_type(res_name)
                # Use canonical res_id format internally
                res_id = f"{chain}-{base_type}-{seq_with_ins}"
                residues[dssr_key] = Residue(res_id=res_id, base_type=base_type)

            residues[dssr_key].add_atom(
                name=atom_name,
                coords=np.array([x, y, z]),
                element=element
            )

    return residues


def kabsch_align_legacy(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Wrapper for kabsch_align that returns (R, t) instead of (R, centroid_Q, centroid_P).

    This maintains backward compatibility with the original interface.
    """
    R, centroid_Q, centroid_P = kabsch_align(P, Q)
    t = centroid_Q - R @ centroid_P
    return R, t


def align_template_to_target(
    template_res1: Dict[str, np.ndarray],
    template_res2: Dict[str, np.ndarray],
    target_res1: Residue,
    target_res2: Residue,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], float]:
    """Align template to target using common ring atoms."""
    ring_atoms = ["C2", "C4", "C5", "C6", "N1", "N3", "N7", "C8", "N9"]

    template_points = []
    target_points = []

    for atom_name in ring_atoms:
        if atom_name in template_res1 and atom_name in target_res1.atoms:
            template_points.append(template_res1[atom_name])
            target_points.append(target_res1.atoms[atom_name].coords)
        if atom_name in template_res2 and atom_name in target_res2.atoms:
            template_points.append(template_res2[atom_name])
            target_points.append(target_res2.atoms[atom_name].coords)

    if len(template_points) < 4:
        raise ValueError(f"Not enough common atoms for alignment: {len(template_points)}")

    template_points = np.array(template_points)
    target_points = np.array(target_points)

    R, t = kabsch_align_legacy(template_points, target_points)

    aligned_res1 = {name: R @ coords + t for name, coords in template_res1.items()}
    aligned_res2 = {name: R @ coords + t for name, coords in template_res2.items()}

    aligned_template = np.array([R @ p + t for p in template_points])
    rmsd = compute_rmsd(aligned_template, target_points)

    return aligned_res1, aligned_res2, rmsd


def get_base_letter(res_name: str) -> str:
    """Convert residue name to single letter.

    Wrapper around normalize_base_type for backward compatibility.
    """
    return normalize_base_type(res_name)


def find_template(
    sequence: str,
    lw_class: str,
    idealized_dir: Path,
    exemplar_dir: Path,
) -> Optional[Path]:
    """Find template PDB file for given sequence and LW class."""
    # Try idealized first
    idealized_path = idealized_dir / lw_class / f"{sequence}.pdb"
    if idealized_path.exists():
        return idealized_path

    idealized_path = idealized_dir / lw_class / f"{sequence[0]}_{sequence[1].lower()}.pdb"
    if idealized_path.exists():
        return idealized_path

    # Try exemplars
    for pattern in [
        f"{sequence[0]}-{sequence[1]}-{lw_class}.pdb",
        f"{sequence[0]}plus{sequence[1]}-{lw_class}.pdb",
        f"{sequence[0].lower()}-{sequence[1]}-{lw_class}.pdb",
    ]:
        exemplar_path = exemplar_dir / pattern
        if exemplar_path.exists():
            return exemplar_path

    return None


def write_aligned_template_pdb(
    aligned_res1: Dict[str, np.ndarray],
    aligned_res2: Dict[str, np.ndarray],
    sequence: str,
    output_path: Path,
) -> None:
    """Write aligned template coordinates to PDB file.

    Wrapper around write_pair_pdb that adds a REMARK line.
    """
    # Use core module to write the PDB
    write_pair_pdb(
        res1_atoms=aligned_res1,
        res2_atoms=aligned_res2,
        res1_name=sequence[0],
        res2_name=sequence[1],
        output_path=output_path,
        chain="T",
    )

    # Prepend REMARK line
    with open(output_path, "r") as f:
        content = f.read()

    with open(output_path, "w") as f:
        f.write(f"REMARK   1 Aligned template for {sequence}\n")
        f.write(content)


def generate_pymol_script(
    pdb_path: Path,
    res1_id: str,
    res2_id: str,
    template_pdb: Path,
    aligned_template_pdb: Path,
    sequence: str,
    lw_class: str,
    rmsd: float,
    output_path: Path,
) -> None:
    """Generate PyMOL visualization script for a single pair."""

    def parse_res_id(res_id: str) -> Tuple[str, str, str]:
        parts = res_id.split(".")
        chain = parts[0]
        rest = parts[1]
        i = 0
        while i < len(rest) and not rest[i].isdigit() and rest[i] != '-':
            i += 1
        resn = rest[:i]
        resi = rest[i:]
        return chain, resn, resi

    chain1, resn1, resi1 = parse_res_id(res1_id)
    chain2, resn2, resi2 = parse_res_id(res2_id)

    script = f'''# PyMOL script for template overlay visualization
# Target pair: {res1_id} - {res2_id} ({sequence}-{lw_class})
# Template RMSD: {rmsd:.3f} Å

load {pdb_path.absolute()}, target
load {aligned_template_pdb.absolute()}, template

hide everything
select target_pair, (target and chain {chain1} and resi {resi1}) or (target and chain {chain2} and resi {resi2})
show sticks, target_pair
color marine, target_pair and elem C
show sticks, template
color salmon, template and elem C
color red, elem O
color blue, elem N

select donors, (target_pair and elem N+O)
select acceptors, (target_pair and elem N+O)
distance hbonds_target, donors, acceptors, 3.5, mode=2
color yellow, hbonds_target
hide labels, hbonds_target

center target_pair
zoom target_pair, 5
set stick_radius, 0.15
bg_color white

print("Template: {sequence}-{lw_class}, RMSD: {rmsd:.3f} Angstroms")
'''

    with open(output_path, "w") as f:
        f.write(script)


def generate_combined_pymol_script(
    pdb_path: Path,
    overlays: List[Dict],
    output_path: Path,
    lw_class: str,
    save_pse: bool = False,
) -> None:
    """Generate a combined PyMOL script showing all overlays."""

    def parse_res_id(res_id):
        parts = res_id.split(".")
        chain = parts[0]
        rest = parts[1]
        j = 0
        while j < len(rest) and not rest[j].isdigit() and rest[j] != '-':
            j += 1
        resi = rest[j:]
        return chain, resi

    # Group by sequence for coloring
    sequences = list(set(o["sequence"] for o in overlays))
    seq_colors = {}
    colors = ["marine", "forest", "purple", "orange", "cyan", "magenta", "yellow", "lime"]
    for i, seq in enumerate(sorted(sequences)):
        seq_colors[seq] = colors[i % len(colors)]

    script = f'''# PyMOL session: {lw_class} base pairs
# PDB: {pdb_path.name}
# Total pairs: {len(overlays)}
# Sequences: {", ".join(sorted(sequences))}

# Load structure
load {pdb_path.absolute()}, target
hide everything
show cartoon, target and polymer.nucleic
set cartoon_transparency, 0.7
color gray80, target

'''

    for i, overlay in enumerate(overlays):
        res1_id = overlay["res1_id"]
        res2_id = overlay["res2_id"]
        sequence = overlay["sequence"]
        rmsd = overlay["rmsd"]
        template_pdb = overlay["template_pdb"]

        chain1, resi1 = parse_res_id(res1_id)
        chain2, resi2 = parse_res_id(res2_id)

        template_name = f"tmpl_{i+1}_{sequence}"
        pair_name = f"pair_{i+1}_{sequence}"
        color = seq_colors[sequence]

        script += f'''
# Pair {i+1}: {res1_id} - {res2_id} ({sequence}) RMSD={rmsd:.3f}
load {template_pdb}, {template_name}
select {pair_name}, (target and chain {chain1} and resi {resi1}) or (target and chain {chain2} and resi {resi2})
show sticks, {pair_name}
show sticks, {template_name}
color {color}, {pair_name} and elem C
color salmon, {template_name} and elem C
'''

    script += '''
# Standard element colors
color red, elem O
color blue, elem N
color white, elem H

# Settings
set stick_radius, 0.12
set stick_ball, 0
bg_color white
set ray_opaque_background, 1
set depth_cue, 0

# Create groups for easy toggling
'''

    # Group by sequence
    for seq in sorted(sequences):
        pairs = [f"pair_{i+1}_{seq}" for i, o in enumerate(overlays) if o["sequence"] == seq]
        tmpls = [f"tmpl_{i+1}_{seq}" for i, o in enumerate(overlays) if o["sequence"] == seq]
        script += f"group {seq}_pairs, {' '.join(pairs)}\n"
        script += f"group {seq}_templates, {' '.join(tmpls)}\n"

    script += '''
group all_pairs, *_pairs
group all_templates, *_templates

# Zoom to show all pairs
zoom all_pairs, 3

print("=" * 60)
print("Loaded {num_pairs} {lw_class} pairs")
print("=" * 60)
print("Toggle visibility:")
print("  disable all_templates  # hide all templates")
print("  enable all_templates   # show all templates")
'''.format(num_pairs=len(overlays), lw_class=lw_class)

    # Add legend
    script += f'''
print("")
print("Sequences and colors:")
'''
    for seq in sorted(sequences):
        count = sum(1 for o in overlays if o["sequence"] == seq)
        script += f'print("  {seq}: {seq_colors[seq]} ({count} pairs)")\n'

    if save_pse:
        pse_path = output_path.with_suffix(".pse")
        script += f'''
# Save session
save {pse_path.absolute()}
print("")
print("Session saved to: {pse_path}")
'''

    with open(output_path, "w") as f:
        f.write(script)


def generate_multi_pdb_pymol_script(
    overlays: List[Dict],
    output_path: Path,
    lw_class: str,
    pdb_dir: Path,
    save_pse: bool = False,
) -> None:
    """Generate PyMOL script for overlays from multiple PDBs.

    Each pair and its template are grouped together. Only pair residues
    are shown (not full PDB structures). Target pairs are colored marine,
    templates are colored salmon.
    """

    def parse_res_id(res_id):
        parts = res_id.split(".")
        chain = parts[0]
        rest = parts[1]
        j = 0
        while j < len(rest) and not rest[j].isdigit() and rest[j] != '-':
            j += 1
        resi = rest[j:]
        return chain, resi

    pdb_ids = sorted(set(o["pdb_id"] for o in overlays))
    sequences = sorted(set(o["sequence"] for o in overlays))

    script = f'''# PyMOL session: {lw_class} base pairs
# Total pairs: {len(overlays)}
# From {len(pdb_ids)} PDBs
# Sequences: {", ".join(sequences)}
# Target pairs: marine | Templates: salmon

'''

    # Load PDBs (will extract just the pairs we need)
    for pdb_id in pdb_ids:
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            pdb_path = pdb_dir / f"{pdb_id.lower()}.pdb"
        script += f"load {pdb_path.absolute()}, {pdb_id}_full\n"

    script += "\nhide everything\n\n"

    # Process each overlay - create pair + template group
    for i, overlay in enumerate(overlays):
        pdb_id = overlay["pdb_id"]
        res1_id = overlay["res1_id"]
        res2_id = overlay["res2_id"]
        sequence = overlay["sequence"]
        rmsd = overlay["rmsd"]
        template_pdb = overlay["template_pdb"]

        chain1, resi1 = parse_res_id(res1_id)
        chain2, resi2 = parse_res_id(res2_id)

        group_name = f"bp_{i+1}_{pdb_id}_{sequence}"
        pair_name = f"target_{i+1}"
        template_name = f"tmpl_{i+1}"

        script += f'''# {i+1}: {pdb_id} {res1_id}-{res2_id} ({sequence}) RMSD={rmsd:.3f}
load {template_pdb}, {template_name}
select {pair_name}, ({pdb_id}_full and chain {chain1} and resi {resi1}) or ({pdb_id}_full and chain {chain2} and resi {resi2})
show sticks, {pair_name}
show sticks, {template_name}
color marine, {pair_name} and elem C
color salmon, {template_name} and elem C
group {group_name}, {pair_name} {template_name}

'''

    script += '''# Element colors
color red, elem O
color blue, elem N

# Settings
set stick_radius, 0.15
bg_color white

# Hide full PDBs (we only need the extracted pairs)
'''
    for pdb_id in pdb_ids:
        script += f"disable {pdb_id}_full\n"

    script += f'''
# Group all
group all_targets, target_*
group all_templates, tmpl_*
group all_pairs, bp_*

zoom all_targets, 3

print("Loaded {len(overlays)} {lw_class} pairs")
print("  Target: marine | Template: salmon")
print("  Toggle: disable/enable all_templates")
'''

    if save_pse:
        pse_path = output_path.with_suffix(".pse")
        script += f'''
save {pse_path.absolute()}
print("Saved: {pse_path}")
'''

    with open(output_path, "w") as f:
        f.write(script)


def load_dssr_pairs(dssr_path: Path) -> List[Dict]:
    """Load pairs from DSSR JSON."""
    with open(dssr_path) as f:
        data = json.load(f)
    return data.get("pairs", [])


def list_lw_types(dssr_path: Path) -> None:
    """List all LW types in a DSSR file with counts."""
    pairs = load_dssr_pairs(dssr_path)

    lw_counts = Counter()
    seq_lw_counts = Counter()

    for pair in pairs:
        lw = pair.get("LW", "")
        bp = pair.get("bp", "")
        if lw:
            lw_counts[lw] += 1
            if bp and "-" in bp:
                parts = bp.split("-")
                if len(parts) == 2:
                    seq = parts[0][-1].upper() + parts[1][-1].upper()
                    seq_lw_counts[(seq, lw)] += 1

    print(f"\nLW types in {dssr_path.name}:")
    print("=" * 50)
    for lw, count in lw_counts.most_common():
        print(f"  {lw}: {count} pairs")

    print(f"\nBy sequence-LW (top 30):")
    print("-" * 50)
    for (seq, lw), count in seq_lw_counts.most_common(30):
        is_standard = seq in STANDARD_WC_PAIRS and lw == "cWW"
        marker = " [standard]" if is_standard else ""
        print(f"  {seq}-{lw}: {count}{marker}")


def run_pymol_save_pse(pml_path: Path) -> bool:
    """Run PyMOL to execute script and save session."""
    try:
        result = subprocess.run(
            ["pymol", "-cq", str(pml_path)],
            capture_output=True,
            text=True,
            timeout=60
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"Warning: Could not run PyMOL automatically: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Overlay base pair templates for PyMOL visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single pair
  %(prog)s 1EHZ A.G1 A.C72 --lw cWW

  # All tSH pairs from DSSR
  %(prog)s 1EHZ --from-dssr --lw tSH

  # Non-standard cWW only (exclude GC, CG, AU, UA)
  %(prog)s 1EHZ --from-dssr --lw cWW --exclude-standard

  # Limit to 100 pairs from one PDB
  %(prog)s 1EHZ --from-dssr --lw cWW --max-pairs 100

  # Scan ALL PDBs for first 100 AG-cWW pairs
  %(prog)s --scan-all --lw cWW --seq AG --max-pairs 100

  # Scan for GA-tSH pairs
  %(prog)s --scan-all --lw tSH --seq GA --max-pairs 50

  # Save as PyMOL session file
  %(prog)s 1EHZ --from-dssr --lw tWH --save-pse

  # List available LW types
  %(prog)s 1EHZ --list-types
"""
    )
    parser.add_argument("pdb_id", nargs="?", help="PDB ID or path to PDB file")
    parser.add_argument("res1", nargs="?", help="Residue 1 ID (e.g., A.G1)")
    parser.add_argument("res2", nargs="?", help="Residue 2 ID (e.g., A.C72)")
    parser.add_argument("--lw", default="cWW", help="LW class (default: cWW)")
    parser.add_argument(
        "--seq", type=str, default=None,
        help="Filter by sequence (e.g., AG, GA, AC)"
    )
    parser.add_argument(
        "--scan-all", action="store_true",
        help="Scan all DSSR files for pairs (use with --max-pairs)"
    )
    parser.add_argument(
        "--from-dssr", action="store_true",
        help="Use DSSR to find pairs"
    )
    parser.add_argument(
        "--exclude-standard", action="store_true",
        help="Exclude standard WC pairs (GC, CG, AU, UA, GU, UG)"
    )
    parser.add_argument(
        "--list-types", action="store_true",
        help="List available LW types and exit"
    )
    parser.add_argument(
        "--save-pse", action="store_true",
        help="Generate PyMOL session file (.pse)"
    )
    parser.add_argument(
        "--pdb-dir", type=Path, default=Path("data/pdb"),
        help="Directory containing PDB files"
    )
    parser.add_argument(
        "--dssr-dir", type=Path, default=Path("data/json_dssr"),
        help="Directory containing DSSR JSON files"
    )
    parser.add_argument(
        "--idealized-dir", type=Path, default=Path("basepair-idealized"),
        help="Directory containing idealized templates"
    )
    parser.add_argument(
        "--exemplar-dir", type=Path, default=Path("basepair-exemplars"),
        help="Directory containing exemplar templates"
    )
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output path"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("."),
        help="Output directory for generated files"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--max-pairs", type=int, default=0,
        help="Maximum number of pairs to process (0 = all)"
    )
    parser.add_argument(
        "--min-rmsd", type=float, default=0.0,
        help="Only include pairs with RMSD above this threshold (e.g., 0.5)"
    )

    args = parser.parse_args()

    # Handle --scan-all mode (scan all DSSR files)
    if args.scan_all:
        if args.max_pairs <= 0:
            print("Error: --scan-all requires --max-pairs to limit results")
            sys.exit(1)

        pairs_to_process = []
        dssr_files = sorted(args.dssr_dir.glob("*.json"))

        if args.verbose:
            print(f"Scanning {len(dssr_files)} DSSR files for {args.lw} pairs...")

        for dssr_path in dssr_files:
            if args.max_pairs > 0 and len(pairs_to_process) >= args.max_pairs:
                break

            pdb_id = dssr_path.stem
            dssr_pairs = load_dssr_pairs(dssr_path)

            for pair in dssr_pairs:
                if args.max_pairs > 0 and len(pairs_to_process) >= args.max_pairs:
                    break

                if pair.get("LW") == args.lw:
                    bp = pair.get("bp", "")
                    if bp and "-" in bp:
                        parts = bp.split("-")
                        if len(parts) == 2:
                            seq = parts[0][-1].upper() + parts[1][-1].upper()

                            if args.exclude_standard and seq in STANDARD_WC_PAIRS:
                                continue

                            # Filter by sequence if specified
                            if args.seq and seq.upper() != args.seq.upper():
                                continue

                            pairs_to_process.append({
                                "pdb_id": pdb_id,
                                "res1": pair["nt1"],
                                "res2": pair["nt2"],
                                "lw": pair["LW"],
                                "bp": bp,
                                "sequence": seq,
                            })

        if not pairs_to_process:
            print(f"No {args.lw} pairs found in any DSSR file")
            sys.exit(0)

        if args.verbose:
            pdb_ids = set(p["pdb_id"] for p in pairs_to_process)
            print(f"Found {len(pairs_to_process)} pairs from {len(pdb_ids)} PDBs")

        # Process pairs from multiple PDBs
        args.output_dir.mkdir(parents=True, exist_ok=True)
        overlays = []
        pdb_cache: Dict[str, Dict[str, Residue]] = {}

        for i, pair_info in enumerate(pairs_to_process):
            pdb_id = pair_info["pdb_id"]
            res1_id = pair_info["res1"]
            res2_id = pair_info["res2"]
            lw_class = pair_info["lw"]
            sequence = pair_info["sequence"]

            # Load PDB if not cached
            if pdb_id not in pdb_cache:
                pdb_path = args.pdb_dir / f"{pdb_id}.pdb"
                if not pdb_path.exists():
                    pdb_path = args.pdb_dir / f"{pdb_id.lower()}.pdb"
                if not pdb_path.exists():
                    if args.verbose:
                        print(f"Warning: PDB not found: {pdb_id}")
                    continue
                pdb_cache[pdb_id] = parse_pdb_residues(pdb_path)

            residues = pdb_cache[pdb_id]

            if res1_id not in residues or res2_id not in residues:
                if args.verbose:
                    print(f"Warning: Residue not found in {pdb_id}")
                continue

            target_res1 = residues[res1_id]
            target_res2 = residues[res2_id]

            if args.verbose:
                print(f"Processing {i+1}/{len(pairs_to_process)}: {pdb_id} {res1_id}-{res2_id} ({sequence}-{lw_class})")

            # Find and align template
            template_path = find_template(
                sequence, lw_class,
                args.idealized_dir, args.exemplar_dir
            )

            if not template_path:
                base1 = target_res1.base_type
                base2 = target_res2.base_type
                rev_sequence = base2 + base1
                template_path = find_template(
                    rev_sequence, lw_class,
                    args.idealized_dir, args.exemplar_dir
                )
                if template_path:
                    target_res1, target_res2 = target_res2, target_res1
                    res1_id, res2_id = res2_id, res1_id
                    sequence = rev_sequence

            if not template_path:
                if args.verbose:
                    print(f"  Warning: No template for {sequence}-{lw_class}")
                continue

            template_res1, template_res2 = parse_template_pdb(template_path)

            try:
                aligned_res1, aligned_res2, rmsd = align_template_to_target(
                    template_res1, template_res2,
                    target_res1, target_res2
                )
            except ValueError as e:
                if args.verbose:
                    print(f"  Warning: Alignment failed: {e}")
                continue

            # Filter by min RMSD threshold
            if args.min_rmsd > 0 and rmsd < args.min_rmsd:
                if args.verbose:
                    print(f"  RMSD: {rmsd:.3f} Å (skipped, below {args.min_rmsd})")
                continue

            if args.verbose:
                print(f"  RMSD: {rmsd:.3f} Å")

            safe_name = f"{pdb_id}_{res1_id}_{res2_id}_{lw_class}".replace(".", "_")
            aligned_pdb_path = args.output_dir / f"{safe_name}_template.pdb"
            write_aligned_template_pdb(aligned_res1, aligned_res2, sequence, aligned_pdb_path)

            overlays.append({
                "pdb_id": pdb_id,
                "res1_id": res1_id,
                "res2_id": res2_id,
                "sequence": sequence,
                "lw_class": lw_class,
                "rmsd": rmsd,
                "template_pdb": str(aligned_pdb_path.absolute()),
            })

        if not overlays:
            print("No pairs could be processed")
            sys.exit(1)

        # Generate multi-PDB script
        suffix = ""
        if args.seq:
            suffix += f"_{args.seq.upper()}"
        if args.exclude_standard:
            suffix += "_nonstandard"
        if args.min_rmsd > 0:
            suffix += f"_rmsd{args.min_rmsd}"
        pml_path = args.output if args.output else args.output_dir / f"scan_{args.lw}{suffix}.pml"
        generate_multi_pdb_pymol_script(
            overlays, pml_path, args.lw, args.pdb_dir, save_pse=args.save_pse
        )

        print(f"\nGenerated: {pml_path}")
        print(f"  {len(overlays)} pairs from {len(set(o['pdb_id'] for o in overlays))} PDBs")

        seq_counts = Counter(o["sequence"] for o in overlays)
        print(f"\nBy sequence:")
        for seq, count in seq_counts.most_common():
            rmsds = [o["rmsd"] for o in overlays if o["sequence"] == seq]
            print(f"  {seq}: {count} pairs, avg RMSD={np.mean(rmsds):.3f} Å")

        if args.save_pse:
            print(f"\nGenerating .pse file...")
            if run_pymol_save_pse(pml_path):
                print(f"  Saved: {pml_path.with_suffix('.pse')}")
            else:
                print(f"  Run manually: pymol {pml_path}")
        else:
            print(f"\nRun: pymol {pml_path}")

        sys.exit(0)

    # Single PDB mode - need pdb_id
    if not args.pdb_id:
        print("Error: pdb_id required (or use --scan-all)")
        sys.exit(1)

    # Find PDB file
    if Path(args.pdb_id).exists():
        pdb_path = Path(args.pdb_id)
        pdb_id = pdb_path.stem
    else:
        pdb_id = args.pdb_id.upper()
        pdb_path = args.pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            pdb_path = args.pdb_dir / f"{pdb_id.lower()}.pdb"

    if not pdb_path.exists():
        print(f"Error: PDB file not found: {pdb_path}")
        sys.exit(1)

    # Handle --list-types
    if args.list_types:
        dssr_path = args.dssr_dir / f"{pdb_id}.json"
        if not dssr_path.exists():
            print(f"Error: DSSR file not found: {dssr_path}")
            sys.exit(1)
        list_lw_types(dssr_path)
        sys.exit(0)

    # Parse PDB
    if args.verbose:
        print(f"Loading PDB: {pdb_path}")
    residues = parse_pdb_residues(pdb_path)
    if args.verbose:
        print(f"  Found {len(residues)} residues")

    # Get pairs to process
    pairs_to_process = []

    if args.from_dssr:
        dssr_path = args.dssr_dir / f"{pdb_id}.json"
        if not dssr_path.exists():
            print(f"Error: DSSR file not found: {dssr_path}")
            sys.exit(1)

        dssr_pairs = load_dssr_pairs(dssr_path)
        for pair in dssr_pairs:
            if pair.get("LW") == args.lw:
                bp = pair.get("bp", "")
                if bp and "-" in bp:
                    parts = bp.split("-")
                    if len(parts) == 2:
                        seq = parts[0][-1].upper() + parts[1][-1].upper()

                        # Skip standard pairs if requested
                        if args.exclude_standard and seq in STANDARD_WC_PAIRS:
                            continue

                        pairs_to_process.append({
                            "res1": pair["nt1"],
                            "res2": pair["nt2"],
                            "lw": pair["LW"],
                            "bp": bp,
                            "sequence": seq,
                        })

        if args.verbose:
            excluded = " (excluding standard)" if args.exclude_standard else ""
            print(f"Found {len(pairs_to_process)} {args.lw} pairs from DSSR{excluded}")
    else:
        if not args.res1 or not args.res2:
            print("Error: Must specify res1 and res2, or use --from-dssr")
            sys.exit(1)

        pairs_to_process.append({
            "res1": args.res1,
            "res2": args.res2,
            "lw": args.lw,
        })

    if not pairs_to_process:
        print(f"No pairs found for {args.lw}")
        if args.exclude_standard:
            print("  (Try without --exclude-standard)")
        sys.exit(0)

    # Limit number of pairs if requested
    if args.max_pairs > 0 and len(pairs_to_process) > args.max_pairs:
        if args.verbose:
            print(f"Limiting to {args.max_pairs} pairs (from {len(pairs_to_process)})")
        pairs_to_process = pairs_to_process[:args.max_pairs]

    # Process each pair
    args.output_dir.mkdir(parents=True, exist_ok=True)
    overlays = []

    for i, pair_info in enumerate(pairs_to_process):
        res1_id = pair_info["res1"]
        res2_id = pair_info["res2"]
        lw_class = pair_info["lw"]

        if res1_id not in residues:
            if args.verbose:
                print(f"Warning: Residue {res1_id} not found in PDB")
            continue
        if res2_id not in residues:
            if args.verbose:
                print(f"Warning: Residue {res2_id} not found in PDB")
            continue

        target_res1 = residues[res1_id]
        target_res2 = residues[res2_id]

        base1 = target_res1.base_type
        base2 = target_res2.base_type
        sequence = base1 + base2

        if args.verbose:
            print(f"Processing pair {i+1}: {res1_id} - {res2_id} ({sequence}-{lw_class})")

        # Find template
        template_path = find_template(
            sequence, lw_class,
            args.idealized_dir, args.exemplar_dir
        )

        if not template_path:
            rev_sequence = base2 + base1
            template_path = find_template(
                rev_sequence, lw_class,
                args.idealized_dir, args.exemplar_dir
            )
            if template_path:
                target_res1, target_res2 = target_res2, target_res1
                res1_id, res2_id = res2_id, res1_id
                sequence = rev_sequence

        if not template_path:
            if args.verbose:
                print(f"  Warning: No template found for {sequence}-{lw_class}")
            continue

        template_res1, template_res2 = parse_template_pdb(template_path)

        try:
            aligned_res1, aligned_res2, rmsd = align_template_to_target(
                template_res1, template_res2,
                target_res1, target_res2
            )
        except ValueError as e:
            if args.verbose:
                print(f"  Warning: Alignment failed: {e}")
            continue

        if args.verbose:
            print(f"  RMSD: {rmsd:.3f} Å")

        # Generate output files
        safe_name = f"{pdb_id}_{res1_id}_{res2_id}_{lw_class}".replace(".", "_")
        aligned_pdb_path = args.output_dir / f"{safe_name}_template.pdb"

        write_aligned_template_pdb(aligned_res1, aligned_res2, sequence, aligned_pdb_path)

        overlays.append({
            "res1_id": res1_id,
            "res2_id": res2_id,
            "sequence": sequence,
            "lw_class": lw_class,
            "rmsd": rmsd,
            "template_pdb": str(aligned_pdb_path.absolute()),
        })

    if not overlays:
        print("No pairs could be processed")
        sys.exit(1)

    # Generate combined script
    if args.output:
        pml_path = args.output
    else:
        suffix = "_nonstandard" if args.exclude_standard else ""
        pml_path = args.output_dir / f"{pdb_id}_{args.lw}{suffix}.pml"

    generate_combined_pymol_script(
        pdb_path, overlays, pml_path, args.lw, save_pse=args.save_pse
    )

    print(f"\nGenerated: {pml_path}")
    print(f"  {len(overlays)} pairs processed")

    # Summary by sequence
    seq_counts = Counter(o["sequence"] for o in overlays)
    seq_rmsds = {}
    for o in overlays:
        seq = o["sequence"]
        if seq not in seq_rmsds:
            seq_rmsds[seq] = []
        seq_rmsds[seq].append(o["rmsd"])

    print(f"\nBy sequence:")
    for seq in sorted(seq_counts.keys()):
        rmsds = seq_rmsds[seq]
        avg_rmsd = np.mean(rmsds)
        print(f"  {seq}: {seq_counts[seq]} pairs, avg RMSD={avg_rmsd:.3f} Å")

    if args.save_pse:
        print(f"\nGenerating .pse file...")
        if run_pymol_save_pse(pml_path):
            pse_path = pml_path.with_suffix(".pse")
            print(f"  Saved: {pse_path}")
        else:
            print(f"  Run manually: pymol {pml_path}")
            print(f"  Then in PyMOL: save {pml_path.with_suffix('.pse')}")
    else:
        print(f"\nRun: pymol {pml_path}")


if __name__ == "__main__":
    main()
