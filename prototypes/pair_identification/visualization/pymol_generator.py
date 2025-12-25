"""PyMOL script generator for base pair visualization."""

from pathlib import Path
from typing import List, Tuple
from dataclasses import dataclass


@dataclass
class HBondViz:
    """H-bond visualization data."""
    donor_res: int       # Residue number (1 or 2)
    donor_atom: str
    acceptor_res: int
    acceptor_atom: str
    distance: float
    color: str = "yellow"
    label: str = ""


@dataclass
class TemplateViz:
    """Template visualization data."""
    name: str           # e.g., "cWW", "tHS"
    pdb_path: Path
    rmsd: float
    color: str
    enabled: bool = True


class PyMOLScriptGenerator:
    """Generate PyMOL scripts for structural visualization."""

    # Default colors for different elements
    COLORS = {
        "target": "cyan",
        "cWW": "green",
        "tWW": "orange",
        "cWS": "magenta",
        "tHS": "pink",
        "observed_hbond": "yellow",
        "expected_hbond": "green",
        "missing_hbond": "red",
    }

    def __init__(self, output_dir: Path):
        """Initialize generator with output directory.

        Args:
            output_dir: Directory to write .pml scripts.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_pair_script(
        self,
        target_pdb: Path,
        templates: List[TemplateViz],
        observed_hbonds: List[HBondViz],
        expected_hbonds: List[HBondViz],
        output_name: str,
        title: str = "",
    ) -> Path:
        """Generate PyMOL script comparing target pair to templates.

        Args:
            target_pdb: Path to target pair PDB.
            templates: List of template visualizations.
            observed_hbonds: H-bonds we detected (yellow).
            expected_hbonds: Expected H-bonds for classification (green/red).
            output_name: Output filename (without .pml).
            title: Optional title comment.

        Returns:
            Path to generated .pml script.
        """
        lines = []

        # Header
        lines.append(f"# PyMOL visualization: {output_name}")
        if title:
            lines.append(f"# {title}")
        lines.append("")

        # Load target
        lines.append(f"load {target_pdb.absolute()}, target")
        lines.append("")

        # Load templates
        for tmpl in templates:
            lines.append(f"load {tmpl.pdb_path.absolute()}, {tmpl.name}")
        lines.append("")

        # Hide everything, then show sticks
        lines.append("hide everything")
        lines.append("")

        # Style target
        lines.append("# Target structure (thick sticks)")
        lines.append("show sticks, target")
        lines.append("set stick_radius, 0.15, target")
        lines.append(f"color {self.COLORS['target']}, target")
        lines.append("")

        # Style templates
        lines.append("# Templates (thin sticks)")
        for tmpl in templates:
            lines.append(f"show sticks, {tmpl.name}")
            lines.append(f"set stick_radius, 0.08, {tmpl.name}")
            lines.append(f"color {tmpl.color}, {tmpl.name}")
            if not tmpl.enabled:
                lines.append(f"disable {tmpl.name}")
        lines.append("")

        # Color by atom type
        lines.append("# Color by atom type")
        lines.append("color nitrogen, name N*")
        lines.append("color red, name O*")
        lines.append(f"color {self.COLORS['target']}, target and name C*")
        lines.append("")

        # Draw observed H-bonds
        if observed_hbonds:
            lines.append("# Observed H-bonds (detected)")
            for hb in observed_hbonds:
                dist_name = f"obs_{hb.donor_atom}_{hb.acceptor_atom}"
                lines.append(
                    f"distance {dist_name}, "
                    f"target and resi {hb.donor_res} and name {hb.donor_atom}, "
                    f"target and resi {hb.acceptor_res} and name {hb.acceptor_atom}, "
                    f"mode=0"
                )
                lines.append(f"color {hb.color}, {dist_name}")
                lines.append(f"set dash_width, 3, {dist_name}")
            lines.append("")

        # Draw expected H-bonds (for comparison)
        if expected_hbonds:
            lines.append("# Expected H-bonds (for cWW classification)")
            for hb in expected_hbonds:
                dist_name = f"exp_{hb.donor_atom}_{hb.acceptor_atom}"
                lines.append(
                    f"distance {dist_name}, "
                    f"target and resi {hb.donor_res} and name {hb.donor_atom}, "
                    f"target and resi {hb.acceptor_res} and name {hb.acceptor_atom}, "
                    f"mode=0"
                )
                lines.append(f"color {hb.color}, {dist_name}")
                lines.append(f"set dash_width, 2, {dist_name}")
                lines.append(f"set dash_gap, 0.3, {dist_name}")
            lines.append("")

        # Labels for key atoms
        lines.append("# Atom labels")
        base_atoms = ["N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6"]
        for atom in base_atoms:
            lines.append(f'label target and name {atom}, "{atom}"')
        lines.append("set label_size, 12")
        lines.append("set label_color, white")
        lines.append("")

        # View settings
        lines.append("# View settings")
        lines.append("center target")
        lines.append("zoom target, 3")
        lines.append("orient target")
        lines.append("bg_color black")
        lines.append("")

        # Legend as comments
        lines.append("# ============================================")
        lines.append("# LEGEND")
        lines.append("# ============================================")
        lines.append(f"# Target = {self.COLORS['target']} (thick sticks)")
        for tmpl in templates:
            status = "" if tmpl.enabled else " [disabled]"
            lines.append(
                f"# {tmpl.name} = {tmpl.color} (RMSD={tmpl.rmsd:.3f}Å){status}"
            )
        lines.append("#")
        lines.append("# Yellow dashed = Observed H-bonds")
        lines.append("# Green dashed = Expected H-bonds (found)")
        lines.append("# Red dashed = Missing H-bonds")
        lines.append("")
        lines.append("# Toggle templates: disable <name> / enable <name>")

        # Write script
        output_path = self.output_dir / f"{output_name}.pml"
        with open(output_path, 'w') as f:
            f.write("\n".join(lines))

        return output_path

    def generate_simple_pair_script(
        self,
        pdb_path: Path,
        res1_num: int,
        res2_num: int,
        hbonds: List[Tuple[str, str, float]],
        output_name: str,
    ) -> Path:
        """Generate simple script to view a pair with H-bonds.

        Args:
            pdb_path: Path to PDB file.
            res1_num: First residue number.
            res2_num: Second residue number.
            hbonds: List of (donor_atom, acceptor_atom, distance) tuples.
            output_name: Output filename.

        Returns:
            Path to generated script.
        """
        lines = [
            f"# Simple pair visualization: {output_name}",
            "",
            f"load {pdb_path.absolute()}",
            "",
            f"select pair, resi {res1_num}+{res2_num}",
            "hide everything",
            "show sticks, pair",
            "set stick_radius, 0.15, pair",
            "color cyan, pair",
            "color nitrogen, name N*",
            "color red, name O*",
            "",
        ]

        # Add H-bond distances
        for donor, acceptor, dist in hbonds:
            # Determine which residue has each atom
            donor_res = res1_num  # Assume first residue is donor
            acceptor_res = res2_num

            dist_name = f"hb_{donor}_{acceptor}"
            lines.append(
                f"distance {dist_name}, "
                f"resi {donor_res} and name {donor}, "
                f"resi {acceptor_res} and name {acceptor}"
            )
            lines.append(f"color yellow, {dist_name}")

        lines.extend([
            "",
            "center pair",
            "zoom pair, 3",
            "bg_color black",
        ])

        output_path = self.output_dir / f"{output_name}.pml"
        with open(output_path, 'w') as f:
            f.write("\n".join(lines))

        return output_path
