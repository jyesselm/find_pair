"""Generate pair templates from DSSR data and idealized PDB files.

This module scans DSSR JSON files to collect statistics for each (sequence, LW class)
combination, and loads idealized PDB templates for atom coordinates. The result is a
comprehensive template registry for pair classification.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import re

import numpy as np


@dataclass
class DSSRPairInstance:
    """A single pair instance extracted from DSSR data.

    Attributes:
        pdb_id: Source PDB identifier.
        nt1_id: Nucleotide 1 ID (e.g., "A.G1").
        nt2_id: Nucleotide 2 ID (e.g., "A.C72").
        sequence: Two-letter sequence (e.g., "GC").
        lw_class: Leontis-Westhof class (e.g., "cWW").
        n1n9_dist: Distance between N1/N9 atoms.
        interbase_angle: Angle between base planes (degrees).
        planarity: Planarity measure (Angstroms).
        hbonds_num: Number of hydrogen bonds.
        hbonds_desc: H-bond description string.
        frame_origin: Pair frame origin.
        frame_x: Pair frame X-axis.
        frame_y: Pair frame Y-axis.
        frame_z: Pair frame Z-axis.
    """

    pdb_id: str
    nt1_id: str
    nt2_id: str
    sequence: str
    lw_class: str
    n1n9_dist: float
    interbase_angle: float
    planarity: float
    hbonds_num: int
    hbonds_desc: str
    frame_origin: Optional[np.ndarray] = None
    frame_x: Optional[np.ndarray] = None
    frame_y: Optional[np.ndarray] = None
    frame_z: Optional[np.ndarray] = None


@dataclass
class TemplateStats:
    """Statistics for a (sequence, LW class) template.

    Attributes:
        sequence: Two-letter sequence (e.g., "GC").
        lw_class: Leontis-Westhof class (e.g., "cWW").
        count: Number of instances observed.
        n1n9_dist_mean: Mean N1/N9 distance.
        n1n9_dist_std: Std of N1/N9 distance.
        interbase_angle_mean: Mean interbase angle.
        interbase_angle_std: Std of interbase angle.
        planarity_mean: Mean planarity.
        planarity_std: Std of planarity.
        hbonds_num_mean: Mean number of H-bonds.
        expected_hbonds: Most common H-bond pattern.
    """

    sequence: str
    lw_class: str
    count: int = 0
    n1n9_dist_mean: float = 0.0
    n1n9_dist_std: float = 0.0
    interbase_angle_mean: float = 0.0
    interbase_angle_std: float = 0.0
    planarity_mean: float = 0.0
    planarity_std: float = 0.0
    hbonds_num_mean: float = 0.0
    expected_hbonds: str = ""


@dataclass
class IdealizedTemplate:
    """An idealized base pair template with coordinates.

    Attributes:
        sequence: Two-letter sequence (e.g., "GC").
        lw_class: Leontis-Westhof class (e.g., "cWW").
        stats: Statistics from DSSR data.
        res1_atoms: Atom coordinates for residue 1.
        res2_atoms: Atom coordinates for residue 2.
        pdb_file: Source PDB file path.
    """

    sequence: str
    lw_class: str
    stats: Optional[TemplateStats] = None
    res1_atoms: Dict[str, np.ndarray] = field(default_factory=dict)
    res2_atoms: Dict[str, np.ndarray] = field(default_factory=dict)
    pdb_file: Optional[str] = None


class TemplateGenerator:
    """Generate templates from DSSR data and idealized PDB files.

    This class scans DSSR JSON files to collect statistics for each pair type,
    and loads idealized PDB templates for atom coordinates.

    Example:
        generator = TemplateGenerator(
            dssr_dir=Path("data/json_dssr"),
            idealized_dir=Path("basepair-idealized")
        )
        generator.collect_pairs_from_dssr()
        generator.load_idealized_templates()
        generator.save_registry(Path("templates.json"))
    """

    def __init__(
        self,
        dssr_dir: Path | str,
        idealized_dir: Optional[Path | str] = None,
    ):
        """Initialize template generator.

        Args:
            dssr_dir: Directory containing DSSR JSON files.
            idealized_dir: Directory containing idealized PDB templates.
        """
        self.dssr_dir = Path(dssr_dir)
        self.idealized_dir = Path(idealized_dir) if idealized_dir else None

        # (sequence, lw_class) -> list of instances
        self.pair_instances: Dict[Tuple[str, str], List[DSSRPairInstance]] = {}

        # (sequence, lw_class) -> template
        self.templates: Dict[Tuple[str, str], IdealizedTemplate] = {}

    def collect_pairs_from_dssr(
        self,
        max_files: Optional[int] = None,
        verbose: bool = False,
    ) -> None:
        """Scan all DSSR JSON files and collect classified pairs.

        Args:
            max_files: Maximum number of files to process (for testing).
            verbose: Print progress messages.
        """
        dssr_files = sorted(self.dssr_dir.glob("*.json"))

        if max_files:
            dssr_files = dssr_files[:max_files]

        if verbose:
            print(f"Processing {len(dssr_files)} DSSR files...")

        for i, dssr_file in enumerate(dssr_files):
            if verbose and i % 500 == 0:
                print(f"  Processing {i}/{len(dssr_files)}: {dssr_file.stem}")

            try:
                self._process_dssr_file(dssr_file)
            except (json.JSONDecodeError, KeyError, TypeError) as e:
                if verbose:
                    print(f"  Warning: Failed to parse {dssr_file.name}: {e}")

        if verbose:
            print(f"Collected {sum(len(v) for v in self.pair_instances.values())} pairs")
            print(f"Unique (sequence, LW) types: {len(self.pair_instances)}")

    def _process_dssr_file(self, dssr_file: Path) -> None:
        """Process a single DSSR JSON file.

        Args:
            dssr_file: Path to DSSR JSON file.
        """
        with open(dssr_file) as f:
            data = json.load(f)

        pdb_id = dssr_file.stem

        if "pairs" not in data:
            return

        for pair in data["pairs"]:
            # Skip pairs without LW classification
            if "LW" not in pair or not pair["LW"]:
                continue

            # Extract sequence from bp field (e.g., "G-C" -> "GC")
            bp = pair.get("bp", "")
            if not bp or "-" not in bp:
                continue

            parts = bp.split("-")
            if len(parts) != 2:
                continue

            # Normalize sequence (single letter, uppercase)
            seq1 = self._normalize_base(parts[0])
            seq2 = self._normalize_base(parts[1])
            if not seq1 or not seq2:
                continue

            sequence = seq1 + seq2
            lw_class = pair["LW"]

            # Extract frame if available
            frame_origin = None
            frame_x = None
            frame_y = None
            frame_z = None

            if "frame" in pair and pair["frame"]:
                frame = pair["frame"]
                if "origin" in frame:
                    frame_origin = np.array(frame["origin"], dtype=np.float64)
                if "x_axis" in frame:
                    frame_x = np.array(frame["x_axis"], dtype=np.float64)
                if "y_axis" in frame:
                    frame_y = np.array(frame["y_axis"], dtype=np.float64)
                if "z_axis" in frame:
                    frame_z = np.array(frame["z_axis"], dtype=np.float64)

            # Create instance
            instance = DSSRPairInstance(
                pdb_id=pdb_id,
                nt1_id=pair.get("nt1", ""),
                nt2_id=pair.get("nt2", ""),
                sequence=sequence,
                lw_class=lw_class,
                n1n9_dist=pair.get("N1N9_dist", 0.0),
                interbase_angle=pair.get("interBase_angle", 0.0),
                planarity=pair.get("planarity", 0.0),
                hbonds_num=pair.get("hbonds_num", 0),
                hbonds_desc=pair.get("hbonds_desc", ""),
                frame_origin=frame_origin,
                frame_x=frame_x,
                frame_y=frame_y,
                frame_z=frame_z,
            )

            # Add to collection
            key = (sequence, lw_class)
            if key not in self.pair_instances:
                self.pair_instances[key] = []
            self.pair_instances[key].append(instance)

    def _normalize_base(self, base: str) -> str:
        """Normalize base name to single uppercase letter.

        Args:
            base: Base name (e.g., "G", "DG", "g").

        Returns:
            Single uppercase letter (A, C, G, U/T) or empty string.
        """
        base = base.strip().upper()

        # Handle DNA prefixes
        if base.startswith("D") and len(base) == 2:
            base = base[1]

        # Map T to U for consistency
        if base == "T":
            base = "U"

        # Validate
        if base in ["A", "C", "G", "U"]:
            return base

        return ""

    def compute_statistics(self, min_count: int = 10) -> None:
        """Compute statistics for each (sequence, LW class) combination.

        Args:
            min_count: Minimum instances required to include template.
        """
        for key, instances in self.pair_instances.items():
            if len(instances) < min_count:
                continue

            sequence, lw_class = key

            # Collect values (filter None values)
            n1n9_dists = [
                inst.n1n9_dist for inst in instances
                if inst.n1n9_dist is not None and inst.n1n9_dist > 0
            ]
            angles = [
                inst.interbase_angle for inst in instances
                if inst.interbase_angle is not None
            ]
            planarities = [
                inst.planarity for inst in instances
                if inst.planarity is not None
            ]
            hbonds_nums = [
                inst.hbonds_num for inst in instances
                if inst.hbonds_num is not None
            ]

            # Find most common H-bond pattern
            hbond_patterns: Dict[str, int] = {}
            for inst in instances:
                if inst.hbonds_desc:
                    hbond_patterns[inst.hbonds_desc] = (
                        hbond_patterns.get(inst.hbonds_desc, 0) + 1
                    )
            expected_hbonds = ""
            if hbond_patterns:
                expected_hbonds = max(hbond_patterns.items(), key=lambda x: x[1])[0]

            # Compute stats
            stats = TemplateStats(
                sequence=sequence,
                lw_class=lw_class,
                count=len(instances),
                n1n9_dist_mean=np.mean(n1n9_dists) if n1n9_dists else 0.0,
                n1n9_dist_std=np.std(n1n9_dists) if len(n1n9_dists) > 1 else 0.0,
                interbase_angle_mean=np.mean(angles) if angles else 0.0,
                interbase_angle_std=np.std(angles) if len(angles) > 1 else 0.0,
                planarity_mean=np.mean(planarities) if planarities else 0.0,
                planarity_std=np.std(planarities) if len(planarities) > 1 else 0.0,
                hbonds_num_mean=np.mean(hbonds_nums) if hbonds_nums else 0.0,
                expected_hbonds=expected_hbonds,
            )

            # Create template
            template = IdealizedTemplate(
                sequence=sequence,
                lw_class=lw_class,
                stats=stats,
            )
            self.templates[key] = template

    def load_idealized_templates(self) -> None:
        """Load atom coordinates from idealized PDB files.

        This populates res1_atoms and res2_atoms for templates that have
        matching idealized PDB files.
        """
        if self.idealized_dir is None:
            return

        if not self.idealized_dir.exists():
            return

        # Iterate over LW class directories
        for lw_dir in self.idealized_dir.iterdir():
            if not lw_dir.is_dir():
                continue

            lw_class = lw_dir.name

            # Process each PDB file
            for pdb_file in lw_dir.glob("*.pdb"):
                # Extract sequence from filename (e.g., "GC.pdb" -> "GC")
                sequence = pdb_file.stem.upper()

                # Handle special cases (e.g., "A_a.pdb" -> "AA")
                if "_" in sequence:
                    parts = sequence.split("_")
                    if len(parts) == 2:
                        sequence = parts[0].upper() + parts[1].upper()

                key = (sequence, lw_class)

                # Load atoms from PDB
                res1_atoms, res2_atoms = self._parse_pdb_atoms(pdb_file)

                if key in self.templates:
                    self.templates[key].res1_atoms = res1_atoms
                    self.templates[key].res2_atoms = res2_atoms
                    self.templates[key].pdb_file = str(pdb_file)
                else:
                    # Create template without stats
                    template = IdealizedTemplate(
                        sequence=sequence,
                        lw_class=lw_class,
                        res1_atoms=res1_atoms,
                        res2_atoms=res2_atoms,
                        pdb_file=str(pdb_file),
                    )
                    self.templates[key] = template

    def _parse_pdb_atoms(
        self, pdb_file: Path
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """Parse atom coordinates from idealized PDB file.

        Args:
            pdb_file: Path to PDB file.

        Returns:
            Tuple of (res1_atoms, res2_atoms) dictionaries mapping
            atom names to coordinates.
        """
        res1_atoms: Dict[str, np.ndarray] = {}
        res2_atoms: Dict[str, np.ndarray] = {}

        with open(pdb_file) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                # Parse PDB ATOM record
                atom_name = line[12:16].strip()
                res_seq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                coords = np.array([x, y, z], dtype=np.float64)

                if res_seq == 1:
                    res1_atoms[atom_name] = coords
                else:
                    res2_atoms[atom_name] = coords

        return res1_atoms, res2_atoms

    def get_template(
        self, sequence: str, lw_class: str
    ) -> Optional[IdealizedTemplate]:
        """Get template for a specific sequence and LW class.

        Args:
            sequence: Two-letter sequence (e.g., "GC").
            lw_class: Leontis-Westhof class (e.g., "cWW").

        Returns:
            IdealizedTemplate or None if not found.
        """
        return self.templates.get((sequence, lw_class))

    def get_all_lw_classes(self) -> List[str]:
        """Get list of all LW classes in templates.

        Returns:
            Sorted list of LW class names.
        """
        classes = set()
        for _, lw_class in self.templates.keys():
            classes.add(lw_class)
        return sorted(classes)

    def get_sequences_for_lw(self, lw_class: str) -> List[str]:
        """Get all sequences available for a specific LW class.

        Args:
            lw_class: Leontis-Westhof class (e.g., "cWW").

        Returns:
            List of sequence strings.
        """
        sequences = []
        for (seq, lw), _ in self.templates.items():
            if lw == lw_class:
                sequences.append(seq)
        return sorted(sequences)

    def to_json(self) -> Dict:
        """Serialize templates to JSON-compatible dictionary.

        Returns:
            Dictionary with all template data.
        """
        templates_json = []

        for key, template in sorted(self.templates.items()):
            template_dict = {
                "sequence": template.sequence,
                "lw_class": template.lw_class,
                "pdb_file": template.pdb_file,
            }

            # Add stats if available
            if template.stats:
                template_dict["stats"] = {
                    "count": template.stats.count,
                    "n1n9_dist_mean": template.stats.n1n9_dist_mean,
                    "n1n9_dist_std": template.stats.n1n9_dist_std,
                    "interbase_angle_mean": template.stats.interbase_angle_mean,
                    "interbase_angle_std": template.stats.interbase_angle_std,
                    "planarity_mean": template.stats.planarity_mean,
                    "planarity_std": template.stats.planarity_std,
                    "hbonds_num_mean": template.stats.hbonds_num_mean,
                    "expected_hbonds": template.stats.expected_hbonds,
                }

            # Add atom coords if available
            if template.res1_atoms:
                template_dict["res1_atoms"] = {
                    k: v.tolist() for k, v in template.res1_atoms.items()
                }
            if template.res2_atoms:
                template_dict["res2_atoms"] = {
                    k: v.tolist() for k, v in template.res2_atoms.items()
                }

            templates_json.append(template_dict)

        return {
            "num_templates": len(templates_json),
            "lw_classes": self.get_all_lw_classes(),
            "templates": templates_json,
        }

    def save_registry(self, path: Path | str) -> None:
        """Save template registry to JSON file.

        Args:
            path: Output file path.
        """
        path = Path(path)
        with open(path, "w") as f:
            json.dump(self.to_json(), f, indent=2)

    @classmethod
    def load_registry(cls, path: Path | str) -> "TemplateGenerator":
        """Load template registry from JSON file.

        Args:
            path: Input file path.

        Returns:
            TemplateGenerator instance with loaded templates.
        """
        path = Path(path)
        with open(path) as f:
            data = json.load(f)

        generator = cls(dssr_dir=Path("."))

        for template_data in data["templates"]:
            sequence = template_data["sequence"]
            lw_class = template_data["lw_class"]

            # Create stats if available
            stats = None
            if "stats" in template_data:
                s = template_data["stats"]
                stats = TemplateStats(
                    sequence=sequence,
                    lw_class=lw_class,
                    count=s["count"],
                    n1n9_dist_mean=s["n1n9_dist_mean"],
                    n1n9_dist_std=s["n1n9_dist_std"],
                    interbase_angle_mean=s["interbase_angle_mean"],
                    interbase_angle_std=s["interbase_angle_std"],
                    planarity_mean=s["planarity_mean"],
                    planarity_std=s["planarity_std"],
                    hbonds_num_mean=s["hbonds_num_mean"],
                    expected_hbonds=s["expected_hbonds"],
                )

            # Create template
            template = IdealizedTemplate(
                sequence=sequence,
                lw_class=lw_class,
                stats=stats,
                pdb_file=template_data.get("pdb_file"),
            )

            # Load atoms
            if "res1_atoms" in template_data:
                template.res1_atoms = {
                    k: np.array(v, dtype=np.float64)
                    for k, v in template_data["res1_atoms"].items()
                }
            if "res2_atoms" in template_data:
                template.res2_atoms = {
                    k: np.array(v, dtype=np.float64)
                    for k, v in template_data["res2_atoms"].items()
                }

            generator.templates[(sequence, lw_class)] = template

        return generator

    def print_summary(self) -> None:
        """Print summary of collected templates."""
        print(f"\nTemplate Summary")
        print("=" * 60)
        print(f"Total templates: {len(self.templates)}")
        print(f"LW classes: {', '.join(self.get_all_lw_classes())}")

        # Count by LW class
        lw_counts: Dict[str, int] = {}
        for _, lw_class in self.templates.keys():
            lw_counts[lw_class] = lw_counts.get(lw_class, 0) + 1

        print(f"\nBy LW class:")
        for lw_class in sorted(lw_counts.keys()):
            count = lw_counts[lw_class]
            print(f"  {lw_class}: {count} templates")

        # Top 10 by instance count
        print(f"\nTop 10 by DSSR instance count:")
        templates_with_stats = [
            (key, t) for key, t in self.templates.items() if t.stats
        ]
        templates_with_stats.sort(key=lambda x: x[1].stats.count, reverse=True)

        for (seq, lw), template in templates_with_stats[:10]:
            print(
                f"  {seq}-{lw}: {template.stats.count:,} instances, "
                f"N1N9={template.stats.n1n9_dist_mean:.2f}±{template.stats.n1n9_dist_std:.2f} Å"
            )


def main():
    """Generate templates from DSSR data."""
    import argparse

    parser = argparse.ArgumentParser(description="Generate pair templates from DSSR")
    parser.add_argument(
        "--dssr-dir",
        type=Path,
        default=Path("data/json_dssr"),
        help="Directory containing DSSR JSON files",
    )
    parser.add_argument(
        "--idealized-dir",
        type=Path,
        default=Path("basepair-idealized"),
        help="Directory containing idealized PDB templates",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("prototypes/pair_identification/data/templates.json"),
        help="Output JSON file",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help="Minimum instances for a template",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Maximum DSSR files to process (for testing)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Generate templates
    generator = TemplateGenerator(
        dssr_dir=args.dssr_dir,
        idealized_dir=args.idealized_dir,
    )

    generator.collect_pairs_from_dssr(max_files=args.max_files, verbose=args.verbose)
    generator.compute_statistics(min_count=args.min_count)
    generator.load_idealized_templates()

    generator.print_summary()

    # Save registry
    generator.save_registry(args.output)
    print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
