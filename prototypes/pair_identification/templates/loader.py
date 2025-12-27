"""Template file discovery and loading."""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


class TemplateLoader:
    """Discovers and loads base pair template PDB files.

    Manages template file discovery from idealized and exemplar directories,
    with caching for efficient repeated access.
    """

    def __init__(
        self,
        idealized_dir: Path = Path("basepair-idealized"),
        exemplar_dir: Path = Path("basepair-exemplars"),
    ):
        """Initialize template loader.

        Args:
            idealized_dir: Directory with idealized templates organized by LW class.
            exemplar_dir: Directory with exemplar templates (flat structure).
        """
        self.idealized_dir = idealized_dir
        self.exemplar_dir = exemplar_dir
        self._cache: Dict[str, Tuple[Dict, Dict]] = {}

    def find_template(self, lw_class: str, sequence: str) -> Optional[Path]:
        """Find template PDB for given LW class and sequence.

        Searches idealized directory first (organized by LW class subdirectories),
        then exemplar directory with various naming conventions.

        Args:
            lw_class: Leontis-Westhof class (e.g., "cWW", "tWH").
            sequence: Two-letter sequence (e.g., "GC", "AU").

        Returns:
            Path to template PDB file, or None if not found.
        """
        idealized_path = self._find_idealized(lw_class, sequence)
        if idealized_path:
            return idealized_path

        return self._find_exemplar(lw_class, sequence)

    def _find_idealized(self, lw_class: str, sequence: str) -> Optional[Path]:
        """Find template in idealized directory."""
        base_path = self.idealized_dir / lw_class

        patterns = [
            f"{sequence}.pdb",
            f"{sequence[0]}_{sequence[1].lower()}.pdb",
        ]

        for pattern in patterns:
            path = base_path / pattern
            if path.exists():
                return path

        return None

    def _find_exemplar(self, lw_class: str, sequence: str) -> Optional[Path]:
        """Find template in exemplar directory."""
        patterns = [
            f"{sequence[0]}-{sequence[1]}-{lw_class}.pdb",
            f"{sequence[0]}plus{sequence[1]}-{lw_class}.pdb",
            f"{sequence[0].lower()}-{sequence[1]}-{lw_class}.pdb",
            f"{sequence}-{lw_class}.pdb",
        ]

        for pattern in patterns:
            path = self.exemplar_dir / pattern
            if path.exists():
                return path

        return None

    def load_template(
        self, template_path: Path
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """Load template PDB and return atom coords for both residues.

        Uses caching to avoid repeated file I/O.

        Args:
            template_path: Path to template PDB file.

        Returns:
            Tuple of (res1_atoms, res2_atoms) dicts mapping atom names to coords.
        """
        cache_key = str(template_path)
        if cache_key in self._cache:
            return self._cache[cache_key]

        res1_atoms, res2_atoms = self._parse_pdb_atoms(template_path)
        self._cache[cache_key] = (res1_atoms, res2_atoms)
        return res1_atoms, res2_atoms

    def _parse_pdb_atoms(
        self, pdb_file: Path
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        """Parse atom coordinates from template PDB file.

        Assumes residue 1 has res_seq=1, residue 2 has res_seq=2.

        Args:
            pdb_file: Path to PDB file.

        Returns:
            Tuple of (res1_atoms, res2_atoms) dictionaries.
        """
        res1_atoms: Dict[str, np.ndarray] = {}
        res2_atoms: Dict[str, np.ndarray] = {}

        with open(pdb_file) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                coords = self._parse_atom_coords(line)
                res_seq = int(line[22:26].strip())

                if res_seq == 1:
                    res1_atoms[coords[0]] = coords[1]
                else:
                    res2_atoms[coords[0]] = coords[1]

        return res1_atoms, res2_atoms

    def _parse_atom_coords(self, line: str) -> Tuple[str, np.ndarray]:
        """Parse single ATOM line to extract name and coordinates."""
        atom_name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords = np.array([x, y, z], dtype=np.float64)
        return atom_name, coords

    def list_available(self) -> Dict[str, List[str]]:
        """List all available templates by LW class.

        Scans both idealized and exemplar directories to build a complete
        inventory of available templates.

        Returns:
            Dict mapping LW class to list of available sequences.
        """
        templates: Dict[str, List[str]] = {}

        if self.idealized_dir.exists():
            templates.update(self._scan_idealized())

        if self.exemplar_dir.exists():
            templates.update(self._scan_exemplars())

        return templates

    def _scan_idealized(self) -> Dict[str, List[str]]:
        """Scan idealized directory structure."""
        templates: Dict[str, List[str]] = {}

        for lw_dir in self.idealized_dir.iterdir():
            if not lw_dir.is_dir():
                continue

            sequences = self._scan_lw_directory(lw_dir)
            if sequences:
                templates[lw_dir.name] = sequences

        return templates

    def _scan_lw_directory(self, lw_dir: Path) -> List[str]:
        """Scan single LW class directory for sequences."""
        sequences = []
        for pdb_file in lw_dir.glob("*.pdb"):
            seq = self._extract_sequence(pdb_file.stem)
            if seq:
                sequences.append(seq)
        return sorted(set(sequences))

    def _scan_exemplars(self) -> Dict[str, List[str]]:
        """Scan exemplar directory (flat structure)."""
        templates: Dict[str, List[str]] = {}

        for pdb_file in self.exemplar_dir.glob("*.pdb"):
            lw_class, seq = self._parse_exemplar_name(pdb_file.stem)
            if lw_class and seq:
                if lw_class not in templates:
                    templates[lw_class] = []
                templates[lw_class].append(seq)

        for lw_class in templates:
            templates[lw_class] = sorted(set(templates[lw_class]))

        return templates

    def _extract_sequence(self, stem: str) -> str:
        """Extract sequence from filename stem (e.g., 'GC' or 'G_c')."""
        stem = stem.upper()

        if len(stem) == 2 and stem.isalpha():
            return stem

        if "_" in stem:
            parts = stem.split("_")
            if len(parts) == 2 and all(p.isalpha() for p in parts):
                return parts[0] + parts[1].upper()

        return ""

    def _parse_exemplar_name(self, stem: str) -> Tuple[str, str]:
        """Parse LW class and sequence from exemplar filename."""
        if "-" not in stem:
            return "", ""

        parts = stem.split("-")
        if len(parts) < 3:
            return "", ""

        lw_class = parts[-1]
        seq_parts = parts[:-1]

        if len(seq_parts) == 2:
            seq = seq_parts[0][0].upper() + seq_parts[1][0].upper()
            return lw_class, seq

        return "", ""

    def clear_cache(self) -> None:
        """Clear the template cache."""
        self._cache.clear()
