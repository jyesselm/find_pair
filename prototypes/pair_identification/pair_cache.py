"""Cache pre-computed frames and validation results for base pairs.

This module loads pre-computed reference frames and atom coordinates, finds all
potential pairs within distance, validates them geometrically, and caches the
results for efficient access.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import json

import numpy as np
from scipy.spatial import KDTree

from prototypes.pair_identification.frame_loader import FrameLoader, ReferenceFrame
from prototypes.pair_identification.geometric_validator import (
    GeometricValidator,
    ValidationResult,
)


@dataclass
class AtomCoords:
    """Coordinates for key atoms of a residue.

    Attributes:
        c1_prime: C1' position (sugar carbon).
        n1_or_n9: N1 (pyrimidines) or N9 (purines) position.
        ring_atoms: Dictionary of ring atom positions.
        hbond_atoms: Dictionary of H-bond donor/acceptor positions.
    """

    c1_prime: Optional[np.ndarray] = None
    n1_or_n9: Optional[np.ndarray] = None
    ring_atoms: Dict[str, np.ndarray] = field(default_factory=dict)
    hbond_atoms: Dict[str, np.ndarray] = field(default_factory=dict)


@dataclass
class CachedPair:
    """Cached data for a potential base pair.

    Attributes:
        res1_id: Residue 1 ID (format: "chain-name-seq[ins]").
        res2_id: Residue 2 ID.
        res1_name: Residue 1 single-letter code.
        res2_name: Residue 2 single-letter code.
        frame1: Reference frame for residue 1.
        frame2: Reference frame for residue 2.
        validation: Geometric validation result.
        hbonds: List of H-bond dicts (to be populated by H-bond finder).
    """

    res1_id: str
    res2_id: str
    res1_name: str
    res2_name: str
    frame1: ReferenceFrame
    frame2: ReferenceFrame
    validation: ValidationResult
    hbonds: List[Dict] = field(default_factory=list)


class PairCache:
    """Cache of validated base pairs with pre-computed frames.

    This class loads reference frames and atom coordinates from JSON output,
    finds all potential pairs within a distance cutoff, validates them
    geometrically, and caches the results.

    Example:
        cache = PairCache("1EHZ", Path("data/json"))
        cache.build_cache(max_distance=15.0)
        valid_pairs = cache.get_valid_pairs()
        print(f"Found {len(valid_pairs)} valid pairs")
        cache.save(Path("cache.json"))
    """

    def __init__(self, pdb_id: str, json_dir: Path | str):
        """Initialize pair cache.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ").
            json_dir: Directory containing JSON output (e.g., "data/json").
        """
        self.pdb_id = pdb_id
        self.json_dir = Path(json_dir)
        self.pairs: List[CachedPair] = []
        self.frames: Dict[str, ReferenceFrame] = {}
        self.atoms: Dict[str, AtomCoords] = {}

    def build_cache(self, max_distance: float = 15.0) -> None:
        """Build cache from pre-computed frames.

        This method:
        1. Loads reference frames from ls_fitting JSON
        2. Loads atom coordinates from pdb_atoms JSON
        3. Builds KDTree on frame origins
        4. Finds pairs within distance and validates them

        Args:
            max_distance: Maximum distance between frame origins (Angstroms).

        Raises:
            FileNotFoundError: If required JSON files do not exist.
            ValueError: If JSON format is invalid.
        """
        # 1. Load frames from ls_fitting JSON
        frame_loader = FrameLoader(self.json_dir)
        self.frames = frame_loader.load_frames(self.pdb_id)

        if not self.frames:
            raise ValueError(f"No frames loaded for {self.pdb_id}")

        # 2. Load atom coordinates from pdb_atoms JSON
        self._load_atom_coords()

        # 3. Build KDTree on frame origins
        res_ids = list(self.frames.keys())
        origins = np.array([self.frames[r].origin for r in res_ids])
        kdtree = KDTree(origins)

        # 4. Find pairs within distance and validate
        validator = GeometricValidator()
        for i, res1_id in enumerate(res_ids):
            # Query neighbors within distance
            neighbors = kdtree.query_ball_point(origins[i], r=max_distance)

            for j in neighbors:
                if j <= i:  # Skip self-pairs and duplicate pairs
                    continue

                res2_id = res_ids[j]

                # Get frames
                frame1 = self.frames[res1_id]
                frame2 = self.frames[res2_id]

                # Get N1/N9 positions
                if res1_id not in self.atoms or res2_id not in self.atoms:
                    continue

                n1n9_pos1 = self.atoms[res1_id].n1_or_n9
                n1n9_pos2 = self.atoms[res2_id].n1_or_n9

                if n1n9_pos1 is None or n1n9_pos2 is None:
                    continue

                # Validate geometry
                validation = validator.validate(frame1, frame2, n1n9_pos1, n1n9_pos2)

                # Extract 1-letter residue names
                res1_name = self._extract_residue_name(res1_id)
                res2_name = self._extract_residue_name(res2_id)

                # Create cached pair
                pair = CachedPair(
                    res1_id=res1_id,
                    res2_id=res2_id,
                    res1_name=res1_name,
                    res2_name=res2_name,
                    frame1=frame1,
                    frame2=frame2,
                    validation=validation,
                )
                self.pairs.append(pair)

    def _load_atom_coords(self) -> None:
        """Load atom coordinates from pdb_atoms JSON.

        Reads from {json_dir}/pdb_atoms/{pdb_id}.json and extracts:
        - C1' positions
        - N1 (pyrimidines) or N9 (purines) positions
        - Ring atoms
        - H-bond donor/acceptor atoms

        Raises:
            FileNotFoundError: If pdb_atoms JSON does not exist.
            ValueError: If JSON format is invalid.
        """
        atoms_file = self.json_dir / "pdb_atoms" / f"{self.pdb_id}.json"

        if not atoms_file.exists():
            raise FileNotFoundError(
                f"Atoms JSON not found: {atoms_file}\n"
                f"Run: ./build/generate_modern_json data/pdb/{self.pdb_id}.pdb "
                f"data/json --stage=atoms"
            )

        with open(atoms_file) as f:
            data = json.load(f)

        if not isinstance(data, list) or len(data) == 0:
            raise ValueError(f"Expected non-empty list in atoms JSON, got {type(data)}")

        if "atoms" not in data[0]:
            raise ValueError("Missing 'atoms' key in JSON data")

        atoms_list = data[0]["atoms"]

        # Group atoms by res_id
        residue_atoms: Dict[str, List[dict]] = {}
        for atom in atoms_list:
            res_id = atom["res_id"]
            if res_id not in residue_atoms:
                residue_atoms[res_id] = []
            residue_atoms[res_id].append(atom)

        # Extract coordinates for each residue
        for res_id, atoms in residue_atoms.items():
            coords = AtomCoords()

            # Get residue name for determining base type
            residue_name = atoms[0]["residue_name"]

            # Determine if purine or pyrimidine
            is_purine = residue_name in ["A", "G", "DA", "DG"]

            for atom in atoms:
                atom_name = atom["atom_name"]
                xyz = np.array(atom["xyz"], dtype=np.float64)

                # C1' position
                if atom_name == "C1'":
                    coords.c1_prime = xyz

                # N1 or N9 position (glycosidic nitrogen)
                if is_purine and atom_name == "N9":
                    coords.n1_or_n9 = xyz
                elif not is_purine and atom_name == "N1":
                    coords.n1_or_n9 = xyz

                # Ring atoms (all carbon and nitrogen in base)
                if atom_name in ["C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"]:
                    coords.ring_atoms[atom_name] = xyz

                # H-bond donors/acceptors (O, N atoms in base)
                if atom_name in [
                    "N1",
                    "N2",
                    "N3",
                    "N4",
                    "N6",
                    "N7",
                    "O2",
                    "O4",
                    "O6",
                ]:
                    coords.hbond_atoms[atom_name] = xyz

            self.atoms[res_id] = coords

    def _extract_residue_name(self, res_id: str) -> str:
        """Extract single-letter residue name from res_id.

        Args:
            res_id: Residue ID (format: "chain-name-seq[ins]").

        Returns:
            Single-letter residue name (e.g., "G", "C").

        Raises:
            ValueError: If res_id format is invalid.
        """
        parts = res_id.split("-")
        if len(parts) < 2:
            raise ValueError(f"Invalid res_id format: {res_id}")

        # Handle multi-letter residue names (e.g., "DA", "DG")
        name = parts[1]
        if name.startswith("D") and len(name) == 2:
            return name[1]  # "DA" -> "A"
        return name

    def get_valid_pairs(self) -> List[CachedPair]:
        """Return pairs where validation.is_valid is True.

        Returns:
            List of CachedPair instances with valid geometry.
        """
        return [pair for pair in self.pairs if pair.validation.is_valid]

    def to_json(self) -> Dict:
        """Serialize cache to JSON-compatible dictionary.

        Returns:
            Dictionary with pdb_id, frames, atoms, and pairs.
        """
        # Convert frames to JSON
        frames_json = {}
        for res_id, frame in self.frames.items():
            frames_json[res_id] = {
                "origin": frame.origin.tolist(),
                "rotation": frame.rotation.tolist(),
                "rmsd_fit": frame.rmsd_fit,
            }

        # Convert atoms to JSON
        atoms_json = {}
        for res_id, coords in self.atoms.items():
            atoms_json[res_id] = {
                "c1_prime": coords.c1_prime.tolist()
                if coords.c1_prime is not None
                else None,
                "n1_or_n9": coords.n1_or_n9.tolist()
                if coords.n1_or_n9 is not None
                else None,
                "ring_atoms": {k: v.tolist() for k, v in coords.ring_atoms.items()},
                "hbond_atoms": {k: v.tolist() for k, v in coords.hbond_atoms.items()},
            }

        # Convert pairs to JSON
        pairs_json = []
        for pair in self.pairs:
            pair_dict = {
                "res1_id": pair.res1_id,
                "res2_id": pair.res2_id,
                "res1_name": pair.res1_name,
                "res2_name": pair.res2_name,
                "frame1": {
                    "origin": pair.frame1.origin.tolist(),
                    "rotation": pair.frame1.rotation.tolist(),
                    "rmsd_fit": pair.frame1.rmsd_fit,
                },
                "frame2": {
                    "origin": pair.frame2.origin.tolist(),
                    "rotation": pair.frame2.rotation.tolist(),
                    "rmsd_fit": pair.frame2.rmsd_fit,
                },
                "validation": {
                    "dorg": pair.validation.dorg,
                    "d_v": pair.validation.d_v,
                    "plane_angle": pair.validation.plane_angle,
                    "dNN": pair.validation.dNN,
                    "dir_x": pair.validation.dir_x,
                    "dir_y": pair.validation.dir_y,
                    "dir_z": pair.validation.dir_z,
                    "quality_score": pair.validation.quality_score,
                    "distance_check": pair.validation.distance_check,
                    "d_v_check": pair.validation.d_v_check,
                    "plane_angle_check": pair.validation.plane_angle_check,
                    "dNN_check": pair.validation.dNN_check,
                    "is_valid": pair.validation.is_valid,
                },
                "hbonds": pair.hbonds,
            }
            pairs_json.append(pair_dict)

        return {
            "pdb_id": self.pdb_id,
            "frames": frames_json,
            "atoms": atoms_json,
            "pairs": pairs_json,
        }

    def save(self, path: Path | str) -> None:
        """Save cache to JSON file.

        Args:
            path: Path to output JSON file.
        """
        path = Path(path)
        with open(path, "w") as f:
            json.dump(self.to_json(), f, indent=2)

    @classmethod
    def load(cls, path: Path | str) -> "PairCache":
        """Load cache from JSON file.

        Args:
            path: Path to input JSON file.

        Returns:
            PairCache instance with loaded data.

        Raises:
            FileNotFoundError: If file does not exist.
            ValueError: If JSON format is invalid.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Cache file not found: {path}")

        with open(path) as f:
            data = json.load(f)

        # Create cache instance
        cache = cls(pdb_id=data["pdb_id"], json_dir=Path("."))

        # Load frames
        for res_id, frame_data in data["frames"].items():
            cache.frames[res_id] = ReferenceFrame(
                origin=np.array(frame_data["origin"], dtype=np.float64),
                rotation=np.array(frame_data["rotation"], dtype=np.float64),
                rmsd_fit=frame_data["rmsd_fit"],
            )

        # Load atoms
        for res_id, atom_data in data["atoms"].items():
            coords = AtomCoords()
            if atom_data["c1_prime"] is not None:
                coords.c1_prime = np.array(atom_data["c1_prime"], dtype=np.float64)
            if atom_data["n1_or_n9"] is not None:
                coords.n1_or_n9 = np.array(atom_data["n1_or_n9"], dtype=np.float64)
            coords.ring_atoms = {
                k: np.array(v, dtype=np.float64)
                for k, v in atom_data["ring_atoms"].items()
            }
            coords.hbond_atoms = {
                k: np.array(v, dtype=np.float64)
                for k, v in atom_data["hbond_atoms"].items()
            }
            cache.atoms[res_id] = coords

        # Load pairs
        for pair_data in data["pairs"]:
            frame1 = ReferenceFrame(
                origin=np.array(pair_data["frame1"]["origin"], dtype=np.float64),
                rotation=np.array(pair_data["frame1"]["rotation"], dtype=np.float64),
                rmsd_fit=pair_data["frame1"]["rmsd_fit"],
            )
            frame2 = ReferenceFrame(
                origin=np.array(pair_data["frame2"]["origin"], dtype=np.float64),
                rotation=np.array(pair_data["frame2"]["rotation"], dtype=np.float64),
                rmsd_fit=pair_data["frame2"]["rmsd_fit"],
            )

            val_data = pair_data["validation"]
            validation = ValidationResult(
                dorg=val_data["dorg"],
                d_v=val_data["d_v"],
                plane_angle=val_data["plane_angle"],
                dNN=val_data["dNN"],
                dir_x=val_data["dir_x"],
                dir_y=val_data["dir_y"],
                dir_z=val_data["dir_z"],
                quality_score=val_data["quality_score"],
                distance_check=val_data["distance_check"],
                d_v_check=val_data["d_v_check"],
                plane_angle_check=val_data["plane_angle_check"],
                dNN_check=val_data["dNN_check"],
                is_valid=val_data["is_valid"],
            )

            pair = CachedPair(
                res1_id=pair_data["res1_id"],
                res2_id=pair_data["res2_id"],
                res1_name=pair_data["res1_name"],
                res2_name=pair_data["res2_name"],
                frame1=frame1,
                frame2=frame2,
                validation=validation,
                hbonds=pair_data.get("hbonds", []),
            )
            cache.pairs.append(pair)

        return cache
