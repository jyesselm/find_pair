"""Load reference frames from modern C++ JSON output.

This module provides utilities to load pre-computed reference frames from the
modern C++ find_pair implementation's JSON output. Frames are computed via
least-squares fitting to standard base templates.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import json

import numpy as np


@dataclass
class ReferenceFrame:
    """Reference frame for a nucleotide base.

    The frame is defined by:
    - origin: Position in 3D space (translation vector)
    - rotation: 3x3 rotation matrix with axes as columns
    - rmsd_fit: RMSD from template fitting (quality metric)

    The z-axis (rotation[:, 2]) is the base plane normal.
    """

    origin: np.ndarray  # 3D position (translation from JSON)
    rotation: np.ndarray  # 3x3 rotation matrix
    rmsd_fit: float  # RMSD from template fitting

    @property
    def x_axis(self) -> np.ndarray:
        """First axis of the reference frame."""
        return self.rotation[:, 0]

    @property
    def y_axis(self) -> np.ndarray:
        """Second axis of the reference frame."""
        return self.rotation[:, 1]

    @property
    def z_axis(self) -> np.ndarray:
        """Third axis of the reference frame (base plane normal)."""
        return self.rotation[:, 2]


class FrameLoader:
    """Load reference frames from modern C++ JSON output.

    Frames are loaded from {json_dir}/ls_fitting/{PDB}.json.
    Each frame is computed by least-squares fitting to standard base templates.

    Example:
        loader = FrameLoader("/path/to/data/json")
        frames = loader.load_frames("1EHZ")
        frame = frames["A-G-1"]  # Access by res_id
    """

    def __init__(self, json_dir: Path | str):
        """Initialize frame loader.

        Args:
            json_dir: Directory containing JSON output (e.g., data/json).
        """
        self.json_dir = Path(json_dir)
        self.frames_dir = self.json_dir / "ls_fitting"

    def load_frames(self, pdb_id: str) -> Dict[str, ReferenceFrame]:
        """Load all reference frames for a PDB structure.

        Args:
            pdb_id: PDB identifier (e.g., "1EHZ").

        Returns:
            Dictionary mapping res_id -> ReferenceFrame.
            res_id format: "chain-name-seq[ins]" (e.g., "A-G-1", "A-C-10A").

        Raises:
            FileNotFoundError: If JSON file does not exist.
            ValueError: If JSON format is invalid or data is missing.
        """
        json_file = self.frames_dir / f"{pdb_id}.json"

        if not json_file.exists():
            raise FileNotFoundError(
                f"Frame JSON not found: {json_file}\n"
                f"Run: ./build/generate_modern_json data/pdb/{pdb_id}.pdb "
                f"data/json --stage=frames"
            )

        with open(json_file) as f:
            data = json.load(f)

        if not isinstance(data, list):
            raise ValueError(f"Expected list of frames, got {type(data)}")

        frames = {}
        for entry in data:
            res_id = self._extract_res_id(entry)
            frame = self._parse_frame(entry)
            frames[res_id] = frame

        return frames

    def _extract_res_id(self, entry: dict) -> str:
        """Extract res_id from JSON entry.

        Args:
            entry: JSON entry with res_id field.

        Returns:
            res_id string (e.g., "A-G-1" or "A-C-10A").

        Raises:
            ValueError: If res_id is missing.
        """
        if "res_id" not in entry:
            raise ValueError(f"Missing res_id in entry: {entry}")
        return entry["res_id"]

    def _parse_frame(self, entry: dict) -> ReferenceFrame:
        """Parse ReferenceFrame from JSON entry.

        Args:
            entry: JSON entry with rotation_matrix, translation, rms_fit.

        Returns:
            ReferenceFrame instance.

        Raises:
            ValueError: If required fields are missing or invalid.
        """
        if "translation" not in entry:
            raise ValueError(f"Missing translation in entry: {entry}")
        if "rotation_matrix" not in entry:
            raise ValueError(f"Missing rotation_matrix in entry: {entry}")
        if "rms_fit" not in entry:
            raise ValueError(f"Missing rms_fit in entry: {entry}")

        # Parse translation (origin)
        translation = np.array(entry["translation"], dtype=np.float64)
        if translation.shape != (3,):
            raise ValueError(
                f"Expected translation shape (3,), got {translation.shape}"
            )

        # Parse rotation matrix
        rotation_list = entry["rotation_matrix"]
        if len(rotation_list) != 3:
            raise ValueError(
                f"Expected 3 rows in rotation_matrix, got {len(rotation_list)}"
            )
        rotation = np.array(rotation_list, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise ValueError(
                f"Expected rotation_matrix shape (3, 3), got {rotation.shape}"
            )

        # Parse RMSD
        rmsd_fit = float(entry["rms_fit"])

        return ReferenceFrame(origin=translation, rotation=rotation, rmsd_fit=rmsd_fit)
