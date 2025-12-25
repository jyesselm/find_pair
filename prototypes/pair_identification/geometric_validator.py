"""Geometric validation for base pairs.

This module computes geometric validation metrics from two ReferenceFrame objects,
matching the C++ BasePairValidator algorithm exactly.

Validation checks:
    1. Origin distance (dorg): Distance between frame origins
    2. Vertical distance (d_v): Projection onto average helix axis
    3. Plane angle: Angle between base plane normals (0-90 degrees)
    4. N1/N9 distance (dNN): Glycosidic nitrogen distance
    5. Quality score: Weighted combination of metrics
"""

from dataclasses import dataclass

import numpy as np

from prototypes.pair_identification.frame_loader import ReferenceFrame


# Validation constants (from C++ BasePairValidator)
MAX_DORG = 15.0
MAX_D_V = 2.5
MAX_PLANE_ANGLE = 65.0
MIN_DNN = 4.5
OVERLAP_THRESHOLD = 0.01
D_V_WEIGHT = 1.5
PLANE_ANGLE_DIVISOR = 180.0


@dataclass
class ValidationResult:
    """Result of geometric validation between two base pairs.

    Attributes:
        dorg: Distance between frame origins (Angstroms).
        d_v: Vertical distance along helix axis (Angstroms).
        plane_angle: Angle between base plane normals (degrees, 0-90).
        dNN: Distance between N1/N9 atoms (Angstroms).
        dir_x: Dot product of x-axes.
        dir_y: Dot product of y-axes.
        dir_z: Dot product of z-axes.
        quality_score: Weighted metric for pair quality.
        distance_check: True if dorg is valid.
        d_v_check: True if d_v is valid.
        plane_angle_check: True if plane_angle is valid.
        dNN_check: True if dNN is valid.
        is_valid: True if all geometric checks pass.
    """

    dorg: float
    d_v: float
    plane_angle: float
    dNN: float
    dir_x: float
    dir_y: float
    dir_z: float
    quality_score: float
    distance_check: bool
    d_v_check: bool
    plane_angle_check: bool
    dNN_check: bool
    is_valid: bool


class GeometricValidator:
    """Validates base pair geometry from reference frames.

    This class implements the geometric validation algorithm from the C++
    BasePairValidator, computing distances, angles, and quality scores.

    Example:
        validator = GeometricValidator()
        result = validator.validate(
            frame1, frame2,
            n1n9_pos1, n1n9_pos2
        )
        if result.is_valid:
            print(f"Quality score: {result.quality_score:.2f}")
    """

    def __init__(
        self,
        max_dorg: float = MAX_DORG,
        max_d_v: float = MAX_D_V,
        max_plane_angle: float = MAX_PLANE_ANGLE,
        min_dNN: float = MIN_DNN,
    ):
        """Initialize geometric validator.

        Args:
            max_dorg: Maximum origin distance threshold.
            max_d_v: Maximum vertical distance threshold.
            max_plane_angle: Maximum plane angle threshold.
            min_dNN: Minimum N1/N9 distance threshold.
        """
        self.max_dorg = max_dorg
        self.max_d_v = max_d_v
        self.max_plane_angle = max_plane_angle
        self.min_dNN = min_dNN

    def validate(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
        n1n9_pos1: np.ndarray,
        n1n9_pos2: np.ndarray,
    ) -> ValidationResult:
        """Validate geometry between two base pair frames.

        Args:
            frame1: Reference frame for first residue.
            frame2: Reference frame for second residue.
            n1n9_pos1: N1 or N9 position for residue 1 (3D array).
            n1n9_pos2: N1 or N9 position for residue 2 (3D array).

        Returns:
            ValidationResult with all computed metrics and checks.
        """
        # 1. Origin distance (matches C++ line 57: dorg = frame1.origin() - frame2.origin())
        dorg_vec = frame1.origin - frame2.origin
        dorg = float(np.linalg.norm(dorg_vec))

        # 2. Direction vectors (matches C++ lines 135-137: dot products of axes)
        dir_x = float(np.dot(frame1.x_axis, frame2.x_axis))
        dir_y = float(np.dot(frame1.y_axis, frame2.y_axis))
        dir_z = float(np.dot(frame1.z_axis, frame2.z_axis))

        # 3. Average z-axis (matches C++ get_bp_zoave, lines 147-170)
        if dir_z > 0:
            zave = frame1.z_axis + frame2.z_axis
        else:
            zave = frame2.z_axis - frame1.z_axis

        zave_len = np.linalg.norm(zave)
        if zave_len > 1e-10:
            zave = zave / zave_len
        else:
            # Fallback if vectors cancel out
            zave = frame1.z_axis

        # 4. Vertical distance (matches C++ line 64: d_v = abs(dorg.dot(zave)))
        d_v = float(abs(np.dot(dorg_vec, zave)))

        # 5. Plane angle (matches C++ z1_z2_angle_in_0_to_90, lines 173-186)
        dot = np.clip(np.dot(frame1.z_axis, frame2.z_axis), -1.0, 1.0)
        angle_rad = np.arccos(abs(dot))
        plane_angle = float(np.degrees(angle_rad))

        # 6. N1/N9 distance (matches C++ lines 70-77)
        dNN = float(np.linalg.norm(n1n9_pos1 - n1n9_pos2))

        # 7. Quality score (matches C++ lines 80-81)
        quality_score = dorg + D_V_WEIGHT * d_v + plane_angle / PLANE_ANGLE_DIVISOR

        # 8. Validation checks (matches C++ lines 84-87)
        distance_check = dorg <= self.max_dorg
        d_v_check = d_v <= self.max_d_v
        plane_angle_check = plane_angle <= self.max_plane_angle
        dNN_check = dNN >= self.min_dNN

        # All geometric checks must pass
        is_valid = distance_check and d_v_check and plane_angle_check and dNN_check

        return ValidationResult(
            dorg=dorg,
            d_v=d_v,
            plane_angle=plane_angle,
            dNN=dNN,
            dir_x=dir_x,
            dir_y=dir_y,
            dir_z=dir_z,
            quality_score=quality_score,
            distance_check=distance_check,
            d_v_check=d_v_check,
            plane_angle_check=plane_angle_check,
            dNN_check=dNN_check,
            is_valid=is_valid,
        )
