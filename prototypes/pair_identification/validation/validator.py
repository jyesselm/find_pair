"""Geometric validator for base pairs."""

import numpy as np

from prototypes.pair_identification.frame_loader import ReferenceFrame
from prototypes.pair_identification.validation.result import ValidationResult
from prototypes.pair_identification.validation.thresholds import ValidationThresholds


class GeometricValidator:
    """Validates base pair geometry from reference frames.

    This class implements the geometric validation algorithm from the C++
    BasePairValidator, computing distances, angles, and quality scores to
    determine if two bases can form a valid pair.

    Attributes:
        thresholds: ValidationThresholds object with configurable limits.

    Example:
        validator = GeometricValidator()
        result = validator.validate(
            frame1, frame2,
            n1n9_pos1, n1n9_pos2
        )
        if result.is_valid:
            print(f"Quality score: {result.quality_score:.2f}")
    """

    def __init__(self, thresholds: ValidationThresholds | None = None):
        """Initialize geometric validator.

        Args:
            thresholds: Optional custom thresholds. If None, uses defaults.

        Example:
            >>> validator = GeometricValidator()
            >>> strict_validator = GeometricValidator(
            ...     ValidationThresholds.strict()
            ... )
        """
        self.thresholds = thresholds if thresholds is not None else ValidationThresholds.default()

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
        metrics = self._compute_metrics(frame1, frame2, n1n9_pos1, n1n9_pos2)
        checks = self._run_checks(metrics)
        return self._build_result(metrics, checks)

    def _compute_metrics(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
        n1n9_pos1: np.ndarray,
        n1n9_pos2: np.ndarray,
    ) -> dict:
        """Compute all geometric metrics."""
        dorg = self._compute_origin_distance(frame1, frame2)
        dir_x, dir_y, dir_z = self._compute_direction_cosines(frame1, frame2)
        zave = self._compute_average_z_axis(frame1, frame2, dir_z)
        d_v = self._compute_vertical_distance(frame1, frame2, zave)
        plane_angle = self._compute_plane_angle(frame1, frame2)
        dNN = self._compute_n1n9_distance(n1n9_pos1, n1n9_pos2)
        quality_score = self.thresholds.compute_quality_score(dorg, d_v, plane_angle)

        return {
            "dorg": dorg, "d_v": d_v, "plane_angle": plane_angle, "dNN": dNN,
            "dir_x": dir_x, "dir_y": dir_y, "dir_z": dir_z,
            "quality_score": quality_score,
        }

    def _run_checks(self, metrics: dict) -> dict:
        """Run validation checks on metrics."""
        return {
            "distance_check": metrics["dorg"] <= self.thresholds.max_dorg,
            "d_v_check": metrics["d_v"] <= self.thresholds.max_d_v,
            "plane_angle_check": metrics["plane_angle"] <= self.thresholds.max_plane_angle,
            "dNN_check": metrics["dNN"] >= self.thresholds.min_dNN,
        }

    def _build_result(self, metrics: dict, checks: dict) -> ValidationResult:
        """Build ValidationResult from metrics and checks."""
        is_valid = all(checks.values())
        rejection_reason = None if is_valid else self._build_rejection_reason(checks)

        return ValidationResult(
            is_valid=is_valid,
            rejection_reason=rejection_reason,
            **metrics,
            **checks,
        )

    def _compute_origin_distance(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
    ) -> float:
        """Compute distance between frame origins.

        Args:
            frame1: First reference frame.
            frame2: Second reference frame.

        Returns:
            Origin-to-origin distance in Angstroms.
        """
        dorg_vec = frame1.origin - frame2.origin
        return float(np.linalg.norm(dorg_vec))

    def _compute_direction_cosines(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
    ) -> tuple[float, float, float]:
        """Compute dot products of frame axes.

        Args:
            frame1: First reference frame.
            frame2: Second reference frame.

        Returns:
            Tuple of (dir_x, dir_y, dir_z) dot products.
        """
        dir_x = float(np.dot(frame1.x_axis, frame2.x_axis))
        dir_y = float(np.dot(frame1.y_axis, frame2.y_axis))
        dir_z = float(np.dot(frame1.z_axis, frame2.z_axis))
        return dir_x, dir_y, dir_z

    def _compute_average_z_axis(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
        dir_z: float,
    ) -> np.ndarray:
        """Compute average helix axis from base plane normals.

        Matches C++ get_bp_zoave algorithm (lines 147-170).

        Args:
            frame1: First reference frame.
            frame2: Second reference frame.
            dir_z: Dot product of z-axes.

        Returns:
            Normalized average z-axis vector.
        """
        if dir_z > 0:
            zave = frame1.z_axis + frame2.z_axis
        else:
            zave = frame2.z_axis - frame1.z_axis

        zave_len = np.linalg.norm(zave)
        if zave_len > 1e-10:
            return zave / zave_len
        return frame1.z_axis

    def _compute_vertical_distance(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
        zave: np.ndarray,
    ) -> float:
        """Compute vertical distance along average helix axis.

        Args:
            frame1: First reference frame.
            frame2: Second reference frame.
            zave: Average z-axis vector.

        Returns:
            Vertical distance in Angstroms.
        """
        dorg_vec = frame1.origin - frame2.origin
        return float(abs(np.dot(dorg_vec, zave)))

    def _compute_plane_angle(
        self,
        frame1: ReferenceFrame,
        frame2: ReferenceFrame,
    ) -> float:
        """Compute angle between base plane normals.

        Matches C++ z1_z2_angle_in_0_to_90 algorithm (lines 173-186).
        Returns angle in range [0, 90] degrees.

        Args:
            frame1: First reference frame.
            frame2: Second reference frame.

        Returns:
            Plane angle in degrees.
        """
        dot = np.clip(np.dot(frame1.z_axis, frame2.z_axis), -1.0, 1.0)
        angle_rad = np.arccos(abs(dot))
        return float(np.degrees(angle_rad))

    def _compute_n1n9_distance(
        self,
        n1n9_pos1: np.ndarray,
        n1n9_pos2: np.ndarray,
    ) -> float:
        """Compute distance between glycosidic nitrogens.

        Args:
            n1n9_pos1: N1 or N9 position for residue 1.
            n1n9_pos2: N1 or N9 position for residue 2.

        Returns:
            N1/N9 distance in Angstroms.
        """
        return float(np.linalg.norm(n1n9_pos1 - n1n9_pos2))

    def _build_rejection_reason(self, checks: dict) -> str:
        """Build human-readable rejection reason from failed checks.

        Args:
            checks: Dictionary of check results.

        Returns:
            Rejection reason string listing failed checks.
        """
        failed = [name.replace("_check", "") for name, passed in checks.items() if not passed]
        return f"Failed: {', '.join(failed)}"
