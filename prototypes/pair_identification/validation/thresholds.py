"""Configurable thresholds for geometric validation."""

from dataclasses import dataclass

from prototypes.pair_identification.core.constants import (
    MAX_DORG,
    MAX_D_V,
    MAX_PLANE_ANGLE,
    MIN_DNN,
    D_V_WEIGHT,
    PLANE_ANGLE_DIVISOR,
)


@dataclass
class ValidationThresholds:
    """Thresholds for geometric validation checks.

    This class encapsulates all configurable parameters used in base pair
    geometric validation, including distance limits, angle limits, and
    quality score weights.

    Attributes:
        max_dorg: Maximum origin distance threshold (Angstroms).
        max_d_v: Maximum vertical distance threshold (Angstroms).
        max_plane_angle: Maximum plane angle threshold (degrees).
        min_dNN: Minimum N1/N9 distance threshold (Angstroms).
        d_v_weight: Weight for d_v in quality score calculation.
        plane_angle_divisor: Divisor for plane_angle in quality score.

    Example:
        thresholds = ValidationThresholds.strict()
        validator = GeometricValidator(thresholds)
    """

    max_dorg: float = MAX_DORG
    max_d_v: float = MAX_D_V
    max_plane_angle: float = MAX_PLANE_ANGLE
    min_dNN: float = MIN_DNN
    d_v_weight: float = D_V_WEIGHT
    plane_angle_divisor: float = PLANE_ANGLE_DIVISOR

    @classmethod
    def default(cls) -> "ValidationThresholds":
        """Create default thresholds matching legacy X3DNA behavior.

        Returns:
            ValidationThresholds with default values.

        Example:
            >>> thresholds = ValidationThresholds.default()
            >>> thresholds.max_dorg
            15.0
        """
        return cls()

    @classmethod
    def strict(cls) -> "ValidationThresholds":
        """Create stricter thresholds for higher quality pair detection.

        Returns:
            ValidationThresholds with tighter limits.

        Example:
            >>> thresholds = ValidationThresholds.strict()
            >>> thresholds.max_dorg
            12.0
        """
        return cls(
            max_dorg=12.0,
            max_d_v=2.0,
            max_plane_angle=45.0,
        )

    @classmethod
    def relaxed(cls) -> "ValidationThresholds":
        """Create more permissive thresholds for flexible structures.

        Returns:
            ValidationThresholds with looser limits.

        Example:
            >>> thresholds = ValidationThresholds.relaxed()
            >>> thresholds.max_dorg
            18.0
        """
        return cls(
            max_dorg=18.0,
            max_d_v=3.0,
            max_plane_angle=75.0,
        )

    def compute_quality_score(
        self,
        dorg: float,
        d_v: float,
        plane_angle: float,
    ) -> float:
        """Compute quality score from geometric metrics.

        Quality score is a weighted combination of origin distance,
        vertical distance, and plane angle. Lower scores indicate
        better geometric fit.

        Args:
            dorg: Origin distance (Angstroms).
            d_v: Vertical distance (Angstroms).
            plane_angle: Plane angle (degrees).

        Returns:
            Quality score (lower is better).

        Example:
            >>> thresholds = ValidationThresholds.default()
            >>> score = thresholds.compute_quality_score(8.5, 1.2, 15.0)
            >>> score
            10.383333333333333
        """
        return dorg + self.d_v_weight * d_v + plane_angle / self.plane_angle_divisor
