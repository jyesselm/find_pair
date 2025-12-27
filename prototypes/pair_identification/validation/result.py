"""Validation result dataclass for base pair geometric checks."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ValidationResult:
    """Result of geometric validation between two base frames.

    This dataclass contains all geometric metrics computed during base pair
    validation, including distances, angles, direction cosines, and pass/fail
    status for each check.

    Attributes:
        dorg: Distance between frame origins (Angstroms).
        d_v: Vertical distance along average helix axis (Angstroms).
        plane_angle: Angle between base plane normals (degrees, 0-90).
        dNN: Distance between glycosidic N1/N9 atoms (Angstroms).
        dir_x: Dot product of x-axes (base long axis alignment).
        dir_y: Dot product of y-axes (base short axis alignment).
        dir_z: Dot product of z-axes (base plane normal alignment).
        quality_score: Weighted metric combining dorg, d_v, and plane_angle.
        distance_check: True if dorg passes threshold.
        d_v_check: True if d_v passes threshold.
        plane_angle_check: True if plane_angle passes threshold.
        dNN_check: True if dNN passes threshold.
        is_valid: True if all geometric checks pass.
        rejection_reason: Optional description of why validation failed.

    Example:
        result = ValidationResult(
            dorg=8.5, d_v=1.2, plane_angle=15.0, dNN=8.8,
            dir_x=0.95, dir_y=0.92, dir_z=0.98,
            quality_score=10.3,
            distance_check=True, d_v_check=True,
            plane_angle_check=True, dNN_check=True,
            is_valid=True
        )
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
    rejection_reason: Optional[str] = None

    def get_failed_checks(self) -> list[str]:
        """Get list of failed validation checks.

        Returns:
            List of check names that failed, empty if all passed.

        Example:
            >>> result = ValidationResult(...)
            >>> result.get_failed_checks()
            ['dorg', 'plane_angle']
        """
        failed = []
        if not self.distance_check:
            failed.append("dorg")
        if not self.d_v_check:
            failed.append("d_v")
        if not self.plane_angle_check:
            failed.append("plane_angle")
        if not self.dNN_check:
            failed.append("dNN")
        return failed
