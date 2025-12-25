"""Quick test to verify geometric_validator.py works correctly."""

import numpy as np

from prototypes.pair_identification.frame_loader import ReferenceFrame
from prototypes.pair_identification.geometric_validator import GeometricValidator


def test_basic_validation():
    """Test geometric validation with simple test frames."""
    # Create two test frames (identical)
    origin1 = np.array([0.0, 0.0, 0.0])
    rotation1 = np.eye(3)  # Identity rotation
    frame1 = ReferenceFrame(origin=origin1, rotation=rotation1, rmsd_fit=0.1)

    # Frame 2 shifted along x-y plane (2 A along x, 1 A along y)
    # This keeps d_v low while providing reasonable dorg
    origin2 = np.array([2.0, 1.0, 0.5])
    rotation2 = np.eye(3)
    frame2 = ReferenceFrame(origin=origin2, rotation=rotation2, rmsd_fit=0.1)

    # N1/N9 positions (typical for Watson-Crick pair)
    n1n9_pos1 = np.array([0.0, 0.0, 0.0])
    n1n9_pos2 = np.array([2.0, 11.5, 0.5])  # ~10.8 Angstrom separation

    # Validate
    validator = GeometricValidator()
    result = validator.validate(frame1, frame2, n1n9_pos1, n1n9_pos2)

    print("Geometric Validation Test:")
    print(f"  dorg: {result.dorg:.2f} A (check: {result.distance_check})")
    print(f"  d_v: {result.d_v:.2f} A (check: {result.d_v_check})")
    print(
        f"  plane_angle: {result.plane_angle:.2f} deg (check: {result.plane_angle_check})"
    )
    print(f"  dNN: {result.dNN:.2f} A (check: {result.dNN_check})")
    print(f"  Quality score: {result.quality_score:.3f}")
    print(
        f"  Direction vectors: x={result.dir_x:.3f}, y={result.dir_y:.3f}, z={result.dir_z:.3f}"
    )
    print(f"  Is valid: {result.is_valid}")

    # Expected values
    expected_dorg = np.sqrt(2.0**2 + 1.0**2 + 0.5**2)  # ~2.29
    assert abs(result.dorg - expected_dorg) < 0.01, (
        f"Expected dorg ~{expected_dorg}, got {result.dorg}"
    )
    assert result.d_v < 1.0, f"Expected d_v < 1.0, got {result.d_v}"
    assert abs(result.plane_angle - 0.0) < 0.01, (
        f"Expected plane_angle ~0, got {result.plane_angle}"
    )
    expected_dNN = np.sqrt(2.0**2 + 11.5**2 + 0.5**2)  # ~11.75
    assert abs(result.dNN - expected_dNN) < 0.1, (
        f"Expected dNN ~{expected_dNN}, got {result.dNN}"
    )
    assert result.dir_x == 1.0, f"Expected dir_x=1.0, got {result.dir_x}"
    assert result.dir_y == 1.0, f"Expected dir_y=1.0, got {result.dir_y}"
    assert result.dir_z == 1.0, f"Expected dir_z=1.0, got {result.dir_z}"
    assert result.is_valid, "Expected validation to pass"

    print("\nTest passed!")


if __name__ == "__main__":
    test_basic_validation()
