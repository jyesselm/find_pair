"""Test script for frame_loader module."""

from pathlib import Path
import numpy as np

from prototypes.pair_identification import FrameLoader


def test_load_frames():
    """Test loading frames from JSON."""
    # Use the data/json directory
    json_dir = Path(__file__).parent.parent.parent / "data" / "json"
    loader = FrameLoader(json_dir)

    # Load frames for 1BNA (B-DNA dodecamer)
    frames = loader.load_frames("1BNA")

    print(f"Loaded {len(frames)} frames for 1BNA")
    print(f"\nAvailable res_ids: {sorted(frames.keys())}")

    # Check a specific frame
    if "A-DC-1" in frames:
        frame = frames["A-DC-1"]
        print("\nFrame for A-DC-1:")
        print(f"  Origin: {frame.origin}")
        print(f"  RMSD fit: {frame.rmsd_fit:.6f}")
        print(f"  X-axis: {frame.x_axis}")
        print(f"  Y-axis: {frame.y_axis}")
        print(f"  Z-axis (normal): {frame.z_axis}")

        # Verify rotation matrix properties
        print("\nRotation matrix properties:")
        print(f"  Shape: {frame.rotation.shape}")
        print(f"  Det(R): {np.linalg.det(frame.rotation):.6f} (should be ~1)")

        # Check orthonormality
        x_dot_y = np.dot(frame.x_axis, frame.y_axis)
        x_dot_z = np.dot(frame.x_axis, frame.z_axis)
        y_dot_z = np.dot(frame.y_axis, frame.z_axis)
        print(f"  X·Y: {x_dot_y:.6f} (should be ~0)")
        print(f"  X·Z: {x_dot_z:.6f} (should be ~0)")
        print(f"  Y·Z: {y_dot_z:.6f} (should be ~0)")

        x_norm = np.linalg.norm(frame.x_axis)
        y_norm = np.linalg.norm(frame.y_axis)
        z_norm = np.linalg.norm(frame.z_axis)
        print(f"  ||X||: {x_norm:.6f} (should be ~1)")
        print(f"  ||Y||: {y_norm:.6f} (should be ~1)")
        print(f"  ||Z||: {z_norm:.6f} (should be ~1)")


def test_error_handling():
    """Test error handling for missing files."""
    json_dir = Path(__file__).parent.parent.parent / "data" / "json"
    loader = FrameLoader(json_dir)

    print("\n\nTesting error handling:")
    try:
        loader.load_frames("NONEXISTENT_PDB")
        print("ERROR: Should have raised FileNotFoundError")
    except FileNotFoundError as e:
        print(f"Correctly raised FileNotFoundError: {e}")


if __name__ == "__main__":
    test_load_frames()
    test_error_handling()
