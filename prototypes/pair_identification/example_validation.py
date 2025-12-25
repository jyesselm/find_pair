"""Example: Validate base pair geometry from C++ JSON output.

This script demonstrates loading reference frames from modern C++ output
and performing geometric validation using the GeometricValidator.

Usage:
    python example_validation.py
"""

from pathlib import Path


from prototypes.pair_identification import FrameLoader, GeometricValidator


def validate_pair_example():
    """Load frames from JSON and validate a pair."""
    # Paths
    json_dir = Path("/Users/jyesselman2/local/code/cpp/find_pair_2/data/json")
    pdb_id = "1EHZ"

    # Load frames
    print(f"Loading frames for {pdb_id}...")
    loader = FrameLoader(json_dir)
    frames = loader.load_frames(pdb_id)
    print(f"Loaded {len(frames)} frames\n")

    # Get two frames to validate
    # 1EHZ has Watson-Crick pairs like A-DG-2 : A-DC-11
    res_ids = list(frames.keys())
    if len(res_ids) < 2:
        print("Not enough frames to validate")
        return

    # Try to find a known pair (or just take first two)
    res_id1 = res_ids[0]
    res_id2 = res_ids[1]

    print(f"Validating pair: {res_id1} <-> {res_id2}")
    frame1 = frames[res_id1]
    frame2 = frames[res_id2]

    # For demonstration, we'll use frame origins as N1/N9 positions
    # (In practice, you'd extract these from atoms JSON)
    n1n9_pos1 = frame1.origin
    n1n9_pos2 = frame2.origin

    # Validate
    validator = GeometricValidator()
    result = validator.validate(frame1, frame2, n1n9_pos1, n1n9_pos2)

    # Print results
    print("\nGeometric Validation Results:")
    print("-" * 60)
    print(
        f"Origin distance (dorg):      {result.dorg:8.3f} A  {'PASS' if result.distance_check else 'FAIL'}"
    )
    print(
        f"Vertical distance (d_v):     {result.d_v:8.3f} A  {'PASS' if result.d_v_check else 'FAIL'}"
    )
    print(
        f"Plane angle:                 {result.plane_angle:8.3f} °  {'PASS' if result.plane_angle_check else 'FAIL'}"
    )
    print(
        f"N1/N9 distance (dNN):        {result.dNN:8.3f} A  {'PASS' if result.dNN_check else 'FAIL'}"
    )
    print(f"Quality score:               {result.quality_score:8.3f}")
    print("-" * 60)
    print("Direction vectors:")
    print(f"  dir_x (x-axes):            {result.dir_x:8.3f}")
    print(f"  dir_y (y-axes):            {result.dir_y:8.3f}")
    print(f"  dir_z (z-axes):            {result.dir_z:8.3f}")
    print("-" * 60)
    print(f"Overall validation:          {'VALID' if result.is_valid else 'INVALID'}")
    print()

    # Try a few more pairs
    print("\nValidating multiple pairs:")
    print("-" * 80)
    print(
        f"{'Res1':<15} {'Res2':<15} {'dorg':>8} {'d_v':>8} {'angle':>8} {'dNN':>8} {'Valid'}"
    )
    print("-" * 80)

    count = 0
    for i, res_id1 in enumerate(res_ids[:10]):
        for res_id2 in res_ids[i + 1 : 10]:
            frame1 = frames[res_id1]
            frame2 = frames[res_id2]
            n1n9_pos1 = frame1.origin
            n1n9_pos2 = frame2.origin

            result = validator.validate(frame1, frame2, n1n9_pos1, n1n9_pos2)

            print(
                f"{res_id1:<15} {res_id2:<15} {result.dorg:8.2f} {result.d_v:8.2f} "
                f"{result.plane_angle:8.2f} {result.dNN:8.2f} {'  YES' if result.is_valid else '   NO'}"
            )

            count += 1
            if count >= 20:
                break
        if count >= 20:
            break


if __name__ == "__main__":
    validate_pair_example()
