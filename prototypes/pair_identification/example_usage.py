"""Example usage of PairCache module."""

from pathlib import Path
from prototypes.pair_identification import PairCache


def main():
    """Demonstrate basic PairCache usage."""
    # Set up paths
    pdb_id = "1EHZ"
    json_dir = Path("data/json")
    output_file = Path("pair_cache_1EHZ.json")

    # Create and build cache
    print(f"Building pair cache for {pdb_id}...")
    cache = PairCache(pdb_id, json_dir)
    cache.build_cache(max_distance=15.0)

    # Show statistics
    total_pairs = len(cache.pairs)
    valid_pairs = cache.get_valid_pairs()
    print(f"\nStatistics:")
    print(f"  Total residues: {len(cache.frames)}")
    print(f"  Total pairs tested: {total_pairs}")
    print(f"  Valid pairs: {len(valid_pairs)} ({100*len(valid_pairs)/total_pairs:.1f}%)")

    # Show validation breakdown
    distance_fails = sum(1 for p in cache.pairs if not p.validation.distance_check)
    d_v_fails = sum(1 for p in cache.pairs if not p.validation.d_v_check)
    plane_fails = sum(1 for p in cache.pairs if not p.validation.plane_angle_check)
    dNN_fails = sum(1 for p in cache.pairs if not p.validation.dNN_check)

    print(f"\nValidation failures:")
    print(f"  Distance (dorg > 15.0): {distance_fails}")
    print(f"  Vertical distance (d_v > 2.5): {d_v_fails}")
    print(f"  Plane angle (> 65.0 deg): {plane_fails}")
    print(f"  N1/N9 distance (dNN < 4.5): {dNN_fails}")

    # Show best quality pairs
    print(f"\nTop 10 pairs by quality score:")
    sorted_pairs = sorted(valid_pairs, key=lambda p: p.validation.quality_score)
    for i, pair in enumerate(sorted_pairs[:10], 1):
        print(f"  {i:2d}. {pair.res1_id:12s} ({pair.res1_name}) - "
              f"{pair.res2_id:12s} ({pair.res2_name}): "
              f"score={pair.validation.quality_score:.2f}")

    # Save cache
    print(f"\nSaving cache to {output_file}...")
    cache.save(output_file)
    print(f"Done! Cache saved to {output_file}")


if __name__ == "__main__":
    main()
