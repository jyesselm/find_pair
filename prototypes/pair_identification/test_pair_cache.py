"""Test script for PairCache module."""

from pathlib import Path
from prototypes.pair_identification.pair_cache import PairCache


def test_pair_cache():
    """Test PairCache on 1EHZ structure."""
    # Set up paths
    pdb_id = "1EHZ"
    json_dir = Path("/Users/jyesselman2/local/code/cpp/find_pair_2/data/json")

    # Create cache
    print(f"Building pair cache for {pdb_id}...")
    cache = PairCache(pdb_id, json_dir)

    # Build cache
    cache.build_cache(max_distance=15.0)

    # Print summary
    print(f"\nCache summary:")
    print(f"  Total residues: {len(cache.frames)}")
    print(f"  Residues with atoms: {len(cache.atoms)}")
    print(f"  Total pairs: {len(cache.pairs)}")

    # Get valid pairs
    valid_pairs = cache.get_valid_pairs()
    print(f"  Valid pairs: {len(valid_pairs)}")

    # Show first few valid pairs
    print(f"\nFirst 5 valid pairs:")
    for i, pair in enumerate(valid_pairs[:5]):
        print(f"  {i+1}. {pair.res1_id} ({pair.res1_name}) - {pair.res2_id} ({pair.res2_name})")
        print(f"     dorg={pair.validation.dorg:.2f}, d_v={pair.validation.d_v:.2f}, "
              f"plane_angle={pair.validation.plane_angle:.2f}, dNN={pair.validation.dNN:.2f}")
        print(f"     quality_score={pair.validation.quality_score:.2f}")

    # Test save/load
    cache_file = Path("/tmp/test_cache.json")
    print(f"\nSaving cache to {cache_file}...")
    cache.save(cache_file)

    print(f"Loading cache from {cache_file}...")
    loaded_cache = PairCache.load(cache_file)

    print(f"\nLoaded cache summary:")
    print(f"  Total residues: {len(loaded_cache.frames)}")
    print(f"  Residues with atoms: {len(loaded_cache.atoms)}")
    print(f"  Total pairs: {len(loaded_cache.pairs)}")
    print(f"  Valid pairs: {len(loaded_cache.get_valid_pairs())}")

    # Verify data integrity
    assert len(cache.pairs) == len(loaded_cache.pairs), "Pair count mismatch"
    assert len(cache.frames) == len(loaded_cache.frames), "Frame count mismatch"
    assert len(cache.atoms) == len(loaded_cache.atoms), "Atom count mismatch"

    print("\nAll tests passed!")


if __name__ == "__main__":
    test_pair_cache()
