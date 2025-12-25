# Pair Identification Prototype

Python prototype for base pair identification algorithms, using pre-computed reference frames from the modern C++ implementation.

## Purpose

This prototype allows experimentation with pair identification strategies without modifying the C++ codebase. It reuses the high-quality reference frames computed by the modern C++ `BaseFrameCalculator`.

## Modules

- `frame_loader.py` - Load reference frames from JSON output
- `geometric_validator.py` - Geometric validation for base pairs
- `pair_cache.py` - Cache pre-computed frames and validation results
- `test_frame_loader.py` - Test script for frame loading
- `test_pair_cache.py` - Test script for pair caching
- `example_usage.py` - Example usage of PairCache

## Usage

### Basic Example

```python
from pathlib import Path
from prototypes.pair_identification import PairCache

# Create cache
cache = PairCache("1EHZ", Path("data/json"))
cache.build_cache(max_distance=15.0)

# Get valid pairs
valid_pairs = cache.get_valid_pairs()
print(f"Found {len(valid_pairs)} valid pairs")

# Access pair data
for pair in valid_pairs[:5]:
    print(f"{pair.res1_id} - {pair.res2_id}")
    print(f"  Quality score: {pair.validation.quality_score:.2f}")
```

### Save and Load Cache

```python
# Save cache to file
cache.save(Path("cache.json"))

# Load from file
loaded_cache = PairCache.load(Path("cache.json"))
```

## Running Examples

```bash
# Set PYTHONPATH
export PYTHONPATH=/Users/jyesselman2/local/code/cpp/find_pair_2:$PYTHONPATH

# Run frame loader test
python prototypes/pair_identification/test_frame_loader.py

# Run pair cache test
python prototypes/pair_identification/test_pair_cache.py

# Run usage example
python prototypes/pair_identification/example_usage.py
```

## Data Structures

### ReferenceFrame

Reference frame for a nucleotide base with origin, rotation matrix, and RMSD.

Properties:
- `origin`: 3D position (translation vector)
- `rotation`: 3x3 orthonormal rotation matrix
- `rmsd_fit`: RMSD from template fitting (quality metric)
- `x_axis`, `y_axis`, `z_axis`: Convenience properties for axes
- Z-axis is the base plane normal

### AtomCoords

Coordinates for key atoms of a residue:
- `c1_prime`: C1' position
- `n1_or_n9`: N1 (pyrimidines) or N9 (purines)
- `ring_atoms`: Dictionary of ring atom positions
- `hbond_atoms`: Dictionary of H-bond donor/acceptor positions

### ValidationResult

Geometric validation metrics:
- `dorg`: Distance between frame origins
- `d_v`: Vertical distance along helix axis
- `plane_angle`: Angle between base plane normals
- `dNN`: Distance between N1/N9 atoms
- `quality_score`: Weighted metric for pair quality
- `is_valid`: True if all geometric checks pass

### CachedPair

Cached data for a potential base pair:
- `res1_id`, `res2_id`: Residue identifiers
- `res1_name`, `res2_name`: Single-letter codes
- `frame1`, `frame2`: Reference frames
- `validation`: Geometric validation result
- `hbonds`: List of H-bond dicts (for future use)

## Algorithm

1. Load pre-computed reference frames from ls_fitting JSON
2. Load atom coordinates from pdb_atoms JSON
3. Build KDTree on frame origins for efficient neighbor search
4. For each residue pair within max_distance:
   - Extract frames and N1/N9 positions
   - Validate geometry using GeometricValidator
   - Cache result as CachedPair
5. Filter valid pairs using `get_valid_pairs()`

## Data Format

Frames are loaded from `data/json/ls_fitting/{PDB}.json`:

```json
[
    {
        "res_id": "A-DC-1",
        "rotation_matrix": [[...], [...], [...]],
        "translation": [x, y, z],
        "rms_fit": 0.020423,
        "residue_name": "DC",
        "chain_id": "A",
        "residue_seq": 1
    }
]
```

Atoms are loaded from `data/json/pdb_atoms/{PDB}.json`:

```json
[
    {
        "atoms": [
            {
                "atom_idx": 1,
                "atom_name": "N9",
                "chain_id": "A",
                "res_id": "A-G-1",
                "residue_name": "G",
                "residue_seq": 1,
                "xyz": [51.628, 45.992, 53.798]
            }
        ]
    }
]
```

## Requirements

- NumPy
- SciPy (for KDTree)

## Future Extensions

- H-bond detection and scoring
- Pair classification (Watson-Crick, Hoogsteen, etc.)
- Mutual best-pair selection
- Integration with hbond_optimizer prototype
