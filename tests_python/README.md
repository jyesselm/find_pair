# Python Tests

This directory contains Python tests for validating X3DNA modern C++ implementation against legacy C code.

**For complete testing documentation, see [../docs/TESTING_GUIDE.md](../docs/TESTING_GUIDE.md).**

## Quick Reference

### Running Tests

```bash
# Run all tests
pytest tests_python/ -v

# Run stage validation
pytest tests_python/integration/test_stage_validation.py -v

# Run specific stage
pytest tests_python/integration/test_stage_validation.py -v -k "stage3"

# Limit PDBs tested
pytest tests_python/integration/test_stage_validation.py -v --max-pdbs=50

# Test single PDB
pytest tests_python/integration/test_stage_validation.py -v --pdb=1EHZ
```

### CLI Mode

The stage validation script can also run directly:

```bash
# Run stage 1 (atoms) on 100 PDBs
python tests_python/integration/test_stage_validation.py 1 --max-pdbs 100

# Run all frame stages with verbose output
python tests_python/integration/test_stage_validation.py frames -v

# Run all stages on a single PDB
python tests_python/integration/test_stage_validation.py all --pdb 1EHZ
```

## Directory Structure

```
tests_python/
├── conftest.py          # Shared pytest fixtures
├── README.md            # This file
├── integration/         # Integration tests (require executables + data)
│   ├── test_stage_validation.py   # Main stage validation (RECOMMENDED)
│   ├── test_atoms_batch.py        # Atom comparison tests
│   ├── test_frames_batch.py       # Frame calculation tests
│   ├── test_ls_fitting_batch.py   # LS fitting tests
│   └── test_residue_indices_batch.py
└── unit/                # Unit tests (no external dependencies)
```

## Validation Stages

| Stage | ID | Description |
|-------|-----|-------------|
| 1 | `pdb_atoms` | Atom parsing |
| 2 | `residue_indices` | Residue index mapping |
| 3 | `base_frame_calc` | Base frame calculation |
| 4 | `ls_fitting` | Least squares fitting |
| 5 | `frame_calc` | Reference frame calculation |
| 6 | `pair_validation` | Pair validation |
| 7 | `distance_checks` | Distance measurements |
| 8 | `hbond_list` | H-bond detection |
| 9 | `base_pair` | Base pair records |
| 10 | `find_bestpair_selection` | Final pair selection (PRIMARY) |
| 11 | `bpstep_params` | Step parameters |
| 12 | `helical_params` | Helical parameters |

## Stage Groups

| Group | Stages |
|-------|--------|
| `atoms` | 1 |
| `residue` | 2 |
| `frames` | 3, 4, 5 |
| `pairs` | 6, 7, 9, 10 |
| `hbonds` | 8 |
| `steps` | 11, 12 |
| `all` | 1-12 |

## Prerequisites

1. **Build executables:**
   ```bash
   # Modern
   cd build && cmake .. && make -j
   
   # Legacy
   cd org/build && cmake .. && make -j
   ```

2. **Install pytest:**
   ```bash
   pip install -e ".[dev]"
   ```

3. **Verify environment:**
   ```bash
   fp2-validate info
   ```

## Test Output

Results are saved to `data/validation_results/`:
- `stage{N}_{id}_pytest.json` - pytest results
- `stage{N}_{id}_cli.json` - CLI results

## See Also

- [TESTING_GUIDE.md](../docs/TESTING_GUIDE.md) - Complete testing documentation
- [JSON_DATA_TYPES_AND_COMPARISONS.md](../docs/JSON_DATA_TYPES_AND_COMPARISONS.md) - JSON schema reference
- [x3dna_json_compare/](../x3dna_json_compare/) - Comparison library
