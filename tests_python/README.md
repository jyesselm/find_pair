# Python Tests

This directory contains Python tests for X3DNA JSON generation and comparison.

**All Python tests use pytest.** This is enforced by project rules - see `.cursorrules`.

## Structure

- `integration/` - Integration tests that require executables and data files
- `unit/` - Unit tests for utility functions
- `conftest.py` - Shared pytest fixtures and configuration
- `README.md` - This file

## Running Tests

### Run all tests
```bash
pytest tests_python/
```

### Run only integration tests
```bash
pytest tests_python/integration/
```

### Run a specific test file
```bash
pytest tests_python/integration/test_atoms_batch.py
```

### Run tests with verbose output
```bash
pytest tests_python/ -v
```

### Run tests in parallel (if pytest-xdist is installed)
```bash
pytest tests_python/ -n auto
```

## Test Scripts

The integration tests can also be run as standalone scripts:

```bash
# Test atoms JSON generation
python tests_python/integration/test_atoms_batch.py

# Test residue indices
python tests_python/integration/test_residue_indices_batch.py

# Test ls_fitting (single PDB)
python tests_python/integration/test_ls_fitting.py <PDB_ID>

# Test ls_fitting (batch - all fast PDBs)
python tests_python/integration/test_ls_fitting_batch.py

# Test ls_fitting (batch - limited)
python tests_python/integration/test_ls_fitting_batch.py --max-pdbs 100

# Test frames
python tests_python/integration/test_frames_batch.py
```

### Monitoring Batch Tests

For long-running batch tests, use the monitoring scripts:

```bash
# Monitor ls_fitting batch test
python scripts/monitor_ls_fitting_batch.py

# Watch log file
tail -f data/validation_results/ls_fitting_batch.log
```

## Installation

Install pytest and development dependencies:

```bash
# From project root
pip install -e ".[dev]"
```

This installs:
- `pytest>=7.0` - Testing framework
- `pytest-timeout>=2.0` - Timeout support
- `black>=22.0` - Code formatting (optional)
- `flake8>=5.0` - Linting (optional)

## Requirements

Tests require:
- Python 3.8+
- pytest (install with `pip install -e ".[dev]"`)
- Built executables in `build/` and `org/build/bin/`
- PDB files in `data/pdb/`
- Legacy JSON files in `data/json_legacy/` (optional, will be generated if missing)

## Test Markers

Tests are marked with pytest markers:
- `@pytest.mark.integration` - Integration tests
- `@pytest.mark.slow` - Tests that take a long time
- `@pytest.mark.requires_legacy` - Requires legacy executable
- `@pytest.mark.requires_modern` - Requires modern executable

Run only fast tests:
```bash
pytest tests_python/ -m "not slow"
```

