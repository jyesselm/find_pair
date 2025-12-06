# Testing Guide

**Last Updated**: December 6, 2025  
**Purpose**: Centralized guide for all testing and JSON comparison workflows

---

## Quick Start: Unified CLI (`fp2-validate`)

The **primary tool** for all validation is the `fp2-validate` CLI. This unified tool replaces all the individual validation scripts.

### Installation

```bash
# Install from project root
pip install -e .
```

### Basic Usage

```bash
# Validate everything (all stages, all fast PDBs)
fp2-validate

# Validate specific stages
fp2-validate frames               # Frames only
fp2-validate hbonds               # H-bonds only  
fp2-validate atoms                # Atoms only
fp2-validate pairs                # Pairs only
fp2-validate steps                # Step parameters only

# Validate specific PDB(s)
fp2-validate --pdb 1EHZ -v        # Single PDB, verbose
fp2-validate --pdb 1EHZ --pdb 1BNA --pdb 1H4S  # Multiple PDBs

# Use test sets (predefined collections)
fp2-validate --test-set 10        # Quick 10-PDB test
fp2-validate --test-set 100       # 100-PDB test
fp2-validate --test-set 1000      # Large 1000-PDB test

# Limit number of PDBs
fp2-validate --max 50             # First 50 PDBs

# Debug mode - stop at first failure
fp2-validate --stop-on-first      # Great for debugging

# Document differences to file
fp2-validate --diff               # Saves to data/validation_results/

# Quiet mode (for CI/scripts)
fp2-validate --quiet              # Exit code only: 0=pass, 1=fail

# Show environment info
fp2-validate info                 # Shows executables, PDB counts, test sets

# List available PDBs
fp2-validate list-pdbs            # All fast PDBs
fp2-validate list-pdbs --test-set 50  # PDBs in test set 50
```

### Command Reference

| Command | Description |
|---------|-------------|
| `fp2-validate` | Validate all stages, all fast PDBs |
| `fp2-validate validate [STAGES]` | Validate specific stages |
| `fp2-validate atoms` | Validate atom records |
| `fp2-validate frames` | Validate frame calculations |
| `fp2-validate hbonds` | Validate H-bond lists |
| `fp2-validate pairs` | Validate base pairs |
| `fp2-validate steps` | Validate step parameters |
| `fp2-validate info` | Show environment info |
| `fp2-validate list-pdbs` | List available PDBs |

### Options Reference

| Option | Description |
|--------|-------------|
| `--pdb, -p TEXT` | Specific PDB(s) to validate |
| `--max, -n INT` | Maximum number of PDBs to process |
| `--test-set [10\|50\|100\|500\|1000]` | Use a predefined test set |
| `--workers, -w INT` | Number of parallel workers (default: CPU cores - 1) |
| `--quiet, -q` | Suppress output, exit code only |
| `--verbose, -v` | Show per-PDB results |
| `--stop-on-first, -s` | Stop at first failure |
| `--diff` | Document differences to file |
| `--diff-file PATH` | Custom differences output file |
| `--checkpoint PATH` | Save progress to checkpoint file for resume |
| `--resume` | Skip already-passed PDBs (requires --checkpoint) |
| `--clean-on-match` | Delete modern JSON files that match legacy |

### Common Workflows

```bash
# Quick validation during development
fp2-validate frames --pdb 1EHZ -v

# Full test before commit
fp2-validate --test-set 100

# Debug a failing PDB
fp2-validate --pdb FAILING_PDB --stop-on-first -v

# CI pipeline
fp2-validate --test-set 100 --quiet

# Generate difference report
fp2-validate --test-set 100 --diff

# Long validation run with checkpoint (can resume if interrupted)
fp2-validate validate --test-set 1000 --checkpoint run.json

# Resume from checkpoint
fp2-validate validate --test-set 1000 --checkpoint run.json --resume

# Validate and clean up matched files (save disk space)
fp2-validate pairs --test-set 1000 --clean-on-match
```

---

## Alternative: Python Scripts

For specialized analysis, the `compare_json.py` script is still available:

```bash
# Compare all available PDBs (atoms and frames)
python3 scripts/compare_json.py compare

# Compare specific PDB file(s)
python3 scripts/compare_json.py compare 1H4S
python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

# Compare only atoms
python3 scripts/compare_json.py atoms 1H4S

# Compare only frames
python3 scripts/compare_json.py frames 1H4S

# Compare ring atoms
python3 scripts/compare_json.py ring-atoms 1H4S

# Use test sets (predefined collections of PDBs)
python3 scripts/compare_json.py compare --test-set 100

# Save report to file
python3 scripts/compare_json.py compare --output report.md

# Show only files with differences
python3 scripts/compare_json.py compare --diff-only

# Verbose output with detailed differences
python3 scripts/compare_json.py compare --verbose
```

### What Gets Compared

The `compare_json.py` script compares the following JSON record types:

1. **`pdb_atoms`** - Atom records from PDB parsing
2. **`base_frame_calc`** - Base frame calculation results
3. **`frame_calc`** / **`ref_frame`** - Reference frame calculations
4. **`base_pair`** - Base pair records
5. **`pair_validation`** - Pair validation results
6. **`distance_checks`** - Distance and geometric checks
7. **`hbond_list`** - Hydrogen bond lists
8. **`find_bestpair_selection`** - Selected base pairs
9. **`bpstep_params`** - Step parameters (Shift, Slide, Rise, Tilt, Roll, Twist) ⭐ **NEW**
10. **`helical_params`** - Helical parameters (x_displacement, y_displacement, rise, inclination, tip, twist) ⭐ **NEW**

See [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) for detailed information about each record type and what is checked.

### Comparison Library

The comparison functionality is provided by the `x3dna_json_compare` Python package, which includes:

- **`JsonComparator`** - Main comparison engine
- **`PdbFileReader`** - Efficient PDB file reader
- **`JsonValidator`** - JSON validation utilities

See [x3dna_json_compare/README.md](../x3dna_json_compare/README.md) for API documentation.

---

## Test Types

### 1. Python Integration Tests (pytest)

**Location**: `tests_python/integration/`

**Purpose**: Test JSON generation and comparison workflows

**Framework**: pytest

**Run**:
```bash
# All Python tests
pytest tests_python/

# Only integration tests
pytest tests_python/integration/

# Specific test file
pytest tests_python/integration/test_atoms_batch.py

# With verbose output
pytest tests_python/ -v

# Exclude slow tests
pytest tests_python/ -m "not slow"
```

**Or run as standalone scripts:**
```bash
# Test atoms JSON generation
python tests_python/integration/test_atoms_batch.py

# Test residue indices
python tests_python/integration/test_residue_indices_batch.py

# Test ls_fitting
python tests_python/integration/test_ls_fitting.py

# Test frames
python tests_python/integration/test_frames_batch.py
```

**Key Test Files**:
- `test_atoms_batch.py` - Batch testing of atoms JSON generation and comparison
- `test_residue_indices_batch.py` - Batch testing of residue indices JSON
- `test_ls_fitting.py` - LS fitting JSON generation and comparison
- `test_frames_batch.py` - Batch testing of frames JSON (base_frame_calc, frame_calc, ls_fitting)

**Requirements**:
- Python 3.8+
- pytest (install with `pip install -e ".[dev]"`)
- Built executables in `build/` and `org/build/bin/`
- PDB files in `data/pdb/`
- Legacy JSON files in `data/json_legacy/` (optional, will be generated if missing)

**Shared Fixtures** (from `tests_python/conftest.py`):
- `project_root_path` - Project root directory
- `executables` - Legacy and modern executables
- `temp_output_dir` - Temporary output directory for test files
- `valid_pdbs_fast` - List of valid PDB IDs from `valid_pdbs_fast.json`

**Test Markers**:
- `@pytest.mark.integration` - Integration tests
- `@pytest.mark.slow` - Tests that take a long time
- `@pytest.mark.requires_legacy` - Requires legacy executable
- `@pytest.mark.requires_modern` - Requires modern executable

See [tests_python/README.md](../tests_python/README.md) for complete pytest documentation.

### 2. C++ Unit Tests

**Location**: `tests/unit/`

**Purpose**: Test individual classes and functions in isolation

**Framework**: Google Test

**Run**:
```bash
make test
# or
cd build && ctest
```

**Coverage**: Core classes, geometry, I/O, algorithms

**Key Test Files**:
- `test_atom.cpp` - Atom class tests
- `test_residue.cpp` - Residue class tests
- `test_chain.cpp` - Chain class tests
- `test_structure.cpp` - Structure class tests
- `test_pdb_parser.cpp` - PDB parsing tests
- `test_json_writer.cpp` - JSON writing tests
- `test_base_pair_finder.cpp` - Base pair finding tests
- `test_base_pair_validator.cpp` - Validation tests

### 3. C++ Integration Tests

**Location**: `tests/integration/`

**Purpose**: Test interactions between multiple components

**Framework**: Google Test

**Run**:
```bash
make test
# or
cd build && ctest -R integration
```

**Key Test Files**:
- `test_single_pdb.cpp` - Single PDB processing
- `test_json_generation.cpp` - JSON generation workflow
- `test_io_integration.cpp` - I/O integration
- `test_base_pair_integration.cpp` - Base pair finding integration
- `test_frame_calculation_legacy.cpp` - Frame calculation comparison

### 4. JSON Regression Tests

**Location**: Python scripts using `compare_json.py`

**Purpose**: Ensure output matches legacy JSON files exactly

**Run**:
```bash
# Compare all PDBs
python3 scripts/compare_json.py compare

# Compare test set
python3 scripts/compare_json.py compare --test-set 100

# Compare specific PDBs
python3 scripts/compare_json.py compare 1H4S 2BNA
```

**Tolerance Values**:
- Coordinates: 1e-6
- Distances: 1e-6
- Angles: 1e-6
- Rotation matrices: < 0.0001 per element
- RMS fit: < 0.001

### 5. Legacy Test Tools

**Location**: `org/src/test_*.c`

**Purpose**: Isolated component testing of legacy code

**Build**:
```bash
cd org/build
cmake ..
make
```

**Tools**:
- `test_hbond_initial` - Test initial H-bond detection
- `test_hbond_conflict` - Test conflict resolution
- `test_hbond_validation` - Test H-bond validation

**Usage**:
```bash
cd org/build
./bin/test_hbond_initial ../data/pdb/3G8T.pdb 946 947
./bin/test_hbond_conflict ../data/pdb/3G8T.pdb 946 947
./bin/test_hbond_validation ../data/pdb/3G8T.pdb 946 947
```

See [LEGACY_TEST_TOOLS.md](LEGACY_TEST_TOOLS.md) for detailed documentation.

---

## JSON Comparison Workflow

### Step 1: Generate Modern JSON

```bash
# Generate JSON for a single PDB
./build/generate_modern_json data/pdb/1H4S.pdb data/json/1H4S.json

# Or use the rebuild script
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only

# Generate .inp for downstream tools if needed
./build/find_pair_app data/pdb/1H4S.pdb output.inp
```

### Step 2: Compare with Legacy

```bash
# Quick comparison (atoms, frames, and steps - frames checked first automatically)
python3 scripts/compare_json.py compare 1H4S

# Detailed comparison with report
python3 scripts/compare_json.py compare 1H4S --verbose --output report.md

# Compare specific record types
python3 scripts/compare_json.py atoms 1H4S
python3 scripts/compare_json.py frames 1H4S
python3 scripts/compare_json.py ring-atoms 1H4S

# Compare step parameters (frames are automatically verified first)
python3 scripts/compare_json.py steps 1H4S

# Note: Step parameter files are in split directories:
#   - Modern: data/json/bpstep_params/<PDB_ID>.json
#   - Legacy: data/json_legacy/bpstep_params/<PDB_ID>.json
```

**Important**: When comparing step parameters, frames are automatically checked first since step parameters depend on frames. If frames don't match, step parameter comparisons may be unreliable.

### Step 3: Analyze Differences

If differences are found:

1. **Check the report** - The comparison report shows:
   - Summary statistics
   - Files with differences
   - Detailed mismatch information (with `--verbose`)

2. **Use specialized analysis scripts** (if needed):
   ```bash
   # Analyze H-bond differences
   python3 scripts/compare_hbond_detection.py 1H4S data/json_legacy data/json
   
   # Analyze ring atom matching
   python3 scripts/compare_ring_atoms.py 1H4S
   ```

3. **Debug specific issues**:
   ```bash
   # Debug H-bond mismatches
   python3 scripts/debug_hbond_mismatch.py 1H4S <base_i> <base_j>
   
   # Debug validation mismatches
   python3 scripts/debug_validation_mismatch.py 1H4S <base_i> <base_j>
   ```

### Step 4: Verify Fixes

After making changes:

1. **Regenerate JSON**:
   ```bash
   python3 scripts/rebuild_json.py regenerate 1H4S --modern-only
   ```

2. **Re-run comparison**:
   ```bash
   python3 scripts/compare_json.py compare 1H4S --verbose
   ```

3. **Check for improvements**:
   - Match rate should increase
   - Differences should decrease

---

## Step Parameter Generation

### Using analyze_app

The `analyze_app` can work with both modern and legacy input files:

**Modern input files** (residue indices):
```bash
# Generate input file with modern find_pair
./build/find_pair_app data/pdb/1H4S.pdb /tmp/1H4S.inp

# Generate step parameters
./build/analyze_app /tmp/1H4S.inp
```

**Legacy input files** (atom indices):
```bash
# Use legacy-generated input file directly
# Atom indices are automatically converted to residue indices
./build/analyze_app 1H4S.inp
```

**Note**: Modern `analyze_app` automatically detects and converts atom indices to residue indices when processing legacy input files. You'll see a message like "Converted N atom indices to residue indices" if conversion occurs.

**See**: [ATOM_INDEX_CONVERSION_FIX.md](ATOM_INDEX_CONVERSION_FIX.md) for technical details.

---

## Test Data

### PDB Files

**Location**: `data/pdb/`

**Test Sets**: `data/test_sets/`
- `test_set_100.json` - 100 PDBs for comprehensive testing
- `test_set_10.json` - 10 PDBs for quick testing

**Generate test sets**:
```bash
python3 scripts/compare_json.py generate-test-sets
```

### JSON Files

**Directory Structure**:
- `data/json/` - All modern JSON files (main location)
- `data/json/atoms/` - Atom comparison data (modern)
- `data/json_legacy/` - All legacy JSON files (reference)
- `data/json_legacy/atoms/` - Atom comparison data (legacy)

**Legacy JSON**: `data/json_legacy/`
- Reference output from original code
- Used for regression testing

**Modern JSON**: `data/json/`
- Output from modernized code
- Compared against legacy JSON

**See [DATA_STRUCTURE.md](DATA_STRUCTURE.md) for complete data directory documentation.**

### Valid PDBs List

**Location**: `data/valid_pdbs.json`

**Purpose**: List of PDBs that have valid JSON files for both legacy and modern

---

## Python Tools Reference

### Utility Tools (in `tools/`)

| Tool | Purpose | Usage |
|------|---------|-------|
| `tools/download_pdbs.py` | Download PDB files from RCSB | `python tools/download_pdbs.py --test-set 100` |
| `tools/find_slow_pdbs.py` | Identify PDBs that take >30s to process | `python tools/find_slow_pdbs.py` |

**Download PDBs:**
```bash
# Download PDBs from test set
python tools/download_pdbs.py --test-set 100

# Download specific PDB IDs
python tools/download_pdbs.py 1H4S 2BNA 3DNA

# Download from file containing PDB IDs
python tools/download_pdbs.py --from-file pdb_list.txt
```

**Find Slow PDBs:**
```bash
# Identify PDBs that take more than 30 seconds to generate legacy output
python tools/find_slow_pdbs.py

# Output: data/slow_pdbs.json
# These PDBs can be excluded from batch testing to avoid timeouts
```

## Comparison Tools Reference

### Primary Tools (Use These)

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `scripts/compare_json.py` | **Main comparison tool** | All JSON comparisons |
| `scripts/rebuild_json.py` | Regenerate JSON files | When regenerating test data |
| `x3dna_json_compare` | Comparison library | When writing custom comparison scripts |

### Specialized Analysis Tools (Use When Needed)

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `scripts/compare_hbond_detection.py` | H-bond specific comparison | Debugging H-bond differences |
| `scripts/compare_ring_atoms.py` | Ring atom matching | Debugging frame calculation differences |
| `scripts/debug_hbond_mismatch.py` | Debug H-bond mismatches | Investigating specific H-bond issues |
| `scripts/debug_validation_mismatch.py` | Debug validation mismatches | Investigating validation differences |

### Residue Indexing Tools (For Comparison with Legacy)

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `build/test_residue_matching_by_pdb_props` | Test PDB properties matching | Verifying legacy index derivation |
| `build/fix_residue_indices_from_json` | Fix indices standalone tool | Historical debugging (rarely needed now) |
| `build/check_residue_indices` | Check for duplicate indices | Verifying residue mapping |
| `build/debug_bp_type_id_step_params` | Debug bp_type_id calculation | Debugging bp_type_id issues (auto-fixes indices) |

> **Note**: Residue and atom indices are now assigned directly from the PDB scan order, so the old `--fix-indices` tooling is no longer required during comparisons. The helpers above remain available for forensic debugging only.

### Step Parameter Tools

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `build/analyze_app` | Generate step parameters from input file | After find_pair generates input file |
| `python3 scripts/compare_json.py steps` | Compare step parameters | Validating step parameter calculations |
| **Note**: `analyze_app` automatically handles both modern (residue indices) and legacy (atom indices) input files |

### Reference Frames Comparison Tools

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `build/find_pair_app --legacy-inp=FILE` | Generate ref_frames_modern.dat with legacy ordering | When comparing ref_frames.dat with legacy output |
| `python3 scripts/compare_ref_frames.py` | Compare ref_frames.dat files | Validating ref_frames output matches legacy |
| **Note**: For exact legacy matching, use `--legacy-inp` option to ensure correct strand ordering |

### C++ Comparison Tools (Advanced)

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `tools/compare_hbond_stages.cpp` | Compare H-bond stages | Deep debugging of H-bond algorithm |
| `tools/compare_initial_hbonds.cpp` | Compare initial H-bonds | Debugging initial H-bond detection |
| `tools/compare_atom_selection.cpp` | Compare atom selection | Debugging atom selection differences |

---

## Redundant/Deprecated Tools

The following tools are **redundant** and should **not be used** for new work. They are kept for historical reference but are superseded by `compare_json.py`:

### Deprecated Python Scripts

- ❌ `scripts/quick_compare.py` - Use `compare_json.py compare` instead
- ❌ `scripts/compare_all_pdbs.py` - Use `compare_json.py compare` instead
- ❌ `scripts/detailed_analysis.py` - Use `compare_json.py compare --verbose` instead
- ❌ `scripts/analyze_all_pdbs.py` - Use `compare_json.py compare` instead
- ❌ `scripts/analyze_pair_differences.py` - Use `compare_json.py compare --verbose` instead

### Analysis Scripts (Use Sparingly)

These scripts provide specialized analysis but most functionality is available in `compare_json.py`:

- ⚠️ `scripts/analyze_hbond_differences.py` - Specialized H-bond analysis
- ⚠️ `scripts/analyze_missing_pairs.py` - Specialized pair analysis
- ⚠️ `scripts/investigate_missing_pairs.py` - Specialized investigation

**Recommendation**: Use `compare_json.py` first. Only use specialized scripts if you need functionality not available in the main tool.

---

## Running Tests

### All Tests

```bash
# Run all C++ tests
make test

# Run with verbose output
make test-verbose
```

### Specific Test Categories

```bash
# Run only unit tests
cd build && ctest -R unit

# Run only integration tests
cd build && ctest -R integration

# Run specific test
cd build && ./tests/unit/test_atom
```

### JSON Comparison Tests

```bash
# Compare all PDBs
python3 scripts/compare_json.py compare

# Compare test set
python3 scripts/compare_json.py compare --test-set 100

# Compare specific PDBs
python3 scripts/compare_json.py compare 1H4S 2BNA
```

---

## Test Maintenance

### Adding New Tests

1. **Unit Tests**: Add to `tests/unit/`
   - Follow naming: `test_<class_name>.cpp`
   - Use Google Test framework
   - Add to `tests/unit/CMakeLists.txt`

2. **Integration Tests**: Add to `tests/integration/`
   - Test component interactions
   - Use test fixtures from `test_fixtures.hpp`
   - Add to `tests/integration/CMakeLists.txt`

3. **JSON Regression Tests**: Use `compare_json.py`
   - No code changes needed
   - Just run comparison after changes

### Updating Test Data

1. **Regenerate Legacy JSON** (if needed):
   ```bash
   python3 scripts/rebuild_json.py regenerate --legacy-only
   ```

2. **Regenerate Modern JSON**:
   ```bash
   python3 scripts/rebuild_json.py regenerate --modern-only
   ```

3. **Validate JSON Files**:
   ```bash
   python3 scripts/rebuild_json.py validate
   ```

### Cleaning Test Data

```bash
# Clean invalid/empty JSON files (dry-run)
python3 scripts/rebuild_json.py clean

# Actually remove invalid files
python3 scripts/rebuild_json.py clean --execute
```

---

## Troubleshooting

### Tests Fail

1. **Check test data**: Ensure PDB and JSON files exist
2. **Check tolerance values**: May need adjustment for floating-point differences
3. **Check comparison logic**: Verify comparison functions are correct
4. **Debug step-by-step**: Use verbose output to see intermediate values

### JSON Comparison Issues

1. **Missing files**: Ensure both legacy and modern JSON files exist
2. **Invalid JSON**: Run `python3 scripts/rebuild_json.py validate`
3. **Threading issues**: Use `--threads 1` for single-threaded execution

### Residue Indexing Issues

1. **Wrong residues matched**: Regenerate modern JSON to ensure the latest parser assigned indices from PDB order.

2. **Residue indices don't match legacy**: 
   - Check if legacy JSON file exists: `data/json_legacy/base_frame_calc/<PDB_ID>.json`
   - Use `build/test_residue_matching_by_pdb_props` to verify matching
   - Use `build/check_residue_indices` to check for duplicate indices

3. **Frame origins don't match**: 
   - May be due to residue indexing mismatch in legacy data
   - Use `build/debug_frame_json` to compare frame calculations

### Step Parameter Issues

1. **"Calculated 0 step parameters"**:
   - Check if input file has base pairs
   - Verify frames are calculated (should happen automatically)
   - If using legacy input file, atom indices should be automatically converted
   - Check for conversion message: "Converted N atom indices to residue indices"

2. **Step parameters don't match legacy**:
   - Verify frames match first (step parameters depend on frames)
   - Check if base pair selection matches (different pairs = different parameters)
   - Use `python3 scripts/compare_json.py steps <PDB_ID> --verbose` for detailed comparison

### Reference Frames Issues

1. **ref_frames.dat doesn't match legacy**:
   - Use `--legacy-inp` option when generating modern ref_frames:
     ```bash
     ./build/find_pair_app --legacy-inp org/1H4S.inp data/pdb/1H4S.pdb output.inp
     ```
   - This ensures correct strand ordering (strand 2 first, strand 1 second)
   - Compare using: `python3 scripts/compare_ref_frames.py org/ref_frames.dat ref_frames_modern.dat --verbose`

2. **Y/Z axes differ by 180° rotation**:
   - This is expected if legacy ordering is not used
   - Legacy uses variable ordering based on scan order
   - Modern uses consistent ordering (smaller-residue-first)
   - Use `--legacy-inp` option for exact matching

3. **Origin and X-axis match but Y/Z don't**:
   - Usually indicates strand ordering issue
   - Verify legacy .inp file is correct
   - Check that `--legacy-inp` option is being used

### Performance Issues

1. **Use test sets**: Compare smaller test sets for faster iteration
2. **Use parallel processing**: Default uses CPU cores - 1 (leaves one core free)
3. **Profile slow tests**: Use profiling tools to identify bottlenecks

---

## Best Practices

1. **Use `compare_json.py`** for all JSON comparisons
2. **Use `--legacy-inp` option** when comparing ref_frames.dat with legacy
3. **Run tests regularly** during development
4. **Compare with legacy** before committing changes
5. **Use test sets** for quick validation
6. **Document test failures** and fixes
7. **Keep test data up to date**
8. **Use verbose output** when debugging
9. **Save reports** for important comparisons
10. **Verify residue matching** if frame origins or dorg don't match

---

## Related Documentation

- [COMPARISON_WORKFLOW.md](COMPARISON_WORKFLOW.md) ⭐ **NEW** - Frame-first comparison workflow
- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - Detailed JSON record types and comparison checks
- [LEGACY_TEST_TOOLS.md](LEGACY_TEST_TOOLS.md) - Legacy test tools documentation
- [ATOM_INDEX_CONVERSION_FIX.md](ATOM_INDEX_CONVERSION_FIX.md) ⭐ **NEW** - Atom index conversion for legacy input files
- [REF_FRAMES_COMPARISON_FIX.md](REF_FRAMES_COMPARISON_FIX.md) ⭐ **NEW** - Ref frames comparison and strand ordering fix
- [RESIDUE_INDEXING_COMPLETE.md](RESIDUE_INDEXING_COMPLETE.md) - Residue indexing solution details
- [scripts/README.md](../scripts/README.md) - Scripts reference guide
- [x3dna_json_compare/README.md](../x3dna_json_compare/README.md) - Comparison library API
- [BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md) - Build instructions
- [QUICK_START.md](QUICK_START.md) - Quick start guide

---

*This guide consolidates all testing and comparison workflows. For questions or issues, refer to the related documentation or check the script help: `python3 scripts/compare_json.py --help`*

