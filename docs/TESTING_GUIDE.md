# Testing Guide

**Last Updated**: 2025-01-XX  
**Purpose**: Centralized guide for all testing and JSON comparison workflows

---

## Quick Start: JSON Comparison

The **primary tool** for JSON comparisons is `scripts/compare_json.py`. This is the unified, comprehensive comparison tool that should be used for all JSON comparison needs.

### Basic Usage

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

See [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) for detailed information about each record type and what is checked.

### Comparison Library

The comparison functionality is provided by the `x3dna_json_compare` Python package, which includes:

- **`JsonComparator`** - Main comparison engine with caching
- **`PdbFileReader`** - Efficient PDB file reader
- **`JsonValidator`** - JSON validation utilities
- **`ComparisonCache`** - Result caching for performance

See [x3dna_json_compare/README.md](../x3dna_json_compare/README.md) for API documentation.

---

## Test Types

### 1. C++ Unit Tests

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

### 2. C++ Integration Tests

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

### 3. JSON Regression Tests

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

### 4. Legacy Test Tools

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
3. **Cache issues**: Use `--no-cache` flag to bypass caching
4. **Threading issues**: Use `--threads 1` for single-threaded execution

### Performance Issues

1. **Use caching**: Comparison results are cached by default
2. **Use test sets**: Compare smaller test sets for faster iteration
3. **Use parallel processing**: Default uses all CPU cores
4. **Profile slow tests**: Use profiling tools to identify bottlenecks

---

## Best Practices

1. **Use `compare_json.py`** for all JSON comparisons
2. **Run tests regularly** during development
3. **Compare with legacy** before committing changes
4. **Use test sets** for quick validation
5. **Document test failures** and fixes
6. **Keep test data up to date**
7. **Use verbose output** when debugging
8. **Save reports** for important comparisons

---

## Related Documentation

- [COMPARISON_WORKFLOW.md](COMPARISON_WORKFLOW.md) ⭐ **NEW** - Frame-first comparison workflow
- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - Detailed JSON record types and comparison checks
- [LEGACY_TEST_TOOLS.md](LEGACY_TEST_TOOLS.md) - Legacy test tools documentation
- [scripts/README.md](../scripts/README.md) - Scripts reference guide
- [x3dna_json_compare/README.md](../x3dna_json_compare/README.md) - Comparison library API
- [BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md) - Build instructions
- [QUICK_START.md](QUICK_START.md) - Quick start guide

---

*This guide consolidates all testing and comparison workflows. For questions or issues, refer to the related documentation or check the script help: `python3 scripts/compare_json.py --help`*

