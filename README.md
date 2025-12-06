# X3DNA Modernized - Modern C++ Library

A modern C++ rewrite of the X3DNA v2.4 codebase with strong object-oriented design, comprehensive testing, and full JSON serialization support.

## Status

✅ **100% Match Achieved on Test Set!** - All 10 PDBs in test set show perfect matches across all comparison types.

**Current Status**: Residue index matching complete. Running full validation. See [docs/VALIDATION_PROGRESS.md](docs/VALIDATION_PROGRESS.md) for detailed status and next steps.

---

## Quick Start: Validation CLI

The primary tool for all validation is the `fp2-validate` CLI:

```bash
# Install
pip install -e .

# Run validation
fp2-validate                      # All stages, all 3602 fast PDBs
fp2-validate --pdb 1EHZ -v        # Single PDB, verbose output
fp2-validate --test-set 100       # 100-PDB test set
fp2-validate frames --max 50      # Frames only, first 50 PDBs
fp2-validate --stop-on-first      # Debug mode - stop at first failure

# Show info
fp2-validate info                 # Environment & available test sets
fp2-validate list-pdbs            # List all available PDBs
```

See [docs/TESTING_GUIDE.md](docs/TESTING_GUIDE.md) for complete CLI reference.

---

## Getting Started

### Initial Setup

1. **Install Dependencies**
   ```bash
   # Install Python dependencies (including pytest)
   pip install -e ".[dev]"
   ```

2. **Build Both Legacy and Modern Code**
   ```bash
   # Build modern code (Release mode)
   make release

   # Build legacy code (Release mode)
   make org-release
   ```

3. **Download Test PDBs (Optional)**
   ```bash
   # Download PDBs from test set
   python tools/download_pdbs.py --test-set 100

   # Or download specific PDBs
   python tools/download_pdbs.py 1H4S 2BNA 3DNA
   ```

### Running Tests by Stage

Validation is done in **10 stages**, each testing a specific aspect of the code. Work through stages sequentially, validating each before moving to the next.

#### Stage 1: Atom Parsing ✅ **COMPLETE**

**Status**: ✅ Validated - 3602/3602 PDBs (100% pass rate)

**What it tests**: PDB file parsing and atom extraction

**Run tests**:
```bash
# Using pytest (recommended)
pytest tests_python/integration/test_atoms_batch.py

# Or run as standalone script
python tests_python/integration/test_atoms_batch.py

# Or compare specific PDB
python3 scripts/compare_json.py atoms 1H4S
```

**Generate JSON** (if needed):
```bash
# Generate atoms JSON for a PDB
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=atoms
```

---

#### Stage 2: Residue Indices ✅ **COMPLETE**

**Status**: ✅ Validated

**What it tests**: Mapping of residues to atom ranges (seidx)

**Run tests**:
```bash
# Using pytest
pytest tests_python/integration/test_residue_indices_batch.py

# Or run as standalone script
python tests_python/integration/test_residue_indices_batch.py
```

**Generate JSON**:
```bash
# Generate residue_indices JSON
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=residue_indices
```

---

#### Stage 3: LS Fitting ✅ **COMPLETE**

**Status**: ✅ Validated

**What it tests**: Least-squares fitting for frame calculation

**Run tests**:
```bash
# Using pytest
pytest tests_python/integration/test_ls_fitting.py

# Or run as standalone script
python tests_python/integration/test_ls_fitting.py
```

**Generate JSON**:
```bash
# Generate ls_fitting JSON
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=ls_fitting
```

---

#### Stage 4: Reference Frames ✅ **COMPLETE**

**Status**: ✅ Validated

**What it tests**: Reference frame calculations (base_frame_calc, frame_calc)

**Run tests**:
```bash
# Using pytest
pytest tests_python/integration/test_frames_batch.py

# Or run as standalone script
python tests_python/integration/test_frames_batch.py

# Or compare specific PDB
python3 scripts/compare_json.py frames 1H4S
```

**Generate JSON**:
```bash
# Generate frames JSON (base_frame_calc + frame_calc, but not ls_fitting)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=frames
```

---

#### Stage 5: Distance Checks ⏳ **TODO**

**Status**: ⏳ Not yet validated

**What it tests**: Geometric distance and angle calculations (dorg, dNN, plane_angle, d_v, overlap_area)

**Run tests** (when implemented):
```bash
# Compare distance checks
python3 scripts/compare_json.py distance-checks 1H4S
```

**Generate JSON**:
```bash
# Generate distance_checks JSON (when stage is implemented)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=distances
```

---

#### Stage 6: Hydrogen Bonds ⏳ **TODO**

**Status**: ⏳ Not yet validated

**What it tests**: H-bond detection and conflict resolution

**Run tests** (when implemented):
```bash
# Compare H-bonds
python3 scripts/compare_json.py hbonds 1H4S
```

**Generate JSON**:
```bash
# Generate hbond_list JSON (when stage is implemented)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=hbonds
```

---

#### Stage 7: Pair Validation ⏳ **TODO**

**Status**: ⏳ Not yet validated

**What it tests**: Validation logic (which pairs pass/fail geometric checks)

**Run tests** (when implemented):
```bash
# Compare pair validation
python3 scripts/compare_json.py validation 1H4S
```

**Generate JSON**:
```bash
# Generate pair_validation JSON (when stage is implemented)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=validation
```

---

#### Stage 8: Base Pair Selection ⏳ **TODO** ⭐ **CRITICAL**

**Status**: ⏳ Not yet validated

**What it tests**: Final selected base pairs (THE PRIMARY OUTPUT)

**Run tests** (when implemented):
```bash
# Compare base pairs
python3 scripts/compare_json.py pairs 1H4S
```

**Generate JSON**:
```bash
# Generate base_pair and find_bestpair_selection JSON (when stage is implemented)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=selection
```

---

#### Stage 9: Step Parameters ⏳ **TODO**

**Status**: ⏳ Not yet validated

**What it tests**: Step parameter calculations (Shift, Slide, Rise, Tilt, Roll, Twist)

**Run tests**:
```bash
# Compare step parameters
python3 scripts/compare_json.py steps 1H4S
```

**Generate JSON**:
```bash
# Generate bpstep_params JSON
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=steps
```

---

#### Stage 10: Helical Parameters ⏳ **TODO**

**Status**: ⏳ Not yet validated

**What it tests**: Helical parameter calculations

**Run tests** (when implemented):
```bash
# Compare helical parameters
python3 scripts/compare_json.py helical 1H4S
```

**Generate JSON**:
```bash
# Generate helical_params JSON (when stage is implemented)
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=helical
```

---

### Quick Test Workflow

**For a single PDB**:
```bash
# 1. Generate modern JSON for current stage
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=<stage_name>

# 2. Compare with legacy
python3 scripts/compare_json.py <stage_command> 1H4S --verbose

# 3. Fix any issues and repeat
```

**For batch testing**:
```bash
# Run pytest tests (recommended)
pytest tests_python/integration/test_<stage>_batch.py

# Or run standalone test script
python tests_python/integration/test_<stage>_batch.py
```

**For all available PDBs**:
```bash
# Compare all PDBs for a specific stage
python3 scripts/compare_json.py <stage_command> --test-set 100
```

### Stage Status Summary

| Stage | Name | Status | Test Command |
|-------|------|--------|--------------|
| 1 | Atoms | ✅ Complete | `pytest tests_python/integration/test_atoms_batch.py` |
| 2 | Residue Indices | ✅ Complete | `pytest tests_python/integration/test_residue_indices_batch.py` |
| 3 | LS Fitting | ✅ Complete | `pytest tests_python/integration/test_ls_fitting.py` |
| 4 | Frames | ✅ Complete | `pytest tests_python/integration/test_frames_batch.py` |
| 5 | Distance Checks | ⏳ TODO | `python3 scripts/compare_json.py distance-checks <PDB>` |
| 6 | H-Bonds | ⏳ TODO | `python3 scripts/compare_json.py hbonds <PDB>` |
| 7 | Pair Validation | ⏳ TODO | `python3 scripts/compare_json.py validation <PDB>` |
| 8 | Base Pair Selection | ⏳ TODO ⭐ | `python3 scripts/compare_json.py pairs <PDB>` |
| 9 | Step Parameters | ⏳ TODO | `python3 scripts/compare_json.py steps <PDB>` |
| 10 | Helical Parameters | ⏳ TODO | `python3 scripts/compare_json.py helical <PDB>` |

See [docs/STAGED_VALIDATION_PLAN.md](docs/STAGED_VALIDATION_PLAN.md) for detailed stage information.

---

## Generating JSON Output Files

**Note**: For staged validation, only generate JSON for the stage you're currently working on. See "Running Tests by Stage" above for stage-specific commands.

### Legacy JSON Output (Reference)

The `org/` directory contains the original X3DNA code. Legacy output serves as the **reference** that modern code must match.

**Location:** `data/json_legacy/`

#### Method 1: Using Python Script (Recommended)

```bash
# Single PDB
python3 scripts/rebuild_json.py regenerate 1H4S --legacy-only

# Multiple PDB
python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA --legacy-only

# Using test set
python3 scripts/rebuild_json.py regenerate --test-set 100 --legacy-only

# All available PDBs
python3 scripts/rebuild_json.py regenerate --legacy-only
```

#### Method 2: Direct Executable

```bash
# From project root
cd org
./build/bin/find_pair_analyze ../data/pdb/1H4S.pdb
```

**Note:** Run from the `org/` directory as it uses relative paths for templates.

---

### Modern JSON Output

Modern code output is compared against legacy output to verify correctness.

**Location:** `data/json/`

**Important**: Use `--stage=<stage_name>` to generate only the JSON types needed for the current validation stage.

#### Method 1: Using Direct Executable (Recommended for Staged Validation)

```bash
# Generate specific stage JSON
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=atoms
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=residue_indices
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=ls_fitting
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=frames
./build/generate_modern_json data/pdb/1H4S.pdb data/json/ --stage=all  # All stages

# Available stages: atoms, residue_indices, ls_fitting, frames, distances, hbonds, 
#                   validation, selection, steps, helical, all
```

#### Method 2: Using Python Script

```bash
# Single PDB (all stages)
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only

# Multiple PDBs
python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA --modern-only

# Using test set
python3 scripts/rebuild_json.py regenerate --test-set 100 --modern-only
```

#### Method 3: Using find_pair_app

```bash
# Generate .inp file and ref_frames_modern.dat
./build/find_pair_app data/pdb/1H4S.pdb output.inp
```

This creates:
- `output.inp` - Input file for analyze_app
- `ref_frames_modern.dat` - Reference frames in legacy-compatible format

---

### Regenerate Both Legacy and Modern

```bash
# Single PDB
python3 scripts/rebuild_json.py regenerate 1H4S

# Multiple PDBs
python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA

# Using test set
python3 scripts/rebuild_json.py regenerate --test-set 100

# All available PDBs
python3 scripts/rebuild_json.py regenerate
```

---

## Comparing Outputs

The primary comparison tool is `scripts/compare_json.py`. It compares modern JSON against legacy JSON to identify differences.

**For complete testing documentation, see [docs/TESTING_GUIDE.md](docs/TESTING_GUIDE.md)**

**For staged validation, use stage-specific comparison commands** (see "Running Tests by Stage" above).

### Basic Comparison

```bash
# Compare all available PDBs
python3 scripts/compare_json.py compare

# Compare specific PDB(s)
python3 scripts/compare_json.py compare 1H4S
python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

# Compare using test sets (10, 50, 100, 500, or 1000 PDBs)
python3 scripts/compare_json.py compare --test-set 100
```

### Comparison Options

```bash
# Verbose output with detailed differences
python3 scripts/compare_json.py compare 1H4S --verbose

# Save report to file
python3 scripts/compare_json.py compare --output report.md

# Show only files with differences
python3 scripts/compare_json.py compare --diff-only

# Auto-regenerate missing JSON files
python3 scripts/compare_json.py compare --regenerate

# Custom thread count
python3 scripts/compare_json.py compare --threads 4
```

### Compare Specific Record Types

```bash
python3 scripts/compare_json.py atoms 1H4S       # Atoms only
python3 scripts/compare_json.py frames 1H4S      # Reference frames only
python3 scripts/compare_json.py ring-atoms 1H4S  # Ring atoms only
python3 scripts/compare_json.py steps 1H4S       # Step parameters only
```

### What Gets Compared

| Record Type | Description |
|-------------|-------------|
| `pdb_atoms` | Atom records from PDB parsing |
| `base_frame_calc` | Base frame calculation results |
| `frame_calc` / `ref_frame` | Reference frame calculations |
| `base_pair` | Base pair records |
| `pair_validation` | Pair validation results |
| `distance_checks` | Distance and geometric checks |
| `hbond_list` | Hydrogen bond lists |
| `find_bestpair_selection` | Selected base pairs |
| `bpstep_params` | Step parameters (Shift, Slide, Rise, Tilt, Roll, Twist) |
| `helical_params` | Helical parameters |

See [docs/JSON_DATA_TYPES_AND_COMPARISONS.md](docs/JSON_DATA_TYPES_AND_COMPARISONS.md) for detailed information.

---

## Complete Workflow Example

```bash
# 1. Build both codebases
make release
make org-release

# 2. Generate JSON for a new PDB
cp /path/to/new.pdb data/pdb/
python3 scripts/rebuild_json.py regenerate new --legacy-only
python3 scripts/rebuild_json.py regenerate new --modern-only

# 3. Compare outputs
python3 scripts/compare_json.py compare new --verbose

# 4. If differences found, fix modern code and repeat
#    ... edit src/x3dna/*.cpp ...
make release
python3 scripts/rebuild_json.py regenerate new --modern-only
python3 scripts/compare_json.py compare new --verbose

# 5. Once matching, validate against full test set
python3 scripts/compare_json.py compare --test-set 100
```

---

## Executables Reference

### Modern Executables (`./build/`)

| Executable | Purpose | Usage |
|------------|---------|-------|
| `find_pair_app` | Find base pairs (modern find_pair) | `./build/find_pair_app [--legacy-inp=FILE] data/pdb/1H4S.pdb output.inp` |
| `analyze_app` | Calculate step parameters | `./build/analyze_app output.inp` |
| `generate_modern_json` | Generate JSON for comparison | `./build/generate_modern_json data/pdb/1H4S.pdb data/json/` |

### Legacy Executables (`org/build/bin/`)

| Executable | Purpose | Usage |
|------------|---------|-------|
| `find_pair_original` | Original find_pair | `./build/bin/find_pair_original ../data/pdb/1H4S.pdb` |
| `analyze_original` | Original analyze | `./build/bin/analyze_original 1H4S.inp` |
| `find_pair_analyze` | Combined find_pair + JSON output | `./build/bin/find_pair_analyze ../data/pdb/1H4S.pdb` |

### Comparison

| Legacy | Modern | Output |
|--------|--------|--------|
| `find_pair_original` | `find_pair_app` | `.inp` file, `ref_frames.dat` |
| `analyze_original` | `analyze_app` | Step parameters |
| `find_pair_analyze` | `generate_modern_json` | JSON files |

Modern parsing now assigns legacy-style residue and atom indices directly from the PDB scan order, so no extra flags are required when comparing outputs. If you previously relied on `--fix-indices`, simply regenerate with the current tools—the indices will already match the legacy JSON.

### The `--legacy-inp` Option for Exact Frame Matching

When generating `ref_frames_modern.dat`, the Y/Z axis orientation depends on residue ordering within pairs. Legacy's `find_pair` uses a variable ordering based on scan order, while modern consistently uses smaller-residue-first.

**To match legacy's frame orientation exactly:**

```bash
# First generate legacy .inp and ref_frames.dat
cd org && ./build/bin/find_pair_original ../data/pdb/1H4S.pdb && cd ..

# Then generate modern with legacy ordering
./build/find_pair_app --legacy-inp=1H4S.inp data/pdb/1H4S.pdb output.inp
```

This parses the legacy `.inp` file to determine which residue was "first" in each pair, then uses the same ordering for frame calculations. This achieves **100% match** on all axes (origin, x, y, z).

**Without `--legacy-inp`**: Origin and X-axis will match 100%, but Y/Z axes may be flipped (180° rotation around X) for ~35% of pairs.

> **Note**: The base pair local parameters (Shear, Stretch, Stagger, Buckle, Propeller, Opening) calculated by `analyze` may also be affected by the Y/Z frame orientation difference. These parameters are not yet compared in JSON format. If exact matching is needed, use `--legacy-inp` to ensure consistent frame orientation.

---

## Building the Code

### Modern Code

**Using Makefile (recommended):**

```bash
# Build in Release mode (optimized, default)
make release

# Build in Debug mode (with debug symbols)
make debug

# Build in RelWithDebInfo mode (optimized with debug symbols)
make relwithdebinfo

# Clean build artifacts
make clean

# Remove entire build directory
make clean-all

# Rebuild from scratch
make rebuild
```

**Or manually with CMake:**

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

**Build Targets:**
- `make release` - Optimized build (`-O3`, `-DNDEBUG`)
- `make debug` - Debug build (`-g`, `-O0`)
- `make relwithdebinfo` - Optimized with debug symbols (`-O2`, `-g`)
- `make test` - Run all tests
- `make test-verbose` - Run tests with verbose output
- `make format` - Format all C++ source files
- `make format-check` - Check code formatting
- `make compile_commands` - Generate `compile_commands.json` for IDE support

---

### Legacy Code

**Using Makefile (recommended):**

```bash
# Build in Release mode (default, optimized with -O3)
make org-release

# Build in Debug mode (with debug symbols, -g, -O0)
make org-debug

# Build in RelWithDebInfo mode (optimized with debug symbols, -O2, -g)
make org-relwithdebinfo

# Clean build artifacts
make org-clean

# Remove entire build directory
make org-clean-all
```

**Or manually with CMake:**

```bash
cd org
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

**Executables created in `org/build/bin/`:**
- `find_pair_analyze` - Main executable that processes PDB files and generates JSON output
- `find_pair_original` - Original find_pair executable
- `analyze_original` - Original analyze executable

---

## Running Tests

### C++ Tests

```bash
# Run all C++ tests
make test

# Run tests with verbose output
make test-verbose

# Or using ctest directly
cd build
ctest --output-on-failure

# Run specific test category
cd build && ctest -R unit        # Unit tests only
cd build && ctest -R integration # Integration tests only
```

### Python Tests (pytest)

**All Python tests use pytest and are located in `tests_python/`.**

```bash
# Install pytest (if not already installed)
pip install -e ".[dev]"

# Run all Python tests
pytest tests_python/

# Run only integration tests
pytest tests_python/integration/

# Run a specific test file
pytest tests_python/integration/test_atoms_batch.py

# Run with verbose output
pytest tests_python/ -v

# Run only fast tests (exclude slow markers)
pytest tests_python/ -m "not slow"
```

**Test scripts can also be run standalone:**
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

See [tests_python/README.md](tests_python/README.md) for detailed pytest documentation.

### Test Types

1. **C++ Unit Tests** - Test individual classes and functions (`tests/unit/`)
2. **C++ Integration Tests** - Test component interactions (`tests/integration/`)
3. **Python Integration Tests** - JSON generation and comparison tests (`tests_python/integration/`)
4. **JSON Regression Tests** - Compare modern JSON with legacy JSON (primary test method)

The primary testing method is JSON regression testing, which ensures modern code output exactly matches legacy code output.

## Development Workflow

### Standard Development Cycle

1. **Edit modern code** (`src/x3dna/` or `include/x3dna/`)
2. **Build**: `make release`
3. **Generate modern JSON**: `python3 scripts/rebuild_json.py regenerate <PDB> --modern-only`
4. **Compare with legacy**: `python3 scripts/compare_json.py compare <PDB>`
5. **Fix differences** → Repeat until 100% match

### Key Principles

- **Minimize file creation** - Edit existing files rather than creating new ones
- **Focus on matching legacy output** - The goal is 100% JSON match, not perfect architecture
- **Use regression testing** - Primary test method is JSON comparison
- **Legacy code uses 1-based indexing** - Modern code uses 0-based indexing (handled automatically)

See [docs/SIMPLIFIED_DEVELOPMENT.md](docs/SIMPLIFIED_DEVELOPMENT.md) for detailed development guidelines.

---

## JSON Tools Reference

### `compare_json.py` - Compare JSON Files

The main tool for comparing legacy and modern JSON outputs.

```bash
# Basic usage
python3 scripts/compare_json.py compare 1H4S

# Compare multiple PDBs
python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

# Compare all available PDBs
python3 scripts/compare_json.py compare

# Compare test set
python3 scripts/compare_json.py compare --test-set 100

# Options
python3 scripts/compare_json.py compare --verbose      # Detailed differences
python3 scripts/compare_json.py compare --diff-only    # Only show differences
python3 scripts/compare_json.py compare --output report.md  # Save to file
python3 scripts/compare_json.py compare --regenerate   # Auto-regenerate missing files

# Compare specific record types
python3 scripts/compare_json.py atoms 1H4S       # Atoms
python3 scripts/compare_json.py frames 1H4S      # Reference frames
python3 scripts/compare_json.py ring-atoms 1H4S  # Ring atoms
python3 scripts/compare_json.py steps 1H4S       # Step parameters

# Generate test sets
python3 scripts/compare_json.py generate-test-sets
python3 scripts/compare_json.py generate-test-sets --force  # Overwrite existing
```

### `rebuild_json.py` - Regenerate, Validate, and Clean JSON

Unified tool for managing JSON files.

```bash
# Regenerate JSON files
python3 scripts/rebuild_json.py regenerate 1H4S              # Single PDB (both)
python3 scripts/rebuild_json.py regenerate 1H4S --legacy-only   # Legacy only
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only   # Modern only
python3 scripts/rebuild_json.py regenerate --test-set 100    # Test set

# Validate JSON files
python3 scripts/rebuild_json.py validate
python3 scripts/rebuild_json.py validate --legacy-dir data/json_legacy --modern-dir data/json

# Clean invalid/empty JSON files
python3 scripts/rebuild_json.py clean           # Dry run
python3 scripts/rebuild_json.py clean --execute # Actually remove files
```

For more detailed documentation, see `scripts/README.md` and [docs/TESTING_GUIDE.md](docs/TESTING_GUIDE.md).

## Project Structure

```
find_pair_2/
├── apps/                   # Main applications
│   ├── find_pair_app.cpp   # Modern find_pair implementation
│   └── analyze_app.cpp     # Modern analyze implementation
├── include/x3dna/          # Public headers (modern code)
├── src/x3dna/              # Implementation files (modern code)
├── org/                    # Legacy/original X3DNA code (reference)
│   ├── src/                # Legacy C source files
│   ├── include/            # Legacy headers
│   └── build/              # Legacy build directory
├── tools/                  # Standalone utility tools
│   └── generate_modern_json.cpp  # Main JSON generation tool
├── tests/                  # Test suite
│   ├── unit/              # Unit tests
│   └── integration/       # Integration tests
├── scripts/                # Essential Python scripts
│   ├── compare_json.py    # Main JSON comparison tool
│   ├── rebuild_json.py    # JSON regeneration tool
│   ├── download_pdbs.py   # PDB download utility
│   └── cluster/           # Cluster computing scripts
├── resources/              # Runtime resources
│   ├── templates/         # Base frame templates (Atomic_*.pdb)
│   ├── config/            # Configuration files
│   └── test_sets/         # Test set definitions
├── data/                   # Generated data and inputs
│   ├── pdb/               # PDB structure files (inputs)
│   ├── json/              # Modern JSON output
│   └── json_legacy/       # Legacy JSON output (reference)
├── docs/                   # Documentation
│   ├── README.md          # Documentation index
│   ├── QUICK_START.md     # Getting started guide
│   ├── TESTING_GUIDE.md   # Testing and comparison guide
│   ├── CODE_FLOW.md       # Architecture and code flow
│   ├── legacy/            # Legacy code documentation
│   ├── modernization/     # Modernization plan and stages
│   └── archive/           # Archived debugging/completed docs
└── x3dna_json_compare/     # Python comparison library
    └── *.py               # Comparison modules
```

## Directory Organization

### Essential Directories

**Source Code:**
- `apps/` - Main application entry points
- `include/x3dna/`, `src/x3dna/` - Modern C++ library
- `org/` - Legacy reference code (do not edit except for JSON output)

**Resources (Runtime):**
- `resources/templates/` - Base frame templates (Atomic_*.pdb)
- `resources/config/` - Configuration files
- `resources/test_sets/` - Test set definitions

**Data (Generated):**
- `data/pdb/` - Input PDB files
- `data/json/` - Modern JSON output
- `data/json_legacy/` - Legacy JSON output (reference)

**Scripts & Tools:**
- `scripts/` - Essential Python tools (compare_json.py, rebuild_json.py, test_utils.py)
- `tools/` - Utility programs (C++ tools and Python tools like download_pdbs.py, find_slow_pdbs.py)
- `x3dna_json_compare/` - Python comparison library

**Testing:**
- `tests/unit/` - C++ unit tests
- `tests/integration/` - C++ integration tests
- `tests_python/` - Python tests (pytest)
  - `tests_python/integration/` - Integration tests for JSON generation/comparison
  - `tests_python/unit/` - Unit tests for utility functions
  - `tests_python/conftest.py` - Shared pytest fixtures

**Documentation:**
- `docs/` - All documentation (see [docs/README.md](docs/README.md))

### Archive Directories

Completed debugging work and temporary files are archived in:
- `docs/archive/` - Archived documentation
- `scripts/archive/` - Archived scripts
- `tools/archive/` - Archived debugging tools
- `data/archive/` - Archived result files

## Documentation

### Primary Documents

- **[100_PERCENT_MATCH_PLAN.md](docs/100_PERCENT_MATCH_PLAN.md)** ⭐ - Current status and progress toward 100% match
- **[TESTING_GUIDE.md](docs/TESTING_GUIDE.md)** - Complete testing and JSON comparison guide
- **[BUILD_INSTRUCTIONS.md](docs/BUILD_INSTRUCTIONS.md)** - Detailed build instructions
- **[QUICK_START.md](docs/QUICK_START.md)** - Quick start guide for new developers
- **[SIMPLIFIED_DEVELOPMENT.md](docs/SIMPLIFIED_DEVELOPMENT.md)** - Development workflow and best practices

### Reference Documents

- **[JSON_DATA_TYPES_AND_COMPARISONS.md](docs/JSON_DATA_TYPES_AND_COMPARISONS.md)** - JSON structure and comparison methodology
- **[ALGORITHM_CRITICAL_GUIDE.md](docs/ALGORITHM_CRITICAL_GUIDE.md)** - Critical algorithm details and edge cases
- **[DEBUGGING_STRATEGIES.md](docs/DEBUGGING_STRATEGIES.md)** - Debugging strategies and tools
- **[MODERNIZATION_PLAN.md](docs/MODERNIZATION_PLAN.md)** - Overall modernization strategy
- **[OOP_CLASS_HIERARCHY.md](docs/OOP_CLASS_HIERARCHY.md)** - Class hierarchy design

### Legacy Code Documentation

- **[legacy/](docs/legacy/)** - Complete legacy code documentation
  - [00_INDEX.md](docs/legacy/00_INDEX.md) - Navigation guide
  - [01_ARCHITECTURE.md](docs/legacy/01_ARCHITECTURE.md) - Program structure
  - [02_DATA_STRUCTURES.md](docs/legacy/02_DATA_STRUCTURES.md) - Data organization
  - [03_CORE_FUNCTIONS.md](docs/legacy/03_CORE_FUNCTIONS.md) - Critical functions
  - [04_ALGORITHMS.md](docs/legacy/04_ALGORITHMS.md) - Algorithm details
  - [10_JSON_STRUCTURE.md](docs/legacy/10_JSON_STRUCTURE.md) - Legacy JSON format

See [docs/README.md](docs/README.md) for complete documentation index.

## Dependencies

- **C++17 or later** - Required for modern C++ features
- **CMake 3.15+** - Build system
- **Python 3.6+** - For JSON comparison and rebuilding scripts
- **nlohmann/json** - Header-only JSON library (fetched automatically via CMake)
- **Google Test** - Testing framework (fetched automatically via CMake)

### Optional Dependencies

- **Ninja** - Faster build system (auto-detected if available)
- **clang-format** - Code formatting (optional, for `make format`)

## Important Notes

### Indexing Differences

- **Legacy code**: Uses **1-based indexing** (arrays start at index 1)
- **Modern code**: Uses **0-based indexing** (standard C++ convention)
- **Comparison tools**: Automatically handle index conversion

### Code Organization Rules

- **Legacy code** (`org/`): Do not edit except to add JSON output and debug statements
- **Modern code** (`src/x3dna/`, `include/x3dna/`): Main development area
- **Functions**: Keep under 50 lines when possible
- **Indentation**: Maximum 3 levels of indentation per function
- **File creation**: Minimize new files - build on existing files when possible

### JSON File Locations

- **Legacy JSON**: `data/json_legacy/` - Reference output from original code
- **Modern JSON**: `data/json/` - Output from modernized code
- **Comparison**: Modern code must exactly match legacy code for all calculations

## Troubleshooting

### Build Issues

```bash
# Clean and rebuild
make clean-all
make release

# Check build system info
make info
```

### Test Failures

1. Ensure both legacy and modern JSON files exist
2. Regenerate JSON files if needed: `python3 scripts/rebuild_json.py regenerate <PDB>`
3. Check tolerance values in comparison tools
4. Use verbose output: `python3 scripts/compare_json.py compare <PDB> --verbose`

### JSON Comparison Issues

1. **Missing files**: Regenerate with `python3 scripts/rebuild_json.py regenerate <PDB>`
2. **Invalid JSON**: Validate with `python3 scripts/rebuild_json.py validate`
3. **Cache issues**: Use `--no-cache` flag in comparison tools
4. **Threading issues**: Use `--threads 1` for single-threaded execution

## License

[To be determined]

## Contributing

See `docs/modernization/` for the detailed implementation plan and `docs/SIMPLIFIED_DEVELOPMENT.md` for development guidelines.

