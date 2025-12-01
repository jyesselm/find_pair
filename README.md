# X3DNA Modernized - Modern C++ Library

A modern C++ rewrite of the X3DNA v2.4 codebase with strong object-oriented design, comprehensive testing, and full JSON serialization support.

## Status

✅ **100% Match Achieved on Test Set!** - All 10 PDBs in test set show perfect matches across all comparison types.

**Current Status**: Modern code matches legacy code output exactly. See [docs/100_PERCENT_MATCH_PLAN.md](docs/100_PERCENT_MATCH_PLAN.md) for detailed status.

---

## Quick Start

### 1. Build Both Legacy and Modern Code

```bash
# Build modern code (Release mode)
make release

# Build legacy code (Release mode)
make org-release
```

### 2. Generate JSON Output Files

```bash
# Generate legacy JSON (reference output)
python3 scripts/rebuild_json.py regenerate 1H4S --legacy-only

# Generate modern JSON (to be compared against legacy)
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only

# Or generate both at once
python3 scripts/rebuild_json.py regenerate 1H4S
```

### 3. Compare Outputs

```bash
# Compare a specific PDB
python3 scripts/compare_json.py compare 1H4S

# Compare with verbose output
python3 scripts/compare_json.py compare 1H4S --verbose
```

---

## Generating Output Files

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

#### Method 1: Using Python Script (Recommended)

```bash
# Single PDB
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only

# Multiple PDBs
python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA --modern-only

# Using test set
python3 scripts/rebuild_json.py regenerate --test-set 100 --modern-only

# All available PDBs
python3 scripts/rebuild_json.py regenerate --modern-only
```

#### Method 2: Direct Executable

```bash
# From project root
./build/generate_modern_json data/pdb/1H4S.pdb data/json/
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
| `find_pair_app` | Find base pairs (modern find_pair) | `./build/find_pair_app [--fix-indices] [--legacy-inp=FILE] data/pdb/1H4S.pdb output.inp` |
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

### Important: The `--fix-indices` Option

**When comparing with legacy output, always use `--fix-indices`:**

```bash
./build/find_pair_app --fix-indices data/pdb/1H4S.pdb output.inp
```

**Why?** The modern and legacy code assign different indices to residues (off by 1). The `--fix-indices` option:
- Reads the legacy JSON file (`data/json_legacy/base_frame_calc/<PDB>.json`)
- Maps residues by their PDB properties (chain, seq_num, insertion code)
- Assigns matching legacy indices to modern residues

**Without `--fix-indices`**: Base pairs will have different residue indices and appear non-matching even when they're the same physical pairs.

### The `--legacy-inp` Option for Exact Frame Matching

When generating `ref_frames_modern.dat`, the Y/Z axis orientation depends on residue ordering within pairs. Legacy's `find_pair` uses a variable ordering based on scan order, while modern consistently uses smaller-residue-first.

**To match legacy's frame orientation exactly:**

```bash
# First generate legacy .inp and ref_frames.dat
cd org && ./build/bin/find_pair_original ../data/pdb/1H4S.pdb && cd ..

# Then generate modern with legacy ordering
./build/find_pair_app --fix-indices --legacy-inp=1H4S.inp data/pdb/1H4S.pdb output.inp
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

### Test Types

1. **C++ Unit Tests** - Test individual classes and functions (`tests/unit/`)
2. **C++ Integration Tests** - Test component interactions (`tests/integration/`)
3. **JSON Regression Tests** - Compare modern JSON with legacy JSON (primary test method)

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
├── include/x3dna/          # Public headers (modern code)
├── src/x3dna/              # Implementation files (modern code)
├── org/                    # Legacy/original X3DNA code (reference)
│   ├── src/                # Legacy C source files
│   ├── include/            # Legacy headers
│   └── build/              # Legacy build directory
├── tools/                  # Standalone tools and utilities
├── tests/                   # Test suite
│   ├── unit/               # Unit tests
│   ├── integration/        # Integration tests
│   └── regression/         # Regression tests
├── scripts/                 # Python scripts for JSON comparison/rebuilding
├── data/                    # Data files
│   ├── pdb/                # PDB structure files
│   ├── templates/          # Template files for base frame calculation
│   ├── json/               # Modern JSON output (segmented by record type)
│   ├── json_legacy/        # Legacy JSON output (reference for comparison)
│   └── test_sets/          # Test set definitions
└── docs/                    # Documentation
    ├── legacy/             # Legacy code documentation
    └── modernization/      # Stage-by-stage modernization plan
```

## Directory Structure

### Resources Directory (Required Runtime Resources)
- **`resources/templates/`** - Template files (Atomic_*.pdb) for base frame calculation
- **`resources/config/`** - Configuration files (special_residues.json)
- **`resources/test_sets/`** - Test set definitions (test_set_10.json, test_set_100.json, etc.)

### Data Directory (Test Inputs and Generated Outputs)
- **`data/pdb/`** - PDB structure files (test inputs)
- **`data/json/`** - Modern JSON output files (segmented by record type, generated)
- **`data/json_legacy/`** - Legacy JSON output files (reference for comparison, generated)

See [docs/DATA_STRUCTURE.md](docs/DATA_STRUCTURE.md) for complete data directory documentation.

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

