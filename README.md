# X3DNA Modernized - Modern C++ Library

A modern C++ rewrite of the X3DNA v2.4 codebase with strong object-oriented design, comprehensive testing, and full JSON serialization support.

## Status

✅ **100% Match Achieved on Test Set!** - All 10 PDBs in test set show perfect matches across all comparison types.

**Current Status**: Modern code matches legacy code output exactly. See [docs/100_PERCENT_MATCH_PLAN.md](docs/100_PERCENT_MATCH_PLAN.md) for detailed status.

## Testing & JSON Comparison

**For complete testing documentation, see [docs/TESTING_GUIDE.md](docs/TESTING_GUIDE.md)**

### Quick Comparison

```bash
# Compare all available PDBs
python3 scripts/compare_json.py compare

# Compare specific PDB(s)
python3 scripts/compare_json.py compare 1H4S
python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

# Compare test set (10, 50, 100, 500, or 1000 PDBs)
python3 scripts/compare_json.py compare --test-set 100

# Compare specific record types
python3 scripts/compare_json.py atoms 1H4S      # Atoms only
python3 scripts/compare_json.py frames 1H4S       # Frames only
python3 scripts/compare_json.py ring-atoms 1H4S  # Ring atoms only

# Verbose output with detailed differences
python3 scripts/compare_json.py compare 1H4S --verbose

# Save report to file
python3 scripts/compare_json.py compare --output report.md

# Show only files with differences
python3 scripts/compare_json.py compare --diff-only
```

### Test Types

1. **C++ Unit Tests** - Test individual classes and functions (`tests/unit/`)
2. **C++ Integration Tests** - Test component interactions (`tests/integration/`)
3. **JSON Regression Tests** - Compare modern JSON with legacy JSON (primary test method)

The primary testing method is JSON regression testing, which ensures modern code output exactly matches legacy code output.

## Building and Running Legacy Code

The `org/` directory contains the original X3DNA code with JSON output support. This is the reference implementation that modern code must match exactly.

### Building Legacy Code

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

### Running Legacy Code

To process a PDB file and generate legacy JSON output:

```bash
# From the org/ directory
cd org
./build/bin/find_pair_analyze ../data/pdb/1H4S.pdb
```

The JSON output will be written to `data/json_legacy/` in segmented format (one file per record type).

**Note:** The executable should be run from the `org/` directory (or with the PDB path relative to the project root), as it may use relative paths for templates or other resources.

### Generating Legacy JSON

```bash
# Generate legacy JSON for a single PDB
cd org
./build/bin/find_pair_analyze ../data/pdb/1H4S.pdb

# Or use the rebuild script
python3 scripts/rebuild_json.py regenerate 1H4S --legacy-only

# Generate for multiple PDBs
python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA --legacy-only
```

**Important:** Legacy code uses 1-based indexing, while modern code uses 0-based indexing. This is handled automatically in the comparison tools.

## Quick Start

### Building Modern Code

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

### Running Tests

```bash
# Run all tests
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

### Generating Modern JSON

```bash
# Generate JSON for a single PDB
./build/generate_modern_json data/pdb/1H4S.pdb data/json/

# Or use the rebuild script
python3 scripts/rebuild_json.py regenerate 1H4S --modern-only
```

## Development Workflow

### Standard Development Cycle

1. **Edit modern code** (`src/x3dna/` or `include/x3dna/`)
2. **Build**: `make release`
3. **Generate modern JSON**: `./build/generate_modern_json data/pdb/<PDB>.pdb data/json/`
4. **Compare with legacy**: `python3 scripts/compare_json.py compare <PDB>`
5. **Fix differences** → Repeat until 100% match

### Key Principles

- **Minimize file creation** - Edit existing files rather than creating new ones
- **Focus on matching legacy output** - The goal is 100% JSON match, not perfect architecture
- **Use regression testing** - Primary test method is JSON comparison
- **Legacy code uses 1-based indexing** - Modern code uses 0-based indexing (handled automatically)

See [docs/SIMPLIFIED_DEVELOPMENT.md](docs/SIMPLIFIED_DEVELOPMENT.md) for detailed development guidelines.

## JSON Comparison and Rebuilding Tools

The project includes two main Python scripts for working with JSON files:

### `compare_json.py` - Compare JSON Files

Compare legacy and modern JSON files to identify differences.

**Basic Usage:**
```bash
# Compare all available PDBs
python3 scripts/compare_json.py compare

# Compare specific PDB file(s)
python3 scripts/compare_json.py compare 1H4S
python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

# Compare specific record types
python3 scripts/compare_json.py atoms 1H4S      # Atoms only
python3 scripts/compare_json.py frames 1H4S     # Frames only
python3 scripts/compare_json.py ring-atoms 1H4S # Ring atoms only

# Compare test set
python3 scripts/compare_json.py compare --test-set 100
```

**Options:**
```bash
# Verbose output with detailed differences
python3 scripts/compare_json.py compare --verbose

# Save report to file
python3 scripts/compare_json.py compare --output report.md

# Show only files with differences
python3 scripts/compare_json.py compare --diff-only

# Regenerate missing JSON files automatically
python3 scripts/compare_json.py compare --regenerate

# Custom thread count
python3 scripts/compare_json.py compare --threads 4
```

**Test Sets:**
```bash
# Generate test sets of different sizes (10, 50, 100, 500, 1000)
python3 scripts/compare_json.py generate-test-sets

# Regenerate test sets (overwrite existing)
python3 scripts/compare_json.py generate-test-sets --force
```

### `rebuild_json.py` - Rebuild JSON Files

Regenerate, validate, and clean JSON files (both legacy and modern).

**Regenerate JSON Files:**
```bash
# Regenerate all JSON files (legacy and modern)
python3 scripts/rebuild_json.py regenerate

# Regenerate only legacy JSON files
python3 scripts/rebuild_json.py regenerate --legacy-only

# Regenerate only modern JSON files
python3 scripts/rebuild_json.py regenerate --modern-only

# Regenerate specific PDB file(s)
python3 scripts/rebuild_json.py regenerate 1H4S
python3 scripts/rebuild_json.py regenerate 1H4S 2BNA 3DNA

# Use test set
python3 scripts/rebuild_json.py regenerate --test-set 100
```

**Validate JSON Files:**
```bash
# Validate all existing JSON files
python3 scripts/rebuild_json.py validate

# Validate specific directories
python3 scripts/rebuild_json.py validate --legacy-dir data/json_legacy --modern-dir data/json
```

**Clean Invalid JSON Files:**
```bash
# Dry run (shows what would be removed)
python3 scripts/rebuild_json.py clean

# Actually remove invalid/empty files
python3 scripts/rebuild_json.py clean --execute
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

## Data Directory Structure

- **`data/pdb/`** - PDB structure files
- **`data/templates/`** - Template files (Atomic_*.pdb) for base frame calculation
- **`data/json/`** - Modern JSON output files (segmented by record type)
- **`data/json_legacy/`** - Legacy JSON output files (reference for comparison)
- **`data/test_sets/`** - Test set definitions (test_set_10.json, test_set_100.json, etc.)

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

