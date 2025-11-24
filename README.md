# X3DNA Modernized - Modern C++ Library

A modern C++ rewrite of the X3DNA v2.4 codebase with strong object-oriented design, comprehensive testing, and full JSON serialization support.

## Status

🚧 **In Development** - Stage 0: Project Setup

## Compiling and Running the Org Code

The `org/` directory contains the original X3DNA code with JSON output support. To compile and run it:

### Compilation

**Using Makefile (recommended):**

```bash
# Build in Release mode (default, optimized)
make org-release

# Build in Debug mode
make org-debug

# Build in RelWithDebInfo mode (optimized with debug symbols)
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
make
```

This will create executables in `org/build/bin/`:
- `find_pair_analyze` - Main executable that processes PDB files and generates JSON output
- `find_pair_original` - Original find_pair executable
- `analyze_original` - Original analyze executable

### Running

To process a PDB file and generate JSON output:

```bash
# From the org/ directory
cd org
./build/bin/find_pair_analyze ../data/pdb/1H4S.pdb
```

The JSON output will be written to `data/json_legacy/1H4S.json`.

**Note:** The executable should be run from the `org/` directory (or with the PDB path relative to the project root), as it may use relative paths for templates or other resources.

### Build Types

You can specify different build types using the Makefile:

- **Release** (default, optimized with `-O3`): `make org-release`
- **Debug** (with debug symbols, `-g`, `-O0`): `make org-debug`
- **RelWithDebInfo** (optimized with debug symbols, `-O2`, `-g`): `make org-relwithdebinfo`

Or manually with CMake:

- **Release**:
  ```bash
  cmake .. -DCMAKE_BUILD_TYPE=Release
  ```

- **Debug**:
  ```bash
  cmake .. -DCMAKE_BUILD_TYPE=Debug
  ```

- **RelWithDebInfo**:
  ```bash
  cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
  ```

## Quick Start

### Building

```bash
mkdir build && cd build
cmake ..
make
```

Or use the build script:

```bash
./scripts/build.sh Release
```

### Running Tests

```bash
cd build
ctest
```

## JSON Comparison and Rebuilding

The project includes two main Python scripts for working with JSON files:

### `compare_json.py` - Compare JSON Files

Compare legacy and modern JSON files to identify differences.

**Basic Usage:**
```bash
# Compare all available files (atoms and frames)
python3 scripts/compare_json.py compare

# Compare specific PDB file(s)
python3 scripts/compare_json.py compare 1H4S
python3 scripts/compare_json.py compare 1H4S 2BNA 3DNA

# Compare only atoms
python3 scripts/compare_json.py atoms 1H4S

# Compare only frames
python3 scripts/compare_json.py frames 1H4S

# Compare ring atom matching
python3 scripts/compare_json.py ring-atoms 1H4S

# List available PDB files
python3 scripts/compare_json.py list
```

**Options:**
```bash
# Use legacy mode JSON files
python3 scripts/compare_json.py compare --legacy-mode

# Save report to file
python3 scripts/compare_json.py compare --output report.md

# Verbose output with detailed differences
python3 scripts/compare_json.py compare --verbose

# Show only files with differences
python3 scripts/compare_json.py compare --diff-only

# Custom thread count
python3 scripts/compare_json.py compare --threads 4

# Regenerate missing JSON files automatically
python3 scripts/compare_json.py compare --regenerate

# Use a test set (10, 50, 100, 500, or 1000 PDBs)
python3 scripts/compare_json.py compare --test-set 10
python3 scripts/compare_json.py compare --test-set 50
```

**Test Sets:**
```bash
# Generate test sets of different sizes (10, 50, 100, 500, 1000)
python3 scripts/compare_json.py generate-test-sets

# Regenerate test sets (overwrite existing)
python3 scripts/compare_json.py generate-test-sets --force

# Use custom random seed for reproducibility
python3 scripts/compare_json.py generate-test-sets --seed 12345
```

Test sets are saved in `data/test_sets/` and can be reused for consistent testing across different runs.

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

# Use legacy mode for modern JSON generation
python3 scripts/rebuild_json.py regenerate --legacy-mode

# Custom thread count
python3 scripts/rebuild_json.py regenerate --threads 4

# Use a test set (10, 50, 100, 500, or 1000 PDBs)
python3 scripts/rebuild_json.py regenerate --test-set 10
python3 scripts/rebuild_json.py regenerate --test-set 50
```

**Test Sets:**
```bash
# Generate test sets of different sizes (10, 50, 100, 500, 1000)
python3 scripts/rebuild_json.py generate-test-sets

# Regenerate test sets (overwrite existing)
python3 scripts/rebuild_json.py generate-test-sets --force

# Use custom random seed for reproducibility
python3 scripts/rebuild_json.py generate-test-sets --seed 12345
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

# Clean only legacy or modern
python3 scripts/rebuild_json.py clean --legacy-only --execute
python3 scripts/rebuild_json.py clean --modern-only --execute
```

**Reorganize JSON Files:**
```bash
# Reorganize JSON file (array to grouped format)
python3 scripts/rebuild_json.py reorganize data/json/1H4S.json
```

For more detailed documentation, see `scripts/README.md`.

## Project Structure

```
find_pair_2/
├── include/x3dna/     # Public headers
├── src/x3dna/         # Implementation files
├── apps/              # Application executables
├── tests/             # Test suite
│   ├── unit/          # Unit tests
│   ├── integration/   # Integration tests
│   └── regression/   # Regression tests
└── docs/              # Documentation
    └── modernization/ # Stage-by-stage plan
```

## Documentation

- **Modernization Plan**: `docs/MODERNIZATION_PLAN.md`
- **OOP Class Hierarchy**: `docs/OOP_CLASS_HIERARCHY.md`
- **Stage-by-Stage Plan**: `docs/modernization/`
- **Integration Testing**: `docs/modernization/INTEGRATION_TESTING.md`

## Dependencies

- C++17 or later
- CMake 3.15+
- nlohmann/json (header-only, fetched automatically)
- Google Test (for testing, fetched automatically)

## License

[To be determined]

## Contributing

See `docs/modernization/` for the detailed implementation plan.

