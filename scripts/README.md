# Scripts Directory

## Overview

This directory contains scripts organized into:
- **Build/Maintenance scripts**: For building, regenerating, and maintaining JSON files
- **Debug scripts**: For debugging and analyzing differences (in `debug/` subdirectory)
- **lib/**: Reusable Python modules for JSON comparison, validation, and utilities

## Directory Structure

```
scripts/
├── lib/                    # Reusable Python modules
│   ├── json_comparison.py  # JSON comparison engine
│   ├── json_validator.py   # JSON validation utilities
│   ├── atom_comparison.py  # Atom-level comparison
│   ├── frame_comparison.py # Frame calculation comparison
│   ├── pdb_utils.py        # PDB file utilities
│   ├── result_cache.py     # Caching system
│   ├── parallel_executor.py # Parallel processing
│   └── models.py           # Data models
├── debug/                  # Debug and analysis scripts
│   ├── analyze_frame_failures.py
│   ├── compare_problematic_pdbs.py
│   ├── debug_specific_failure.py
│   ├── investigate_batch.sh
│   └── investigate_pdb.sh
├── build.sh                # Build script
├── format_code.sh          # Code formatting script
├── refresh_all_json.py     # Refresh all JSON files with validation
├── regenerate_legacy_json.py  # Regenerate legacy JSON files
├── regenerate_problematic_pdbs.py  # Regenerate problematic PDBs
├── clean_json_legacy.py    # Clean corrupted JSON files
├── remove_empty_json.py    # Remove empty JSON files
└── validate_removed_atoms.py  # Validate removed atom filtering
```

## Build/Maintenance Scripts

### refresh_all_json.py

Refreshes all JSON files (both legacy and modern) with validation. Automatically removes corrupted or empty JSON files instead of saving them.

**Usage:**
```bash
# Refresh all JSON files (legacy and modern)
python scripts/refresh_all_json.py

# Validate existing JSON files only (don't regenerate)
python scripts/refresh_all_json.py --validate-only

# Refresh only legacy JSON files
python scripts/refresh_all_json.py --legacy-only

# Refresh only modern JSON files
python scripts/refresh_all_json.py --modern-only
```

**Features:**
- Automatic validation using `lib.json_validator`
- Automatic cleanup of invalid/corrupted/empty JSON files
- Parallel processing for fast execution
- Progress reporting and summary statistics

### regenerate_legacy_json.py

Regenerates legacy JSON files for problematic PDBs using the original C code.

**Usage:**
```bash
python scripts/regenerate_legacy_json.py
```

**Features:**
- Uses `lib.json_validator` for validation
- Automatically removes invalid files
- Parallel processing

### clean_json_legacy.py

Removes corrupted or empty JSON files from json_legacy directory.

**Usage:**
```bash
# Dry run (default)
python scripts/clean_json_legacy.py

# Actually remove files
python scripts/clean_json_legacy.py --execute
```

**Features:**
- Uses `lib.json_validator` for validation
- Dry-run mode by default
- Reports summary of files to remove

### remove_empty_json.py

Removes empty JSON files from a directory.

**Usage:**
```bash
# Dry run (default)
python scripts/remove_empty_json.py

# Actually remove files
python scripts/remove_empty_json.py --execute

# Non-recursive search
python scripts/remove_empty_json.py --no-recursive
```

**Features:**
- Uses `lib.json_validator` for validation
- Recursive search by default
- Dry-run mode by default

### validate_removed_atoms.py

Validates that the modern parser correctly filters removed atoms.

**Usage:**
```bash
python scripts/validate_removed_atoms.py <pdb_id>
```

## Debug Scripts

Debug scripts are located in the `debug/` subdirectory. These scripts are for analyzing differences and debugging issues.

### debug/analyze_frame_failures.py

Analyzes frame calculation failures from detailed JSON reports.

**Usage:**
```bash
python scripts/debug/analyze_frame_failures.py
```

### debug/debug_specific_failure.py

Debugs a specific failing residue to understand atom matching differences.

**Usage:**
```bash
python scripts/debug/debug_specific_failure.py <pdb_name> <chain_id> <seq_num>
```

### debug/compare_problematic_pdbs.py

Compares problematic PDBs between legacy and modern code, generating detailed reports.

**Usage:**
```bash
python scripts/debug/compare_problematic_pdbs.py
```

## Reusable Library Modules (lib/)

All scripts use reusable modules from `lib/` to avoid code duplication:

### lib/json_validator.py

Provides JSON validation utilities:
- `JsonValidator.validate_file()` - Validate a JSON file
- `JsonValidator.validate_data()` - Validate JSON data structure
- `JsonValidator.remove_invalid_file()` - Remove invalid files
- `JsonValidator.validate_and_clean()` - Validate and optionally remove

### lib/json_comparison.py

Main JSON comparison engine:
- `JsonComparator` - Compare legacy and modern JSON files
- Caching support
- Parallel processing

### lib/atom_comparison.py

Atom-level comparison utilities:
- `compare_atoms()` - Compare atom records

### lib/frame_comparison.py

Frame calculation comparison:
- `compare_frames()` - Compare frame calculation records

### lib/pdb_utils.py

PDB file utilities:
- `PdbFileReader` - Read and parse PDB files
- Get PDB lines for atoms

### lib/parallel_executor.py

Parallel processing utilities:
- `ParallelExecutor` - Execute tasks in parallel
- `print_progress()` - Progress reporting

### lib/result_cache.py

Caching system for comparison results:
- `ComparisonCache` - Cache comparison results

### lib/models.py

Data models for comparison results:
- `ComparisonResult` - Complete comparison result
- `AtomComparison` - Atom-level comparison result
- `FrameComparison` - Frame calculation comparison result

## Code Reuse

All scripts are refactored to use the reusable `lib/` modules:

- **JSON validation**: All scripts use `lib.json_validator.JsonValidator`
- **JSON comparison**: Debug scripts use `lib.json_comparison.JsonComparator`
- **PDB utilities**: Scripts use `lib.pdb_utils.PdbFileReader`
- **Parallel processing**: Scripts use `lib.parallel_executor.ParallelExecutor`
- **Data models**: Scripts use `lib.models` for structured data

This ensures:
- Consistent validation logic across all scripts
- No code duplication
- Easy maintenance and updates
- Better error handling

## Shell Scripts

### build.sh

Build script for the project.

### format_code.sh

Code formatting script using clang-format.

## Best Practices

1. **Always use lib modules**: Don't duplicate validation or comparison logic
2. **Use validation**: Always validate JSON files before using them
3. **Dry-run first**: Use dry-run mode for destructive operations
4. **Check return codes**: Scripts return appropriate exit codes
5. **Progress reporting**: Long-running scripts show progress
