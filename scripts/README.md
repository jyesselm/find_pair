# Scripts Reference Guide

Complete reference for all scripts in the `scripts/` directory.

---

## JSON Comparison

### `compare_json.py` - Unified JSON Comparison Tool

**Purpose:** Multithreaded tool for comparing JSON files between legacy and modern implementations with detailed human-readable reports.

**Usage:**
```bash
# Compare all available files
python3 scripts/compare_json.py

# Compare specific PDB file(s)
python3 scripts/compare_json.py 1H4S
python3 scripts/compare_json.py 1H4S 2BNA 3DNA

# Use legacy mode JSON files
python3 scripts/compare_json.py --legacy-mode

# Save report to file
python3 scripts/compare_json.py --output report.md

# Verbose output with detailed differences
python3 scripts/compare_json.py --verbose

# Show only files with differences
python3 scripts/compare_json.py --diff-only

# Custom thread count
python3 scripts/compare_json.py --threads 4
```

**Features:**
- ✅ Multithreaded processing (uses all CPU cores by default)
- ✅ Detailed human-readable reports
- ✅ Supports single file or batch comparison
- ✅ Legacy mode support
- ✅ Verbose output option
- ✅ Save reports to file
- ✅ Filter to show only differences

**Report Includes:**
- Summary: Total files, matches, differences, errors
- Frame Calculation Statistics: Total residues, missing residues, mismatches
- Atom Matching Statistics: Exact matches, set matches, real differences
- RMS Fit Differences: Average and maximum RMS differences
- Base Type Statistics: Breakdown by base type (A, G, C, T, U)
- Files with Most Differences: Top 20 files sorted by difference count
- Detailed Differences: (when using `--verbose` flag)

---

## Maintenance Scripts

### `refresh_all_json.py` - Refresh All JSON Files

**Purpose:** Refreshes all JSON files (both legacy and modern) with validation. Automatically removes corrupted or empty JSON files.

**Usage:**
```bash
# Refresh all JSON files (legacy and modern)
python3 scripts/refresh_all_json.py

# Validate existing JSON files only (don't regenerate)
python3 scripts/refresh_all_json.py --validate-only

# Refresh only legacy JSON files
python3 scripts/refresh_all_json.py --legacy-only

# Refresh only modern JSON files
python3 scripts/refresh_all_json.py --modern-only
```

**Features:**
- Automatic validation
- Automatic cleanup of invalid/corrupted/empty JSON files
- Parallel processing for fast execution
- Progress reporting and summary statistics

---

### `regenerate_legacy_json.py` - Regenerate Legacy JSON Files

**Purpose:** Regenerates legacy JSON files for problematic PDBs using the original C code.

**Usage:**
```bash
python3 scripts/regenerate_legacy_json.py
```

**Features:**
- Uses original C code for regeneration
- Automatic validation
- Automatically removes invalid files
- Parallel processing

---

### `regenerate_problematic_pdbs.py` - Regenerate Problematic PDBs

**Purpose:** Regenerate modern JSON files only for PDBs with known differences. Uses diff lists from docs/ directory.

**Usage:**
```bash
python3 scripts/regenerate_problematic_pdbs.py
```

---

### `verify_all_legacy_mode.py` - Comprehensive Legacy Mode Verification

**Purpose:** Comprehensive legacy mode verification for all PDB files. Regenerates JSON with --legacy flag and compares with legacy JSON.

**Usage:**
```bash
python3 scripts/verify_all_legacy_mode.py
```

**Features:**
- Regenerates JSON with --legacy flag
- Compares with legacy JSON
- Parallel processing
- Detailed statistics

---

### `clean_json_legacy.py` - Clean Legacy JSON Files

**Purpose:** Removes corrupted or empty JSON files from json_legacy directory.

**Usage:**
```bash
# Dry run (default)
python3 scripts/clean_json_legacy.py

# Actually remove files
python3 scripts/clean_json_legacy.py --execute
```

**Features:**
- Uses validation for detection
- Dry-run mode by default
- Reports summary of files to remove

---

### `remove_empty_json.py` - Remove Empty JSON Files

**Purpose:** Removes empty JSON files from a directory.

**Usage:**
```bash
# Dry run (default)
python3 scripts/remove_empty_json.py

# Actually remove files
python3 scripts/remove_empty_json.py --execute

# Non-recursive search
python3 scripts/remove_empty_json.py --no-recursive
```

**Features:**
- Recursive search by default
- Dry-run mode by default
- Uses validation for detection

---

## Analysis Scripts

### `analyze_real_issues.py` - Analyze Real Atom Set Differences

**Purpose:** Analyze real atom set differences and significant RMS differences. Focuses on understanding what's wrong and if legacy mode fixes it.

**Usage:**
```bash
python3 scripts/analyze_real_issues.py
```

---

### `analyze_base_type_differences.py` - Analyze Base Type Differences

**Purpose:** Analyze differences for a specific base type. Extracts and analyzes residues of a specific base type that have differences between legacy and modern JSON files.

**Usage:**
```bash
# Analyze A (adenine) differences
python3 scripts/analyze_base_type_differences.py A

# Limit examples shown
python3 scripts/analyze_base_type_differences.py A --limit 50

# Verbose output
python3 scripts/analyze_base_type_differences.py a --verbose
```

**Features:**
- Identifies common difference patterns
- Shows top patterns by frequency
- Provides detailed examples
- Saves JSON report for further analysis

---

### `analyze_case_sensitivity.py` - Analyze Case Sensitivity Differences

**Purpose:** Investigate why lowercase bases (especially 'a') have much worse match rates than uppercase bases. Compares uppercase vs lowercase base handling.

**Usage:**
```bash
# Analyze all case sensitivity issues
python3 scripts/analyze_case_sensitivity.py

# Focus on specific base type
python3 scripts/analyze_case_sensitivity.py --base a

# Verbose output
python3 scripts/analyze_case_sensitivity.py --verbose
```

**Features:**
- Compares uppercase vs lowercase match rates
- Identifies residue name patterns
- Shows examples of problematic lowercase residues
- Analyzes missing residues

---

### `analyze_high_rms_differences.py` - Analyze High RMS Differences

**Purpose:** Extract and analyze residues with RMS differences above a threshold to understand if these are correlated with atom set differences.

**Usage:**
```bash
# Analyze RMS differences >= 0.1 Å
python3 scripts/analyze_high_rms_differences.py

# Custom threshold
python3 scripts/analyze_high_rms_differences.py --threshold 1.0

# Limit examples
python3 scripts/analyze_high_rms_differences.py --threshold 0.1 --limit 50
```

**Features:**
- Identifies residues with high RMS differences
- Correlates RMS differences with atom set differences
- Shows distribution of RMS differences
- Base type statistics

---

### `analyze_ring_fitting_differences.py` - Analyze Ring Fitting Differences

**Purpose:** Comprehensive analysis of ring fitting differences between legacy and modern code. Analyzes JSON files and PDB files to understand why different atoms are matched during ring fitting.

**Usage:**
```bash
python3 scripts/analyze_ring_fitting_differences.py <pdb_id>
```

---

### `analyze_all_ring_fitting_differences.py` - Analyze All Ring Fitting Differences

**Purpose:** Analyze ring fitting differences for all PDBs with differences. Processes all PDBs and generates comprehensive reports.

**Usage:**
```bash
python3 scripts/analyze_all_ring_fitting_differences.py
```

---

### `diagnose_modern_frame_calculation.py` - Diagnose Frame Calculation Issues

**Purpose:** Diagnose why modern code isn't calculating frames. Checks if atoms are being parsed correctly, if residue types are being detected correctly, etc.

**Usage:**
```bash
python3 scripts/diagnose_modern_frame_calculation.py <pdb_id>
```

---

### `test_residue_name_detection.py` - Test Residue Name Detection

**Purpose:** Test residue name detection - compare legacy vs modern logic.

**Usage:**
```bash
python3 scripts/test_residue_name_detection.py
```

---

## Debug Scripts

Located in `scripts/debug/` subdirectory.

### `debug/analyze_frame_failures.py` - Analyze Frame Calculation Failures

**Purpose:** Analyze frame calculation failures from detailed JSON reports. Focuses on actual numerical differences, not just atom list ordering.

**Usage:**
```bash
python3 scripts/debug/analyze_frame_failures.py
```

---

### `debug/debug_specific_failure.py` - Debug Specific Failure

**Purpose:** Debug a specific failing residue to understand atom matching differences.

**Usage:**
```bash
python3 scripts/debug/debug_specific_failure.py <pdb_name> <chain_id> <seq_num>
```

---

### `debug/compare_problematic_pdbs.py` - Compare Problematic PDBs

**Purpose:** Compare problematic PDBs between legacy and modern code, generating detailed reports with actual PDB line content for differences.

**Usage:**
```bash
python3 scripts/debug/compare_problematic_pdbs.py
```

---

### `debug/analyze_missing_frames.py` - Analyze Missing Frames

**Purpose:** Analyze why modern JSON is missing frame calculations. Investigates specific PDBs to understand why residues that have frame calculations in legacy JSON don't appear in modern JSON.

**Usage:**
```bash
python3 scripts/debug/analyze_missing_frames.py
```

---

## Validation Scripts

### `validate_removed_atoms.py` - Validate Removed Atoms

**Purpose:** Validation script that uses removed atom information from legacy JSON to verify that the modern parser correctly filters out atoms.

**Usage:**
```bash
python3 scripts/validate_removed_atoms.py <pdb_id>
```

---

## Utility Scripts

### `quick_test.py` - Quick Test Script

**Purpose:** Quick testing utility.

**Usage:**
```bash
python3 scripts/quick_test.py
```

---

### `rerun_check.py` - Rerun Verification

**Purpose:** Rerun verification using subprocess with clean environment.

**Usage:**
```bash
python3 scripts/rerun_check.py
```

---

## Shell Scripts

### `build.sh` - Build Script

Build script for the project.

**Usage:**
```bash
./scripts/build.sh
```

---

### `format_code.sh` - Code Formatting Script

Code formatting script using clang-format.

**Usage:**
```bash
./scripts/format_code.sh
```

---

### `debug/investigate_batch.sh` - Investigate Batch

Shell script for batch investigation.

**Usage:**
```bash
./scripts/debug/investigate_batch.sh
```

---

### `debug/investigate_pdb.sh` - Investigate PDB

Shell script for PDB investigation.

**Usage:**
```bash
./scripts/debug/investigate_pdb.sh <pdb_id>
```

---

### `regenerate_problematic_pdbs.sh` - Regenerate Problematic PDBs (Shell)

Shell script version for regenerating problematic PDBs.

**Usage:**
```bash
./scripts/regenerate_problematic_pdbs.sh
```

---

### `run_verification.sh` - Run Verification

Shell script for running verification.

**Usage:**
```bash
./scripts/run_verification.sh
```

---

### `rerun_verification.sh` - Rerun Verification

Shell script for rerunning verification.

**Usage:**
```bash
./scripts/rerun_verification.sh
```

---

## Library Modules

The scripts use reusable modules from `x3dna_json_compare` package:

- **`json_comparison.py`**: Main JSON comparison engine with caching and parallel processing
- **`json_validator.py`**: JSON validation utilities
- **`atom_comparison.py`**: Atom-level comparison
- **`frame_comparison.py`**: Frame calculation comparison
- **`pdb_utils.py`**: PDB file utilities
- **`result_cache.py`**: Caching system for comparison results
- **`parallel_executor.py`**: Parallel processing utilities
- **`models.py`**: Data models for comparison results

---

## Best Practices

1. **Always use the unified comparison script**: Use `compare_json.py` for all JSON comparisons
2. **Use validation**: Always validate JSON files before using them
3. **Dry-run first**: Use dry-run mode for destructive operations
4. **Check return codes**: Scripts return appropriate exit codes
5. **Progress reporting**: Long-running scripts show progress
6. **Multithreading**: Use multithreaded scripts for batch operations

---

## Quick Reference

| Task | Script |
|------|--------|
| Compare JSON files | `compare_json.py` |
| Analyze base type differences | `analyze_base_type_differences.py` |
| Analyze case sensitivity | `analyze_case_sensitivity.py` |
| Analyze high RMS differences | `analyze_high_rms_differences.py` |
| Refresh all JSON files | `refresh_all_json.py` |
| Regenerate legacy JSON | `regenerate_legacy_json.py` |
| Clean corrupted JSON | `clean_json_legacy.py` |
| Remove empty JSON | `remove_empty_json.py` |
| Validate removed atoms | `validate_removed_atoms.py` |
| Debug specific residue | `debug/debug_specific_failure.py` |
| Analyze frame failures | `debug/analyze_frame_failures.py` |
| Analyze ring fitting | `analyze_ring_fitting_differences.py` |

---

## Fix Plan

For a comprehensive plan to fix JSON comparison differences, see:
- **`FIX_PLAN.md`** (in project root) - Detailed 4-phase plan to address differences

The plan covers:
- Phase 1: Investigation & Analysis
- Phase 2: Root Cause Fixes
- Phase 3: Validation & Testing
- Phase 4: Documentation & Cleanup

---

*Last updated: 2025-11-22*
