# Scripts Reference Guide

Complete reference for all scripts in the `scripts/` directory.

**⚠️ Important**: For complete testing and JSON comparison documentation, see **[docs/TESTING_GUIDE.md](../docs/TESTING_GUIDE.md)**.

**For information about deprecated/redundant scripts, see [docs/REDUNDANT_SCRIPTS.md](../docs/REDUNDANT_SCRIPTS.md)**.

---

## Main Scripts

### `compare_json.py` - Unified JSON Comparison Tool

**Purpose:** Comprehensive tool for comparing JSON files between legacy and modern implementations with detailed human-readable reports.

**Usage:**
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

# Compare ring atoms
python3 scripts/compare_json.py ring-atoms 1H4S
python3 scripts/compare_json.py ring-atoms  # All files

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

# List available PDB files
python3 scripts/compare_json.py list
```

**Features:**
- ✅ Multithreaded processing (uses all CPU cores by default)
- ✅ Detailed human-readable reports
- ✅ Supports single file or batch comparison
- ✅ Legacy mode support
- ✅ Verbose output option
- ✅ Save reports to file
- ✅ Filter to show only differences
- ✅ Atom comparison
- ✅ Frame calculation comparison
- ✅ Ring atom matching comparison

**Report Includes:**
- Summary: Total files, matches, differences, errors
- Frame Calculation Statistics: Total residues, missing residues, mismatches
- Atom Matching Statistics: Exact matches, set matches, real differences
- RMS Fit Differences: Average and maximum RMS differences
- Base Type Statistics: Breakdown by base type (A, G, C, T, U)
- Files with Most Differences: Top 20 files sorted by difference count
- Detailed Differences: (when using `--verbose` flag)

---

### `rebuild_json.py` - Unified JSON Rebuilding Tool

**Purpose:** Comprehensive tool for rebuilding, regenerating, validating, and cleaning JSON files (both legacy and modern).

**Usage:**
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

# Validate existing JSON files
python3 scripts/rebuild_json.py validate

# Clean invalid/empty JSON files (dry-run by default)
python3 scripts/rebuild_json.py clean

# Actually remove invalid files
python3 scripts/rebuild_json.py clean --execute

# Clean only legacy or modern
python3 scripts/rebuild_json.py clean --legacy-only --execute
python3 scripts/rebuild_json.py clean --modern-only --execute

# Reorganize JSON file (array to grouped format)
python3 scripts/rebuild_json.py reorganize data/json/1H4S.json

# Custom thread count
python3 scripts/rebuild_json.py regenerate --threads 4
```

**Features:**
- ✅ Regenerate legacy JSON files (from org code)
- ✅ Regenerate modern JSON files (from new code)
- ✅ Validate existing JSON files
- ✅ Clean corrupted/empty JSON files
- ✅ Reorganize JSON structure
- ✅ Parallel processing for fast execution
- ✅ Automatic validation
- ✅ Progress reporting and summary statistics

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

1. **Use the main scripts**: Use `compare_json.py` for all JSON comparisons and `rebuild_json.py` for all JSON rebuilding operations
2. **See TESTING_GUIDE.md**: For complete testing workflow, see [docs/TESTING_GUIDE.md](../docs/TESTING_GUIDE.md)
3. **Avoid deprecated scripts**: See [docs/REDUNDANT_SCRIPTS.md](../docs/REDUNDANT_SCRIPTS.md) for deprecated scripts and migration guide
2. **Use validation**: Always validate JSON files before using them
3. **Dry-run first**: Use dry-run mode for destructive operations (clean command)
4. **Check return codes**: Scripts return appropriate exit codes
5. **Progress reporting**: Long-running scripts show progress
6. **Multithreading**: Use multithreaded scripts for batch operations

---

## Quick Reference

| Task | Script |
|------|--------|
| Compare JSON files | `compare_json.py compare` |
| Compare atoms only | `compare_json.py atoms` |
| Compare frames only | `compare_json.py frames` |
| Compare ring atoms | `compare_json.py ring-atoms` |
| Regenerate all JSON | `rebuild_json.py regenerate` |
| Regenerate legacy JSON | `rebuild_json.py regenerate --legacy-only` |
| Regenerate modern JSON | `rebuild_json.py regenerate --modern-only` |
| Validate JSON files | `rebuild_json.py validate` |
| Clean invalid JSON | `rebuild_json.py clean --execute` |
| Reorganize JSON | `rebuild_json.py reorganize <file>` |
| Analyze base type differences | `analyze_base_type_differences.py` |
| Analyze case sensitivity | `analyze_case_sensitivity.py` |
| Analyze high RMS differences | `analyze_high_rms_differences.py` |
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
