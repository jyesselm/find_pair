# Scripts Directory

This directory contains Python scripts for comparing, analyzing, and managing JSON data from legacy and modern codebases.

## Core Tools

### Comparison
- **`compare_json.py`** - Main comparison tool for comparing legacy and modern JSON outputs
  ```bash
  python3 scripts/compare_json.py --pdb-id 3CF5
  python3 scripts/compare_json.py --test-set 100
  ```

### JSON Management
- **`rebuild_json.py`** - Unified tool for regenerating, validating, and cleaning JSON files
  ```bash
  # Regenerate all JSON
  python3 scripts/rebuild_json.py regenerate
  
  # Regenerate only modern JSON
  python3 scripts/rebuild_json.py regenerate --modern-only --test-set 100
  
  # Regenerate only legacy JSON
  python3 scripts/rebuild_json.py regenerate --legacy-only --test-set 100
  
  # Validate JSON files
  python3 scripts/rebuild_json.py validate
  
  # Clean invalid JSON files
  python3 scripts/rebuild_json.py clean --execute
  ```

## Analysis Tools

### Difference Analysis
- **`analyze_find_bestpair_differences.py`** - Analyze differences in find_bestpair_selection records
  ```bash
  python3 scripts/analyze_find_bestpair_differences.py --test-set 100
  python3 scripts/analyze_find_bestpair_differences.py 3CF5
  ```

### Investigation Tools
- **`investigate_specific_pairs.py`** - Investigate why specific pairs differ
  ```bash
  python3 scripts/investigate_specific_pairs.py 3CME 3844 3865
  ```

- **`investigate_bp_type_id_differences.py`** - Debug bp_type_id calculation differences

- **`investigate_hbond_mismatches.py`** - Debug hydrogen bond mismatches

## Testing Tools

- **`test_multiple_pdbs.py`** - Test multiple PDBs with various filters
  ```bash
  python3 scripts/test_multiple_pdbs.py --filter small
  python3 scripts/test_multiple_pdbs.py --filter large
  ```

## Cleanup Tools

- **`cleanup_optional_json.py`** - Remove optional JSON record types to save space
  ```bash
  python3 scripts/cleanup_optional_json.py --execute
  ```

- **`cleanup_json_data.py`** - General JSON cleanup utilities

## Utility Tools

- **`check_legacy_json_progress.py`** - Check progress of legacy JSON generation

- **`parse_legacy_debug_output.py`** - Parse debug output from legacy code

- **`compare_step_params_for_pairs.py`** - Compare step parameters for specific pairs

## Specialized Tools

- **`generate_hbond_list_from_json.py`** - Generate hydrogen bond lists from JSON

- **`create_legacy_hbond_list.py`** - Create legacy hydrogen bond lists

- **`generate_and_compare_residue_ordering_batch.py`** - Compare residue ordering across PDBs

## Deprecated Scripts

The following scripts have been removed and their functionality consolidated:

- `deep_analyze_find_bestpair_differences.py` → Use `analyze_find_bestpair_differences.py`
- `analyze_remaining_differences.py` → Use `analyze_find_bestpair_differences.py`
- `analyze_mismatched_pairs.py` → Use `compare_json.py`
- `generate_legacy_json_batch.py` → Use `rebuild_json.py regenerate --legacy-only`
- `generate_missing_modern_json.py` → Use `rebuild_json.py regenerate --modern-only`
- `test_tie_breaking_fix.py` → Obsolete (fix complete)
- `test_all_small_pdbs.py` → Use `test_multiple_pdbs.py`
- `test_pdbs_by_size.py` → Use `test_multiple_pdbs.py`
- `compare_existing_legacy_pdbs.py` → Use `compare_json.py`
- `cleanup_large_legacy_json.py` → One-time cleanup, no longer needed

See `CLEANUP_PLAN.md` for details on the cleanup and migration guide.
