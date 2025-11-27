# Script Cleanup Plan

## Scripts to Remove (Redundant/Obsolete)

### 1. Analysis Scripts
- `deep_analyze_find_bestpair_differences.py` - Functionality merged into `analyze_find_bestpair_differences.py`
- `analyze_remaining_differences.py` - Obsolete (we fixed the 0.5% differences)
- `analyze_mismatched_pairs.py` - Functionality available in `compare_json.py`

### 2. Generation Scripts (Replaced by `rebuild_json.py`)
- `generate_legacy_json_batch.py` - Use `rebuild_json.py regenerate --legacy-only`
- `generate_missing_modern_json.py` - Use `rebuild_json.py regenerate --modern-only`
- `check_and_generate_json.py` - Functionality in `rebuild_json.py`

### 3. Test Scripts
- `test_tie_breaking_fix.py` - Obsolete (fix is complete and tested)
- `test_all_small_pdbs.py` - Consolidate into `test_multiple_pdbs.py`
- `test_pdbs_by_size.py` - Consolidate into `test_multiple_pdbs.py`

### 4. Comparison Scripts
- `compare_existing_legacy_pdbs.py` - Functionality in `compare_json.py`

### 5. Cleanup Scripts
- `cleanup_large_legacy_json.py` - One-time cleanup, no longer needed

## Scripts to Keep

### Core Tools
- `compare_json.py` - Main comparison tool
- `rebuild_json.py` - Unified JSON generation/rebuilding tool

### Analysis Tools
- `analyze_find_bestpair_differences.py` - Main analysis tool for differences

### Investigation Tools
- `investigate_specific_pairs.py` - Debug specific pairs
- `investigate_bp_type_id_differences.py` - Debug bp_type_id issues
- `investigate_hbond_mismatches.py` - Debug hbond issues

### Testing Tools
- `test_multiple_pdbs.py` - General testing framework

### Utility Tools
- `cleanup_optional_json.py` - Cleanup optional JSON files
- `cleanup_json_data.py` - General JSON cleanup
- `check_legacy_json_progress.py` - Check generation progress
- `parse_legacy_debug_output.py` - Parse debug output
- `compare_step_params_for_pairs.py` - Compare step parameters

### Specialized Tools
- `generate_hbond_list_from_json.py` - Generate hbond lists
- `create_legacy_hbond_list.py` - Create legacy hbond lists
- `generate_and_compare_residue_ordering_batch.py` - Residue ordering comparison

## Migration Guide

### Replacing Removed Scripts

**Instead of `generate_legacy_json_batch.py`:**
```bash
python3 scripts/rebuild_json.py regenerate --legacy-only --test-set 100
```

**Instead of `generate_missing_modern_json.py`:**
```bash
python3 scripts/rebuild_json.py regenerate --modern-only --test-set 100
```

**Instead of `test_all_small_pdbs.py` or `test_pdbs_by_size.py`:**
```bash
python3 scripts/test_multiple_pdbs.py --filter small
python3 scripts/test_multiple_pdbs.py --filter large
```

**Instead of `analyze_remaining_differences.py`:**
```bash
python3 scripts/analyze_find_bestpair_differences.py --test-set 100
```

