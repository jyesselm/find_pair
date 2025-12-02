# Comparison Index Fixes Summary

**Date**: 2025-01-XX  
**Status**: ✅ Completed

---

## Problem

Comparison tools were not consistently using `legacy_residue_idx`, leading to comparisons of different residues between legacy and modern code.

---

## Fixes Applied

### 1. `scripts/compare_frames.py`
- ✅ Changed to prefer `legacy_residue_idx` over `residue_idx`
- ✅ Added warning when `residue_idx` is used instead
- ✅ Updated documentation to specify legacy indices are required

### 2. `tools/compare_quality_scores.cpp`
- ✅ Changed to prefer `base_i`/`base_j` (always legacy indices) over `residue1_idx`/`residue2_idx`
- ✅ Added comments explaining why `base_i`/`base_j` are preferred
- ✅ Updated usage message to specify legacy indices

### 3. Documentation
- ✅ Created `COMPARISON_INDEX_RULES.md` with comprehensive rules
- ✅ Created this summary document

---

## Rules Established

1. **Always prefer `legacy_residue_idx`** over `residue_idx` in frame records
2. **Always prefer `base_i`/`base_j`** over `residue1_idx`/`residue2_idx` in pair records
3. **Document all command-line arguments** to specify legacy indices are required
4. **Add warnings** when non-legacy indices are used

---

## Files Still Needing Review

The following files should be reviewed to ensure they use legacy indices:

### Python Scripts
- [ ] `scripts/compare_base_pairs.py`
- [ ] `scripts/compare_validation_geometry.py`
- [ ] `scripts/compare_best_partner.py`
- [ ] `scripts/compare_mutual_best.py`
- [ ] `scripts/compare_iteration.py`
- [ ] `x3dna_json_compare/*.py` - All comparison modules

### C++ Tools
- [ ] `tools/compare_hbond_detection.cpp`
- [ ] `tools/compare_validation_discrepancy.cpp`
- [ ] `tools/compare_hbond_stages.cpp`
- [ ] `tools/compare_initial_hbonds.cpp`
- [ ] All other `tools/compare_*.cpp` files

---

## Testing

Before using any comparison tool:
1. Verify modern JSON was generated with `--fix-indices`
2. Check that JSON records contain `legacy_residue_idx` or `base_i`/`base_j`
3. Verify tool documentation states it uses legacy indices

---

## Reference

See `COMPARISON_INDEX_RULES.md` for complete rules and patterns.

