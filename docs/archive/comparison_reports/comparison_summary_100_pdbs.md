# Comparison Summary: 100 PDBs

**Date**: 2025-11-27  
**Test Set**: `data/test_sets/test_set_100.json`

## Executive Summary

Out of 100 PDBs in the test set:
- **29 PDBs** have both legacy and modern JSON files available for comparison
- **71 PDBs** are missing legacy JSON files (need to be generated)
- **0 PDBs** are missing modern JSON files

## Comparison Results (29 PDBs with Both Files)

### Overall Status
- **3 PDBs** show perfect matches across all comparison types
- **26 PDBs** show some differences in specific categories

### Detailed Statistics

#### ✅ Perfect Match Categories
1. **FIND_BESTPAIR Selection** (Actual Selected Pairs)
   - Legacy: 677 pairs
   - Modern: 677 pairs
   - **100% Match** - All selected pairs match exactly

2. **Residue Indices**
   - **100% Match** - All residue indices match perfectly across 26 PDBs

#### ⚠️ Categories with Differences

1. **Base Pairs** (Pairs Passing All Checks)
   - Legacy: 1,513 pairs
   - Modern: 800 pairs
   - Common: 781 pairs
   - Missing in modern: 732 pairs
   - Extra in modern: 19 pairs
   - Mismatched: 702 pairs
   - **Status**: Significant differences in base pair detection/validation

2. **Pair Validation** (All Validated Pairs)
   - Legacy: 5,164 validations
   - Modern: 66,014 validations
   - Missing in modern: 2,173 pairs
   - Extra in modern: 63,023 pairs
   - Mismatched: 518 validations
   - **Status**: Large discrepancy - modern code validates many more pairs

3. **H-Bond Lists**
   - Legacy: 800 H-bond lists
   - Modern: 800 H-bond lists
   - Common pairs: 800
   - Mismatched pairs: 404
   - **Status**: All pairs have H-bond lists, but 404 pairs have mismatched H-bond data

## PDBs with Both Files (29 total)

These PDBs were successfully compared:
- 1F7Y, 1Q96, 1VBY, 1VC0, 1ZX7
- 2B8R, 2F4S, 2PXE, 2YGH
- 3C5D, 3KNC
- 485D
- 4KYY, 4MCF
- 5DHC, 5L00, 5W7O
- 6HC5, 6N5O, 6N5P, 6N5Q, 6ZWU
- 7D7X, 7L0Z, 7Q82, 7Y2P
- 8CLM, 8U5Z, 8VFS

## PDBs Missing Legacy JSON (71 total)

These PDBs need legacy JSON files to be generated:
- 1QU3, 1T0K
- 2ATW, 2IZ8, 2O5J, 2QEX, 2ZH2
- 3AVY, 3CF5, 3CME, 3FOZ, 3DIZ, 3G78, 3G8T, 3IWN, 3RZD
- 4AL5, 4C8Y, 4IFD, 4JV5, 4KZE, 4M2Z, 4Q5V
- 5AXM, 5IWA, 5JXS, 5MSF, 5O7H, 5ON3, 5OOM, 5T16, 5UJ2, 5VZ8, 5XN0, 5XH7
- 6CAP, 6CAQ, 6CVQ, 6G5I, 6HLQ, 6IFU, 6J6G, 6K32, 6LTU, 6RJA, 6T0U, 6ZXH
- 7K9E, 7MKO, 7R0E, 7S36, 7S3B, 7TO2, 7UXA, 7XHT, 7XUE, 7Y7P, 7Y8Y, 7Z4H
- 8D5L, 8DH1, 8I7N, 8J1J, 8J3R, 8J9G, 8OKD, 8T2T, 8VMA, 8WAW
- 9G9K, 9JA1

## Key Findings

### Positive Results ✅
1. **Selected pairs match perfectly** - The final output (find_bestpair_selection) matches 100% between legacy and modern code
2. **Residue indices match perfectly** - All residue ordering and indexing is correct
3. **H-bond lists exist for all pairs** - No missing H-bond data

### Areas Needing Attention ⚠️
1. **Base pair detection differences** - Modern code finds fewer base pairs (800 vs 1,513), suggesting stricter validation criteria
2. **Pair validation differences** - Modern code validates many more pairs (66,014 vs 5,164), suggesting different validation thresholds
3. **H-bond data mismatches** - 404 pairs have different H-bond data, which may affect quality scores

## Recommendations

1. **Generate legacy JSON files** for the 71 missing PDBs to enable full comparison:
   ```bash
   python3 scripts/compare_json.py compare --test-set 100 --regenerate
   ```

2. **Investigate base pair differences** - The large discrepancy (732 missing pairs) suggests validation criteria may be too strict in modern code

3. **Investigate pair validation differences** - The huge difference (63,023 extra validations) suggests validation thresholds may be too lenient in modern code

4. **Review H-bond mismatches** - 404 pairs with mismatched H-bond data need investigation to ensure quality score calculations match

## Full Report

See `comparison_report_100_pdbs.md` for detailed comparison results.

