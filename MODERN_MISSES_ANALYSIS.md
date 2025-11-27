# Analysis: Does Modern Code Miss Pairs?

## Summary

**Yes, but it's complicated.** Modern code can miss pairs that legacy found, but in most cases this appears to be due to:

1. **Different residue indexing schemes** - The pairs may be the same but indexed differently
2. **Complete set differences** - Some PDBs show NO overlap between legacy and modern pairs, suggesting different interpretation of the structure
3. **Validation threshold differences** - Modern may have stricter validation in some cases

## Statistics

Based on analysis of cases where both found pairs but sets differ:

- **7 cases** where modern missed pairs that legacy found
- **Total pairs missed**: 220 (avg 31.4 per case)
- **Total pairs modern found extra**: 242 (avg 34.6 per case)

## Cases Where Modern Missed Pairs

| PDB  | Legacy | Modern | Missed | Extra | Notes |
|------|--------|--------|--------|-------|-------|
| 1ASY | 49     | 48     | 49     | 48    | NO overlap - completely different sets |
| 1ASZ | 47     | 48     | 47     | 48    | NO overlap - completely different sets |
| 1EFW | 46     | 49     | 46     | 49    | NO overlap - completely different sets |
| 1E8O | 21     | 20     | 21     | 20    | NO overlap - completely different sets |
| 1B23 | 20     | 26     | 20     | 26    | NO overlap - completely different sets |
| 1C0A | 19     | 25     | 19     | 25    | NO overlap - completely different sets |
| 1D4R | 18     | 26     | 18     | 26    | NO overlap - completely different sets |

## Key Finding: No Overlap Cases

**All 7 cases show NO overlap** between legacy and modern pair sets. This suggests:

1. **Residue indexing differences** - Same pairs, different numbering
2. **Chain interpretation differences** - Different chain assignments
3. **Structure parsing differences** - Different interpretation of the PDB file

## Analysis

### What This Means

The fact that there's NO overlap in these cases strongly suggests:

1. **Same structure, different indexing**: Legacy and modern may be identifying the same base pairs but using different residue indices
2. **Different chain assignments**: One may assign residues to different chains
3. **PDB parsing differences**: How the PDB is parsed may differ

### When Modern Actually Misses

In cases where:
- Legacy finds pairs and modern finds NONE
- There IS overlap but modern misses some specific pairs
- The structures are parsed the same way

Then we can say modern truly "missed" pairs.

### Patterns Observed

1. **Perfect matches (26%)**: When counts match, pairs match exactly ✅
2. **Modern finds more (30%)**: Modern finds pairs legacy doesn't ✅
3. **No overlap cases (14%)**: Different interpretation of same structure ⚠️
4. **Legacy finds none (30%)**: Modern finds pairs legacy doesn't ✅

## Recommendations

1. **Investigate no-overlap cases first**
   - Compare residue indexing between legacy and modern
   - Check chain assignment differences
   - Verify PDB parsing

2. **For true misses (when there IS overlap)**
   - Check validation thresholds
   - Compare hydrogen bond detection
   - Verify frame calculations

3. **Overall assessment**
   - Modern code is more comprehensive (finds more pairs overall)
   - When pairs match, parameters match perfectly
   - Need to investigate indexing/parsing differences

## Next Steps

1. Compare residue indices in no-overlap cases to see if they represent the same pairs
2. Check if chain assignments differ
3. Investigate PDB parsing differences
4. For true misses (with overlap), investigate validation differences

