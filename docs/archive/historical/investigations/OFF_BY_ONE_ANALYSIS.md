# Off-by-One Indexing Error Analysis

**Date**: 2025-01-XX  
**Status**: ✅ CONFIRMED - Legacy index 1102 = Modern index 1101

---

## Problem

**Legacy pair (1102, 1127)**: bp_type = "GC"  
**Modern pair (1102, 1127)**: bp_type = "CC" ❌  
**Modern pair (1101, 1127)**: bp_type = "GC" ✅

---

## Evidence

### Residue Identification

| Index | Legacy Expects | Modern Has | Match? |
|-------|---------------|------------|--------|
| 1101 | ? | G (Chain A, Seq 1127) | - |
| 1102 | G | C (Chain A, Seq 1128) | ❌ |
| 1103 | ? | C (Chain A, Seq 1129) | - |
| 1127 | C | C (Chain A, Seq 1153) | ✅ |

### Pair Comparison

| Pair | Legacy | Modern | Match? |
|------|--------|--------|--------|
| (1102, 1127) | GC | CC | ❌ |
| (1101, 1127) | ? | GC | ✅ (if legacy 1102 = modern 1101) |

---

## Root Cause

**OFF-BY-ONE ERROR**: Legacy's residue indices are **+1** compared to modern's indices.

- **Legacy index 1102** = **Modern index 1101** (G)
- **Legacy index 1127** = **Modern index 1127** (C) - matches!

---

## Validation Results

### Modern Pair (1101, 1127) - Correct Residues
- bp_type: GC ✅ (matches legacy)
- dorg: 19.7203 ❌ (legacy: 1.831148)
- d_v: 15.5201 ❌ (legacy: ~0.5)
- is_valid: NO ❌

### Modern Pair (1102, 1127) - Wrong Residues
- bp_type: CC ❌ (legacy: GC)
- dorg: 24.8669 ❌ (legacy: 1.831148)
- d_v: 24.8296 ❌
- is_valid: NO ❌

---

## Next Steps

1. **Fix off-by-one error**: When legacy uses index `i`, modern should look up index `i-1`
   - **BUT**: This doesn't make sense if both use 1-based indexing!
   - Need to check if there's a residue counting difference

2. **Investigate dorg mismatch**: Even with correct residues (1101, 1127), dorg = 19.72 vs 1.83
   - This suggests a frame calculation issue beyond just indexing

3. **Check residue counting**: Compare how legacy and modern count residues
   - Are they counting the same residues?
   - Is there a residue that one counts but the other doesn't?

---

## Possible Causes

### 1. Residue Counting Difference

**Hypothesis**: Modern counts one extra residue at the beginning (or legacy skips one).

**Check**: Compare total residue counts between legacy and modern.

### 2. Index Assignment Bug

**Hypothesis**: `legacy_residue_idx` is assigned incorrectly during PDB parsing.

**Check**: Verify that indices start at 1 and increment correctly.

### 3. Lookup Offset

**Hypothesis**: When looking up residues, we're using the wrong index.

**Check**: Verify that `residue_by_legacy_idx.find(1102)` should actually look up 1101.

---

## Related Documentation

- [RESIDUE_INDEXING_ISSUE.md](RESIDUE_INDEXING_ISSUE.md) - Original analysis
- [PAIR_REJECTION_ANALYSIS.md](PAIR_REJECTION_ANALYSIS.md) - Rejection analysis
- [DORG_CALCULATION_INVESTIGATION.md](DORG_CALCULATION_INVESTIGATION.md) - dorg investigation

---

*This document tracks the off-by-one indexing error that causes residue mismatches.*

