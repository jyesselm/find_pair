# Residue Indexing Issue: Pair (1102, 1127) in 6CAQ

**Date**: 2025-01-XX  
**Status**: Root cause identified - residue indexing mismatch

---

## Problem Summary

**Pair**: (1102, 1127) in 6CAQ  
**Legacy Status**: ✅ Selected, bp_type = "GC"  
**Modern Status**: ❌ Rejected, both residues identified as "C"

**Key Mismatch**: Legacy shows bp_type "GC" but modern shows both residues as "C"

---

## Evidence

### Legacy base_pair JSON
- bp_type: "GC" (one G, one C)
- org_i: [231.855856, 166.560636, 11.292042]
- org_j: [230.663151, 165.173878, 11.20568]
- dorg: 1.831148

### Modern Validation
- Residue 1102: C (Chain A, Seq 1128)
- Residue 1127: C (Chain A, Seq 1153)
- Both have C-specific atoms (O2, N4) ✓
- Frame origins: [242.602, 170.023, 14.3657] and [219.334, 165.341, 6.94874]
- dorg: 24.8669

---

## Root Cause

**OFF-BY-ONE INDEXING ERROR**

Legacy's residue indices are **offset by +1** compared to modern's indices.

### Evidence

1. **Legacy pair (1102, 1127)**: bp_type = "GC"
2. **Modern pair (1102, 1127)**: Both residues are "C" ❌
3. **Modern pair (1101, 1127)**: bp_type = "GC" ✅

**Conclusion**: Legacy's index 1102 = Modern's index 1101

### Why This Happens

**OFF-BY-ONE ERROR**: When looking up residues, modern is using indices that are **+1** compared to what legacy uses.

- Legacy's residue 1102 = Modern's residue 1101 (G)
- Legacy's residue 1127 = Modern's residue 1127 (C) - matches!

This suggests a systematic +1 offset in how residues are indexed or looked up.

---

## Possible Causes

### 1. Residue Ordering Difference

**Hypothesis**: Legacy and modern count residues in different orders.

**Check**: Compare residue ordering between legacy and modern for 6CAQ.

### 2. Residue Index Assignment Bug

**Hypothesis**: `legacy_residue_idx` is not assigned correctly during PDB parsing.

**Check**: Verify that `atom.set_legacy_residue_idx()` is called correctly.

### 3. Residue Grouping Difference

**Hypothesis**: Legacy and modern group atoms into residues differently.

**Check**: Compare how residues are grouped by (ResName, ChainID, ResSeq, insertion).

---

## Investigation Steps

### 1. Compare Residue Ordering

```bash
# Generate modern residue ordering
./build/generate_residue_ordering_json data/pdb/6CAQ.pdb data/residue_ordering/6CAQ_modern.json

# Generate legacy residue ordering
./org/build/bin/generate_residue_ordering_json data/pdb/6CAQ.pdb data/residue_ordering/6CAQ_legacy.json

# Compare
./build/compare_residue_ordering data/residue_ordering/6CAQ_modern.json data/residue_ordering/6CAQ_legacy.json
```

### 2. Find Legacy's Actual Residues

Check what residues legacy actually uses for indices 1102 and 1127:
- Look at legacy's base_frame_calc JSON
- Check legacy's pdb_atoms JSON
- Compare with modern's residue identification

### 3. Verify Residue Index Assignment

Check if `legacy_residue_idx` is assigned correctly:
- Verify PDB parser assigns indices sequentially
- Check that all residues (including non-nucleotides) are counted
- Verify indices match legacy's counting

---

## Expected Fix

Once the off-by-one error is fixed:
1. When legacy uses index 1102, modern should look up index **1101** (not 1102)
2. bp_type should be "GC" (matches legacy) ✅ Already confirmed with (1101, 1127)
3. Frame origins should match legacy (still investigating dorg = 19.72 vs 1.83)
4. dorg should be ~1.83 (currently 19.72 for pair 1101, 1127)
5. Pair should pass validation

**Note**: Even with correct residues (1101, 1127), dorg is still 19.72 vs legacy's 1.83, suggesting an additional frame calculation issue.

---

## Related Documentation

- [RESIDUE_ORDERING_COMPARISON.md](RESIDUE_ORDERING_COMPARISON.md) - Residue ordering guide
- [PAIR_REJECTION_ANALYSIS.md](PAIR_REJECTION_ANALYSIS.md) - Rejection analysis
- [DORG_CALCULATION_INVESTIGATION.md](DORG_CALCULATION_INVESTIGATION.md) - dorg investigation

---

*This document tracks the residue indexing mismatch that causes pair (1102, 1127) to be rejected.*

