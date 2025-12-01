# Critical Finding: Residue Index Mismatch in 1TTT

**Date**: 2025-01-XX  
**Status**: ⚠️ **CRITICAL** - Modern JSON not generated with --fix-indices

---

## Summary

**Root Cause**: Modern JSON for 1TTT was **NOT generated with `--fix-indices`**, causing residue indices to not match legacy. This means we're comparing **completely different pairs**, not the same residues!

---

## Evidence

### Pair (162, 177) Investigation

**Legacy** (from `base_frame_calc`):
- Residue 162: **2MG** (modified guanine), base_type: g, chain F, seq 10
- Residue 177: **C** (cytosine), base_type: C, chain F, seq 25

**Modern** (from `base_frame_calc`):
- Residue 162: **C** (cytosine), base_type: C, chain F, seq 11
- Residue 177: **M2G** (modified guanine), base_type: ?, chain F, seq 26

**The residues are COMPLETELY DIFFERENT!**

### Impact

1. **Pair (162, 177)**: 
   - Legacy: Comparing 2MG (seq 10) with C (seq 25) ✅ Validates
   - Modern: Comparing C (seq 11) with M2G (seq 26) ❌ Fails validation
   - **These are different pairs!**

2. **All other differences**: Likely also due to index mismatches

---

## Solution

### Regenerate Modern JSON with --fix-indices

```bash
# Step 1: Regenerate modern JSON with legacy indices
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb /tmp/1TTT.inp

# Step 2: Verify indices match
python3 scripts/compare_json.py compare 1TTT --record-type residue_indices

# Step 3: Re-compare pairs
python3 scripts/compare_json.py compare 1TTT --verbose
```

### Expected Result

After regenerating with `--fix-indices`:
- Residue indices should match legacy exactly
- Pair (162, 177) should compare the same residues in both implementations
- Validation differences should be real differences, not index mismatches

---

## Verification Checklist

Before investigating any differences, verify:

- [ ] Modern JSON was generated with `--fix-indices`
- [ ] Residue indices match between legacy and modern `base_frame_calc`
- [ ] Pairs in `find_bestpair_selection` refer to same residues
- [ ] Validation records compare same residue pairs

---

## Action Items

1. **IMMEDIATE**: Regenerate modern JSON for all mismatched PDBs with `--fix-indices`:
   - 1TTT
   - 9CF3
   - 1TN1, 1TN2
   - 3F2T
   - 5V0O

2. **Verify**: Check that residue indices match after regeneration

3. **Re-investigate**: After fixing indices, re-investigate differences - they may disappear or be different

---

## Related Documentation

- [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md) - Complete guide on using legacy indices
- [1TTT_INVESTIGATION.md](1TTT_INVESTIGATION.md) - Detailed 1TTT investigation
- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - --fix-indices option documentation

---

*This finding explains why we're seeing validation differences - we're not comparing the same pairs!*

