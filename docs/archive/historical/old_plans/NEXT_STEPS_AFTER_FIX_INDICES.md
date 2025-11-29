# Next Steps After --fix-indices Implementation

**Date**: 2025-01-XX  
**Status**: Ready for testing and next phase

**⚠️ IMPORTANT**: See [ACTION_PLAN_NEXT_STEPS.md](ACTION_PLAN_NEXT_STEPS.md) for the **immediate action plan** with specific commands and steps.

---

## What We've Accomplished

✅ **Residue Indexing Fix Complete**
- Root cause identified and fixed (PdbParser)
- `--fix-indices` option implemented
- PDB properties matching approach working
- Debug tools updated
- Documentation complete

✅ **Verified Results**
- Pair (1102, 1127) in 6CAQ: Correctly identified
- dorg: 1.83115 (matches legacy 1.831148)
- bp_type_id: 2 (matches legacy)

---

## Immediate Next Steps

### 1. Test --fix-indices with Full Workflow

**Goal**: Verify the option works end-to-end and improves matching

**Steps**:
```bash
# Test with 6CAQ
./build/find_pair_app --fix-indices data/pdb/6CAQ.pdb /tmp/6CAQ_fixed.inp

# Compare with legacy
python3 scripts/compare_json.py compare 6CAQ --verbose

# Or use test script
python3 scripts/test_fix_indices_option.py 6CAQ
```

**Expected Improvements**:
- More pairs should match legacy
- dorg calculations should match better
- bp_type_id assignments should match better

### 2. Address Phase 2: Pairs Not Found in Validation

**Issue**: 9 pairs in 6CAQ not found in validation

**Pairs to Investigate**:
- (495, 498)
- (501, 506)
- (939, 1378)
- (1029, 1184)
- (1380, 1473)
- (1382, 1470)
- (1385, 1467)
- (1489, 1492)

**Questions**:
- Are these pairs being skipped during validation?
- Why are they not being validated?
- Do they pass geometric checks?
- Are frames available for these residues?

**Tools to Use**:
- `build/debug_bp_type_id_step_params` - Debug specific pairs
- `build/check_residue_indices` - Verify residue mapping
- `scripts/compare_json.py` - Compare validation records

### 3. Test on Multiple PDBs

**Goal**: Measure overall improvement with --fix-indices

**Test Set**: Use 10-PDB or 100-set test

**Steps**:
```bash
# Generate JSON with --fix-indices for test set
# (Need to integrate into generate_modern_json or create wrapper)

# Compare results
python3 scripts/compare_json.py compare --test-set 10 --verbose
```

**Metrics to Track**:
- Match rate improvement
- Number of pairs now matching
- dorg calculation improvements
- bp_type_id assignment improvements

---

## Remaining Issues (After Testing)

### Phase 3: Quality Score Differences
- 4 missing pairs, 5 extra pairs in find_bestpair selection
- May be resolved by fixing residue indexing
- Need to verify after --fix-indices testing

### Phase 4: Residue Recognition
- 3KNC: Only 16/66 residues recognized
- 5UJ2: Residue 2 missing
- May need PDB parsing fixes

### Phase 5: Base Pair Validation
- 1227 missing base pairs
- 107 extra base pairs
- Many may be resolved by fixing root causes

---

## Testing Strategy

### Step 1: Single PDB Test (6CAQ)
1. Run with `--fix-indices`
2. Compare JSON output with legacy
3. Verify improvements
4. Document results

### Step 2: Small Test Set (10 PDBs)
1. Run with `--fix-indices` for all 10
2. Compare results
3. Measure overall improvement
4. Identify remaining issues

### Step 3: Full Test Set (100 PDBs)
1. Run with `--fix-indices` for all 100
2. Compare results
3. Measure final match rate
4. Document remaining differences

---

## Integration Options

### Option A: Make --fix-indices Default When Legacy JSON Available
- Auto-detect legacy JSON
- Fix indices automatically
- No user intervention needed

### Option B: Keep as Optional Flag
- User explicitly enables when needed
- Clear separation of concerns
- Easier to debug

### Option C: Hybrid Approach
- Auto-fix when legacy JSON is available AND comparing mode is enabled
- Keep as optional for production use

**Recommendation**: Option B for now, consider Option C later

---

## Success Criteria

### Short Term (After Testing)
- ✅ `--fix-indices` option works correctly
- ✅ Pair (1102, 1127) matches legacy
- ✅ dorg and bp_type_id calculations match
- ✅ Test script validates improvements

### Medium Term (After Phase 2)
- ✅ 9 missing pairs in 6CAQ resolved or explained
- ✅ Match rate improved on 10-PDB test set
- ✅ Quality score differences reduced

### Long Term (100% Match)
- ✅ 100% match on 100-set test
- ✅ All pairs match legacy
- ✅ All parameters match legacy

---

## Related Documentation

- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - Usage guide
- [RESIDUE_INDEXING_COMPLETE.md](RESIDUE_INDEXING_COMPLETE.md) - Solution details
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows

---

*Next steps after completing the --fix-indices implementation.*

