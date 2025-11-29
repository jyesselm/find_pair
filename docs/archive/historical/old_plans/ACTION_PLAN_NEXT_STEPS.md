# Action Plan: Next Steps

**Last Updated**: 2025-01-XX  
**Status**: Ready to proceed with testing and Phase 2 investigation

---

## ✅ What's Been Completed

### Residue Indexing Solution (COMPLETE)
- ✅ Root cause identified and fixed (PdbParser grouping)
- ✅ `--fix-indices` option implemented
- ✅ PDB properties matching approach working
- ✅ All tools built and tested
- ✅ Comprehensive documentation

### Verified Results
- ✅ Pair (1102, 1127) in 6CAQ: Correctly identified
- ✅ dorg: 1.83115 (matches legacy 1.831148)
- ✅ bp_type_id: 2 (matches legacy)

---

## 🎯 Immediate Next Steps

### Step 1: Test --fix-indices with Full Workflow

**Goal**: Verify the option works end-to-end and improves matching

**Commands**:
```bash
# Test with 6CAQ
./build/find_pair_app --fix-indices data/pdb/6CAQ.pdb /tmp/6CAQ_fixed.inp

# Or use test script
python3 scripts/test_fix_indices_option.py 6CAQ

# Compare results
python3 scripts/compare_json.py compare 6CAQ --verbose
```

**Expected Outcomes**:
- More pairs should match legacy
- dorg calculations should match better
- bp_type_id assignments should match better
- Overall match rate should improve

**Time Estimate**: 30 minutes

---

### Step 2: Investigate 9 Missing Pairs in 6CAQ

**Goal**: Understand why these pairs are not found in validation

**Pairs to Investigate**:
- (495, 498)
- (501, 506)
- (939, 1378)
- (1029, 1184)
- (1380, 1473)
- (1382, 1470)
- (1385, 1467)
- (1489, 1492)

**Tool to Use**:
```bash
# Investigate each pair
./build/investigate_missing_pairs data/pdb/6CAQ.pdb 495 498 data/json_legacy/base_frame_calc/6CAQ.json
./build/investigate_missing_pairs data/pdb/6CAQ.pdb 501 506 data/json_legacy/base_frame_calc/6CAQ.json
# ... repeat for all pairs
```

**What to Check**:
1. Are both residues recognized as nucleotides?
2. Do both residues have frames?
3. Do pairs pass validation when tested directly?
4. Are pairs being skipped during Phase 1 iteration?

**Documentation**: Record findings in `docs/6CAQ_MISSING_PAIRS_INVESTIGATION.md`

**Time Estimate**: 1-2 hours

---

### Step 3: Test on Multiple PDBs

**Goal**: Measure overall improvement with --fix-indices

**Commands**:
```bash
# Generate JSON with --fix-indices for test set
# (May need to integrate into generate_modern_json or create wrapper)

# Compare results
python3 scripts/compare_json.py compare --test-set 10 --verbose
```

**Metrics to Track**:
- Match rate improvement
- Number of pairs now matching
- dorg calculation improvements
- bp_type_id assignment improvements

**Time Estimate**: 1 hour

---

## 📋 Priority Order

### Priority 1: Test --fix-indices (HIGH)
**Why**: Verify the solution works and measure improvements  
**Effort**: Low (30 min)  
**Impact**: High (validates entire solution)

### Priority 2: Investigate Missing Pairs (HIGH)
**Why**: These are blocking 100% match on 6CAQ  
**Effort**: Medium (1-2 hours)  
**Impact**: High (resolves 9 missing pairs)

### Priority 3: Test on Multiple PDBs (MEDIUM)
**Why**: Measure overall improvement  
**Effort**: Low (1 hour)  
**Impact**: Medium (validates across structures)

---

## 🔧 Tools Available

### 1. investigate_missing_pairs
**Purpose**: Debug why specific pairs are missing from validation  
**Usage**:
```bash
./build/investigate_missing_pairs <pdb_file> <idx1> <idx2> [legacy_json]
```

### 2. test_fix_indices_option.py
**Purpose**: Test --fix-indices option with full workflow  
**Usage**:
```bash
python3 scripts/test_fix_indices_option.py <pdb_id> [legacy_json_file]
```

### 3. debug_bp_type_id_step_params
**Purpose**: Debug bp_type_id calculation (auto-fixes indices)  
**Usage**:
```bash
./build/debug_bp_type_id_step_params <pdb_file> <idx1> <idx2> <pdb_id> [legacy_json]
```

### 4. compare_json.py
**Purpose**: Compare modern and legacy JSON outputs  
**Usage**:
```bash
python3 scripts/compare_json.py compare <pdb_id> --verbose
```

---

## 📝 Documentation to Create/Update

### After Step 1 (Test --fix-indices)
- [ ] Create `docs/TEST_RESULTS_FIX_INDICES.md` with results
- [ ] Update `100_PERCENT_MATCH_STATUS.md` with improvements

### After Step 2 (Investigate Missing Pairs)
- [ ] Create `docs/6CAQ_MISSING_PAIRS_INVESTIGATION.md` with findings
- [ ] Document root causes and fixes
- [ ] Update `100_PERCENT_MATCH_STATUS.md` with resolution

### After Step 3 (Test on Multiple PDBs)
- [ ] Create `docs/FIX_INDICES_IMPROVEMENTS.md` with metrics
- [ ] Update overall match rate statistics

---

## 🎯 Success Criteria

### Short Term (After Step 1)
- ✅ `--fix-indices` option works correctly
- ✅ Pair (1102, 1127) matches legacy
- ✅ Test script validates improvements

### Medium Term (After Step 2)
- ✅ 9 missing pairs in 6CAQ resolved or explained
- ✅ Root causes documented
- ✅ Fixes implemented (if needed)

### Long Term (After Step 3)
- ✅ Match rate improved on 10-PDB test set
- ✅ Overall improvements measured and documented
- ✅ Ready for Phase 3 (base_pair record generation)

---

## 🚀 Quick Start

**Want to get started immediately?**

1. **Test --fix-indices**:
   ```bash
   python3 scripts/test_fix_indices_option.py 6CAQ
   ```

2. **Investigate first missing pair**:
   ```bash
   ./build/investigate_missing_pairs data/pdb/6CAQ.pdb 495 498 data/json_legacy/base_frame_calc/6CAQ.json
   ```

3. **Compare results**:
   ```bash
   python3 scripts/compare_json.py compare 6CAQ --verbose
   ```

---

## Related Documentation

- [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md) - Detailed usage guide
- [RESIDUE_INDEXING_COMPLETE.md](RESIDUE_INDEXING_COMPLETE.md) - Solution details
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall status
- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing workflows
- [NEXT_STEPS_AFTER_FIX_INDICES.md](NEXT_STEPS_AFTER_FIX_INDICES.md) - Detailed next steps

---

*Action plan for immediate next steps after completing --fix-indices implementation.*

