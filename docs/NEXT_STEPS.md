# Next Steps: Path to 100% Legacy Accuracy

**Created**: December 2, 2025  
**Status**: Residue Index Matching ✅ COMPLETE  
**Current Task**: Full PDB Validation 🔄 IN PROGRESS  
**Goal**: 100% accuracy with legacy indexes across all validation types

---

## What We've Completed ✅

### ✅ Residue Index Matching (100% COMPLETE)
Every residue in the modern code now correctly maps to its legacy index by matching PDB properties (chain ID, sequence number, insertion code).

**Key Achievement**: No more index mismatches - all residues align perfectly with legacy.

### ✅ LS_Fitting Validation (98.7% COMPLETE - December 3, 2025)
Comprehensive validation of ls_fitting calculations across all 3,602 fast PDBs.

**Results**:
- Perfect match: 860 PDBs (24%)
- FP differences only: 2,693 PDBs (75%)  
- Count mismatches: 47 PDBs (1.3%)
- **Success rate: 98.7%**

**Bugs Fixed**:
1. ✅ Purine detection bug (70U misclassified - used wrong template)
2. ✅ Deduplication key bug (insertion codes ignored)
3. ✅ RMSD threshold strategy (warped vs structural variants)

**Next**: Add 12 modified bases to structural_variants whitelist → 100%

---

## What's Next: Immediate Actions

### 1. Complete Current Validation Run 🔄 ACTIVE
**Status**: Currently running (Batch 3-4 of ~18)  
**What to do**: Let it finish, then analyze results

**Current Stats**:
- ✅ ~150+ PDBs passing
- ❌ 1 failure (7EH2 - count mismatch)
- ⏱️ 13 timeouts (large structures)
- ⏭️ ~1000+ skipped (no legacy JSON)

**Action**: Wait for completion, no intervention needed right now.

---

### 2. Investigate the 7EH2 Failure ⚠️ HIGH PRIORITY
**Issue**: Count mismatch detected  
**When**: After current validation completes

**Investigation Commands**:
```bash
# See detailed comparison
python3 scripts/compare_json.py compare 7EH2 --verbose

# Check specific record types
python3 scripts/compare_json.py atoms 7EH2
python3 scripts/compare_json.py frames 7EH2
python3 scripts/compare_json.py ring-atoms 7EH2

# Regenerate if needed
python3 scripts/rebuild_json.py regenerate 7EH2
```

**Goal**: Understand why counts don't match, fix the issue.

---

### 3. Handle Timeout Cases ⚠️ MEDIUM PRIORITY
**Issue**: 13 PDBs timeout (>120s processing time)

**Timeout PDBs**:
- 6ZMT, 6ZLW, 6ZN5, 6ZOJ, 5XYM, 5NRG
- 6ZV6, 6ZUO, 6ZXG, 6ZXH, 6ZXD, 6ZXE
- 7ASM

**Options**:
1. **Accept them as valid** (they're just slow, likely very large structures)
2. **Increase timeout** threshold in validation
3. **Optimize performance** (if needed for production)

**Recommendation**: Check structure sizes first, then decide if optimization is needed.

---

### 4. Validate All 10 Record Types ⏳ NEXT PHASE

Once failures are resolved, validate that ALL record types match 100%:

#### Phase 1: Core (Already Validated ✅)
- ✅ `pdb_atoms` - Atom records
- ✅ `base_frame_calc` - Frame metadata
- ✅ `frame_calc` / `ref_frame` - Reference frames

#### Phase 2: Base Pairs (Currently Testing 🔄)
- 🔄 `base_pair` - Base pair records
- 🔄 `pair_validation` - Validation results
- 🔄 `distance_checks` - Geometric checks
- 🔄 `hbond_list` - H-bond detection
- 🔄 `find_bestpair_selection` - **PRIMARY OUTPUT** (must be 100%)

#### Phase 3: Step Parameters (Pending ⏳)
- ⏳ `bpstep_params` - Step parameters
- ⏳ `helical_params` - Helical parameters

**Validation Command**:
```bash
# Test all record types on test set
python3 scripts/compare_json.py compare --test-set 100 --verbose
```

---

### 5. Run Comprehensive Test Set Validation ⏳ FINAL VALIDATION

**Goal**: Prove 100% accuracy on test set

**Test Progression**:
```bash
# Quick test (10 PDBs)
python3 scripts/compare_json.py compare --test-set 10

# Standard test (100 PDBs)
python3 scripts/compare_json.py compare --test-set 100

# Extended test (500 PDBs)
python3 scripts/compare_json.py compare --test-set 500

# Full validation (all available PDBs)
python3 scripts/compare_json.py compare
```

**Success Criteria**: 100% pass rate on all record types for all PDBs in test set.

---

## Success Milestones

### Short Term (This Week)
- [ ] Current validation run completes
- [ ] Investigate and fix 7EH2 failure
- [ ] Analyze timeout cases
- [ ] Document results in VALIDATION_PROGRESS.md

### Medium Term (Next Week)
- [ ] 100% pass rate on test_set_100
- [ ] All 10 record types validated
- [ ] Step parameters match legacy
- [ ] No unexplained failures

### Long Term (Production Ready)
- [ ] All available PDBs validated (excluding timeouts if acceptable)
- [ ] Comprehensive documentation complete
- [ ] Performance acceptable (<120s for reasonable structures)
- [ ] CI/CD validation pipeline set up

---

## Key Commands Reference

### Current Validation
```bash
# Let current run complete (already running)
# Check terminal for progress

# After completion, check results
cat data/index_validation_status.csv  # or wherever results are saved
```

### Quick Investigation
```bash
# Investigate a failure
python3 scripts/compare_json.py compare <PDB_ID> --verbose

# Regenerate modern JSON
python3 scripts/rebuild_json.py regenerate <PDB_ID> --modern-only

# Regenerate both
python3 scripts/rebuild_json.py regenerate <PDB_ID>
```

### Test Set Validation
```bash
# Quick test
python3 scripts/compare_json.py compare --test-set 10

# Standard test
python3 scripts/compare_json.py compare --test-set 100

# Full test
python3 scripts/compare_json.py compare
```

### Specific Record Type Testing
```bash
# Test specific components
python3 scripts/compare_json.py atoms <PDB_ID>
python3 scripts/compare_json.py frames <PDB_ID>
python3 scripts/compare_json.py steps <PDB_ID>
```

---

## Tracking Progress

**Progress Document**: `docs/VALIDATION_PROGRESS.md`
- Update after each major milestone
- Document all failures and fixes
- Track statistics and match rates

**Update After**:
- Current validation run completes
- Any failure is fixed
- Test set validation runs
- Any code changes that could affect matching

---

## Questions to Answer

As validation progresses, we need to answer:

1. **What causes timeouts?**
   - Structure size?
   - Algorithm complexity?
   - Optimization opportunities?

2. **What's the 7EH2 issue?**
   - Count mismatch - which count?
   - Residues? Atoms? Frames? Pairs?
   - Quick fix or fundamental issue?

3. **Do step parameters match?**
   - Not yet tested comprehensively
   - Critical for 100% accuracy claim

4. **Are we production ready?**
   - After all validations pass
   - Performance acceptable?
   - Documentation complete?

---

## Immediate Action Plan (Right Now)

### TODAY
1. ✅ **Created tracking documents** (VALIDATION_PROGRESS.md, NEXT_STEPS.md)
2. 🔄 **Let validation run complete** (in progress)
3. ⏳ **Monitor for completion** (check terminal periodically)

### TOMORROW (After Validation Completes)
1. **Analyze results**
   - Count passes, failures, timeouts
   - Update VALIDATION_PROGRESS.md
2. **Investigate 7EH2**
   - Run diagnostic commands
   - Fix if possible
3. **Test fix**
   - Regenerate JSON
   - Re-run comparison
   - Verify fix works

### THIS WEEK
1. **Resolve all failures**
2. **Run test_set_100 validation**
3. **Document results**
4. **Plan step parameter validation**

---

## Related Documents

- **[VALIDATION_PROGRESS.md](VALIDATION_PROGRESS.md)** - Detailed progress tracker
- **[TESTING_GUIDE.md](TESTING_GUIDE.md)** - Complete testing documentation
- **[JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)** - What gets compared
- **[README.md](../README.md)** - Project overview

---

**Bottom Line**: You've completed residue indexing (huge milestone! ✅). Now we're validating that EVERYTHING matches legacy across all PDBs. Current run is in progress - let it finish, then we investigate any failures and move to comprehensive test set validation.

The path is clear: finish current validation → fix failures → test sets → 100% accuracy achieved! 🎯

