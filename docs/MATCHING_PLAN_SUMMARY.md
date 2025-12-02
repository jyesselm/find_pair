# Matching Plan - Implementation Complete ✅

**Date**: 2025-01-XX  
**Status**: ✅ **ALL PHASES COMPLETE**

---

## Summary

The complete step-by-step debugging infrastructure has been implemented and tested. You can now systematically debug any differences between legacy and modern code by:

1. **Extracting minimal test cases** with just 2 consecutive base pairs
2. **Generating step-by-step JSON** showing every decision point
3. **Comparing at each step** to find exact divergence points
4. **Fixing systematically** based on precise identification

---

## ✅ What's Been Implemented

### 1. Minimal Test Case Extraction ✅
**Script**: `scripts/extract_minimal_pairs.py`
- Extracts minimal PDB fragments with N consecutive base pairs
- Uses legacy JSON to identify which pairs to extract
- Automatically maps residue indices

### 2. Legacy Code JSON Output ✅
**Functions Added**:
- `json_writer_record_best_partner_candidates()` - All candidates considered
- `json_writer_record_mutual_best_decision()` - Mutual best decisions
- `json_writer_record_iteration_state()` - State after each iteration

**Code Modified**:
- `best_pair()` - Collects and records all candidates
- `find_bestpair()` - Records iterations and mutual best decisions
- Fixed iteration_state to record all pairs per iteration

### 3. Modern Code JSON Output ✅
**Methods Added**:
- `JsonWriter::record_best_partner_candidates()`
- `JsonWriter::record_mutual_best_decision()`
- `JsonWriter::record_iteration_state()`

**Code Modified**:
- `BasePairFinder::find_best_partner()` - Records all candidates
- `BasePairFinder::find_best_pairs()` - Records iterations and mutual best decisions

### 4. Comparison Tools ✅
**Scripts Created**:
- `scripts/compare_best_partner.py` - Compare best partner finding
- `scripts/compare_iteration.py` - Compare iteration states
- `scripts/compare_mutual_best.py` - Compare mutual best decisions
- `scripts/debug_minimal_case.sh` - Automated workflow

---

## Quick Start

### Automated Workflow (Easiest)

```bash
# Debug a minimal case with 2 consecutive base pairs
./scripts/debug_minimal_case.sh <PDB_ID> 2
```

**Example**:
```bash
./scripts/debug_minimal_case.sh 1AQ4 2
```

**Output**:
- Extracts minimal test case
- Generates legacy and modern JSON
- Compares all step-by-step results
- Shows summary of matches/differences

### Manual Workflow (More Control)

See [DEBUGGING_WORKFLOW.md](DEBUGGING_WORKFLOW.md) for detailed manual steps.

---

## Test Results

### Test Case 1: 1AQ4_minimal_pairs_1_2.pdb
- ✅ Best Partner: MATCH
- ✅ Best Score: MATCH
- ✅ Iteration States: MATCH
- ✅ Mutual Best Decisions: MATCH
- **Result**: 100% match ✅

### Test Case 2: 1H4S_minimal_pairs_1_2.pdb
- ✅ Best Partner: MATCH
- ✅ Best Score: MATCH
- ✅ Iteration States: MATCH
- ✅ Mutual Best Decisions: MATCH
- **Result**: 100% match ✅

---

## How to Use for Debugging

### When You Find a Difference in Full PDB:

1. **Extract minimal case**:
   ```bash
   python3 scripts/extract_minimal_pairs.py \
       --pdb data/pdb/<PDB_WITH_DIFF>.pdb \
       --legacy-json data/json_legacy/find_bestpair_selection/<PDB_WITH_DIFF>.json \
       --output-dir data/pdb/minimal \
       --num-pairs 2
   ```

2. **Run automated debugging**:
   ```bash
   ./scripts/debug_minimal_case.sh <PDB_WITH_DIFF> 2
   ```

3. **Analyze differences**:
   - If iterations differ → Check which iteration first shows difference
   - If mutual best differs → Check which pairs have different decisions
   - If candidates differ → Check which candidates are considered differently

4. **Fix and verify**:
   - Make fix
   - Re-run debugging workflow
   - Verify fix resolves difference

---

## Files Created

### Scripts
- `scripts/extract_minimal_pairs.py` - Extract minimal test cases
- `scripts/compare_best_partner.py` - Compare best partner finding
- `scripts/compare_iteration.py` - Compare iteration states
- `scripts/compare_mutual_best.py` - Compare mutual best decisions
- `scripts/debug_minimal_case.sh` - Automated workflow

### Documentation
- `docs/MATCHING_PLAN.md` - This plan (complete)
- `docs/STEP_BY_STEP_COMPLETE.md` - Implementation details
- `docs/DEBUGGING_WORKFLOW.md` - Usage guide

---

## Key Benefits

1. **Precision**: Identify exact step where legacy and modern diverge
2. **Isolation**: Minimal test cases (2 pairs) make debugging manageable
3. **Visibility**: See all intermediate states and decisions
4. **Systematic**: Debug step-by-step, not guess-and-check
5. **Verification**: Verify fixes at each step level

---

## Next: Use for Actual Debugging

The infrastructure is ready. When you encounter a PDB with differences:

1. Use `./scripts/debug_minimal_case.sh <PDB_ID> 2` to extract and compare
2. Review the comparison output to see where differences occur
3. Use the step-by-step JSON to understand why differences happen
4. Make targeted fixes based on precise identification
5. Re-test to verify fix

---

*The step-by-step debugging infrastructure is complete and ready for systematic debugging of any differences between legacy and modern code!*

