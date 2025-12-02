# Step-by-Step Debugging - Implementation Complete ✅

**Date**: 2025-01-XX  
**Status**: ✅ **FULLY IMPLEMENTED AND TESTED**

---

## Summary

Complete step-by-step debugging infrastructure has been implemented and tested. Both legacy and modern code now output detailed JSON at every step, allowing precise identification of where differences occur.

**Test Results**: ✅ **100% Match** on minimal test case `1AQ4_minimal_pairs_1_2.pdb`

---

## ✅ What's Been Implemented

### 1. Minimal Test Case Extraction
- **Script**: `scripts/extract_minimal_pairs.py`
- **Status**: ✅ Working
- **Usage**: Extracts minimal PDB fragments with N consecutive base pairs

### 2. Legacy Code JSON Output
- **Functions Added**:
  - `json_writer_record_best_partner_candidates()`
  - `json_writer_record_mutual_best_decision()`
  - `json_writer_record_iteration_state()`
- **Code Modified**:
  - `best_pair()` - Records all candidates
  - `find_bestpair()` - Records iterations and mutual best decisions
  - `json_writer_record_iteration_state()` - Fixed to record all pairs per iteration
- **Status**: ✅ Working

### 3. Modern Code JSON Output
- **Methods Added** to `JsonWriter`:
  - `record_best_partner_candidates()`
  - `record_mutual_best_decision()`
  - `record_iteration_state()`
- **Code Modified**:
  - `BasePairFinder::find_best_partner()` - Records all candidates
  - `BasePairFinder::find_best_pairs()` - Records iterations and mutual best decisions
- **Status**: ✅ Working

### 4. Comparison Tools
- **Scripts Created**:
  - `scripts/compare_best_partner.py` - Compare best partner finding
  - `scripts/compare_iteration.py` - Compare iteration states
  - `scripts/compare_mutual_best.py` - Compare mutual best decisions
- **Status**: ✅ All working and tested

---

## Test Results: 1AQ4_minimal_pairs_1_2.pdb

### ✅ Best Partner Candidates
- **Best Partner**: Legacy=4, Modern=4 ✅ **MATCH**
- **Best Score**: Legacy=-1.983001, Modern=-1.983001 ✅ **MATCH**
- **Candidates**: All 4 candidates recorded in both ✅ **MATCH**

### ✅ Mutual Best Decisions
- **Total decisions**: Legacy=2, Modern=2 ✅ **MATCH**
- **All decisions match**: ✅ **MATCH**

### ✅ Iteration States
- **Iteration 1**: Both find pairs [(1,4), (2,3)] ✅ **MATCH**
- **Iteration 2**: Both show no new pairs ✅ **MATCH**
- **All iterations match**: ✅ **MATCH**

---

## Complete Workflow

### Step 1: Extract Minimal Test Case
```bash
python3 scripts/extract_minimal_pairs.py \
    --pdb data/pdb/1AQ4.pdb \
    --legacy-json data/json_legacy/find_bestpair_selection/1AQ4.json \
    --output-dir data/pdb/minimal \
    --num-pairs 2
```

### Step 2: Generate Legacy JSON
```bash
./org/build/bin/find_pair_original data/pdb/minimal/1AQ4_minimal_pairs_1_2.pdb
```

**Creates**:
- `data/json_legacy/best_partner_candidates/1AQ4_minimal_pairs_1_2.json`
- `data/json_legacy/mutual_best_decisions/1AQ4_minimal_pairs_1_2.json`
- `data/json_legacy/iteration_states/1AQ4_minimal_pairs_1_2.json`

### Step 3: Generate Modern JSON
```bash
./build/generate_modern_json data/pdb/minimal/1AQ4_minimal_pairs_1_2.pdb \
    data/json/1AQ4_minimal_pairs_1_2.json --fix-indices
```

**Creates**:
- `data/json/1AQ4_minimal_pairs_1_2.json/best_partner_candidates/1AQ4_minimal_pairs_1_2.json`
- `data/json/1AQ4_minimal_pairs_1_2.json/mutual_best_decisions/1AQ4_minimal_pairs_1_2.json`
- `data/json/1AQ4_minimal_pairs_1_2.json/iteration_states/1AQ4_minimal_pairs_1_2.json`

### Step 4: Compare Step-by-Step

**Compare best partner finding**:
```bash
python3 scripts/compare_best_partner.py 1AQ4_minimal_pairs_1_2 1 \
    --legacy-dir data/json_legacy \
    --modern-dir "data/json/1AQ4_minimal_pairs_1_2.json" \
    --verbose
```

**Compare iteration states**:
```bash
python3 scripts/compare_iteration.py 1AQ4_minimal_pairs_1_2 \
    --legacy-dir data/json_legacy \
    --modern-dir "data/json/1AQ4_minimal_pairs_1_2.json" \
    --verbose
```

**Compare mutual best decisions**:
```bash
python3 scripts/compare_mutual_best.py 1AQ4_minimal_pairs_1_2 \
    --legacy-dir data/json_legacy \
    --modern-dir "data/json/1AQ4_minimal_pairs_1_2.json" \
    --verbose
```

---

## JSON Output Structure

### Best Partner Candidates
```json
{
  "type": "best_partner_candidates",
  "res_i": 1,
  "num_candidates": 4,
  "best_partner": 4,
  "best_score": -1.983001,
  "candidates": [
    {
      "res_j": 1,
      "is_eligible": 0,
      "score": 1e18,
      "bp_type_id": 0,
      "is_best": 0
    },
    ...
  ]
}
```

### Mutual Best Decisions
```json
{
  "type": "mutual_best_decision",
  "res1": 1,
  "res2": 4,
  "best_partner_for_res1": 4,
  "best_partner_for_res2": 1,
  "is_mutual": 1,
  "was_selected": 1
}
```

### Iteration States
```json
{
  "type": "iteration_states",
  "iteration_num": 1,
  "num_matched": 4,
  "num_total": 4,
  "pairs_found_in_iteration": [[1, 4], [2, 3]],
  "matched_residues": [1, 2, 3, 4]
}
```

---

## Files Modified/Created

### Legacy Code
- `org/include/json_writer.h` - Added 3 function declarations
- `org/src/json_writer.c` - Implemented 3 functions, fixed iteration_state recording
- `org/src/find_pair.c` - Modified `best_pair()` and `find_bestpair()`, added JSON init

### Modern Code
- `include/x3dna/io/json_writer.hpp` - Added 3 method declarations
- `src/x3dna/io/json_writer.cpp` - Implemented 3 methods, added to type_to_dir map
- `src/x3dna/algorithms/base_pair_finder.cpp` - Modified to record step-by-step data

### Scripts
- `scripts/extract_minimal_pairs.py` - Extract minimal test cases
- `scripts/compare_best_partner.py` - Compare best partner finding
- `scripts/compare_iteration.py` - Compare iteration states
- `scripts/compare_mutual_best.py` - Compare mutual best decisions

---

## Benefits

1. **Precision**: Identify exact step where legacy and modern diverge
2. **Isolation**: Minimal test cases make debugging manageable
3. **Visibility**: See all intermediate states and decisions
4. **Systematic**: Debug step-by-step, not guess-and-check
5. **Verification**: Verify fixes at each step level

---

## Next Steps

With this infrastructure in place, you can now:

1. **Debug systematically**: Use minimal test cases to isolate issues
2. **Find exact divergence**: Compare at each step to see where differences occur
3. **Verify fixes**: Confirm fixes work at each step level
4. **Scale up**: Once minimal cases work, test with full PDBs

---

## Example: Debugging a Difference

When you find a difference between legacy and modern:

1. **Extract minimal case**:
   ```bash
   python3 scripts/extract_minimal_pairs.py --pdb <PDB> --num-pairs 2
   ```

2. **Generate step-by-step JSON** for both legacy and modern

3. **Compare each step**:
   ```bash
   # Check best partner finding
   python3 scripts/compare_best_partner.py <PDB> <res_i> --verbose
   
   # Check iteration states
   python3 scripts/compare_iteration.py <PDB> --verbose
   
   # Check mutual best decisions
   python3 scripts/compare_mutual_best.py <PDB> --verbose
   ```

4. **Identify divergence point**: The comparison tools will show exactly which step differs

5. **Fix and verify**: Make fix, regenerate, compare again

---

*The step-by-step debugging infrastructure is complete and ready for systematic debugging of any differences between legacy and modern code!*

