# Step-by-Step Debugging Implementation

**Created**: 2025-01-XX  
**Status**: ✅ Legacy JSON output complete, comparison tools ready

---

## Overview

We've implemented step-by-step debugging capabilities to identify exactly where legacy and modern code diverge. This allows us to break down the base pair finding process into smaller, debuggable pieces.

---

## What's Been Implemented

### ✅ Phase 1: Minimal Test Cases
- **Script**: `scripts/extract_minimal_pairs.py`
- **Purpose**: Extract minimal PDB fragments with exactly N consecutive base pairs
- **Status**: Working and tested
- **Usage**:
  ```bash
  python3 scripts/extract_minimal_pairs.py \
      --pdb data/pdb/1AQ4.pdb \
      --legacy-json data/json_legacy/find_bestpair_selection/1AQ4.json \
      --output-dir data/pdb/minimal \
      --num-pairs 2
  ```

### ✅ Phase 2: Legacy JSON Output
- **Functions Added**:
  - `json_writer_record_best_partner_candidates()` - Records all candidates considered
  - `json_writer_record_mutual_best_decision()` - Records mutual best decisions
  - `json_writer_record_iteration_state()` - Records state after each iteration
- **Code Modified**:
  - `best_pair()` - Now collects and records all candidates
  - `find_bestpair()` - Now records iteration states and mutual best decisions
  - `find_pair_main()` - Now initializes JSON writer
- **Status**: Working and tested

### ✅ Phase 3: Comparison Tools
- **Scripts Created**:
  - `scripts/compare_best_partner.py` - Compare best partner finding
  - `scripts/compare_iteration.py` - Compare iteration states
- **Status**: Ready (waiting for modern code JSON output)

---

## JSON Output Files

When legacy code runs, it creates:

### `data/json_legacy/best_partner_candidates/<PDB_ID>.json`
Records all candidates considered when finding the best partner for each residue.

**Structure**:
```json
[
  {
    "type": "best_partner_candidates",
    "res_i": 1,
    "num_candidates": 4,
    "best_partner": 4,
    "best_score": -1.983001,
    "candidates": [
      {
        "res_j": 2,
        "is_eligible": 1,
        "score": 1e+18,
        "bp_type_id": 0,
        "is_best": 0
      },
      ...
    ]
  }
]
```

### `data/json_legacy/mutual_best_decisions/<PDB_ID>.json`
Records mutual best match decisions.

**Structure**:
```json
[
  {
    "type": "mutual_best_decision",
    "res1": 1,
    "res2": 4,
    "best_partner_for_res1": 4,
    "best_partner_for_res2": 1,
    "is_mutual": 1,
    "was_selected": 1
  }
]
```

### `data/json_legacy/iteration_states/<PDB_ID>.json`
Records state after each iteration of `find_bestpair`.

**Structure**:
```json
[
  {
    "type": "iteration_state",
    "iteration_num": 1,
    "num_matched": 4,
    "num_total": 4,
    "pairs_found_in_iteration": [[1, 4], [2, 3]],
    "matched_residues": [1, 2, 3, 4]
  }
]
```

---

## Usage

### 1. Extract Minimal Test Case
```bash
python3 scripts/extract_minimal_pairs.py \
    --pdb data/pdb/1AQ4.pdb \
    --legacy-json data/json_legacy/find_bestpair_selection/1AQ4.json \
    --output-dir data/pdb/minimal \
    --num-pairs 2
```

### 2. Run Legacy Code (generates JSON)
```bash
./org/build/bin/find_pair_original data/pdb/minimal/1AQ4_minimal_pairs_1_2.pdb
```

This creates JSON files in:
- `data/json_legacy/best_partner_candidates/`
- `data/json_legacy/mutual_best_decisions/`
- `data/json_legacy/iteration_states/`

### 3. Compare Best Partner Finding
```bash
python3 scripts/compare_best_partner.py <PDB_ID> <res_i> [--verbose]
```

Example:
```bash
python3 scripts/compare_best_partner.py 1AQ4_minimal_pairs_1_2 1 --verbose
```

### 4. Compare Iterations
```bash
python3 scripts/compare_iteration.py <PDB_ID> [<iteration_num>] [--verbose]
```

Example:
```bash
python3 scripts/compare_iteration.py 1AQ4_minimal_pairs_1_2 --verbose
```

---

## Next Steps

### ⏳ Add JSON Output to Modern Code

The modern code needs to output the same JSON files so we can compare:

1. **Add JSON writer functions** to modern code (or adapt existing ones)
2. **Modify `BasePairFinder::find_best_partner()`** to record candidates
3. **Modify `BasePairFinder::find_best_pairs()`** to record iterations
4. **Test** with minimal PDBs

### ⏳ Use for Debugging

Once modern code outputs JSON:

1. Run both legacy and modern on minimal PDB
2. Use comparison tools to find exact point of divergence
3. Fix the difference
4. Verify fix with comparison tools

---

## Example Workflow

```bash
# 1. Extract minimal test case
python3 scripts/extract_minimal_pairs.py \
    --pdb data/pdb/1AQ4.pdb \
    --legacy-json data/json_legacy/find_bestpair_selection/1AQ4.json \
    --output-dir data/pdb/minimal \
    --num-pairs 2

# 2. Run legacy
./org/build/bin/find_pair_original data/pdb/minimal/1AQ4_minimal_pairs_1_2.pdb

# 3. Run modern (once JSON output is added)
./build/find_pair_app --fix-indices data/pdb/minimal/1AQ4_minimal_pairs_1_2.pdb /tmp/test.inp

# 4. Compare
python3 scripts/compare_best_partner.py 1AQ4_minimal_pairs_1_2 1 --verbose
python3 scripts/compare_iteration.py 1AQ4_minimal_pairs_1_2 --verbose
```

---

## Files Modified/Created

### Legacy Code
- `org/include/json_writer.h` - Added 3 new function declarations
- `org/src/json_writer.c` - Implemented 3 new JSON writer functions
- `org/src/find_pair.c` - Modified `best_pair()` and `find_bestpair()`, added JSON init

### Scripts
- `scripts/extract_minimal_pairs.py` - Extract minimal test cases
- `scripts/compare_best_partner.py` - Compare best partner finding
- `scripts/compare_iteration.py` - Compare iteration states

---

## Benefits

1. **Isolation**: Minimal test cases make debugging easier
2. **Visibility**: See exactly what happens at each step
3. **Comparison**: Direct comparison of legacy vs modern at each step
4. **Precision**: Identify exact point of divergence
5. **Verification**: Verify fixes work at each step level

---

*This debugging infrastructure provides the foundation for systematically identifying and fixing differences between legacy and modern code.*

