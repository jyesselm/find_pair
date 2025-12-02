# Plan to Match Legacy and Modern Results Exactly

**Created**: 2025-01-XX  
**Purpose**: Systematic plan to achieve exact matching between legacy and modern code results

---

## ✅ MAJOR MILESTONE: RMSD-Based Residue Recognition Implemented

**Status**: RESOLVED - Modern code now matches legacy for residue recognition  
**Date**: Implementation complete

### Achievement

- ✅ Implemented RMSD-based residue recognition exactly matching legacy behavior
- ✅ Fixed 1TTT residue 16 issue (extra pair (16, 59) now correctly rejected)
- ✅ Perfect match: 230 residues, 92 pairs (matches legacy exactly)

### Documentation

- **Implementation details**: See `docs/RMSD_RESIDUE_RECOGNITION.md`
- **Specific issue**: See `docs/DEBUG_1TTT_RESIDUE_16.md` for complete resolution

### Key Implementation

- Added `check_nt_type_by_rmsd()` function using standard ring geometry
- Integrated RMSD check in `BaseFrameCalculator::calculate_frame_impl()`
- Rejects residues with RMSD > 0.2618 Å (NT_CUTOFF)
- Applies to all modified nucleotides not in NT_LIST

---

## Overview

This plan focuses on three main strategies:
1. **Create minimal test cases** - Extract small PDB fragments with just 2 consecutive base pairs
2. **Break legacy code into smaller testable units** - Isolate each step for independent testing
3. **Subdivide steps into finer-grained operations** - Make each calculation step independently debuggable

---

## Strategy 1: Create Minimal Test Cases with 2 Base Pairs

### Goal
Create minimal PDB files containing exactly 2 consecutive base pairs to isolate issues and simplify debugging.

### Approach

#### Step 1.1: Extract Base Pair Information from Legacy Output

**Tool**: Create a script that analyzes legacy JSON to identify consecutive base pairs.

```bash
# New script: scripts/extract_minimal_pairs.py
# - Reads legacy find_bestpair_selection JSON
# - Identifies consecutive base pairs (i, i+1 in base pair order)
# - Extracts residue ranges needed from PDB
# - Generates minimal PDB files with just those residues
```

**Output**: 
- Minimal PDB files like `data/pdb/minimal/1H4S_pair_1_2.pdb`
- Each contains exactly 2 consecutive base pairs plus necessary surrounding atoms

#### Step 1.2: Use Legacy Output to Determine Residue Ranges

**How it works**:
1. Load `data/json_legacy/find_bestpair_selection/<PDB_ID>.json`
2. For each consecutive pair: `(res1_i, res1_j)` and `(res2_i, res2_j)`
3. Extract all residues between `min(res1_i, res1_j)` and `max(res2_i, res2_j)` from PDB
4. Include all atoms for those residues
5. Save as minimal PDB

**Example**:
```python
# Legacy output shows:
# Pair 1: (162, 177)
# Pair 2: (163, 176)  # consecutive in base pair order

# Extract residues 162-177 from PDB
# Create minimal PDB with just these residues
```

#### Step 1.3: Generate Test Cases

**Tool**: `scripts/create_minimal_pair_test.py`

```bash
# Create minimal test cases from a PDB
python3 scripts/create_minimal_pair_test.py \
    --pdb data/pdb/1H4S.pdb \
    --legacy-json data/json_legacy/find_bestpair_selection/1H4S.json \
    --output-dir data/pdb/minimal \
    --max-pairs 2
```

**Output**:
- `1H4S_minimal_pair_1_2.pdb` - Contains base pairs 1 and 2
- `1H4S_minimal_pair_2_3.pdb` - Contains base pairs 2 and 3
- etc.

#### Step 1.4: Test Minimal Cases

```bash
# Test each minimal case
for min_pdb in data/pdb/minimal/*.pdb; do
    # Generate legacy output
    cd org && ./build/bin/find_pair_original ../$min_pdb
    
    # Generate modern output
    ./build/find_pair_app --fix-indices $min_pdb /tmp/test.inp
    
    # Compare
    python3 scripts/compare_json.py compare $(basename $min_pdb .pdb)
done
```

**Benefits**:
- Smaller files = faster iteration
- Easier to debug (only 2 pairs to track)
- Can identify which specific pair combinations cause issues
- Can verify each step of the algorithm independently

---

## Strategy 2: Break Legacy Code into Smaller Testable Units

### Goal
Isolate each major step of the legacy algorithm so it can be tested independently.

### Current Legacy Code Structure

The main function `find_bestpair()` in `org/src/find_pair.c` (line 1580) performs:

1. **Initialization** - Set up arrays and structures
2. **Best pair finding loop** - Iterative mutual best match
   - Calls `best_pair()` for each residue
   - Checks mutual best match
   - Records selected pairs
3. **Pair validation** - Inside `best_pair()` (calls `check_pair()`)

### Breakdown Plan

#### Step 2.1: Create Standalone Test Functions

**File**: `org/src/test_find_bestpair_steps.c`

Extract each major step into a testable function:

```c
// Step 1: Test pair validation only
int test_validate_pair(int i, int j, ...);

// Step 2: Test best partner finding (single iteration)
int test_find_best_partner(int i, ...);

// Step 3: Test mutual best match check
int test_check_mutual_best(int i, int j, ...);

// Step 4: Test single iteration of find_bestpair
int test_single_iteration(...);

// Step 5: Test full find_bestpair with debug output
int test_find_bestpair_full(...);
```

#### Step 2.2: Add JSON Output for Each Step

Modify legacy code to output intermediate results:

```c
// In best_pair() - output candidate pairs
json_writer_record_best_pair_candidate(i, best_j, score, is_valid);

// In mutual best check - output decision
json_writer_record_mutual_best_decision(i, j, is_mutual);

// After each iteration - output current state
json_writer_record_iteration_state(iteration_num, num_matched, pairs);
```

**Output files**:
- `data/json_legacy/best_pair_candidates/<PDB_ID>.json` - All candidates considered
- `data/json_legacy/mutual_best_decisions/<PDB_ID>.json` - Mutual best decisions
- `data/json_legacy/iteration_states/<PDB_ID>.json` - State after each iteration

#### Step 2.3: Compare Step-by-Step

Create comparison tools for each step:

```bash
# Compare pair validation
./build/compare_pair_validation <PDB_ID> <res1> <res2>

# Compare best partner finding
./build/compare_best_partner <PDB_ID> <res1>

# Compare mutual best check
./build/compare_mutual_best <PDB_ID> <res1> <res2>

# Compare full iteration
./build/compare_iteration <PDB_ID> <iteration_num>
```

**Benefits**:
- Identify exactly which step differs
- Debug each step in isolation
- Verify fixes at each step before proceeding

---

## Strategy 3: Subdivide Steps into Finer-Grained Operations

### Goal
Break down each major step into the smallest possible sub-operations, each independently testable.

### Current Steps Breakdown

#### Step 3.1: Pair Validation Breakdown

`check_pair()` performs multiple checks. Break down:

```c
// Sub-step 1: Distance check
int check_distance(int i, int j, ...);

// Sub-step 2: Hydrogen bond detection
int check_hbonds(int i, int j, ...);

// Sub-step 3: Coplanarity check
int check_coplanarity(int i, int j, ...);

// Sub-step 4: Quality score calculation
double calculate_quality_score(int i, int j, ...);

// Sub-step 5: Final validation decision
int make_validation_decision(int i, int j, ...);
```

**JSON output for each**:
```json
{
  "type": "pair_validation_substeps",
  "res1": 162,
  "res2": 177,
  "distance_check": { "passed": true, "distance": 10.5, "threshold": 20.0 },
  "hbond_detection": { "num_hbonds": 2, "hbonds": [...] },
  "coplanarity_check": { "passed": true, "angle": 12.3 },
  "quality_score": 15.6,
  "final_decision": "valid"
}
```

#### Step 3.2: Best Partner Finding Breakdown

Break `best_pair()` into sub-steps:

```c
// Sub-step 1: Initialize candidate list
void initialize_candidates(int i, ...);

// Sub-step 2: Iterate through possible partners
for (each candidate j) {
    // Sub-step 2a: Check if j is eligible
    if (!is_eligible_partner(i, j, ...)) continue;
    
    // Sub-step 2b: Get validation result (from Phase 1)
    ValidationResult result = get_validation_result(i, j);
    
    // Sub-step 2c: Calculate/retrieve quality score
    double score = result.quality_score;
    
    // Sub-step 2d: Update best candidate
    if (score < best_score) {
        best_score = score;
        best_j = j;
    }
}

// Sub-step 3: Return best partner
return best_j;
```

**JSON output**:
```json
{
  "type": "best_partner_finding",
  "res": 162,
  "candidates_considered": [
    { "res_j": 177, "eligible": true, "score": 15.6, "is_best": true },
    { "res_j": 178, "eligible": false, "reason": "already_matched" },
    ...
  ],
  "best_partner": 177,
  "best_score": 15.6
}
```

#### Step 3.3: Mutual Best Check Breakdown

Break mutual best check into sub-steps:

```c
// Sub-step 1: Find best partner for i
int best_j = find_best_partner(i, ...);

// Sub-step 2: Find best partner for j
int best_i = find_best_partner(j, ...);

// Sub-step 3: Check if mutual
bool is_mutual = (best_j == j && best_i == i);

// Sub-step 4: Make decision
if (is_mutual) {
    select_pair(i, j);
}
```

**JSON output**:
```json
{
  "type": "mutual_best_check",
  "res1": 162,
  "res2": 177,
  "best_partner_for_res1": 177,
  "best_partner_for_res2": 162,
  "is_mutual": true,
  "decision": "select"
}
```

---

## Implementation Plan

### Phase 1: Create Minimal Test Cases (Priority: High)

**Tasks**:
1. ✅ Create `scripts/extract_minimal_pairs.py`
   - Parse legacy JSON to find consecutive pairs
   - Extract residue ranges from PDB
   - Generate minimal PDB files

2. ✅ Create `scripts/create_minimal_pair_test.py`
   - User-friendly script to create test cases
   - Supports batch creation from multiple PDBs

3. ✅ Test with 2-3 PDBs initially
   - Verify minimal PDBs work correctly
   - Ensure legacy/modern both process them

**Estimated time**: 2-3 hours

### Phase 2: Add JSON Output to Legacy Code (Priority: High)

**Tasks**:
1. ✅ Add JSON output in `best_pair()`
   - Record all candidate pairs considered
   - Record best partner found

2. ✅ Add JSON output in mutual best check
   - Record decisions made

3. ✅ Add JSON output for iterations
   - Record state after each iteration

4. ✅ Create comparison scripts
   - Compare candidate lists
   - Compare best partner decisions
   - Compare iteration states

**Estimated time**: 3-4 hours

### Phase 3: Break Down Steps into Sub-operations (Priority: Medium)

**Tasks**:
1. ✅ Add JSON output for validation sub-steps
   - Distance check
   - H-bond detection
   - Coplanarity
   - Quality score

2. ✅ Add JSON output for best partner sub-steps
   - Eligibility checks
   - Score retrieval
   - Best candidate tracking

3. ✅ Create detailed comparison tools
   - Compare each sub-step independently
   - Identify exact point of divergence

**Estimated time**: 4-5 hours

### Phase 4: Test and Debug (Priority: High)

**Tasks**:
1. ✅ Run minimal test cases
   - Generate legacy and modern JSON
   - Compare at each level (sub-step, step, full)

2. ✅ Identify differences
   - Use comparison tools to find mismatches
   - Document each difference found

3. ✅ Fix differences
   - One at a time, starting with sub-steps
   - Verify fix doesn't break other cases

4. ✅ Validate with full PDBs
   - Once minimal cases pass, test full PDBs

**Estimated time**: Ongoing, iterative

---

## Tools to Create

### 1. `scripts/extract_minimal_pairs.py`

```python
"""
Extract minimal PDB fragments containing exactly N consecutive base pairs.

Usage:
    python3 scripts/extract_minimal_pairs.py \
        --pdb data/pdb/1H4S.pdb \
        --legacy-json data/json_legacy/find_bestpair_selection/1H4S.json \
        --output-dir data/pdb/minimal \
        --num-pairs 2
"""
```

**Features**:
- Read legacy find_bestpair_selection JSON
- Identify consecutive base pairs
- Extract residue ranges
- Generate minimal PDB files

### 2. `tools/compare_pair_validation.cpp`

```bash
# Compare validation results for a single pair
./build/compare_pair_validation <PDB_ID> <res1> <res2> --verbose
```

**Output**:
- Shows each sub-step (distance, hbonds, coplanarity, score)
- Highlights differences
- Shows intermediate values

### 3. `tools/compare_best_partner.cpp`

```bash
# Compare best partner finding for a single residue
./build/compare_best_partner <PDB_ID> <res1> --verbose
```

**Output**:
- Shows all candidates considered
- Shows eligibility decisions
- Shows scores for each candidate
- Highlights which candidate was selected and why

### 4. `tools/compare_iteration.cpp`

```bash
# Compare a single iteration of find_bestpair
./build/compare_iteration <PDB_ID> <iteration_num> --verbose
```

**Output**:
- Shows pairs found in this iteration
- Shows current matched state
- Compares legacy vs modern

---

## Workflow for Debugging

### When a Difference is Found:

1. **Identify the minimal case**:
   ```bash
   # Find smallest PDB with the issue
   python3 scripts/extract_minimal_pairs.py --pdb <PDB> --num-pairs 2
   ```

2. **Compare at the step level**:
   ```bash
   # See which step differs
   python3 scripts/compare_json.py compare <PDB> --verbose
   ```

3. **Compare at the sub-step level**:
   ```bash
   # For a specific pair that differs
   ./build/compare_pair_validation <PDB> <res1> <res2> --verbose
   ```

4. **Compare candidate lists**:
   ```bash
   # See what candidates were considered
   ./build/compare_best_partner <PDB> <res1> --verbose
   ```

5. **Fix and verify**:
   ```bash
   # Make fix, regenerate, compare again
   ./build/find_pair_app --fix-indices <PDB> /tmp/test.inp
   python3 scripts/compare_json.py compare <PDB> --verbose
   ```

---

## Success Criteria

✅ **Minimal test cases work**:
- Can extract 2-pair PDBs from any PDB
- Both legacy and modern process them
- Comparison works correctly

✅ **Step-level comparison works**:
- Can identify which step differs
- Can see intermediate states
- Can verify fixes at step level

✅ **Sub-step-level comparison works**:
- Can identify exact sub-operation that differs
- Can see all intermediate values
- Can verify fixes at sub-step level

✅ **Results match exactly**:
- All minimal test cases pass
- All full PDBs pass
- 100% match on all comparisons

---

## Implementation Status

### ✅ Phase 1: Minimal Test Cases - COMPLETE
- ✅ Created `scripts/extract_minimal_pairs.py`
- ✅ Tested successfully with multiple PDBs
- ✅ Can extract 2-pair minimal PDBs automatically

### ✅ Phase 2: JSON Output to Legacy Code - COMPLETE
- ✅ Added JSON writer functions for step-by-step debugging
- ✅ Modified `best_pair()` to record all candidates
- ✅ Modified `find_bestpair()` to record iterations and mutual best decisions
- ✅ Fixed iteration_state recording to show all pairs per iteration

### ✅ Phase 3: JSON Output to Modern Code - COMPLETE
- ✅ Added JSON writer methods to `JsonWriter` class
- ✅ Modified `BasePairFinder::find_best_partner()` to record all candidates
- ✅ Modified `BasePairFinder::find_best_pairs()` to record iterations and mutual best decisions
- ✅ All JSON output working correctly

### ✅ Phase 4: Comparison Tools - COMPLETE
- ✅ Created `scripts/compare_best_partner.py`
- ✅ Created `scripts/compare_iteration.py`
- ✅ Created `scripts/compare_mutual_best.py`
- ✅ Created `scripts/debug_minimal_case.sh` - Automated workflow script
- ✅ All tools tested and working

### ✅ Phase 5: Testing - COMPLETE
- ✅ Tested with 1AQ4_minimal_pairs_1_2.pdb - 100% match
- ✅ Tested with 1H4S_minimal_pairs_1_2.pdb - 100% match (on important fields)
- ✅ Verified all comparison tools work correctly

---

## Next Steps: Using the Infrastructure

Now that the infrastructure is complete, use it to debug differences:

1. **For any PDB with differences**: Use `./scripts/debug_minimal_case.sh <PDB_ID> 2`
2. **Identify divergence point**: Compare step-by-step to find where differences occur
3. **Fix systematically**: Make fixes at the appropriate level (candidate selection, iteration logic, etc.)
4. **Verify fixes**: Re-run debugging workflow to confirm fix works

---

## Related Documentation

- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing and comparison workflows
- [REF_FRAMES_DEBUGGING_APPROACH.md](REF_FRAMES_DEBUGGING_APPROACH.md) - Ref frames debugging
- [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md) - Legacy indices usage
- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - JSON structure

---

*This plan provides a systematic approach to achieving exact matching through incremental testing and debugging.*


---

## Implementation Progress

### ✅ Phase 1: Minimal Test Cases - COMPLETE
- Created `scripts/extract_minimal_pairs.py` 
- Tested successfully with 1AQ4.pdb
- Can extract 2-pair minimal PDBs automatically

### ✅ Phase 2: JSON Output to Legacy Code - COMPLETE
- Added JSON writer functions:
  - `json_writer_record_best_partner_candidates()` - Records all candidates considered
  - `json_writer_record_mutual_best_decision()` - Records mutual best decisions
  - `json_writer_record_iteration_state()` - Records state after each iteration
- Modified `best_pair()` to collect and record all candidates
- Modified `find_bestpair()` to record iteration states and mutual best decisions

**JSON Output Files** (when legacy code runs):
- `data/json_legacy/best_partner_candidates/<PDB_ID>.json` - All candidates for each residue
- `data/json_legacy/mutual_best_decisions/<PDB_ID>.json` - Mutual best decisions
- `data/json_legacy/iteration_states/<PDB_ID>.json` - State after each iteration

### ⏳ Phase 3: Create Comparison Tools - IN PROGRESS
Next: Create comparison tools to analyze the JSON output

