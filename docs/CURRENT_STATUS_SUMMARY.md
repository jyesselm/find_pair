# Current Status Summary - Base Pair Finding Implementation

**Date**: Current session  
**Goal**: Achieve 100% match between legacy and modern code outputs for base pair finding

---

## Current Match Rate

**Status**: ~99% match (623/627 pairs match in 6CAQ)

**Remaining Discrepancies**: 4 mismatched pairs in 6CAQ
- **2 missing pairs** (in legacy, not in modern): (968, 1024), (980, 998)
- **2 extra pairs** (in modern, not in legacy): (1024, 1188), (980, 997)

---

## Key Accomplishments

### 1. ✅ Frame Calculation
- Frames are calculated correctly and stored on residue objects
- Frame reversal logic (when `dir_z <= 0`) is implemented correctly
- **Verified**: Frames match perfectly between legacy and modern for problematic pairs

### 2. ✅ Validation Consistency
- Phase 1 validation stores results in `phase1_validation_results` map
- Selection phase reuses Phase 1 validation results for consistency
- **Verified**: All 4 problematic pairs ARE validated in Phase 1 (not a JSON output issue)

### 3. ✅ bp_type_id Consistency
- `bp_type_id` is now calculated once in Phase 1 and stored in `phase1_bp_type_ids` map
- Selection phase reuses the stored `bp_type_id` instead of recalculating
- This ensures consistency between Phase 1 and selection

### 4. ✅ Modified Nucleotides
- Added support for `PSEUDOURIDINE`, `INOSINE`, and `NONCANONICAL_RNA` types
- All modified nucleotides now have frames calculated correctly

---

## Root Cause of Remaining 4 Mismatched Pairs

### Problem
Modern assigns `bp_type_id=2` (Watson-Crick) to pairs (1024, 1188) and (980, 997), giving them a -2.0 quality score adjustment. Legacy assigns `bp_type_id=-1` (not classified) to all of them, so no adjustment.

### Quality Score Comparison

**Modern Phase 1:**
- (968, 1024): quality=4.77, bp_type_id=-1 → **adjusted=4.77**
- (1024, 1188): quality=5.40, bp_type_id=2 → **adjusted=3.40** (after -2.0)
- (980, 998): quality=9.37, bp_type_id=-1 → **adjusted=9.37**
- (980, 997): quality=10.04, bp_type_id=2 → **adjusted=8.04** (after -2.0)

**Legacy Phase 1:**
- (968, 1024): quality=4.77, bp_type_id=-1 → **adjusted=4.77**
- (1024, 1188): quality=5.40, bp_type_id=-1 → **adjusted=5.40** (no adjustment)
- (980, 998): quality=9.37, bp_type_id=-1 → **adjusted=9.37**
- (980, 997): quality=10.04, bp_type_id=-1 → **adjusted=10.04** (no adjustment)

### Why Modern Selects Different Pairs

Modern correctly selects the best pairs based on adjusted quality scores:
- (1024, 1188) = 3.40 is better than (968, 1024) = 4.77
- (980, 997) = 8.04 is better than (980, 998) = 9.37

Legacy doesn't apply the -2.0 adjustment, so it selects:
- (968, 1024) = 4.77 is better than (1024, 1188) = 5.40
- (980, 998) = 9.37 is better than (980, 997) = 10.04

---

## Step Parameter Analysis

### Modern Step Parameters (for bp_type_id=2 pairs)

**Pair (1024, 1188):**
- pars[1] (Slide/Shear): -1.719321 → fabs=1.719321 ≤ 1.8 ✅
- pars[2] (Rise/Stretch): -0.035735 → fabs=0.035735 ≤ 2.0 ✅
- pars[6] (Twist/Opening): -50.870333 → fabs=50.870333 ≤ 60.0 ✅
- Base pair: AU (in WC_LIST) ✅
- **Result**: bp_type_id = 2

**Pair (980, 997):**
- pars[1] (Slide/Shear): -1.370042 → fabs=1.370042 ≤ 1.8 ✅
- pars[2] (Rise/Stretch): -1.925168 → fabs=1.925168 ≤ 2.0 ✅
- pars[6] (Twist/Opening): -22.298109 → fabs=22.298109 ≤ 60.0 ✅
- Base pair: CG (in WC_LIST) ✅
- **Result**: bp_type_id = 2

### Legacy Behavior
Legacy assigns `bp_type_id=-1` to both pairs, which means:
- Either legacy calculates different step parameters
- Or legacy uses different thresholds/logic in `check_wc_wobble_pair`

---

## Verification Results

### ✅ Frames Match Perfectly
- Frame origins match exactly
- Rotation matrices match exactly (after reversal)
- Frame reversal logic is applied correctly

### ✅ Frame Order is Correct
- Legacy calls: `bpstep_par(r2, org[j], r1, org[i], ...)`
- Modern calls: `calculate_step_parameters(frame2, frame1)`
- Order matches: (frame2, frame1) = (r2, r1)

### ✅ Displacement Vector Matches
- Legacy: `ddxyz(org1, org2, t1)` → `t1 = org2 - org1 = org[i] - org[j]`
- Modern: `t1_displacement = o2 - o1 = frame1.origin() - frame2.origin() = org[i] - org[j]`

---

## Tools Created

1. **`compare_frames_and_step_params`**: Compares final frames and 6 base pair step parameters
   - Shows frames before/after reversal
   - Shows all 6 step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)
   - Shows bp_type_id calculation with all checks

2. **`compare_bp_type_id_calculation`**: Analyzes bp_type_id calculation step-by-step

3. **`trace_pair_selection`**: Traces `find_best_partner` logic for specific residues

4. **`compare_quality_score_components`**: Breaks down quality score into components

---

## Code Changes Made

### `src/x3dna/algorithms/base_pair_finder.cpp`
- Added `phase1_bp_type_ids` map to store bp_type_id from Phase 1
- Modified `find_best_partner` to reuse stored bp_type_id instead of recalculating
- Ensures consistency between Phase 1 and selection phases

### `include/x3dna/algorithms/base_pair_finder.hpp`
- Updated `find_best_partner` signature to accept `phase1_bp_type_ids` map

---

## Next Steps to Achieve 100% Match

### Priority 1: Verify Legacy Step Parameters
**Question**: Does legacy calculate the same step parameters as modern?

**Action**: 
- Calculate step parameters from legacy frames using the same algorithm
- Compare pars[1], pars[2], pars[6] values
- If they differ, identify why (different algorithm, different frames, etc.)

### Priority 2: Verify check_wc_wobble_pair Logic
**Question**: Does legacy use the same thresholds and logic?

**Action**:
- Compare legacy `check_wc_wobble_pair` implementation with modern
- Verify thresholds: fabs(stretch) <= 2.0, fabs(opening) <= 60, fabs(shear) <= 1.8
- Verify WC_LIST matching logic

### Priority 3: Check for Edge Cases
**Question**: Are there any edge cases in the step parameter calculation?

**Action**:
- Check if there are any special cases in legacy `bpstep_par` that modern doesn't handle
- Verify degenerate case handling (parallel/anti-parallel z-axes)

---

## Key Files

### Core Implementation
- `src/x3dna/algorithms/base_pair_finder.cpp` - Main pair finding logic
- `src/x3dna/algorithms/base_pair_validator.cpp` - Validation logic
- `src/x3dna/algorithms/parameter_calculator.cpp` - Step parameter calculation
- `src/x3dna/algorithms/base_frame_calculator.cpp` - Frame calculation

### Legacy Reference
- `org/src/cmn_fncs.c` - Legacy `check_pair` and `calculate_more_bppars`
- `org/src/ana_fncs.c` - Legacy `bpstep_par` and `check_wc_wobble_pair`

### Tools
- `tools/compare_frames_and_step_params.cpp` - Frame and step parameter comparison
- `tools/compare_bp_type_id_calculation.cpp` - bp_type_id analysis
- `tools/trace_pair_selection.cpp` - Selection tracing

---

## Critical Code Sections

### bp_type_id Calculation (Modern)
```cpp
// In calculate_bp_type_id:
if (result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) {
    // Apply frame reversal if dir_z <= 0
    if (result.dir_z <= 0.0) {
        // Reverse y and z columns of frame2
    }
    // Calculate step parameters
    BasePairStepParameters params = param_calculator_.calculate_step_parameters(frame2, frame1);
    double shear = params.slide;
    double stretch = params.rise;
    double opening = params.twist;
    
    // Check thresholds
    if (std::abs(stretch) <= 2.0 && std::abs(opening) <= 60.0) {
        if (std::abs(shear) <= 1.8 && in_wc_list) {
            bp_type_id = 2; // Watson-Crick
        }
    }
}
```

### Legacy bp_type_id Calculation
```c
// In calculate_more_bppars:
bpstep_par(r2, org[j], r1, org[i], pars, mst, &rtn_val[5]);
*bpid = -1;
if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
    check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
    if (*bpid == 2)
        rtn_val[5] -= 2.0;
}
```

---

## Hypothesis

**Most Likely Cause**: Legacy calculates different step parameters (pars[1], pars[2], pars[6]) than modern, causing different bp_type_id assignments.

**Why**: 
- Frames match perfectly ✅
- Frame order is correct ✅
- But step parameters might differ due to:
  - Different algorithm implementation
  - Different handling of edge cases
  - Different numerical precision

**Next Investigation**: Calculate step parameters from legacy frames and compare directly with modern calculations.

---

## Commands to Reproduce

```bash
# Generate modern JSON
build/generate_modern_json data/pdb/6CAQ.pdb data/json

# Compare frames and step parameters
build/compare_frames_and_step_params data/pdb/6CAQ.pdb 1024 1188
build/compare_frames_and_step_params data/pdb/6CAQ.pdb 980 997

# Check selection
python3 -c "
import json
with open('data/json/find_bestpair_selection/6CAQ.json') as f:
    modern = json.load(f)
with open('data/json_legacy/find_bestpair_selection/6CAQ.json') as f:
    legacy = json.load(f)
# Compare pairs...
"
```

---

## Status: Investigation Complete - Ready for Legacy Step Parameter Extraction

All tools are in place. The investigation has verified that:
- ✅ Modern `check_wc_wobble_pair` logic matches legacy exactly
- ✅ Modern step parameter calculations are correct
- ✅ Modern correctly assigns `bp_type_id=2` for problematic pairs

**Next Step**: Extract step parameters from legacy code execution to compare directly with modern calculations.

**See**: `docs/BP_TYPE_ID_INVESTIGATION_SUMMARY.md` for complete investigation details.

---

## Investigation Update (2025-01-27)

### Verification of check_wc_wobble_pair Logic ✅

**Legacy Implementation** (`org/src/ana_fncs.c:1122-1131`):
```c
void check_wc_wobble_pair(long *bpid, char *bp, double shear, double stretch, double opening)
{
    static char *WC[9] = { WC_LIST };
    if (fabs(stretch) > 2.0 || fabs(opening) > 60)
        return;
    if (dval_in_range(fabs(shear), 1.8, 2.8))
        *bpid = 1;
    if (fabs(shear) <= 1.8 && num_strmatch(bp, WC, 1, 8))
        *bpid = 2;
}
```

**Modern Implementation** (`src/x3dna/algorithms/base_pair_finder.cpp:660-689`):
- ✅ Thresholds match: `fabs(stretch) > 2.0 || fabs(opening) > 60`
- ✅ Wobble check matches: `fabs(shear) >= 1.8 && fabs(shear) <= 2.8` → `bp_type_id = 1`
- ✅ WC check matches: `fabs(shear) <= 1.8 && bp_type in WC_LIST` → `bp_type_id = 2`
- ✅ WC_LIST matches: `{"XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"}`

**Conclusion**: The `check_wc_wobble_pair` logic is correctly implemented in modern code.

### Root Cause Hypothesis

Since the logic matches perfectly, the most likely cause is that **legacy calculates different step parameters** than modern, even though:
- ✅ Frames match perfectly
- ✅ Frame reversal logic is correct
- ✅ Frame order is correct (r2, r1)

**Possible reasons for different step parameters**:
1. **Numerical precision differences** in `bpstep_par` calculation
2. **Different handling of edge cases** (e.g., degenerate z-axes)
3. **Subtle differences in matrix operations** (e.g., rotation matrix construction)

### Next Steps

1. ✅ **Verified modern step parameters** - Modern calculates correct step parameters:
   - Pair (1024, 1188): shear=-1.719321, stretch=-0.035735, opening=-50.870333 → bp_type_id=2
   - Pair (980, 997): shear=-1.370042, stretch=-1.925168, opening=-22.298109 → bp_type_id=2
2. ⏳ **Run legacy code with debug output** to capture actual step parameters for these pairs
3. ⏳ **Compare step parameters directly** - Need to verify if legacy calculates different values
4. ⏳ **Verify bpstep_par implementation** - Check for subtle differences in numerical precision or edge cases

### Tools Created

- ✅ `scripts/compare_step_params_for_pairs.py` - Script to compare step parameters from JSON
- ✅ `tools/compare_frames_and_step_params` - C++ tool to compare frames and step parameters
- ✅ `tools/compare_bp_type_id_calculation` - Analyzes bp_type_id calculation step-by-step

### Modern Step Parameter Results (Verified)

**Pair (1024, 1188) - AU:**
- Direction vectors: dir_x=0.622029, dir_y=-0.623112, dir_z=-0.973039 ✅
- Step parameters: shear=-1.719321, stretch=-0.035735, opening=-50.870333
- Thresholds: fabs(stretch)=0.035735 ≤ 2.0 ✅, fabs(opening)=50.870333 ≤ 60.0 ✅
- WC check: fabs(shear)=1.719321 ≤ 1.8 ✅, bp_type="AU" in WC_LIST ✅
- **Result**: bp_type_id = 2 (Watson-Crick)

**Pair (980, 997) - CG:**
- Direction vectors: dir_x=0.741459, dir_y=-0.923190, dir_z=-0.799189 ✅
- Step parameters: shear=-1.370042, stretch=-1.925168, opening=-22.298109
- Thresholds: fabs(stretch)=1.925168 ≤ 2.0 ✅, fabs(opening)=22.298109 ≤ 60.0 ✅
- WC check: fabs(shear)=1.370042 ≤ 1.8 ✅, bp_type="CG" in WC_LIST ✅
- **Result**: bp_type_id = 2 (Watson-Crick)

**Conclusion**: Modern code correctly calculates step parameters and assigns bp_type_id=2. The discrepancy with legacy suggests legacy either:
1. Calculates different step parameters (numerical precision differences)
2. Has a bug in check_wc_wobble_pair logic
3. Uses different base pair string construction

