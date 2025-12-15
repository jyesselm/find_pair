# Step Parameter Investigation

## Current Status (December 15, 2024)

**Pass Rate**: 80% (80/100 PDBs)

**Key Achievement**: bp_idx ordering now matches legacy **100%** (99/99 PDBs with step data)

**Remaining Issue**: 20 PDBs have matching bp_idx pairs but different step parameter VALUES

## What's Working

- ✅ Helix organization (pair ordering) - 100% match
- ✅ bp_idx assignments - 100% match
- ✅ HelixOrganizer five2three algorithm - produces correct pair ordering
- ✅ end_stack_xang threshold (110° → 125°) - matches legacy

## What Needs Investigation

The 20 failing PDBs have identical bp_idx pairs but different calculated parameter values. This means:
1. The pairs being compared are correct
2. The step ordering is correct
3. **The parameter calculation itself differs**

### Root Cause Candidates

1. **strand_swapped differences** - Different swap flags cause different frame selection
2. **ParameterCalculator implementation** - Midstep frame or parameter formulas may differ
3. **Frame selection logic** - How swap affects which frame (org_i vs org_j) is used

### Failing PDBs Analysis (20 total)

| Category | Count | Examples | Issue |
|----------|-------|----------|-------|
| Single step mismatch | 3 | 1Y27, 3GAO, 4FEN | Only last step (29→30) differs |
| Few step mismatches | 8 | 3UCU, 5CCX, 6MWN, 6XKN, 8EXA, 8UBD | 1-2 steps differ |
| Many step mismatches | 9 | 1TTT, 2DU3, 2EEW, 5FJ1, 5Y85, 6ICZ, 7YGA, 7YGB, 8U5Z, 8RUJ, 8Z1P | Multiple steps differ |

### Next Steps to Fix

#### Phase 1: Compare strand_swapped Arrays (Priority: HIGH)

The strand_swapped flag determines which frame is used for step calculations. If modern and legacy have different swap decisions, parameters will differ.

**Files to investigate:**
- `src/x3dna/algorithms/helix_organizer.cpp` - `ensure_five_to_three()` function
- `org/src/find_pair.c` - `five2three()` function (lines 1404-1509)

**How to debug:**
```bash
# Add JSON output for strand_swapped in legacy
# Compare with modern helix_organization JSON

# Example command to compare 1TTT:
./build/generate_modern_json data/pdb/1TTT.pdb data/json --stage=steps
# Check data/json/helix_organization/1TTT.json for strand_swapped values
```

**Key functions in five2three that affect swap:**
1. `first_step()` - Initial swap assignment based on O3'-P linkage
2. `wc_bporien()` - WC pair z-direction alignment check
3. `check_o3dist()` - O3' distance patterns
4. `check_schain()` - Sugar chain connectivity
5. `check_others()` - Frame vector alignment
6. `chain1dir()` - Strand1 reverse linkage check
7. `check_direction()` - Global helix direction analysis
8. `check_strand2()` - Additional strand corrections

#### Phase 2: Investigate Cross-Helix Steps (Priority: MEDIUM)

The 3 PDBs with single-step failures (1Y27, 3GAO, 4FEN) all fail on the LAST step (29→30). This step connects pairs from different helices.

**Questions to answer:**
- Does legacy calculate step params across helix boundaries?
- Is the last step always a cross-helix step?
- Should cross-helix steps use different frame selection?

**Files to check:**
- `tools/generate_modern_json.cpp` - Step calculation loop (lines ~290-340)
- Check if `helix_breaks` vector is being used correctly

#### Phase 3: Debug ParameterCalculator (Priority: LOW if swap is the issue)

If strand_swapped matches but params still differ, investigate the calculation itself.

**Files:**
- `src/x3dna/analyze/ParameterCalculator.cpp`
- `org/analyze.c`

**Key calculations:**
- Midstep frame computation
- Step parameter formulas (Shift, Slide, Rise, Tilt, Roll, Twist)

---

### Detailed Failure Analysis

#### Single-Step Failures (1Y27, 3GAO, 4FEN)

All have:
- 30 base pairs, 29 steps
- Steps 1-28 match perfectly
- Step 29→30 differs completely (different by 10-200+ units)

**Hypothesis**: The last step connects bp_idx 29 (helix boundary) to bp_idx 30 (isolated pair or different helix). Legacy may skip this step or calculate it differently.

**Test**: Compare legacy's bpstep_params JSON - does it even include step 29→30?

#### Multi-Step Failures (1TTT, 2DU3, etc.)

These have multiple steps with different values. Likely caused by strand_swapped differences that propagate through the helix.

**Pattern to look for**: Are the differing steps in specific helix segments? Check if all steps in helix 2+ are inverted.

---

## Historical Progress

**Previous Progress Made (December 2024)**:
1. Fixed base_pair recording to include ALL valid pairs (35 pairs for 1EHZ, not just 30 selected pairs)
2. Fixed bp_idx to use 1-based indexing (matching legacy)
3. Implemented HelixOrganizer for pair ordering (backbone connectivity order)
4. Added strand_swapped support (using org_j instead of org_i when strand is swapped)
5. **Fixed strand_swapped indexing bug** (was using helix position `i` instead of pair index `idx1`/`idx2`)
6. **Fixed check_direction to flip swapped values and reverse helix order** (matching legacy behavior)
7. **Added second check_direction call after check_strand2** (legacy line 1361)
8. **Fixed locate_helices to match legacy end_list traversal algorithm**
9. **Fixed calculate_context neighbor swapping logic** (legacy lines 931-941)
10. **Fixed end_stack_xang threshold** (110° → 125°) - matches legacy END_STACK_XANG
11. **Fixed backbone extraction** - uses legacy_residue_idx from atoms

**Bug Fix Details (December 14, 2024)**:

### Fix 1: strand_swapped indexing
The `strand_swapped` vector in `HelixOrganizer` is indexed by **original pair index**, not by helix position. The code in `generate_modern_json.cpp` was incorrectly using:
```cpp
// WRONG: indexed by helix position
bool swap1 = helix_order.strand_swapped[i];
bool swap2 = helix_order.strand_swapped[i + 1];
```
Fixed to:
```cpp
// CORRECT: indexed by original pair index
bool swap1 = helix_order.strand_swapped[idx1];
bool swap2 = helix_order.strand_swapped[idx2];
```
This fix improved step validation from 50% to 57%.

### Fix 2: check_direction helix flipping
The legacy `check_direction` function not only counts backbone linkage directions but also **modifies swapped values and helix order** under certain conditions. The modern implementation was only computing counts.

Legacy behavior (lines 1310-1322 in find_pair.c):
- If anti-parallel (strand1 forward, strand2 reverse) and `i1_first > j2_last`:
  - Flip ALL swapped values in the helix
  - Reverse the helix order
- If parallel (strand1 forward, strand2 forward) and `i1_first > j1_first`:
  - Flip ALL swapped values in the helix

Added this logic to modern `check_direction`, which improved from 57% to 68%.

**Remaining Issue (December 2024)**: 19% of structures still fail (19/100). These are NOT precision issues - all failures have large algorithmic differences (130-354 units in twist/tilt/roll).

### Remaining 19 Failures Analysis (December 14, 2024)

The 19 failing PDBs fall into several categories:

1. **Mixed direction helices** (e.g., 1TTT):
   - Helix has both forward and reverse linkages on strand1
   - Legacy check_strand2's mixed branch applies different swap logic
   - Modern detects mixed direction but may apply swaps differently

2. **tRNA-like junction patterns** (e.g., 1Y27, 2EEW, 3GAO):
   - All have 9-pair helices with swap difference at bp_idx 8
   - Last pair in helix has different swap in legacy vs modern
   - Related to endpoint handling in five2three

3. **Small helices with junctions** (e.g., 5CCX):
   - 3-pair helix with swap difference at position 2
   - Likely caused by first_step or first pass check differences

4. **Complex scattered swap differences** (e.g., 2DU3):
   - Multiple pairs in helix have different swaps
   - May be caused by subtle differences in check_* function implementations

**Failing PDBs**: 1TTT, 1Y27, 2DU3, 2EEW, 3GAO, 4FEN, 5CCX, 5FJ1, 5Y85, 6ICZ, 6MWN, 6XKN, 7YGA, 7YGB, 8EXA, 8RUJ, 8U5Z, 8UBD, 8Z1P

### Root Cause: Helix Organization Differences

Deep investigation of 1TTT revealed the fundamental issue:

**Legacy step 14** uses pairs at very different positions:
- Pair 24 (base_i=19, base_j=56)
- Pair 40 (base_i=54, base_j=58)

**Modern step 14** uses:
- Pair 14 (base_i=9, base_j=23)
- Pair 16 (base_i=10, base_j=45)

The `HelixOrganizer` produces a completely different pair ordering than legacy's `bp_context`/`locate_helix`/`five2three` algorithm. For 1TTT:
- Modern pair_order: [0,1,2,3,4,5,6, **24,25,26,27,28,29**, **13,15,14**, ...]
- Legacy organizes pairs in a different helix-based order

### Why Position Matching Doesn't Help

The step comparison attempts mst_org position matching, but:
- Position matching threshold is 0.01Å (tight)
- Different pairs produce mst_org positions 30+ Å apart
- Cannot match steps that use different underlying pairs

### Characteristics of 21 Failing Structures

All failures involve complex structures where legacy's helix organization differs from modern's:
- Multiple pairs involving the same base (e.g., base 8 in 4 different pairs)
- Non-WC pairs (AA, GG types)
- Unusual backbone connectivity patterns

### Possible Future Fixes

1. **Match legacy five2three exactly** - Very complex, many edge cases
2. **Compare by residue pairs** - Match steps using the same base combinations regardless of bp_idx
3. **Accept 79% pass rate** - Document that complex structures may have different helix organization

## Residue-Pair Comparison Results (December 2024)

A new comparison method was implemented to match steps by their underlying residue pairs
rather than by bp_idx. This reveals the true accuracy of step parameter calculation.

### Method
Instead of comparing step N vs step N, match steps that compute parameters for the
same pair of base pairs: `((base_i1, base_j1), (base_i2, base_j2))`.

### Results on 100-PDB Test Set

| Category | Count | Description |
|----------|-------|-------------|
| Perfect (100%) | 66 | All steps with matching residue pairs have identical parameters |
| High (>=90%) | 5 | Most steps match, minor differences |
| Medium (50-89%) | 11 | Significant strand swap differences |
| Low (<50%) | 16 | Major strand swap/frame selection differences |
| Errors | 2 | JSON loading errors |

**Overall: 2720/3812 steps match (71.4%)** when comparing same residue pairs.

### Key Insights

1. **66% of structures have PERFECT step calculation** - When the same residue pairs
   are used, the parameters match exactly, proving the calculation algorithm is correct.

2. **Strand swap differences explain remaining failures** - Even when legacy and modern
   use the same residue pairs, the strand swap (which frame to use: org_i vs org_j)
   can differ, causing different parameters.

3. **Two types of failures:**
   - **Helix ordering** (bp_idx differs): Same pairs used, but at different helix positions
   - **Strand swap** (bp_idx same, params differ): Same pairs, same position, but different
     frame selection due to five2three swap logic

### Low Match PDBs (<50%)

These structures have significant strand swap differences:
- 6BJX: 0% (complete swap inversion)
- 5CCX: 3.7%
- 4O26: 12.1%
- 7R6K: 16.0%
- 7S4V: 18.9%
- And 11 others...

### Conclusion

The step parameter CALCULATION is correct. The differences are in:
1. **Helix ordering** (which pairs go at which positions) - affects bp_idx matching
2. **Strand swap logic** (which frame to use for each pair) - affects parameter values

Both are due to the complexity of legacy's `five2three` algorithm which has many
edge cases for non-canonical structures.

## Historical Root Cause Analysis (December 2024)

**Original Problem**: Step parameter mismatches were caused by **different base pair lists** between legacy and modern code. Legacy includes additional base pairs (e.g., triplet interactions, non-WC pairs) that are interleaved in the list, causing bp_idx values to become misaligned.

## Root Cause Analysis (December 2024)

### The Real Problem: Base Pair List Differences

The step parameters are compared by `bp_idx` (base pair index), but legacy and modern have different base pair lists:

**1EHZ Example:**
- Legacy: 35 base pairs
- Modern: 30 base pairs
- 5 extra pairs in legacy: (9,23) AA, (10,45) GG, (22,46) GG, (25,45) CG, (33,35) UA

### How This Causes Step Mismatches

The extra pairs are **interleaved**, not appended:

```
bp_idx | Legacy pair    | Modern pair    | Match?
   1   | ( 1, 72)       | ( 1, 72)       | ✓
   ...
   8   | ( 8, 14)       | ( 8, 14)       | ✓
   9   | *( 9, 23) AA   | (10, 25)       | ✗ <- First divergence
  10   | (10, 25)       | (11, 24)       | ✗
  11   | *(10, 45) GG   | (12, 23)       | ✗
  12   | *(25, 45) CG   | (13, 22)       | ✗
  13   | (11, 24)       | (15, 48)       | ✗
  14   | (12, 23)       | (16, 59)       | ✗
```

When step (13, 14) is compared:
- **Legacy**: Step between (11, 24) → (12, 23)
- **Modern**: Step between (15, 48) → (16, 59)

These are **completely different base pairs**, so the parameters must differ!

### Characteristics of Extra Legacy Pairs

All 5 extra pairs have `dir_xyz[2] = 0.0` (perpendicular z-axes), indicating:
- Non-Watson-Crick base pairs (AA, GG)
- Triplet interactions (residue 45 paired with both 10 and 25)
- Bifurcated pairs

### Why Modern Misses These Pairs

The modern pair selection uses a greedy mutual-best algorithm that may not include:
1. Non-mutual-best pairs
2. Tertiary interactions (triplets)
3. Pairs that fail certain quality criteria

## Solutions

### Option 1: Fix Base Pair Selection (Recommended)
Make modern code select the same base pairs as legacy by:
- Including triplet interactions
- Relaxing quality criteria for certain pair types
- Matching legacy's secondary pair detection logic

### Option 2: Compare by Residue Pairs
Instead of comparing step(bp_idx1, bp_idx2), compare by the actual residue pairs:
- Build a map: (res_i1, res_j1) → (res_i2, res_j2) → step params
- Match steps that connect the same residue pairs
- This is more robust but requires changes to the validation code

### Option 3: Use Position-Based Matching (Current Fallback)
The current code attempts mst_org position matching but fails because:
- Different pairs have different mst_org positions
- The 0.01Å threshold is too tight
- Only works when pairs are in different order, not when pairs are different

## Previous Finding: Validation Bug Fixed

A bug was discovered and fixed in `x3dna_json_compare/step_comparison.py` that was silently masking step parameter mismatches.

### The Bug

The comparison code at lines 185-196 had flawed logic:
1. When bp_idx matching found mismatches, it tried position-based matching (using midstep origin coordinates)
2. Position matching often couldn't find matches (only ~50% matched within 0.01A threshold)
3. When position matching failed to find matches, it returned an empty list
4. The code used whichever method produced "fewer errors"
5. Empty list (0 mismatches) < any real mismatches, so real errors were silently dropped

### The Fix

Only prefer position matching if it matched ALL the common keys:
```python
# OLD (buggy):
if len(position_mismatches) < len(result.mismatched_steps):
    result.mismatched_steps = position_mismatches

# NEW (fixed):
if position_match_count >= len(common_keys) and len(position_mismatches) < len(result.mismatched_steps):
    result.mismatched_steps = position_mismatches
```

## Example: 1EHZ Step Mismatches

After the bug fix, 1EHZ shows 17 mismatched steps out of 29 total. Example mismatch at bp_idx (13, 14):

| Parameter | Legacy | Modern | Difference |
|-----------|--------|--------|------------|
| Shift | 3.59 | -0.72 | 4.31 |
| Slide | -0.80 | 1.25 | 2.05 |
| Rise | 3.23 | 3.22 | 0.01 |
| Tilt | 8.86 | varies | large |
| Roll | 16.38 | varies | large |
| Twist | 50.04 | varies | large |

The midstep frame origins also differ significantly:
- Legacy mst_org: [73.81, 66.95, 37.34]
- Modern midstep_frame.org: [75.82, 65.76, 38.94]
- Distance: ~2.83 Angstrom

## Areas to Investigate

### 1. Midstep Frame Calculation
The midstep frame origin (mst_org) differs between legacy and modern. This is the reference frame used for step parameter calculation. Differences here cascade to all parameters.

Files to investigate:
- `src/x3dna/analyze/ParameterCalculator.cpp` - Modern implementation
- `org/analyze.c` - Legacy implementation
- Look for `mst_org` or midstep frame calculation

### 2. Step Parameter Formula
The six step parameters (Shift, Slide, Rise, Tilt, Roll, Twist) are calculated from the midstep frame. The formulas may differ.

### 3. Pair Ordering for Steps
Steps are calculated between consecutive pairs. The ordering might differ:
- Legacy may use "five-to-three" ordering
- Modern may use sequential bp_idx ordering
- This could cause sign inversions or completely wrong pairs being compared

### 4. Helix Boundary Handling
At helix boundaries or non-consecutive base pairs, the step calculation logic may differ.

## Validation Commands

```bash
# Run step validation
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate steps --test-set 100

# Verbose comparison for specific PDB
X3DNA=/Users/jyesselman2/local/installs/x3dna python -m x3dna_json_compare.cli validate steps --pdb 1EHZ -v

# Compare .par files directly
python tools/compare_par_files.py bp_step.par bp_step_legacy.par
```

## Files Modified

- `x3dna_json_compare/step_comparison.py` - Fixed position matching bug
- `CLAUDE.md` - Updated pass rates to reflect reality

## Next Steps

1. Compare midstep frame calculation between legacy and modern
2. Verify step parameter formulas match
3. Check pair ordering logic (five-to-three vs sequential)
4. Fix modern implementation to match legacy output

## Deep Investigation (December 14, 2024)

### Analysis of 1EHZ Failure

Steps 1-13 match perfectly, but step 14+ diverges completely. The midstep origins jump to different locations:

| Step | Legacy mst_org | Modern mst_org | Match |
|------|----------------|----------------|-------|
| (13,14) | (73.8, 67.0, 37.3) | (73.8, 67.0, 37.3) | YES |
| (14,15) | (72.4, 48.4, 17.6) | (75.1, 61.0, 34.4) | NO |
| (15,16) | (73.4, 31.3, 0.5) | (75.4, 54.1, 31.0) | NO |

The legacy step (14,15) shows a ~20Å jump in z-coordinate - this is a **helix break**. The modern code stays more connected.

### Legacy bp_context / locate_helix Algorithm

The legacy algorithm uses these key data structures:
- `bp_order[i][1]`: -1 for middle pairs (have 2 opposite-side neighbors), 0 for endpoint-like pairs
- `bp_order[i][2]`, `bp_order[i][3]`: neighbor indices
- `end_list`: list of helix endpoints

In `locate_helix`, the critical check is:
```c
if (!bp_order[k][1]) {
    // k is endpoint-like, add at most one more, then break
    break;
}
```

This causes helices to terminate at endpoint-like pairs, creating shorter chains.

### Modern vs Legacy Helix Ordering

The modern `HelixOrganizer::locate_helices()` chains pairs based on:
1. Starting from endpoints
2. Following backbone-connected neighbors
3. Continuing until no unvisited neighbors

The legacy code terminates chains earlier based on `bp_order[k][1]`, leading to different helix boundaries and different pair orderings.

### Attempted Fixes (Reverted)

1. **direction_flag approach**: Added `direction_flag` field to PairContext to match legacy `bp_order[k][1]`. Implementation caused regression from 57% to 27%.

2. **first_step rewrite**: Rewrote `first_step()` to exactly match legacy. Caused regression from 57% to 25%.

The issue is more complex than individual function fixes - the entire helix construction and ordering algorithm differs fundamentally.

### Path Forward

1. **Option A: Detailed tracing** - Add logging to both legacy and modern to trace exact pair ordering decisions
2. **Option B: Match output directly** - Read legacy helix ordering from output files and use that directly
3. **Option C: Accept differences** - The 57% pass rate may represent inherent algorithmic differences that don't affect biological interpretation

### Key Legacy Functions

Located in `org/src/find_pair.c`:
- `bp_context()` (lines 885-1005): Sets neighbor info and direction flags
- `locate_helix()` (lines 1007-1069): Chains pairs into helices
- `five2three()` (lines 1384-1463): Ensures 5'→3' strand direction
- `first_step()` (lines 1082-1109): Initial strand assignment
- `wc_bporien()` (lines 1150-1167): WC pair orientation check
- `check_direction()` (lines 1246-1287): Strand direction counting
- `check_strand2()` (lines 1289-1342): Additional strand corrections

## Detailed Trace Analysis (December 14, 2024)

Added tracing to both legacy and modern code to capture exact decision points. Key findings:

### Finding 1: Pair Ordering Matches Exactly

The `locate_helices` algorithm in modern code produces **identical pair ordering** to legacy:

**Legacy bp_idx (1-based):**
```
1,2,3,4,5,6,7,25,26,27,28,29,30,15, 14,13,8,12,11,10,9,17,18,19,20,21,22,23,24, 16
```

**Modern pair_order (0-based, converted):**
```
1,2,3,4,5,6,7,25,26,27,28,29,30,15, 14,13,8,12,11,10,9,17,18,19,20,21,22,23,24, 16
```

Both produce 3 helices with the same boundaries.

### Finding 2: strand_swapped Values Differ Significantly

The `five2three`/`ensure_five_to_three` algorithm produces **different strand swap assignments**:

| Pair (1-based) | Legacy swapped | Modern swapped | Match |
|----------------|----------------|----------------|-------|
| 1-7            | 0              | 0              | ✓     |
| 8              | 1              | 0              | ✗     |
| 9-14           | 0              | 1              | ✗     |
| 15             | 1              | 1              | ✓     |
| 16             | 0              | 0              | ✓     |
| 17-24          | 1              | 0              | ✗     |
| 25-30          | 0              | 0              | ✓     |

**Critical observation**: The swap values are **inverted for entire helix segments** (pairs 8-14 and 17-24).

### Root Cause: five2three Algorithm Differences

The modern `ensure_five_to_three()` function in `HelixOrganizer` produces different swap decisions than legacy. Specifically:

1. **first_step() behavior differs**: Legacy's `first_step()` looks at just the first step of each helix, while modern analyzes the whole helix for direction voting.

2. **Swap propagation differs**: The chain of `wc_bporien()`, `check_o3dist()`, `check_schain()`, `check_others()` calls may produce different results due to subtle implementation differences.

3. **Helix-wise vs pair-wise**: Modern processes each helix independently, but the swap decision logic may not perfectly match legacy.

### Impact on Step Parameters

When strand is swapped, the step calculation uses `frame2` instead of `frame1` (org_j instead of org_i). Inverting the swap flag causes:
- Different midstep frame calculation
- Different step parameters (often completely different values, not just sign changes)

For step (14,15) in 1EHZ:
- Legacy: shift=15.4, rise=-41.4, twist=-133.8
- Modern: shift=-6.4, rise=-11.4, twist=-121.7

### Next Steps

1. **Deep dive into first_step()**: The initial swap assignment may cascade incorrectly
2. **Compare individual check functions**: `wc_bporien`, `check_o3dist`, etc.
3. **Consider reading legacy swap values**: Use JSON output from legacy to inform modern
4. **Test with simpler structures**: DNA duplexes may have simpler swap patterns

### Tracing Code Locations

Tracing was added to:
- `org/src/find_pair.c`: `locate_helix()` function - prints bp_order table and chaining decisions
- `org/src/find_pair.c`: `five2three()` function - prints final swapped array
- `src/x3dna/algorithms/helix_organizer.cpp`: `locate_helices()` - prints context table and chaining
- `tools/generate_modern_json.cpp`: Step calculation - prints strand_swapped values

To run with tracing:
```bash
# Legacy tracing
X3DNA=/Users/jyesselman2/local/installs/x3dna ./org/build/bin/find_pair_original data/pdb/1EHZ.pdb /dev/null 2>&1 | grep -E '\[LOCATE_HELIX|\[LEGACY'

# Modern tracing
./build/generate_modern_json data/pdb/1EHZ.pdb data/json_test --stage=steps 2>&1 | grep -E '\[MODERN|\[STEP_CALC'
```

## Fix: Pair Origin Calculation (December 14, 2024)

**Result**: Improved pass rate from 68% to 79% (+11 percentage points)

### Problem

The `get_pair_origin()` function in `helix_organizer.cpp` was returning `frame1.origin()` (the origin of the first base frame). However, legacy uses the **averaged origin** of both base frames for pair distance calculations.

This caused incorrect neighbor detection in `calculate_context()`, leading to different helix ordering and ultimately different step parameters.

### Solution

Updated `get_pair_origin()` to return the average of both frame origins:

```cpp
geometry::Vector3D HelixOrganizer::get_pair_origin(const core::BasePair& pair) const {
    // Legacy uses the average of both base origins (morg) for pair distance calculations
    // See refs_right_left() in cmn_fncs.c: morg[j] = (org[base1][j] + org[base2][j]) / 2
    if (!pair.frame1().has_value() && !pair.frame2().has_value()) {
        return geometry::Vector3D(0, 0, 0);
    }
    if (!pair.frame1().has_value()) {
        return pair.frame2().value().origin();
    }
    if (!pair.frame2().has_value()) {
        return pair.frame1().value().origin();
    }
    // Return average of both origins
    auto o1 = pair.frame1().value().origin();
    auto o2 = pair.frame2().value().origin();
    return geometry::Vector3D(
        (o1.x() + o2.x()) / 2.0,
        (o1.y() + o2.y()) / 2.0,
        (o1.z() + o2.z()) / 2.0
    );
}
```

### Legacy Reference

In legacy `cmn_fncs.c`, the `refs_right_left()` function (lines 288-306) computes `morg` as:
```c
for (i = 1; i <= inum_base; i++) {
    morg[j] += org[ik][j];  // Sum of origins for both bases
}
morg[i] /= inum_base;  // Average of origins
```

## Remaining 21 Failures Analysis

The 21 failing PDBs fall into two categories:

### 1. Floating-Point Precision (~10 PDBs)

Differences at the 1e-6 tolerance boundary. These are essentially matches with minor rounding differences:
```
tilt: -10.627700 vs -10.627702 (diff: 2.000000e-06)
roll: 13.141097 vs 13.141100 (diff: 3.000000e-06)
```

### 2. Helix Organization Differences (~11 PDBs)

Complex structures where legacy and modern produce different helix orders for helices 2+. Example from 8Z1P:

**Legacy Helix 2**: `13 8 12 11 10 9 16 17 18 19 20` (11 pairs, starts from pair 13)
**Modern Helix 2**: `9 10 11 12 8 13` (6 pairs, starts from pair 9)

This causes:
- Different pair sequences in step calculations
- Step indices (14-15, 15-16, etc.) refer to different pair combinations
- Large parameter differences (shift: -19.17 vs -3.04)

### Root Cause of Helix Order Differences

The endpoint selection order differs:
1. Modern `find_endpoints()` returns endpoints in pair index order (1, 9, 13, ...)
2. Legacy processes endpoints in discovery order from `bp_context()`
3. When endpoint 9 is processed before endpoint 13, the traversal direction reverses

### Potential Solutions

1. **Match endpoint ordering**: Modify `find_endpoints()` to match legacy's discovery order
2. **Relax tolerance**: Increase tolerance to 1e-5 to pass floating-point edge cases
3. **Use original pair indices**: Output step parameters with original pair indices instead of helix positions
4. **Accept differences**: Document that complex multi-helix structures may have different orderings

## Modern C++ Refactoring (December 14, 2024)

Refactored `helix_organizer.hpp` to use modern C++ idioms:

### 1. LinkDirection enum class

Replaced magic numbers (-1, 0, 1) with strongly-typed enum:
```cpp
enum class LinkDirection {
    None = 0,      ///< No O3'-P linkage detected
    Forward = 1,   ///< i → j linkage (5'→3' direction)
    Reverse = -1   ///< j → i linkage (reverse direction)
};
```

### 2. StrandResidues struct

Replaced output parameters with return struct:
```cpp
struct StrandResidues {
    size_t strand1;  ///< Residue index on strand 1 (1-based)
    size_t strand2;  ///< Residue index on strand 2 (1-based)
};

// Changed from:
void get_ij(pair, swapped, &n1, &n2);
// To:
StrandResidues get_strand_residues(pair, swapped);
```

### 3. Updated function signatures

```cpp
// Old: int is_linked(...)
// New: [[nodiscard]] LinkDirection check_linkage(...)
```

These changes improve type safety, readability, and maintainability while keeping the algorithm behavior identical.

## Helix Organization Debug Tool (December 2024)

Added JSON output to capture helix organization decisions from both legacy and modern code.

### New Files

- `org/src/json_writer.c`: Added `json_writer_record_helix_organization()`
- `src/x3dna/io/json_writer.cpp`: Added `record_helix_organization()`
- `scripts/compare_helix_org.py`: Compare helix org for a single PDB
- `scripts/compare_all_helix_org.py`: Compare helix org across all PDBs

### Comparison Results (21 Failing PDBs)

| PDB | Total | Helix% | Swap% | Status |
|-----|-------|--------|-------|--------|
| 6BJX | 134 | 7.5% | 92.5% | HELIX_DIFF |
| 6ICZ | 123 | 11.4% | 93.5% | HELIX_DIFF |
| 7YGB | 174 | 30.5% | 80.5% | HELIX_DIFF |
| 1Y27 | 30 | 33.3% | 66.7% | HELIX_DIFF |
| 3GAO | 30 | 33.3% | 66.7% | HELIX_DIFF |
| 4FEN | 30 | 33.3% | 66.7% | HELIX_DIFF |
| 8EXA | 61 | 42.6% | 98.4% | HELIX_DIFF |
| 7YGA | 173 | 45.1% | 90.8% | HELIX_DIFF |
| 5Y85 | 60 | 48.3% | 60.0% | HELIX_DIFF |
| 5FJ1 | 82 | 52.4% | 95.1% | OK |
| 5CCX | 28 | 53.6% | 96.4% | OK |
| 2EEW | 30 | 60.0% | 36.7% | SWAP_DIFF |
| 7YG8 | 162 | 60.5% | 98.1% | OK |
| 8RUJ | 164 | 62.2% | 99.4% | OK |
| 8U5Z | 39 | 74.4% | 84.6% | SWAP_DIFF |
| 8Z1P | 26 | 80.8% | 80.8% | SWAP_DIFF |
| 6XKN | 35 | 97.1% | 97.1% | OK |
| 1TTT | 92 | 100.0% | 93.5% | OK |
| 2DU3 | 28 | 100.0% | 57.1% | SWAP_DIFF |
| 6MWN | 76 | 100.0% | 98.7% | OK |
| 8UBD | 57 | 100.0% | 98.2% | OK |

**Summary:**
- **Helix match**: 4 perfect (100%), 5 high (>=90%)
- **Swap match**: 0 perfect (100%), 12 high (>=90%)

### Key Findings

1. **Two distinct failure modes:**
   - **HELIX_DIFF**: Modern assigns pairs to different helices than legacy (9 PDBs)
   - **SWAP_DIFF**: Same helix organization, but different strand swap decisions (5 PDBs)
   - **OK**: >= 90% helix and >= 90% swap match (7 PDBs)

2. **Helix organization is the primary issue:**
   - Only 5/21 failing PDBs have >= 90% helix match
   - Legacy's `bp_context` + `locate_helix` uses different criteria for neighbor detection

3. **Swap decisions are mostly correct when helix org matches:**
   - When helix org matches (>= 90%), swap accuracy is 93-99%
   - This validates the five2three swap algorithm implementation

### Next Steps (Plan for Matching Legacy)

**Phase 2**: Fix helix organization (locate_helices algorithm)
- Debug: Output legacy bp_order and end_list, compare with modern context
- Align neighbor detection criteria (z-direction check, distance thresholds)

**Phase 3-7**: Fix five2three swap algorithm (after helix org is fixed)
- first_step(), wc_bporien, check_o3dist, check_schain, check_others
- check_direction(), check_strand2()
- Second pass for WC check
