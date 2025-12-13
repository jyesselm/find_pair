# Stage 11/12 Step & Helical Parameter Investigation

## Current Status

**Pass Rate: 96.1%** (3374/3512 pass, 138 fail, 90 skip)

This is up from ~63% before fixes were applied.

## What Was Fixed

### 1. Frame Selection (Fixed)
- **Problem**: Was calculating middle frames from pair.frame1() and pair.frame2()
- **Solution**: Legacy uses `frame1` (org_i, orien_i) directly
- **Verification**: `legacy mst_org = (bp1.org_i + bp2.org_i) / 2` - confirmed exact match

### 2. Pair Ordering (Partially Fixed)
- **Problem**: Was iterating through base pairs in sequential order (1,2,3,...)
- **Solution**: Now uses `HelixOrganizer.pair_order` for backbone connectivity order
- **Issue**: HelixOrganizer doesn't exactly match legacy's five2three algorithm

## Root Cause of Remaining 138 Failures

The legacy code uses a complex "five2three" algorithm that traces through helices following backbone connectivity in a specific pattern:

1. Forward on strand 1 (5'→3')
2. At helix boundary, jump to **opposite end** of strand 2
3. Backward on strand 2 (3'→5')
4. At next boundary, connect back to strand 1
5. Continue forward on strand 1

**Example from 1NYI** (19 base pairs, 2 helix sections):

| Step | Legacy Order | Modern Order | Match |
|------|--------------|--------------|-------|
| 1-6  | d1p1→d1p7    | same         | ✅    |
| 7    | d1p7→**d2p19** (jumps to opposite end!) | d1p7→d1p8 | ❌ |
| 8-11 | d2p19→d2p15 (backward) | different | ❌ |
| 12   | d1p8→d2p15 | different | ❌ |
| 13-18| d1p8→d1p14 | same order | ✅ |

## Files to Read

### Core Implementation Files
```
src/x3dna/algorithms/helix_organizer.cpp  - Current helix ordering algorithm
include/x3dna/algorithms/helix_organizer.hpp - HelixOrdering struct definition
tools/generate_modern_json.cpp - Step param generation (lines 255-307)
```

### Validation Files
```
tests_python/validation/comparators/params.py - Step/helical param comparison
tests_python/validation/runner.py - Validation runner with base_pair loading
scripts/validate_all_stages_11_12.py - Batch validation script
```

### Test Data
```
data/json_legacy/bpstep_params/1NYI.json - Example legacy step params
data/json_test/bpstep_params/1NYI.json - Modern step params for comparison
data/json_legacy/base_pair/1NYI.json - Legacy base pairs (38 entries = 19×2 duplexes)
```

## Key Investigation Points

### 1. Understand Legacy five2three Algorithm
The legacy code is in `org/src/analyze.c`. Key functions:
- `five2three()` - Orients strands in 5'→3' direction
- `locate_helix()` - Chains base pairs into helices
- `bp_context()` - Identifies helix endpoints via backbone

Legacy stores pairs interleaved: [d1p1, d2p1, d1p2, d2p2, ...]
- Even indices (0,2,4...) = duplex 1
- Odd indices (1,3,5...) = duplex 2

### 2. HelixOrganizer Improvements Needed
Current `helix_organizer.cpp` uses geometric proximity to chain pairs. It needs:
- Proper backbone O3'-P linkage checking
- Helix boundary detection
- Strand alternation at boundaries (jump to opposite end of complementary strand)

### 3. Debugging Commands

```bash
# Generate for specific PDB
./build/generate_modern_json data/pdb/1NYI.pdb /tmp/test_1NYI --stage=steps

# Compare step params
python3 -c "
import json
with open('/tmp/test_1NYI/bpstep_params/1NYI.json') as f: mod = json.load(f)
with open('data/json_legacy/bpstep_params/1NYI.json') as f: leg = json.load(f)
for i, (l, m) in enumerate(zip(leg[:18], mod)):
    print(f'Step {i+1}: leg_shift={l[\"params\"][\"Shift\"]:.2f} mod_shift={m[\"shift\"]:.2f}')
"

# Find what pairs legacy uses for each step
python3 -c "
import json, numpy as np
with open('data/json_legacy/base_pair/1NYI.json') as f: bp = json.load(f)
with open('data/json_legacy/bpstep_params/1NYI.json') as f: steps = json.load(f)
for step in steps[:18]:
    mst = np.array(step['mst_org'])
    for i in range(len(bp)):
        for j in range(len(bp)):
            if i >= j: continue
            avg = (np.array(bp[i]['org_i']) + np.array(bp[j]['org_i'])) / 2
            if np.linalg.norm(avg - mst) < 0.001:
                print(f'Step ({step[\"bp_idx1\"]},{step[\"bp_idx2\"]}): ({bp[i][\"base_i\"]},{bp[i][\"base_j\"]}) + ({bp[j][\"base_i\"]},{bp[j][\"base_j\"]})')
                break
        else: continue
        break
"
```

## Approach to Reach 100%

### Option A: Fix HelixOrganizer (Recommended)
Modify `helix_organizer.cpp` to implement the exact legacy five2three algorithm:

1. **Add backbone connectivity tracking**
   - Check O3'-P linkage between consecutive residues
   - Mark helix breaks where backbone is discontinuous

2. **Implement strand alternation**
   - At helix boundaries, identify which strand to follow
   - Jump to opposite end of complementary strand

3. **Test incrementally**
   - Start with simple multi-helix structures
   - Progress to complex structures like 1NYI, 1GID

### Option B: Match by Residue Pairs (Alternative)
Instead of matching by bp_idx, match steps by the actual residue pairs involved:
- Already partially implemented in `params.py` but disabled when pair counts match
- Would require tracking which pairs are used in each step

### Option C: Accept 96% (Pragmatic)
The 138 failures are complex structures with unusual helix topologies. The current implementation handles 96% of structures correctly.

## Failed PDB Examples to Investigate

```
1GID, 1I2Y, 1N38, 1NYI, 1NUJ, 1NUV, 1Q2R, 1QBP, 1R3O, 1R3E
```

These all have multiple helix sections with complex junctions.

## Summary

The core issue is that modern's HelixOrganizer chains pairs differently than legacy's five2three algorithm at helix boundaries. To reach 100%, need to either:
1. Replicate the exact legacy pair ordering in HelixOrganizer
2. Match steps by actual residue pairs instead of bp_idx
3. Accept the current 96.1% pass rate as sufficient
