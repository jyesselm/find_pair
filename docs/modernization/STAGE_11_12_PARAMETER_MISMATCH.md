# Stage 11/12: Step and Helical Parameter Mismatch Analysis

## Problem Statement

Modern code generates step parameters (Stage 11) and helical parameters (Stage 12) that **do not match** legacy output. Values differ significantly for non-consecutive base pairs in the helix structure.

**Example (1EHZ):**
- Legacy step (7,8): Shift=0.03, Rise=3.23
- Modern step (7,8): Shift=15.1, Rise=8.29

## Root Cause Analysis

### Legacy Data Flow

Legacy step/helical parameters are recorded by **`analyze.c`**, NOT `find_pair.c`:

```
find_pair.c              analyze.c
    |                        |
    v                        v
find base pairs  -->  .inp file  -->  ref_frames()  -->  get_parameters()
    |                                      |                    |
    v                                      v                    v
base_pair.json              orien[ds][num_bp*9]     bpstep_params.json
                            org[ds][num_bp*3]       helical_params.json
```

### Key Insight: Frame Arrays

Legacy `analyze.c` creates **separate frame arrays** indexed by base pair number:

```c
// In analyze.c
org = dmatrix(1, ds, 1, num_bp * 3);     // ds=1 for single strand, ds=2 for duplex
orien = dmatrix(1, ds, 1, num_bp * 9);

// ref_frames() populates these arrays
for (i = 1; i <= ds; i++) {       // For each strand
    for (j = 1; j <= num_bp; j++) { // For each base pair
        rnum = pair_num[i][j];     // Get residue number for strand i of BP j
        // ... calculate frame for residue rnum ...
        // Store at orien[i][(j-1)*9]
    }
}
```

### Frame Indexing

For `refs_i_j(j, j+1, orien[1], org[1], r1, o1, r2, o2)`:
- `r1, o1` = frame at `orien[1][(j-1)*9]` = frame of residue `pair_num[1][j]`
- `r2, o2` = frame at `orien[1][(j)*9]` = frame of residue `pair_num[1][j+1]`

**For 1EHZ (tRNA):**
```
BP7 = (7, 66)   →  pair_num[1][7] = 7    →  orien[1][54..62] = frame of residue 7
BP8 = (8, 14)   →  pair_num[1][8] = 8    →  orien[1][63..71] = frame of residue 8
```

Step (7,8) calculates params between **residue 7 and residue 8**, which ARE consecutive nucleotides (Rise ≈ 3.2 Å is typical).

### What Modern Code Does

Current `generate_modern_json.cpp`:

```cpp
// For each consecutive pair in base_pairs vector
for (size_t i = 0; i + 1 < base_pairs.size(); ++i) {
    const auto& pair1 = base_pairs[i];   // BP7 = (6, 65) 0-indexed → (7, 66) 1-indexed
    const auto& pair2 = base_pairs[i + 1]; // BP8 = (7, 13) 0-indexed → (8, 14) 1-indexed
    
    // Use frame1() from each pair
    ReferenceFrame frame1 = pair1.frame1().value();  // Frame of residue 7
    ReferenceFrame frame2 = pair2.frame1().value();  // Frame of residue 8
    
    auto step_params = param_calc.calculate_step_parameters(frame1, frame2);
}
```

**The frames used SHOULD be identical:**
- Modern `pair1.frame1()` = frame of `base_i` = residue 7
- Legacy `orien[1][54..62]` = frame of `pair_num[1][7]` = residue 7

### The Actual Discrepancy

Comparing frames directly:
```
Legacy BP7 (7,66) orien_i z-axis: [0.849374, -0.454042, -0.269091]
Modern BP7 (7,66) orien_i z-axis: [0.84937365, -0.45404234, -0.26909099]
```
**Frames MATCH!**

Comparing step params:
```
Legacy step (7,8): Rise=3.23, Twist=34.59
Modern step (7,8): Rise=8.29, Twist=157.11
```
**Results DIFFER significantly!**

### Hypothesis: Different Base Pair Ordering

The legacy `analyze.c` may use a **different ordering** of base pairs than modern `find_pairs()` returns.

Legacy workflow:
1. `find_pair` outputs `.inp` file with base pairs
2. `analyze` reads `.inp` and may reorder/filter pairs
3. `ref_frames()` builds frame arrays based on this ordering

Modern workflow:
1. `BasePairFinder::find_pairs()` returns pairs in selection order
2. Step params calculated directly from this vector

## Investigation Steps

### 1. Compare Base Pair Ordering

Check if legacy `.inp` file has same BP ordering as modern `find_bestpair_selection.json`:

```bash
# Generate .inp file from legacy
cd org && ./find_pair data/pdb/1EHZ.pdb 1EHZ.inp

# Compare with modern
cat data/json/find_bestpair_selection/1EHZ.json
```

### 2. Check Legacy ref_frames() Output

Add JSON recording to legacy `ref_frames()` to capture the exact frame arrays used:

```c
// In ana_fncs.c ref_frames()
json_writer_record_ref_frame_arrays(ds, num_bp, orien, org, pair_num);
```

### 3. Verify bpstep_par Input Frames

Add debug output to compare exact matrices passed to `bpstep_par()`:

```c
// In ana_fncs.c get_parameters()
refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
fprintf(stderr, "Step %ld->%ld: o1=[%.4f,%.4f,%.4f] o2=[%.4f,%.4f,%.4f]\n",
        j, j+1, o1[1], o1[2], o1[3], o2[1], o2[2], o2[3]);
```

## Potential Fixes

### Option A: Match Legacy Frame Array Building

Create equivalent frame arrays in modern code:

```cpp
// In generate_modern_json.cpp or new StepAnalyzer class
std::vector<ReferenceFrame> build_strand_frames(
    const std::vector<BasePair>& pairs,
    const Structure& structure,
    int strand  // 1 or 2
) {
    std::vector<ReferenceFrame> frames;
    for (const auto& pair : pairs) {
        // Get residue index for this strand
        size_t res_idx = (strand == 1) ? pair.residue_idx1() : pair.residue_idx2();
        const Residue* res = get_residue_by_idx(structure, res_idx);
        if (res && res->reference_frame()) {
            frames.push_back(*res->reference_frame());
        }
    }
    return frames;
}
```

### Option B: Use AnalyzeProtocol

The existing `AnalyzeProtocol` class should replicate legacy `analyze.c` behavior:

```cpp
// Current implementation in analyze_protocol.cpp
void AnalyzeProtocol::calculate_parameters(core::Structure& structure) {
    for (size_t i = start_idx; i + 1 < base_pairs_.size(); i += step_size_) {
        const auto& pair1 = base_pairs_[i];
        const auto& pair2 = base_pairs_[i + 1];
        
        // Use frame1() from each pair (matching legacy strand 1 frames)
        core::ReferenceFrame frame1 = pair1.frame1().value();
        core::ReferenceFrame frame2 = pair2.frame1().value();
        
        auto step_params = param_calculator_.calculate_step_parameters(frame1, frame2);
    }
}
```

This should be equivalent to legacy. Need to verify frame loading from `.inp` file matches.

### Option C: Direct Frame Array Comparison

Export frame arrays from both legacy and modern, compare element-by-element:

```python
# Compare frame arrays
def compare_frame_arrays(legacy_frames, modern_frames):
    for i, (leg, mod) in enumerate(zip(legacy_frames, modern_frames)):
        diff = np.linalg.norm(np.array(leg) - np.array(mod))
        if diff > 1e-4:
            print(f"Frame {i}: diff = {diff}")
```

## Implementation Plan

### Phase 1: Investigation

1. [ ] Add JSON output of legacy `orien[]` and `org[]` arrays from `ref_frames()`
2. [ ] Add JSON output of legacy `pair_num[][]` arrays
3. [ ] Compare frame arrays between legacy and modern
4. [ ] Identify where frames diverge

### Phase 2: Fix

1. [ ] If frame arrays differ: fix frame building in modern code
2. [ ] If frame arrays match but params differ: fix `bpstep_par` calculation
3. [ ] If ordering differs: ensure modern uses same BP ordering as legacy `.inp`

### Phase 3: Validation

1. [ ] Regenerate Stage 11/12 JSON for all PDBs
2. [ ] Run validation with updated comparators
3. [ ] Target: 99%+ pass rate

## Files Involved

### Legacy
- `org/src/analyze.c` - Main analyze workflow
- `org/src/ana_fncs.c` - `ref_frames()`, `get_parameters()`, `bpstep_par()`
- `org/src/cmn_fncs.c` - `refs_i_j()`, `ref_frame_i()`

### Modern
- `src/x3dna/algorithms/parameter_calculator.cpp` - `bpstep_par_impl()`
- `src/x3dna/protocols/analyze_protocol.cpp` - `calculate_parameters()`
- `tools/generate_modern_json.cpp` - JSON generation

### Tests
- `tests_python/validation/comparators/params.py` - Stage 11/12 comparison

## Debug Commands

```bash
# Run legacy with debug output
cd org && ./find_pair ../data/pdb/1EHZ.pdb 1EHZ.inp
cd org && ./analyze 1EHZ.inp 2>&1 | tee analyze_debug.log

# Generate modern JSON
./build/generate_modern_json data/pdb/1EHZ.pdb data/json --stage=all

# Compare step params
python3 << 'EOF'
import json
with open('data/json_legacy/bpstep_params/1EHZ.json') as f:
    leg = json.load(f)
with open('data/json/bpstep_params/1EHZ.json') as f:
    mod = json.load(f)
for i in range(min(10, len(leg), len(mod))):
    print(f"Step {leg[i]['bp_idx1']}->{leg[i]['bp_idx2']}:")
    print(f"  Legacy Rise: {leg[i]['params']['Rise']:.4f}")
    print(f"  Modern Rise: {mod[i]['rise']:.4f}")
EOF
```

## Next Steps (To Continue Investigation)

### Immediate Next Action

**Check `pair_checking()` in `org/src/cmn_fncs.c:3528`** to understand how `pair_num[][]` is populated from the `.inp` file.

The key question: **Does legacy use a different base pair ordering than modern?**

### Step-by-Step Plan

1. **Trace `pair_num` population:**
   ```bash
   # In org/src/cmn_fncs.c, find pair_checking() at line 3528
   grep -A 50 "void pair_checking" org/src/cmn_fncs.c
   ```

2. **Check if `.inp` file reorders pairs:**
   - Run legacy `find_pair` on 1EHZ to generate `.inp`
   - Compare pair order in `.inp` vs `find_bestpair_selection.json`

3. **Add debug output to legacy `get_parameters()`:**
   ```c
   // In ana_fncs.c around line 2017
   refs_i_j(j, j + 1, orien[i], org[i], r1, o1, r2, o2);
   if (j == 7) {
       fprintf(stderr, "DEBUG Step 7->8: o1=[%.4f,%.4f,%.4f] o2=[%.4f,%.4f,%.4f]\n",
               o1[1], o1[2], o1[3], o2[1], o2[2], o2[3]);
       fprintf(stderr, "DEBUG Step 7->8: pair_num[1][7]=%ld pair_num[1][8]=%ld\n",
               pair_num[1][7], pair_num[1][8]);
   }
   ```

4. **Compare frame origins directly:**
   - Legacy step (7,8) o1 should be origin of residue `pair_num[1][7]`
   - Modern step (7,8) frame1 is origin of `base_pairs[6].frame1()`
   - If these differ, the ordering is different

### Key Code Locations to Examine

| File | Line | Function | Purpose |
|------|------|----------|---------|
| `org/src/cmn_fncs.c` | 3528 | `pair_checking()` | Reads `.inp`, populates `pair_num` |
| `org/src/ana_fncs.c` | 1228 | `ref_frames()` | Builds frame arrays from `pair_num` |
| `org/src/ana_fncs.c` | 2015 | `get_parameters()` | Iterates and calls `bpstep_par` |
| `tools/generate_modern_json.cpp` | 220 | Step loop | Modern step calculation |

### Hypothesis to Test

**Hypothesis:** Legacy `pair_num[1][7]` and `pair_num[1][8]` give residues 7 and 8, but modern `base_pairs[6]` and `base_pairs[7]` might map to different residues due to ordering.

**Test:**
```python
# Run this after adding debug to legacy
# Compare:
# - Legacy pair_num[1][7], pair_num[1][8]
# - Modern base_pairs[6].residue_idx1(), base_pairs[7].residue_idx1()
```

### If Ordering Differs

If the ordering is different, we need to:
1. Understand how legacy orders pairs (from `.inp` file format)
2. Either:
   - Reorder modern `base_pairs` to match legacy order before step calculation
   - Or use `AnalyzeProtocol` which reads `.inp` files properly

### Quick Command to Resume

```bash
# Start here when resuming:
cd /Users/jyesselman2/Dropbox/2_code/cpp/find_pair_2

# Check pair_checking function
grep -A 80 "void pair_checking" org/src/cmn_fncs.c | head -100

# Generate .inp file for 1EHZ
cd org && ./find_pair ../data/pdb/1EHZ.pdb 1EHZ.inp 2>&1

# Look at .inp file format
cat org/1EHZ.inp
```

## References

- Legacy `bpstep_par()`: `org/src/ana_fncs.c:2019`
- Legacy `ref_frames()`: `org/src/ana_fncs.c:1228`
- Legacy `pair_checking()`: `org/src/cmn_fncs.c:3528`
- Modern `bpstep_par_impl()`: `src/x3dna/algorithms/parameter_calculator.cpp:109`
- JSON comparison docs: `docs/JSON_DATA_TYPES_AND_COMPARISONS.md`

