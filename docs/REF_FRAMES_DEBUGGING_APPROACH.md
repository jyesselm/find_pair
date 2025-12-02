# Ref Frames Debugging Approach

**Date**: 2025-01-XX  
**Status**: 🔧 Systematic debugging approach for ref_frames issues

---

## Problem Summary

We've been trying to fix ref_frames comparison for multiple days. The issues are:

1. **Different pair selection**: Legacy and modern select different base pairs
2. **Y/Z axis inversion**: Even when pairs match, Y/Z axes are inverted (180° rotation)
3. **Strand ordering**: The order of frames when calculating mid-frames may be wrong

---

## Breaking Down the Problem

### Step 1: Isolate the Issue

Use the new `debug_ref_frames.py` script to get detailed analysis:

```bash
# Generate both files
./build/find_pair_app data/pdb/100D.pdb test.inp
# (generate legacy ref_frames.dat separately)

# Detailed analysis
python3 scripts/debug_ref_frames.py legacy_ref_frames.dat ref_frames_modern.dat --verbose
```

This will show:
- Which pairs match by position
- Which pairs match by residue
- Which axes match/invert
- Detailed differences for each pair

### Step 2: Test with Minimal Cases

Create a minimal test case with just 1-2 base pairs:

```bash
# Extract first 2 base pairs from a PDB
# Use a very small PDB or extract a fragment
# Compare just those pairs
```

### Step 3: Verify Frame Calculation

For a single matching pair:
1. Check if origin and X-axis match (they should)
2. Check if Y/Z are inverted (180° rotation)
3. If inverted, the frame ordering is wrong

### Step 4: Verify Strand Assignment

For each pair:
1. Check which residue is strand 1 vs strand 2 in legacy .inp
2. Check which residue is strand 1 vs strand 2 in modern code
3. Verify the mapping is correct

---

## Systematic Testing Plan

### Phase 1: Single Pair Test

1. **Create minimal test case**:
   - Use smallest possible PDB with 1 base pair
   - Or extract a single pair from existing PDB
   - Generate both legacy and modern ref_frames

2. **Compare single pair**:
   ```bash
   python3 scripts/debug_ref_frames.py legacy.dat modern.dat --verbose
   ```

3. **Check results**:
   - Origin should match exactly
   - X-axis should match exactly
   - Y/Z should either match or be inverted (not random)

### Phase 2: Pair Selection Verification

1. **Check if pairs match**:
   - Use `--fix-indices` to ensure same residue selection
   - Compare .inp files to verify same pairs selected
   - If pairs differ, that's the root cause

2. **Fix pair selection first**:
   - This is likely the main issue
   - Once pairs match, frame ordering becomes testable

### Phase 3: Frame Ordering Test

Once pairs match:

1. **Test frame ordering logic**:
   - Add debug output to `write_ref_frames()`
   - Log which frame is used as strand 1 vs strand 2
   - Verify against legacy .inp file

2. **Test both orderings**:
   - Try frame1 first, frame2 second
   - Try frame2 first, frame1 second
   - See which matches legacy

### Phase 4: Integration Test

1. **Test with multiple PDBs**:
   - Start with small PDBs (1-5 pairs)
   - Gradually increase complexity
   - Document which cases work/fail

---

## Tools Created

### 1. `debug_ref_frames.py`

Comprehensive debugging script that:
- Parses both legacy and modern ref_frames files
- Compares by position and by residue pair
- Identifies inverted axes
- Shows detailed differences
- Provides summary statistics

**Usage**:
```bash
python3 scripts/debug_ref_frames.py legacy.dat modern.dat --verbose
```

### 2. `test_ref_frames_single_pair.sh`

Automated test script that:
- Generates both legacy and modern ref_frames
- Tests with and without legacy ordering
- Runs comparison analysis
- Provides clear output

**Usage**:
```bash
./scripts/test_ref_frames_single_pair.sh data/pdb/100D.pdb
```

### 3. `create_minimal_test.py`

Creates minimal test cases by extracting PDB fragments:
- Extract specific residue ranges
- Create small test cases (1-2 pairs)
- Easier to debug

**Usage**:
```bash
python3 scripts/create_minimal_test.py data/pdb/100D.pdb minimal.pdb \
    --start-res 1 --end-res 20 --chain A
```

### 2. Enhanced Logging (TODO)

Add debug output to `write_ref_frames()` to show:
- Which residue is strand 1 vs strand 2
- Which frame is used first vs second
- The calculated mid-frame values

---

## Next Steps

1. **Run debug script** on current test case to see exact differences
2. **Create minimal test case** with 1-2 pairs to isolate the issue
3. **Fix pair selection** if that's the root cause
4. **Test frame ordering** once pairs match
5. **Iterate** with progressively larger test cases

---

## Key Insights

1. **Origin and X-axis always match**: This confirms frames are calculated correctly
2. **Y/Z axes are inverted**: This is a frame ordering issue, not calculation
3. **Pairs differ**: This is a pair selection issue, separate from frame ordering
4. **Both issues need fixing**: But we should fix pair selection first, then frame ordering

---

## Recommendations

1. **Don't try to fix everything at once**: Break into smaller pieces
2. **Test with minimal cases first**: 1-2 pairs is enough to verify logic
3. **Use the debug script**: It provides detailed analysis
4. **Fix pair selection first**: This is likely blocking frame ordering tests
5. **Add debug logging**: See exactly what's happening in the code

---

## Related Files

- `scripts/debug_ref_frames.py` - Debugging script
- `src/x3dna/io/input_file_writer.cpp` - Frame writing code
- `docs/REF_FRAMES_COMPARISON_FIX.md` - Original fix documentation

