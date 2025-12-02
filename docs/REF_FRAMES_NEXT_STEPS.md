# Ref Frames - Next Steps

**Created**: 2025-01-XX  
**Purpose**: Clear action plan to fix ref_frames comparison

---

## What We Know

✅ **Frame calculation is CORRECT**
- Origin matches perfectly
- X-axis matches perfectly
- This means the algorithm works!

❌ **Two issues remain**:
1. **Y/Z axes inverted** (180° rotation) - Frame ordering issue
2. **Different pairs selected** - Pair selection issue

---

## Immediate Next Steps (Do These First)

### Step 1: Test Frame Ordering (5 minutes)

The Y/Z axes are inverted, which means we need to swap the frame order.

**Action**: In `src/x3dna/io/input_file_writer.cpp`, line 138, try:

```cpp
// Current (probably wrong):
auto mid_frame = calc.calculate_pair_frame(frame1.value(), frame2.value());

// Try this instead:
auto mid_frame = calc.calculate_pair_frame(frame2.value(), frame1.value());
```

**Test**:
```bash
make -j4
./build/find_pair_app data/pdb/100D.pdb test.inp
python3 scripts/debug_ref_frames.py /tmp/test_ref_frames/ref_frames.dat ref_frames_modern.dat
```

**Expected**: Y/Z axes should match (or at least not be inverted)

### Step 2: Use --fix-indices (10 minutes)

To test frame ordering properly, we need matching pairs.

**Action**:
```bash
# Generate with --fix-indices to match legacy pair selection
./build/find_pair_app --fix-indices data/pdb/100D.pdb test.inp

# Compare
python3 scripts/debug_ref_frames.py legacy_ref_frames.dat ref_frames_modern.dat --verbose
```

**Expected**: Pairs should match, then we can verify frame ordering

### Step 3: Verify Solution (5 minutes)

Once frame ordering is fixed and pairs match:

```bash
# Run automated test
./scripts/test_ref_frames_single_pair.sh data/pdb/100D.pdb

# Should show: Perfect matches: 10/10 (or close)
```

---

## Tools Available

### 1. Debug Script
```bash
python3 scripts/debug_ref_frames.py legacy.dat modern.dat --verbose
```
- Shows exactly what's wrong
- Identifies inverted axes
- Shows pair differences

### 2. Automated Test
```bash
./scripts/test_ref_frames_single_pair.sh data/pdb/100D.pdb
```
- Generates both legacy and modern
- Runs comparison
- Shows results

### 3. Minimal Test Cases
```bash
python3 scripts/create_minimal_test.py data/pdb/100D.pdb minimal.pdb \
    --start-res 1 --end-res 20 --chain A
```
- Create small test cases
- Easier to debug

---

## Recommended Workflow

1. **Run debug script** → See what's wrong
2. **Try frame order swap** → Fix Y/Z inversion
3. **Use --fix-indices** → Get matching pairs
4. **Verify with debug script** → Confirm fix
5. **Test with multiple PDBs** → Ensure it works generally

---

## Key Files

- `scripts/debug_ref_frames.py` - Main debugging tool
- `scripts/test_ref_frames_single_pair.sh` - Automated testing
- `src/x3dna/io/input_file_writer.cpp:138` - Frame ordering code
- `docs/REF_FRAMES_QUICK_START.md` - Quick reference
- `docs/REF_FRAMES_DEBUGGING_APPROACH.md` - Detailed approach

---

## Success Criteria

✅ **Fixed when**:
- Origin matches: ✓ (already working)
- X-axis matches: ✓ (already working)
- Y-axis matches: ✓ (needs frame order fix)
- Z-axis matches: ✓ (needs frame order fix)
- Pairs match: ✓ (needs --fix-indices)

---

## Don't Give Up!

The hard part (frame calculation) is already working. We just need to:
1. Swap frame order (1 line change)
2. Use --fix-indices (1 flag)

Then test and verify. You've got this! 💪

