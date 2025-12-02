# Ref Frames Debugging - Quick Start Guide

**For when you're stuck and need a systematic approach**

---

## The Problem (What We Know)

1. ✅ **Origin and X-axis match perfectly** - Frame calculation is correct
2. ❌ **Y/Z axes are inverted** - Frame ordering issue (180° rotation)
3. ❌ **Pairs are different** - Pair selection issue (different residues)

---

## Quick Test (5 minutes)

```bash
# 1. Run the automated test
./scripts/test_ref_frames_single_pair.sh data/pdb/100D.pdb

# 2. Look at the output - it will show:
#    - Which pairs match/don't match
#    - Which axes are inverted
#    - Detailed differences
```

---

## Step-by-Step Debugging

### Step 1: Verify Frame Calculation (Already Working!)

The fact that origin and X-axis match means:
- ✅ Frame calculation algorithm is correct
- ✅ Mid-frame calculation is correct
- ❌ Only the Y/Z orientation is wrong

**Conclusion**: This is a frame ordering issue, not a calculation issue.

### Step 2: Test Frame Ordering

The Y/Z axes are inverted, which means we're using the wrong frame order.

**Test both orderings**:
```cpp
// Current: frame1 first, frame2 second
auto mid_frame = calc.calculate_pair_frame(frame1.value(), frame2.value());

// Try: frame2 first, frame1 second  
auto mid_frame = calc.calculate_pair_frame(frame2.value(), frame1.value());
```

**Action**: Try swapping the frame order in code and test.

### Step 3: Fix Pair Selection

Pairs are different, so we can't test frame ordering properly.

**Solution**: Use `--fix-indices` to ensure same pair selection:
```bash
./build/find_pair_app --fix-indices data/pdb/100D.pdb test.inp
```

**Then test again** to see if frame ordering works when pairs match.

---

## Recommended Approach

### Phase 1: Fix Frame Ordering (Easier)

1. **Ignore pair selection for now**
2. **Focus on one pair that exists in both**
3. **Test both frame orderings**:
   - Try `calculate_pair_frame(frame1, frame2)`
   - Try `calculate_pair_frame(frame2, frame1)`
4. **See which one matches legacy**

### Phase 2: Fix Pair Selection (Harder)

1. **Use `--fix-indices`** to match legacy pair selection
2. **Verify pairs match** using debug script
3. **Then test frame ordering** with matching pairs

### Phase 3: Integration

1. **Test with multiple PDBs**
2. **Verify both fixes work together**
3. **Document the solution**

---

## Key Insight

**The frame calculation is correct!** 

The only issue is:
1. Which frame goes first vs second (frame ordering)
2. Which pairs are selected (pair selection)

Both are fixable, but we should fix frame ordering first (easier), then pair selection.

---

## Next Action Items

1. ✅ **Use debug script** to see exact differences
2. ⏭️ **Try swapping frame order** in code
3. ⏭️ **Test with `--fix-indices`** to get matching pairs
4. ⏭️ **Verify solution** with multiple test cases

---

## When You're Stuck

1. **Run the debug script** - it shows exactly what's wrong
2. **Test with minimal case** - 1-2 pairs is enough
3. **Try both frame orderings** - one will match
4. **Use `--fix-indices`** - ensures same pairs
5. **Document what works** - so we don't repeat mistakes

