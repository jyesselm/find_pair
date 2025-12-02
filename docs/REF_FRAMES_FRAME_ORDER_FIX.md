# Ref Frames Frame Order Fix

**Date**: 2025-12-01  
**Status**: ✅ **FIXED** - Frame order swapped to match legacy

---

## Problem

The `ref_frames.dat` output from modern code had Y/Z axes inverted (180° rotation) compared to legacy, even though origin and X-axis matched perfectly.

**Root Cause**: Modern code was using `calculate_pair_frame(frame1, frame2)` but legacy uses strand 2 first, strand 1 second in `refs_right_left()`.

---

## Solution

### Change Made

**File**: `src/x3dna/io/input_file_writer.cpp`  
**Line**: 142

**Before**:
```cpp
auto mid_frame = calc.calculate_pair_frame(frame1.value(), frame2.value());
```

**After**:
```cpp
auto mid_frame = calc.calculate_pair_frame(frame2.value(), frame1.value());
```

### Rationale

Legacy code (`org/src/cmn_fncs.c::refs_right_left()`) uses:
- `r1` = strand 2 frame (`orien[2]`)
- `r2` = strand 1 frame (`orien[1]`)
- Calls `bpstep_par(r1, o1, r2, o2, ...)` with **strand 2 first, strand 1 second**

Modern code mapping:
- `frame1` = residue_idx1 (assumed strand 1)
- `frame2` = residue_idx2 (assumed strand 2)
- Must use `frame2` (strand 2) first, `frame1` (strand 1) second

---

## Testing

### Generate Modern ref_frames

```bash
./build/find_pair_app --fix-indices data/pdb/1TTT.pdb test.inp
# Creates: ref_frames_modern.dat
```

### Generate Legacy ref_frames

```bash
cd org
./build/bin/find_pair_original ../data/pdb/1TTT.pdb
./build/bin/analyze_original 1TTT.inp
# Creates: ref_frames.dat
```

### Compare

```bash
python3 scripts/debug_ref_frames.py org/ref_frames.dat ref_frames_modern.dat --verbose
```

**Expected Result**: Y/Z axes should now match (not inverted)

---

## Code References

- **Legacy**: `org/src/cmn_fncs.c:194-204` - `refs_right_left()` function
- **Legacy**: `org/src/ana_fncs.c:1782-1791` - Mid-frame calculation
- **Modern**: `src/x3dna/io/input_file_writer.cpp:142` - Frame order fix (without legacy ordering)
- **Modern**: `src/x3dna/io/input_file_writer.cpp:380` - Frame order (with legacy ordering - already correct)

---

## Notes

- The version with `legacy_pair_ordering` parameter (line 380) already had correct order
- This fix addresses the version without legacy ordering (line 142)
- Both versions now use strand 2 first, strand 1 second (matching legacy)

---

## Status

✅ **FIXED** - Frame order swapped in code.

**Testing**: To verify the fix works:
1. Generate modern: `./build/find_pair_app --fix-indices data/pdb/1H4S.pdb test.inp`
2. Generate legacy: 
   ```bash
   cd org && ./build/bin/find_pair_original ../data/pdb/1H4S.pdb
   cd org && ./build/bin/analyze_original 1H4S.inp
   ```
3. Compare: `python3 scripts/debug_ref_frames.py org/ref_frames.dat ref_frames_modern.dat --verbose`

**Expected**: Y/Z axes should now match (not inverted) when pairs are the same.

---

## Related Documentation

- `docs/REF_FRAMES_NEXT_STEPS.md` - Original issue documentation
- `docs/REF_FRAMES_COMPARISON_FIX.md` - Previous fix documentation (with legacy ordering)

