# Ref Frames Comparison Fix

**Date**: 2025-01-XX  
**Status**: 🔧 Fixing ref_frames output to match legacy

---

## Problem

The `ref_frames.dat` output from modern code doesn't match legacy output. The frames are calculated correctly, but the ordering of frames when calculating mid-frames may differ.

## Legacy Behavior

### Legacy `ref_frames()` Function
- Calculates frames for each residue in base pairs
- Stores frames in `orien[1][j]` (strand 1) and `orien[2][j]` (strand 2)
- For each base pair position `j`, calculates frame for residue `pair_num[1][j]` (strand 1) and `pair_num[2][j]` (strand 2)

### Legacy `print_analyze_ref_frames()` Function
- Outputs mid-frames calculated from residue frames
- Uses `refs_right_left()` to get frames:
  - `r1` = strand 2 frame (`orien[2]`)
  - `r2` = strand 1 frame (`orien[1]`)
- Calls `bpstep_par(r1, o1, r2, o2, ...)` with **strand 2 first, strand 1 second**
- Stores result in `bp_org` and `bp_orien`
- Outputs to `ref_frames.dat`

### Key Code (org/src/cmn_fncs.c:194-204)
```c
void refs_right_left(long bnum, double **orien, double **org, double **r1, double *o1,
                     double **r2, double *o2)
{
    long ioffset3, ioffset9;
    ioffset3 = (bnum - 1) * 3;
    ioffset9 = (bnum - 1) * 9;
    cpxyz(org[2] + ioffset3, o1);  // r1 = strand 2 (right)
    cpxyz(org[1] + ioffset3, o2);  // r2 = strand 1 (left)
    orien2mst(orien[2], ioffset9, r1);
    orien2mst(orien[1], ioffset9, r2);
}
```

### Legacy Mid-Frame Calculation (org/src/ana_fncs.c:1782-1791)
```c
refs_right_left(i, orien, org, r1, o1, r2, o2);
bpstep_par(r1, o1, r2, o2, bp_par[i], mfi, mfoi);
// ... store mfoi and mfi in bp_org and bp_orien
cpxyz(mfoi, bp_org + (i - 1) * 3);
mst2orien(bp_orien, (i - 1) * 9, mfi);
```

**Key Point**: Legacy uses **strand 2 first, strand 1 second** when calculating mid-frames.

## Modern Behavior

### Modern `write_ref_frames()` Function
- Gets frames from `bp.frame1()` and `bp.frame2()`
- `frame1()` corresponds to `residue_idx1()`
- `frame2()` corresponds to `residue_idx2()`
- Calls `calculate_pair_frame(frame1, frame2)` with **residue_idx1 first, residue_idx2 second**

### Issue
- Modern code uses `residue_idx1` first, `residue_idx2` second
- Legacy uses strand 2 first, strand 1 second
- If `residue_idx1` corresponds to strand 1 and `residue_idx2` corresponds to strand 2, the order is reversed!

## Solution

The modern code needs to:
1. Determine which residue is strand 1 and which is strand 2 from legacy .inp file
2. When calculating mid-frame, use strand 2 first, strand 1 second (matching legacy)
3. Handle legacy pair ordering when available

### Implementation

The fix is in `src/x3dna/io/input_file_writer.cpp`:
- When `legacy_pair_ordering` is available, parse it to determine strand assignment:
  - In legacy .inp: `res1` = strand 1, `res2` = strand 2
  - Use strand 2 frame first, strand 1 frame second when calling `calculate_pair_frame`
- Without legacy ordering, use default assumption (residue_idx1 = strand 1, residue_idx2 = strand 2)
- Always call `calculate_pair_frame` with strand 2 first, strand 1 second (matching legacy `refs_right_left`)

### Changes Made

1. **Updated `write_ref_frames()` with legacy ordering** (lines 306-327):
   - Determines which residue is strand 1 vs strand 2 from legacy .inp file
   - Uses strand 2 first, strand 1 second when calculating mid-frame
   - Matches legacy `refs_right_left` behavior

2. **Updated `write_ref_frames()` without legacy ordering** (line 136):
   - Uses default assumption: residue_idx1 = strand 1, residue_idx2 = strand 2
   - Still uses strand 2 first, strand 1 second (frame2, frame1)
   - Note: For exact legacy matching, use version with legacy ordering

### Testing

Use `scripts/compare_ref_frames.py` to compare:
```bash
# Generate legacy ref_frames.dat
cd org && ./build/bin/find_pair_original ../data/pdb/1H4S.pdb
cd ../org && ./build/bin/analyze_original 1H4S.inp
# ref_frames.dat will be in org/ directory

# Generate modern ref_frames_modern.dat (with legacy ordering)
./build/find_pair_app --legacy-inp org/1H4S.inp data/pdb/1H4S.pdb output.inp

# Compare
python3 scripts/compare_ref_frames.py org/ref_frames.dat ref_frames_modern.dat --verbose
```

### Status

✅ **FIXED** - Modern code now correctly uses strand 2 first, strand 1 second when calculating mid-frames, matching legacy behavior.

---

## Related Files

- `org/src/ana_fncs.c:1782-1791` - Legacy mid-frame calculation
- `org/src/cmn_fncs.c:194-204` - Legacy `refs_right_left()` function
- `src/x3dna/io/input_file_writer.cpp:76-171` - Modern `write_ref_frames()` (without legacy ordering)
- `src/x3dna/io/input_file_writer.cpp:253-362` - Modern `write_ref_frames()` (with legacy ordering)
- `scripts/compare_ref_frames.py` - Comparison script

