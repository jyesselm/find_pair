# Comparison Workflow - Frame-First Approach

**Last Updated**: 2025-11-27  
**Purpose**: Document the frame-first comparison approach

---

## Overview

Step parameters are calculated **FROM** reference frames. Therefore, frames must be verified to match before step parameter comparisons can be trusted.

## Comparison Order

### 1. **Frames First** (Always)

When comparing step parameters, frames are **always** checked first, even if frame comparison wasn't explicitly requested.

**Why**: Step parameters depend on frames:
- Step parameters are calculated between consecutive base pairs
- Each base pair has reference frames for both residues
- Step parameters use these frames to calculate geometric relationships

**If frames don't match**:
- Step parameter differences are expected
- Step parameter comparison results may be unreliable
- **Fix frames first** before trusting step parameter comparisons

### 2. **Step Parameters** (After Frames)

Step parameter comparison only proceeds after frame verification.

**Warnings**: If frames don't match, warnings are added to indicate:
- Frame mismatches detected
- Step parameter results may be unreliable
- Recommendation to fix frames first

---

## Usage

### Compare Everything (Recommended)

```bash
# Compares atoms, frames, and steps (frames checked automatically)
python3 scripts/compare_json.py compare 1H4S
```

### Compare Only Step Parameters

```bash
# Frames are automatically checked first even when comparing steps only
python3 scripts/compare_json.py steps 1H4S
```

**Note**: The `steps` command now automatically enables frame checking to ensure frames are verified before comparing step parameters.

### Compare Only Frames

```bash
# Just verify frames match
python3 scripts/compare_json.py frames 1H4S
```

---

## Warning Messages

### Frame Mismatches Before Step Comparison

If frames don't match, you'll see:

```
WARNING: Frame mismatches detected. Step parameter comparison may be 
unreliable since step parameters are calculated from frames. Fix frame 
differences first before comparing step parameters.

Step parameter comparison proceeding despite frame mismatches 
(5 mismatched frames, 2 missing residues). Results may be unreliable.
```

**Action**: Fix frame calculation differences first, then re-compare step parameters.

---

## Workflow

1. **Verify frames match**: `python3 scripts/compare_json.py frames <PDB>`
2. **If frames match**: Compare step parameters - results are reliable
3. **If frames don't match**: 
   - Fix frame calculation issues first
   - Re-verify frames
   - Then compare step parameters

---

## Implementation Details

### Automatic Frame Checking

The comparison code now:
1. **Always checks frames** when comparing steps (even if `compare_frames=False`)
2. **Warns if frames don't match** before proceeding with step comparison
3. **Still performs step comparison** but marks results as potentially unreliable

### Code Flow

```python
# In json_comparison.py:
# 1. Check if we need frames (for step comparison)
need_frame_check = self.compare_frames or self.compare_steps

# 2. Compare frames FIRST
if need_frame_check:
    frame_comparison = compare_frames(...)
    
    # Check if frames match
    if frame_comparison has mismatches:
        frames_match = False
        add_warning()

# 3. Compare step parameters (with warning if frames don't match)
if self.compare_steps:
    if not frames_match:
        add_warning_about_unreliable_results()
    
    compare_step_parameters(..., frame_comparison=frame_comparison)
```

---

## Best Practices

1. **Always verify frames first** - Use `frames` command before `steps`
2. **Fix frame issues** before trusting step parameter comparisons
3. **Use `compare` command** - It checks everything in the right order
4. **Pay attention to warnings** - They indicate when results may be unreliable

---

## Related Documentation

- [TESTING_GUIDE.md](TESTING_GUIDE.md) - Complete testing guide
- [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md) - JSON structure details
- [ALGORITHM_CRITICAL_GUIDE.md](ALGORITHM_CRITICAL_GUIDE.md) - Algorithm details

