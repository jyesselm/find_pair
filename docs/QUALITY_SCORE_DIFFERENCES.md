# Quality Score Differences Analysis

**Purpose**: Document and investigate quality score differences that cause pair selection mismatches.

**Last Updated**: 2025-01-27

---

## Quality Score Components

The final quality score used for pair selection is calculated as:

```
final_score = base_score + adjust_pairQuality + bp_type_id_adjustment
```

Where:
- **base_score**: Calculated from geometric parameters (dorg, d_v, plane_angle)
- **adjust_pairQuality**: Adjustment based on H-bond quality
- **bp_type_id_adjustment**: -2.0 if `bp_type_id == 2` (Watson-Crick), 0 otherwise

---

## Known Quality Score Difference: 1T0K Pair (491, 492)

### Problem
- **Status**: Extra pair in modern, not in legacy
- **Pair**: (491, 492)
- **Root Cause**: Base quality score difference
- **bp_type_id**: No differences (both -1)

### Analysis Needed

**Components to Compare**:
1. **Base Score Calculation**:
   - `dorg`: Distance between origins
   - `d_v`: Vertical distance (abs(dot(dorg, zave)))
   - `plane_angle`: Angle between z-axes (0-90 degrees)
   - Formula: `base_score = dorg + 2.0 * d_v + plane_angle / 20.0`

2. **adjust_pairQuality**:
   - Based on H-bond quality
   - Formula depends on number of "good" H-bonds

3. **Frame Calculations**:
   - Frame origins (org[i], org[j])
   - Rotation matrices
   - Direction vectors (dir_x, dir_y, dir_z)

4. **H-bond Detection**:
   - Number of H-bonds detected
   - H-bond validation (good_hbatoms)
   - H-bond distances

---

## Investigation Tools

### 1. Compare Quality Score Components
```bash
build/compare_quality_score_components data/pdb/1T0K.pdb 491 492
```

### 2. Compare Frames
```bash
build/compare_frames_and_step_params data/pdb/1T0K.pdb 491 492
```

### 3. Check Validation Results
```python
# Load pair_validation JSON and compare components
python3 -c "
import json
# Load and compare validation records for pair (491, 492)
"
```

---

## Potential Root Causes

### 1. Frame Calculation Differences
- **Symptom**: Different frame origins or rotation matrices
- **Impact**: Affects all geometric calculations
- **Check**: Compare frame origins and rotation matrices

### 2. Distance Calculation Differences
- **Symptom**: Different `dorg`, `d_v`, or `dNN` values
- **Impact**: Directly affects base_score
- **Check**: Compare distance calculations step-by-step

### 3. Angle Calculation Differences
- **Symptom**: Different `plane_angle` values
- **Impact**: Affects base_score (plane_angle / 20.0)
- **Check**: Compare angle calculation methods

### 4. H-bond Detection Differences
- **Symptom**: Different number of H-bonds or H-bond quality
- **Impact**: Affects adjust_pairQuality
- **Check**: Compare H-bond detection and validation logic

### 5. Numerical Precision Differences
- **Symptom**: Small differences in all components
- **Impact**: Accumulates to significant difference in final score
- **Check**: Compare floating-point calculations

---

## Next Steps

1. **Run quality score comparison tool** for 1T0K pair (491, 492)
2. **Compare frame calculations** to verify they match
3. **Compare H-bond detection** to verify counts match
4. **Document findings** in this file
5. **Create fix** if root cause is identified

---

## Related Documentation

- `docs/KNOWN_DIFFERENCES_CATALOG.md`: Catalog of all known differences
- `docs/BP_TYPE_ID_BUG_FIX.md`: bp_type_id bug fix (related but different issue)
- `1T0K_mismatched_pairs_analysis.md`: Detailed analysis of 1T0K mismatch

