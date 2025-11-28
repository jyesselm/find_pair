# Reference Frame Comparison Results

## Finding: Reference Frames Differ

For **1A34**, we compared reference frames for matching base pairs and found differences.

### Differences Found

**Pattern**: Sign differences in rotation matrix elements
- `rotation[0][1]`: Opposite signs (legacy=0.526626, modern=-0.526626)
- `rotation[0][2]`: Opposite signs (legacy=0.344861, modern=-0.344861)

This suggests a **coordinate system convention difference**:
- Different axis orientations (right-handed vs left-handed)
- Different frame definitions
- Different axis ordering

### Impact on Parameters

**This explains the parameter differences!**

If reference frames differ, then:
- Step parameters (calculated from frames) will differ
- Helical parameters (calculated from frames) will differ

This is the **root cause** of the parameter differences we observed.

## Next Steps

1. **Investigate coordinate system conventions**
   - Check if legacy uses right-handed vs modern uses left-handed (or vice versa)
   - Verify axis definitions match

2. **Check frame calculation code**
   - Compare `BaseFrameCalculator` implementation with legacy
   - Verify axis ordering and orientation

3. **Determine if difference is expected**
   - Some differences might be acceptable (different conventions)
   - Or might need to align conventions for exact matching

## Status

✅ **Reference frames extracted and compared**  
⚠️ **Frames differ** - Sign differences in rotation matrices  
🔍 **Root cause of parameter differences identified**

