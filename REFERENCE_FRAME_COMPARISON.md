# Reference Frame Comparison Results

## Finding: Reference Frames Match ✅

For **1A34**, we compared reference frames for matching base pairs.

### Frame Convention Issue (Resolved)

**Initial Finding**: Sign differences in rotation matrix elements  
**Root Cause**: Legacy stores frames for both directions of each pair `(i,j)` and `(j,i)`. The frame orientation depends on which residue is "first" in the pair.

**Solution**: When comparing, use the frame where the residue is the "first" (base_i) in the pair. Modern uses this convention.

### Final Result

✅ **ALL REFERENCE FRAMES MATCH** (tolerance: 0.01)

For **1A34**:
- 9 matching base pairs
- 18 reference frames compared (2 per pair)
- All frames match within tolerance

### Impact on Parameter Differences

**Reference frames match**, so parameter differences must come from:
1. **Parameter calculation algorithm** differences (not frame differences)
2. **Numerical precision** or rounding differences
3. **Step definition** differences (which frames to use for which step)

Since frames match, the parameter differences are likely due to:
- Different algorithms for calculating step/helical parameters from frames
- Different handling of edge cases
- Numerical precision differences in intermediate calculations

## Status

✅ **Reference frames match** - After accounting for frame convention  
✅ **Frames are NOT the root cause** of parameter differences  
🔍 **Parameter differences likely due to calculation algorithm differences**

