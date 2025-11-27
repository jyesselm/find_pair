# Parameter Comparison Investigation

## Status: Parser Working, Parameters Differ

### Fixed Issues ✅

1. **`.par` file creation** - Legacy now creates `.par` files correctly
2. **Argument parsing** - `--no-json` no longer treated as input file
3. **Parameter extraction** - Both legacy and modern parameters are extracted

### Current Status

**6V9Q Test Results**:
- Legacy step params: 6 extracted
- Modern step params: 6 extracted
- Legacy helix params: 6 extracted
- Modern helix params: 6 extracted

**BUT**: All parameters show significant differences (>5.0 for many values)

### Differences Found

Step 1 differences (example):
- Shift: legacy=-20.230, modern=-25.720, diff=5.490
- Slide: legacy=-5.835, modern=0.760, diff=6.595
- Rise: legacy=-1.067, modern=0.370, diff=1.437
- Tilt: legacy=-111.237, modern=-102.15, diff=9.087

These are substantial differences, not just rounding errors.

### Possible Causes

1. **Base pair ordering**: Legacy and modern might order base pairs differently
2. **Step definition**: Steps might be defined differently (which bp to which bp)
3. **Frame calculation**: Reference frames might be calculated differently
4. **Parameter calculation algorithm**: Core calculation might differ

### Next Steps

1. Verify base pair ordering matches between legacy and modern
2. Check if step numbering aligns correctly
3. Compare reference frames for matching base pairs
4. Investigate parameter calculation differences

## Investigation Needed

The parameter differences are too large to be rounding errors. We need to investigate:
- Are we comparing the same steps?
- Are the base pairs in the same order?
- Are the reference frames calculated the same way?

