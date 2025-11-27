# Mismatched Pairs Analysis: 6CAP

Generated: 2025-11-26 21:22:49

## Summary

- **Missing Pairs (in legacy but not modern)**: 4
- **Extra Pairs (in modern but not legacy)**: 5

## Missing Pairs (4)

### Pair (495, 498)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 1.010605

**Quality Scores**:
- Base Score: Legacy=9.069235, Modern=8.058630 (diff=1.010605)
- adjust_pairQuality: Legacy=-0.000000, Modern=-3.000000 (diff=3.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=9.069235, Modern=5.058630 (diff=4.010605)

**H-bonds**:
- Good H-bonds: Legacy=0, Modern=1
- Total H-bonds: Legacy=4, Modern=4

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**

### Pair (443, 459)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0

**Quality Scores**:
- Base Score: Legacy=6.723128, Modern=6.723128 (diff=0.000000)
- adjust_pairQuality: Legacy=-0.000000, Modern=-0.000000 (diff=0.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=6.723128, Modern=6.723128 (diff=0.000000)

**H-bonds**:
- Good H-bonds: Legacy=0, Modern=0
- Total H-bonds: Legacy=2, Modern=2

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**

### Pair (499, 512)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0

**Quality Scores**:
- Base Score: Legacy=5.332441, Modern=5.332441 (diff=0.000000)
- adjust_pairQuality: Legacy=-3.000000, Modern=-3.000000 (diff=0.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=2.332441, Modern=2.332441 (diff=0.000000)

**H-bonds**:
- Good H-bonds: Legacy=1, Modern=1
- Total H-bonds: Legacy=2, Modern=2

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**

### Pair (1123, 1125)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0

**Quality Scores**:
- Base Score: Legacy=9.088725, Modern=9.088724 (diff=0.000001)
- adjust_pairQuality: Legacy=-3.000000, Modern=-3.000000 (diff=0.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=6.088725, Modern=6.088724 (diff=0.000001)

**H-bonds**:
- Good H-bonds: Legacy=1, Modern=1
- Total H-bonds: Legacy=1, Modern=1

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**


## Extra Pairs (5)

### Pair (444, 459)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 2.000000

### Pair (495, 512)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 0.704067

### Pair (1063, 1072)

- **Status**: ✅ EXISTS in legacy validation (is_valid=0)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 3.000000

### Pair (1105, 1125)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 4.000000

### Pair (499, 508)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0
