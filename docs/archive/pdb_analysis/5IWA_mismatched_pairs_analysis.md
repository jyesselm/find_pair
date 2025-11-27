# Mismatched Pairs Analysis: 5IWA

Generated: 2025-11-26 21:20:11

## Summary

- **Missing Pairs (in legacy but not modern)**: 3
- **Extra Pairs (in modern but not legacy)**: 4

## Missing Pairs (3)

### Pair (2551, 2553)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: adjust_pairQuality difference: 3.000000 (H-bond count difference: legacy=2, modern=3)

**Quality Scores**:
- Base Score: Legacy=10.950512, Modern=10.950512 (diff=0.000000)
- adjust_pairQuality: Legacy=-6.000000, Modern=-9.000000 (diff=3.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=4.950512, Modern=1.950512 (diff=3.000000)

**H-bonds**:
- Good H-bonds: Legacy=2, Modern=3
- Total H-bonds: Legacy=9, Modern=9

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**

### Pair (3435, 3579)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: adjust_pairQuality difference: 3.000000 (H-bond count difference: legacy=2, modern=3)

**Quality Scores**:
- Base Score: Legacy=5.072680, Modern=5.072680 (diff=0.000000)
- adjust_pairQuality: Legacy=-6.000000, Modern=-9.000000 (diff=3.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=-0.927320, Modern=-3.927320 (diff=3.000000)

**H-bonds**:
- Good H-bonds: Legacy=2, Modern=3
- Total H-bonds: Legacy=5, Modern=5

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**

### Pair (3473, 3489)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0

**Quality Scores**:
- Base Score: Legacy=12.029782, Modern=12.029782 (diff=0.000000)
- adjust_pairQuality: Legacy=-3.000000, Modern=-3.000000 (diff=0.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=9.029782, Modern=9.029782 (diff=0.000000)

**H-bonds**:
- Good H-bonds: Legacy=1, Modern=1
- Total H-bonds: Legacy=2, Modern=2

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**


## Extra Pairs (4)

### Pair (2551, 2554)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 2.000000

### Pair (3435, 3584)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 2.000000

### Pair (3465, 3474)

- **Status**: ✅ EXISTS in legacy validation (is_valid=0)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 3.000000

### Pair (3471, 3473)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 4.000000
