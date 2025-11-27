# Mismatched Pairs Analysis: 7XHT

Generated: 2025-11-26 20:43:22

## Summary

- **Missing Pairs (in legacy but not modern)**: 2
- **Extra Pairs (in modern but not legacy)**: 3

## Missing Pairs (2)

### Pair (52, 126)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0

**Quality Scores**:
- Base Score: Legacy=12.089820, Modern=12.089820 (diff=0.000000)
- adjust_pairQuality: Legacy=-0.000000, Modern=-0.000000 (diff=0.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=12.089820, Modern=12.089820 (diff=0.000000)

**H-bonds**:
- Good H-bonds: Legacy=0, Modern=0
- Total H-bonds: Legacy=1, Modern=1

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**

### Pair (127, 137)

- **Status**: ✅ EXISTS in validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0

**Quality Scores**:
- Base Score: Legacy=8.114037, Modern=8.114037 (diff=0.000000)
- adjust_pairQuality: Legacy=-6.000000, Modern=-6.000000 (diff=0.000000)
- bp_type_id: Legacy=-1, Modern=0
- Final Score: Legacy=2.114037, Modern=2.114037 (diff=0.000000)

**H-bonds**:
- Good H-bonds: Legacy=2, Modern=2
- Total H-bonds: Legacy=4, Modern=4

⚠️ **QUALITY SCORE DIFFERENCE DETECTED**


## Extra Pairs (3)

### Pair (126, 137)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 2.000000

### Pair (127, 136)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: HIGH
- **Root Cause**: Base quality score difference: 2.000000

### Pair (52, 53)

- **Status**: ✅ EXISTS in legacy validation (is_valid=1)
- **Priority**: MEDIUM
- **Root Cause**: bp_type_id difference: legacy=-1, modern=0
