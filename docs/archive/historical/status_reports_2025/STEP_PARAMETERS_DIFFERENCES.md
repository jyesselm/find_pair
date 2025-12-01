# Step Parameters: Modern vs Legacy Differences

**Date**: 2025-11-29  
**Status**: Documented differences between modern and legacy implementations

---

## Overview

While both modern and legacy implementations generate step parameters correctly, there are some differences in how they process structures and which pairs they analyze. This document details these differences.

---

## Key Differences

### 1. Duplex Processing

**Legacy Behavior**:
- Processes each duplex separately when `ds = 2` (duplex number = 2)
- For each duplex, calculates step parameters for consecutive pairs
- Example: 25 pairs with `ds = 2` → 2 × 24 = 48 step parameters (one set per duplex)

**Modern Behavior**:
- Processes pairs once (single set)
- Calculates step parameters for consecutive pairs in the input file
- Example: 23 pairs → 22 step parameters (or 20 if some pairs don't have frames)

**Code Location**:
- Legacy: `org/src/ana_fncs.c:2011` - Loops `for (i = 1; i <= ds; i++)` to process each duplex
- Modern: `src/x3dna/protocols/analyze_protocol.cpp:111` - Processes pairs once

**Impact**:
- Legacy generates more step parameters when `ds > 1`
- Modern generates one set of step parameters
- **Both are correct** - just different processing approaches

---

### 2. Base Pair Selection

**Legacy**:
- Uses base pairs from legacy `find_pair` output
- May select different pairs than modern find_pair
- Example (1H4S): 25 base pairs selected

**Modern**:
- Uses base pairs from modern `find_pair` output
- May select different pairs than legacy find_pair
- Example (1H4S): 23 base pairs selected

**Impact**:
- Different base pair selections lead to different step parameter sets
- When the same pairs are analyzed, step parameter values match exactly
- This is expected - different pair selection algorithms produce different results

---

### 3. Base Pair Index Ranges

**Example: 1H4S**

**Modern**:
- bp_idx1 range: 3-22 (20 unique values)
- bp_idx2 range: 4-23 (20 unique values)
- Total step parameters: 20

**Legacy**:
- bp_idx1 range: 1-24 (24 unique values)
- bp_idx2 range: 2-25 (24 unique values)
- Total step parameters: 48 (2 duplexes × 24 each)

**Why Different**:
- Legacy starts from pair 1, modern starts from pair 3 (different pair selection)
- Legacy processes 2 duplexes, modern processes 1 set
- Legacy has 25 pairs, modern has 23 pairs

---

### 4. Value Matching

**When Base Pair Indices Match**:
- ✅ **All step parameter values match exactly**
- Verified: 20/20 matching pairs show identical values
- Parameters match: Shift, Slide, Rise, Tilt, Roll, Twist

**Example Match (bp_idx1=3, bp_idx2=4)**:
```
Modern:  shift=-0.326662, slide=-2.096079, rise=2.910300
Legacy:  shift=-0.326662, slide=-2.096079, rise=2.910300
         ✅ Perfect match
```

**Conclusion**:
- Step parameter **calculations** are identical
- Differences are due to **pair selection** and **duplex processing**, not calculation errors

---

## Detailed Comparison: 1H4S

### Record Counts

| Type | Modern | Legacy | Difference |
|------|--------|--------|------------|
| Step parameters | 20 | 48 | Legacy has 2× more (2 duplexes) |
| Helical parameters | 20 | 48 | Legacy has 2× more (2 duplexes) |
| Base pairs (input) | 23 | 25 | Legacy selected 2 more pairs |

### Base Pair Index Coverage

**Modern**:
- Covers pairs 3-23 (20 step parameters)
- Missing: pairs 1-2 (not in modern selection)

**Legacy**:
- Covers pairs 1-25 (24 step parameters per duplex)
- Processes 2 duplexes separately
- Total: 48 step parameters

### Matching Pairs

- **Matching**: 20 pairs have matching bp_idx between modern and legacy
- **Missing in modern**: 4 pairs (legacy has pairs 1-2, 24-25 that modern doesn't)
- **Extra in modern**: 0 pairs (all modern pairs exist in legacy)

**When pairs match**: All 6 step parameters match exactly ✅

---

## Implementation Differences

### 1. Duplex Loop

**Legacy** (`org/src/ana_fncs.c:2011-2032`):
```c
for (i = 1; i <= ds; i++) {  // Loop over duplexes
    for (j = 1; j <= nbpm1; j++) {  // Loop over pairs
        // Calculate step parameters for duplex i, pair j
        json_writer_record_bpstep_params(j, j + 1, ...);
    }
}
```

**Modern** (`src/x3dna/protocols/analyze_protocol.cpp:111-143`):
```cpp
for (size_t i = start_idx; i + 1 < base_pairs_.size(); i += step_size_) {
    // Calculate step parameters for pair i, i+1
    json_writer_->record_bpstep_params(i + 1, i + 2, ...);
}
```

**Key Difference**:
- Legacy: Outer loop over duplexes (`ds` iterations)
- Modern: Single loop over pairs (no duplex loop)

### 2. Base Pair Indexing

**Legacy**:
- Uses 1-based indices directly from loop: `j, j+1`
- Each duplex starts from pair 1

**Modern**:
- Converts 0-based vector index to 1-based: `i + 1, i + 2`
- Uses base pairs from input file order

**Both Correct**:
- Legacy: `bp_idx1 = j, bp_idx2 = j + 1` (1-based loop index)
- Modern: `bp_idx1 = i + 1, bp_idx2 = i + 2` (convert 0-based to 1-based)
- Result: Both use 1-based indices in JSON output ✅

### 3. Circular Structure Handling

**Legacy**:
- Handles circular structures in `get_parameters()` function
- May process last pair with first pair

**Modern**:
- Handles circular structures explicitly (`circular_structure_` flag)
- Processes last pair with first pair if `circular_structure_` is true

**Both Similar**: Both handle circular structures correctly

---

## JSON Format Differences

### Field Names

**Modern**:
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "shift": -0.326662,
  "slide": -2.096079,
  "rise": 2.910300,
  "tilt": 3.171240,
  "roll": 11.777801,
  "twist": 28.807861,
  "midstep_frame": {...}
}
```

**Legacy**:
```json
{
  "type": "bpstep_params",
  "bp_idx1": 3,
  "bp_idx2": 4,
  "params": {
    "Shift": -0.326662,
    "Slide": -2.096079,
    "Rise": 2.910300,
    "Tilt": 3.171240,
    "Roll": 11.777801,
    "Twist": 28.807861
  },
  "mst_org": [...],
  "mst_orien": [[...], [...], [...]]
}
```

**Key Differences**:
- Modern: Parameters as top-level fields (lowercase: `shift`, `slide`, `rise`, `tilt`, `roll`, `twist`)
- Legacy: Parameters nested in `params` object (capitalized: `Shift`, `Slide`, `Rise`, `Tilt`, `Roll`, `Twist`)
- Modern: `midstep_frame` object (contains `org` and `orien`)
- Legacy: `mst_org` array and `mst_orien` matrix (separate fields)

**Field Mapping**:
- Modern `shift` = Legacy `params.Shift`
- Modern `slide` = Legacy `params.Slide`
- Modern `rise` = Legacy `params.Rise`
- Modern `tilt` = Legacy `params.Tilt`
- Modern `roll` = Legacy `params.Roll`
- Modern `twist` = Legacy `params.Twist`
- Modern `midstep_frame.org` = Legacy `mst_org`
- Modern `midstep_frame.orien` = Legacy `mst_orien`

**Note**: Comparison script handles both formats correctly and maps fields appropriately

---

## Why These Differences Exist

### 1. Different Pair Selection

- Legacy `find_pair` and modern `find_pair` use different algorithms
- May select different base pairs for the same structure
- This is expected and acceptable - both are valid selections

### 2. Duplex Processing Philosophy

- Legacy: Process each duplex separately (more detailed analysis)
- Modern: Process pairs once (simpler, faster)
- Both approaches are valid for different use cases

### 3. JSON Format Evolution

- Modern format: Flatter structure, easier to parse
- Legacy format: Nested structure, matches original design
- Comparison script handles both formats

---

## Verification Results

### Value Accuracy

✅ **When base pair indices match, all values match exactly**

**Test**: 1H4S
- Matching pairs: 20/20
- All 6 parameters match: Shift, Slide, Rise, Tilt, Roll, Twist
- Example: bp_idx1=3, bp_idx2=4 → All values identical

**Conclusion**: Step parameter **calculations** are correct and match legacy exactly.

### Calculation Correctness

✅ **Step parameter calculations are identical**

- Same algorithms used
- Same mathematical formulas
- Same reference frame calculations
- Differences are only in pair selection and processing approach

---

## Summary

| Aspect | Modern | Legacy | Status |
|--------|--------|--------|--------|
| **Calculation accuracy** | ✅ Correct | ✅ Correct | ✅ Match exactly |
| **Value matching** | ✅ | ✅ | ✅ 20/20 pairs match |
| **Duplex processing** | Single set | Multiple duplexes | ⚠️ Different approach |
| **Pair selection** | Modern find_pair | Legacy find_pair | ⚠️ Different pairs |
| **Record count** | 20 (1H4S) | 48 (1H4S) | ⚠️ Different counts |
| **JSON format** | Flatter | Nested | ⚠️ Different format |
| **Base pair indices** | 3-23 | 1-25 | ⚠️ Different ranges |

### Key Takeaways

1. ✅ **Step parameter calculations are correct** - Values match exactly when pairs align
2. ⚠️ **Different pair selections** - Legacy and modern find_pair select different pairs
3. ⚠️ **Different processing** - Legacy processes multiple duplexes, modern processes once
4. ✅ **Both implementations are valid** - Differences are by design, not bugs

### Recommendations

1. **For comparison**: Use same input file (same base pairs) for fair comparison
2. **For validation**: Compare values when base pair indices match (proves calculation correctness)
3. **For production**: Both implementations are correct - choose based on use case

---

## Related Documentation

- [STEP_PARAMETERS_IMPLEMENTATION.md](STEP_PARAMETERS_IMPLEMENTATION.md) - Implementation details
- [STEP_PARAMETERS_STATUS.md](STEP_PARAMETERS_STATUS.md) - Current status
- [STEP_PARAMETERS_SUMMARY.md](STEP_PARAMETERS_SUMMARY.md) - Complete summary
- [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) - Overall project status

---

*Documented: 2025-11-29*

