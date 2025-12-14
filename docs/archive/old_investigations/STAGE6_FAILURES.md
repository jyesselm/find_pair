# Stage 6 (Pair Validation) Failures Report

## Summary

After multiple fixes to the pair validation logic, Stage 6 now passes for most PDBs.

| Metric | Count |
|--------|-------|
| Total PDBs tested | 3602 |
| Passed | 3573 |
| Failed | 21 |
| Skipped | 8 |
| **Pass Rate** | 99.2% |

## Root Cause Analysis

The quality_score formula is:
```
quality_score = dorg + 2.0 * d_v + plane_angle / 20.0
                + adjust_pair_quality(hbonds)   // h-bond adjustment
                - 2.0 (if bp_type_id == 2)      // WC/wobble adjustment
```

### H-Bond Adjustment Logic (`adjust_pair_quality`)

```cpp
// Count good h-bonds where:
// 1. type != '*' (non-standard h-bonds are skipped)
// 2. distance in [2.5, 3.5] Angstroms

if (num_good_hb >= 2) return -3.0;
else return -num_good_hb;  // Returns 0, -1, or -2
```

### Remaining Error Categories

| Category | Count | Description |
|----------|-------|-------------|
| dNN | 8 | Modified nucleotide N1/N9 atom lookup |
| quality_score diff<3 | 6 | Minor h-bond counting differences |
| quality_score diff>=3 | 5 | H-bond conflict resolution differences |
| is_valid | 2 | Validation result differences |

## Bug Fixes Applied

### Fix 1: H-bond type filtering

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Issue:** Modern was counting both type `'-'` and `' '` h-bonds, but legacy's `hb_info_string` excludes type `' '`.

```cpp
// Only count '-' type h-bonds (skip '*' and ' ')
if (hbond.type != '-') {
    continue;
}
```

### Fix 2: H-bond distance precision

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Issue:** Legacy uses `%4.2f` format in `hb_info_string`, rounding distances to 2 decimals.

```cpp
// Round to 2 decimals to match legacy %4.2f format
double rounded_dist = std::round(hbond.distance * 100.0) / 100.0;
if (rounded_dist >= 2.5 && rounded_dist <= 3.5) {
    num_good_hb++;
}
```

### Fix 3: N1/N9 atom lookup for modified nucleotides

**File:** `src/x3dna/algorithms/base_pair_validator.cpp`

**Issue:** Modern was checking for C8 atom to determine purine, but some modified pyrimidines (e.g., 70U) have C8 as part of their modification.

```cpp
// Use one_letter_code to determine purine/pyrimidine (matches legacy bseq)
char upper_letter = toupper(one_letter);
bool is_purine = (upper_letter == 'A' || upper_letter == 'G' || upper_letter == 'I');
```

### Fix 4: bp_type_id for INOSINE and PSEUDOURIDINE

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Issue:** `get_base_letter_from_type` didn't handle INOSINE and PSEUDOURIDINE, returning '?' instead of 'I' or 'P', which broke WC_LIST matching for C-I pairs.

```cpp
case ResidueType::INOSINE:
    return 'I';
case ResidueType::PSEUDOURIDINE:
    return 'P';
```

## Remaining Issues (~1% failure rate)

The 27 remaining failures are due to:
1. **dNN differences (8):** Complex modified nucleotides where N1/N9 lookup differs
2. **bp_type_id differences (8):** Different WC/wobble pair classification
3. **quality_score differences (11):** Edge cases in h-bond counting

These represent edge cases that may require deeper investigation of the legacy algorithm.
