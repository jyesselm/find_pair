# Stage 6 (Pair Validation) Debugging Notes

This document captures the detailed debugging process for Stage 6 validation, including root causes, fixes applied, and remaining known issues.

## Overview

Stage 6 validates base pair detection by comparing `pair_validation` JSON output between legacy and modern implementations. The comparison checks:
- `is_valid`: Whether the pair passes validation
- `bp_type_id`: Base pair type (-1, 0, 1, 2)
- `calculated_values`: dorg, dNN, plane_angle, d_v, quality_score

## Final Status

| Metric | Count |
|--------|-------|
| Total PDBs tested | 3602 |
| Passed | 3573 |
| Failed | 21 |
| Skipped | 8 |
| **Pass Rate** | **99.2%** |

---

## Bugs Fixed

### Fix 1: H-bond Type Filtering in `adjust_pair_quality`

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Problem:** The `adjust_pair_quality` function was counting wrong h-bond types for quality score adjustment.

**Root Cause Analysis:**
1. Legacy's `get_hbond_ij` creates `hb_info_string` which EXCLUDES type `' '` h-bonds
2. Legacy's `hb_numlist` parses this string, so it only sees h-bonds that weren't `' '`
3. Legacy's `adjust_pairQuality` then skips type `'*'` via `num_list[k][0]`
4. **Net result:** Only type `'-'` h-bonds are counted

**Original Code (incorrect):**
```cpp
if (hbond.type == '*') {
    continue;  // Only skip '*', counts both '-' and ' '
}
```

**Fixed Code:**
```cpp
if (hbond.type != '-') {
    continue;  // Only count '-', skip both '*' and ' '
}
```

**Impact:** Improved pass rate from ~39% to ~88%

---

### Fix 2: H-bond Distance Precision

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Problem:** H-bond distances at the [2.5, 3.5] range boundary were being counted differently.

**Root Cause Analysis:**
1. Legacy uses `sprintf(stmp, " %s%c%s %4.2f", ...)` to format h-bond info
2. The `%4.2f` format rounds distances to 2 decimal places
3. Example: 2.4995 becomes "2.50" in the string
4. `hb_numlist` parses this string, so it sees 2.50 (IN range)
5. Modern used full precision, so 2.4995 was NOT in range

**Fix:**
```cpp
// Round to 2 decimals to match legacy %4.2f format
double rounded_dist = std::round(hbond.distance * 100.0) / 100.0;
if (rounded_dist >= 2.5 && rounded_dist <= 3.5) {
    num_good_hb++;
}
```

**Impact:** Improved pass rate from ~88% to ~96%

---

### Fix 3: N1/N9 Atom Lookup for Modified Nucleotides

**File:** `src/x3dna/algorithms/base_pair_validator.cpp`

**Problem:** Modified pyrimidines with C8 atoms (e.g., 70U) were incorrectly classified as purines.

**Root Cause Analysis:**
1. Modern was checking for C8 atom presence to determine purine
2. Some modified pyrimidines (like 70U) have C8 as part of their modification
3. This caused modern to look for N9 instead of N1
4. Legacy uses `bseq` (one-letter code) to determine purine/pyrimidine

**Original Code (incorrect):**
```cpp
// Check for N7, C8, or N9 to determine purine
auto n7 = residue.find_atom(" N7 ");
auto c8 = residue.find_atom(" C8 ");
auto n9 = residue.find_atom(" N9 ");
if (n7.has_value() || c8.has_value() || n9.has_value()) {
    is_purine = true;
}
```

**Fixed Code:**
```cpp
// Use one_letter_code to determine purine/pyrimidine (matches legacy bseq)
char upper_letter = toupper(one_letter);
bool is_purine = (upper_letter == 'A' || upper_letter == 'G' || upper_letter == 'I');
```

**Impact:** Improved pass rate from ~96% to ~99%

---

### Fix 4: INOSINE and PSEUDOURIDINE in `get_base_letter_from_type`

**File:** `src/x3dna/algorithms/base_pair_finder.cpp`

**Problem:** C-I (Cytosine-Inosine) pairs weren't getting bp_type_id=2 (Watson-Crick).

**Root Cause Analysis:**
1. `calculate_bp_type_id` uses `get_base_letter_from_type` to get base letters
2. Original function returned '?' for INOSINE and PSEUDOURIDINE
3. WC_LIST contains "CI" and "IC" for Inosine-Cytosine pairs
4. With '?' returned, the pair "C?" didn't match WC_LIST

**Fix:**
```cpp
case ResidueType::INOSINE:
    return 'I';
case ResidueType::PSEUDOURIDINE:
    return 'P';
```

**Impact:** Fixed 6 PDBs with bp_type_id differences (1SAQ, 7WV3, 7WV4, 7WV5, 7WVE, 7WVJ)

---

## Remaining Issues (21 failures, 0.8%)

### Category 1: dNN Differences (8 failures)

**PDBs:** 1Q2R, 1Q2S, 2XD0, 2XDD, 4E8R, 7NQ4, 8OMR, 8TXO

**Root Cause:** Modified nucleotides with unusual purine/pyrimidine classification.

**Example: 9DG (9-deazaguanine) in 1Q2R**

9DG has:
- N7 and C8 (purine-like atoms)
- No N9 (replaced by C9)
- Classified as 'u' (modified uracil) in base_type

Behavior difference:
1. **Legacy:** `residue_ident()` detects N7/C8 → classifies as purine → looks for N9 → not found → falls back to C9
2. **Modern:** Uses base_type='u' → classifies as pyrimidine → uses N1

Result:
- Legacy dNN: 38.884121 (using C9)
- Modern dNN: 41.715711 (using N1)
- Difference: 2.83 Å

**Potential Fix:** Modify `find_n1_n9_position` to use legacy's `residue_ident` logic (check for N7/C8 atoms) instead of relying on `one_letter_code`. This would add complexity and risk introducing other issues.

---

### Category 2: Quality Score Differences (11 failures)

**PDBs with diff<3:** 1QCU, 2OIH, 7EQG, 8UPT, 8UPY, 9K2G  
**PDBs with diff>=3:** 3KTW, 4JNX, 4OQ9, 6T3K, 8K29

**Root Cause:** H-bond conflict resolution algorithm differences.

**Technical Details:**

Legacy's `hb_atompair` function uses a multi-phase conflict resolution:

1. **Phase 1:** Iterative conflict detection
   - For each h-bond, find shortest h-bond for same donor AND same acceptor
   - If both point to same h-bond, mark it as conflict (negate distance)
   - Mark all h-bonds sharing atoms with conflict as "matched" (excluded)
   - Restart from beginning

2. **Phase 2:** Calculate linkage types (lkg_type)
   - Used for filtering in validation phase

3. **Phase 3:** Additional conflict marking (if lkg_type != 18 and distance in [hb_lower, hb_dist2])
   - Note: hb_dist2 defaults to 0.0, so this phase is effectively disabled

**Difference:** In some edge cases, the modern conflict resolution marks different h-bonds as conflicts than legacy. This affects which h-bonds get type `'-'` vs type `' '`, which affects the quality_score adjustment.

**Example: 3KTW pair (102, 145)**
- Legacy: 4 h-bonds, all type `' '` → 0 good h-bonds → adjustment = 0.0
- Modern: 3 h-bonds, 2 type `'-'` → 2 good h-bonds → adjustment = -3.0
- quality_score difference: 3.0

---

### Category 3: is_valid Differences (2 failures)

**PDBs:** 4IQS, 8ANE

**Root Cause:** Edge cases at validation thresholds where legacy and modern compute slightly different values, causing different is_valid results.

---

## Key Legacy Code References

### H-bond Processing Flow

```
get_hbond_ij (cmn_fncs.c)
├── Find all potential h-bonds (within_limits)
├── hb_atompair - Conflict resolution (negates distances for conflicts)
├── validate_hbonds - Assigns types ('-', '*', ' ')
│   └── Only processes NEGATIVE distances (conflicts)
│   └── Positive distances get type ' '
├── Creates hb_info_string (excludes type ' ')
└── Records to JSON (ALL h-bonds including type ' ')

hb_numlist (cmn_fncs.c)
├── Parses hb_info_string (only sees non-' ' h-bonds)
├── Uses %4.2f precision for distances
└── Returns num_list for adjust_pairQuality

adjust_pairQuality (ana_fncs.c)
├── Iterates num_list (from hb_numlist)
├── Skips type '*' (num_list[k][0] == 1)
├── Counts h-bonds with distance in [2.5, 3.5]
└── Returns: -3.0 if >=2, -num_good if <2
```

### N1/N9 Atom Selection

```
glyco_N (cmn_fncs.c)
├── isR (from residue_ident via RY array):
│   ├── True: Look for N9
│   └── False: Look for N1 (or C5 for P/p)
├── If primary atom not found:
│   └── Fallback: Find any atom with '9' or '1' in name
└── If still not found:
    └── Fallback: Find closest base atom to C1'

residue_ident (cmn_fncs.c)
├── Checks for ring atoms: C4, N3, C2, N1, C6, C5, N7, C8, N9
├── kr = count of N7, C8, N9 atoms found
├── If kr > 0 AND rmsd <= NT_CUTOFF:
│   └── Return 1 (purine)
└── Else: Return 0 (pyrimidine)
```

### bp_type_id Calculation

```
check_wc_wobble_pair (ana_fncs.c)
├── Calculate shear, stretch, opening from step parameters
├── If |stretch| > 2.0 OR |opening| > 60: Return (no change)
├── If |shear| in [1.8, 2.8]: bp_type_id = 1 (Wobble)
├── If |shear| <= 1.8 AND bp_type in WC_LIST:
│   └── bp_type_id = 2 (Watson-Crick)
└── WC_LIST: "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
```

---

## Testing Commands

### Run Stage 6 validation on specific PDB:
```bash
python3 -B -c "
import sys
sys.path.insert(0, 'tests_python')
from validation.runner import validate_stage
from pathlib import Path
result = validate_stage(6, ['1ASZ'], Path('data/json'), verbose=True)
print(f'Passed: {result.passed}, Failed: {result.failed}')
"
```

### Regenerate modern JSON and validate:
```bash
./build/generate_modern_json data/pdb/1ASZ.pdb data/json --stage=all --quiet
```

### Compare specific pair values:
```bash
python3 -c "
import json
pdb, pair = '1ASZ', (93, 130)
with open(f'data/json_legacy/pair_validation/{pdb}.json') as f:
    leg = json.load(f)
with open(f'data/json/pair_validation/{pdb}.json') as f:
    mod = json.load(f)
for r in leg:
    if (r['base_i'], r['base_j']) == pair:
        print('Legacy:', r['calculated_values'])
for r in mod:
    if (r['base_i'], r['base_j']) == pair:
        print('Modern:', r['calculated_values'])
"
```

---

## Files Modified

| File | Changes |
|------|---------|
| `src/x3dna/algorithms/base_pair_finder.cpp` | H-bond type filtering, distance rounding, INOSINE/PSEUDOURIDINE handling |
| `src/x3dna/algorithms/base_pair_validator.cpp` | N1/N9 lookup using one_letter_code instead of atom presence |
| `docs/STAGE6_FAILURES.md` | Updated with current status |

---

## Recommendations for Future Work

1. **dNN differences:** Consider implementing `residue_ident`-style logic in `find_n1_n9_position` to handle modified nucleotides that have purine atoms but pyrimidine classification.

2. **Quality score differences:** The h-bond conflict resolution algorithm is complex. Further investigation would require tracing through specific examples with debug output.

3. **Validation threshold edge cases:** These are inherent precision differences that may not be worth fixing.

4. **Overall:** At 99.2% pass rate, the remaining issues are edge cases with complex modified nucleotides. The fixes applied cover the vast majority of real-world use cases.

