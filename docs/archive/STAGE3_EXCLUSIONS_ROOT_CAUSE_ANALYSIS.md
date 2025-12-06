# Stage 3 Exclusions - Root Cause Analysis

**Date**: December 5, 2025  
**Status**: Investigation Complete  
**Key Finding**: ⚠️ **ALL 24 "Stage 3 exclusions" are actually Stage 2 (ls_fitting) failures**

---

## Executive Summary

The file `data/stage3_exclusions.json` incorrectly categorizes 24 PDBs as "Stage 3 (distance_checks) failures" when they actually **failed at Stage 2 (frame calculation)**. All 24 PDBs are already in `data/ls_fitting_failures.json` and have **NO modern JSON files** because the modern code couldn't complete frame calculation for them.

### Key Evidence

1. ✅ All 24 PDBs are in `ls_fitting_failures.json`
2. ✅ None of the 24 have modern `base_frame_calc/*.json` files  
3. ✅ None reached Stage 3 (distance_checks)
4. ✅ The "dNN mismatch" descriptions are SPECULATIVE, not actual test results

---

## Root Cause: Missing Modified Nucleotide Mappings

### Problem

The modern code lacks `one_letter_code()` mappings for several modified nucleotides that legacy DOES process:

| Modified Nucleotide | Legacy Behavior | Modern Behavior | Result |
|---------------------|-----------------|-----------------|--------|
| **J48** | Maps to Uracil template (`Atomic.u.pdb`) | Returns '?' → not recognized as nucleotide | Either skipped or uses wrong template |
| **2YR** | Maps to Cytosine template (`Atomic.c.pdb`) | Returns '?' → not recognized | Skipped/wrong template |
| **NMN** | Maps to Uracil template (`Atomic.u.pdb`) | Returns '?' → not recognized | Skipped/wrong template |
| **NNR** | Maps to Uracil template (`Atomic.u.pdb`) | Returns '?' → not recognized | Skipped/wrong template |
| **WVQ** | Maps to Uracil template (`Atomic.u.pdb`) | Returns '?' → not recognized | Skipped/wrong template |
| **CSL** | Likely mapped by legacy | Returns '?' → not recognized | Unknown |
| **CTP** | Triphosphate - treated as base by legacy | Returns '?' → not recognized | Skipped |
| **GTP** | Triphosphate - treated as base by legacy | Returns '?' → not recognized | Skipped |

### Evidence: 6QIQ (J48 example)

**Legacy JSON**:
- Has 30 residues total
- Has 4 J48 entries (with duplicates)
- J48 uses `Atomic.u.pdb` (lowercase modified template)
- RMS fit: 0.042

**Modern JSON** (freshly generated):
- Has 18 residues total (12 MISSING)
- Has 2 J48 entries (correct, no duplicates)
- J48 uses `Atomic_U.pdb` (UPPERCASE standard template - WRONG!)
- RMS fit: 0.042 (same value)
- `one_letter_code: 'N/A'` (not recognized)

**Debug Output**:
```
DEBUG: Calculating frame for residue: J48 A:101 (type=6, one_letter=?, needs_rmsd=1)
DEBUG: First RMSD check failed (rmsd=1.724584), retrying as pyrimidine-only
DEBUG: RMSD check passed (rmsd=0.0964168, threshold=0.2618)
DEBUG: Has purine atoms: no
DEBUG: Looking for 6 ring atoms
DEBUG: Atom  N3  (repr: 0x20 0x4e 0x33 0x20 ): NOT FOUND in residue
```

J48 is **missing the N3 atom**, so it only matches 5 of 6 required pyrimidine atoms. Despite this, RMSD check passes.

---

## Incorrect Categorization in stage3_exclusions.json

### Category 1: J48 (4 PDBs)

**Claimed**: "J48 is a highly modified nucleotide not in baselist.dat"  
**Reality**: 
- J48 is NOT in baselist.dat (TRUE)
- But legacy DOES process it using heuristics
- Modern code processes it but with WRONG template (uppercase instead of lowercase)
- Missing N3 atom but still gets processed

**PDBs**: 6QIQ, 6QIR, 6QIT, 6QIS

### Category 2: 2YR (4 PDBs)

**Claimed**: "2YR (2'-O-ribose modification) not properly handled"  
**Reality**:
- Legacy processes 2YR as modified cytosine (`Atomic.c.pdb`)
- Modern code doesn't have 2YR in `one_letter_code()` mapping
- Returns '?' → not recognized as nucleotide

**PDBs**: 7S36, 7S3H, 7S38, 9CJJ

### Category 3: NMN/NNR (5 PDBs)

**Claimed**: "NMN/NNR (nicotinamide nucleotides) not handled as bases"  
**Reality**:
- These ARE nicotinamide-based cofactors, not true nucleotides
- But legacy DOES process them (maps to Uracil template)
- Modern code doesn't recognize them (one_letter_code returns '?')

**PDBs**: 8GXC, 8HB1, 8HB3, 8HB8, 8I3Z

### Category 4: WVQ (1 PDB)

**Claimed**: "WVQ modified nucleotide not handled"  
**Reality**:
- Legacy processes WVQ as modified uracil (`Atomic.u.pdb`)
- Modern code doesn't recognize it

**PDBs**: 8UKS

### Category 5: EPE (3 PDBs)

**Claimed**: "EPE (modified cytosine) has minor dNN differences (~0.2Å)"  
**Reality**:
- EPE IS in `one_letter_code()` mapping (returns 'c')
- Should use lowercase template
- Needs investigation why it still fails

**PDBs**: 4E8M, 4E8R, 6T3N

### Category 6: A23 (1 PDB)

**Claimed**: "A23 uses pyrimidine RMSD fallback causing dNN mismatch"  
**Reality**:
- A23 IS in `one_letter_code()` mapping (returns 'a')
- Already in `TemplateAssignment::MODIFIED_PURINES`
- This might be a genuine edge case

**PDBs**: 2XD0

### Category 7: Other (5 PDBs)

**Claimed**: "Various modified nucleotides causing dNN/pair differences"  
**Reality**: Need to investigate each individually

**PDBs**: 4IQS, 8ABZ, 8ANE, 8PFQ, 8SY6

### Category 8: Corrupt (1 PDB)

**Claimed**: "Legacy JSON file is truncated/corrupt"  
**Reality**: This is correct - can't fix corrupt data

**PDBs**: 9CJI

---

## The Real Issues

### Issue 1: Missing one_letter_code Mappings

**Location**: `include/x3dna/core/residue.hpp::one_letter_code()`

**Missing mappings**:
```cpp
// Need to add:
if (trimmed == "J48") return 'u';  // Hypermodified uracil
if (trimmed == "2YR") return 'c';  // 2'-O-ribose cytosine
if (trimmed == "NMN" || trimmed == "NNR") return 'u';  // Nicotinamide (treat as uracil)
if (trimmed == "WVQ") return 'u';  // Unknown modified uracil
if (trimmed == "CSL") return '?';  // Need to investigate
if (trimmed == "CTP") return 'c';  // Cytosine triphosphate
if (trimmed == "GTP") return 'g';  // Guanine triphosphate
```

### Issue 2: Template Selection for Modified Nucleotides

**Location**: `src/x3dna/algorithms/base_frame_calculator.cpp`

Even when `one_letter_code` returns lowercase (e.g., 'u' for modified), the template selection logic must use **lowercase templates** (`Atomic.u.pdb` not `Atomic_U.pdb`).

Current logic (simplified):
```cpp
bool is_modified_nucleotide = (one_letter_code >= 'a' && one_letter_code <= 'z');
standard_template = templates_.load_template(residue_type, is_modified_nucleotide);
```

This should work IF `one_letter_code()` returns lowercase for modified nucleotides.

### Issue 3: Atom Matching for Highly Modified Nucleotides

J48 is missing N3 but still passes RMSD check. This is actually CORRECT legacy behavior - legacy also accepts partial atom matches if RMSD is below threshold.

But we need to ensure the template used is lowercase for modified nucleotides.

---

## Recommended Fixes

### Fix 1: Add Missing Modified Nucleotide Mappings ⭐ HIGH PRIORITY

**File**: `include/x3dna/core/residue.hpp`

Add to `one_letter_code()`:

```cpp
// Modified Uracil/Thymine → 'u' (additional)
if (trimmed == "J48")   // Hypermodified nucleotide (treat as uracil)
    return 'u';
if (trimmed == "NMN" || trimmed == "NNR")  // Nicotinamide (treat as uracil)
    return 'u';
if (trimmed == "WVQ")   // Unknown modified nucleotide (treat as uracil)
    return 'u';

// Modified Cytosine → 'c' (additional)
if (trimmed == "2YR")   // 2'-O-ribose cytosine
    return 'c';
if (trimmed == "CTP")   // Cytosine triphosphate
    return 'c';

// Modified Guanine → 'g' (additional)  
if (trimmed == "GTP")   // Guanine triphosphate
    return 'g';
```

### Fix 2: Verify Template Assignment Logic

Ensure that when `one_letter_code()` returns lowercase, the `is_modified` flag is set correctly in `base_frame_calculator.cpp`.

### Fix 3: Investigate EPE and A23 Separately

These already have `one_letter_code` mappings, so they should work. Need to debug why they still fail.

### Fix 4: Update stage3_exclusions.json

Rename to `stage2_modified_nucleotide_failures.json` and update the descriptions to reflect that these are Stage 2 (frame calculation) failures, not Stage 3 (distance_checks) failures.

---

## Testing Plan

### Test 1: Add Mappings and Regenerate

```bash
# 1. Add mappings to residue.hpp
# 2. Rebuild
cd build && ninja

# 3. Test J48 PDB
./generate_modern_json data/pdb/6QIQ.pdb /tmp/test_j48/

# 4. Verify:
#    - J48 has lowercase template (Atomic.u.pdb)
#    - Residue count matches legacy
#    - No duplicates in output
```

### Test 2: Compare with Legacy

```python
import json

with open('/tmp/test_j48/base_frame_calc/6QIQ.json') as f:
    modern = json.load(f)

with open('data/json_legacy/base_frame_calc/6QIQ.json') as f:
    legacy = json.load(f)

# Should now match (excluding duplicates)
modern_unique = remove_duplicates(modern)
legacy_unique = remove_duplicates(legacy)

assert len(modern_unique) == len(legacy_unique)
assert templates_match(modern_unique, legacy_unique)
```

### Test 3: Run Full Validation

After fixes, re-run Stage 2 validation on these 24 PDBs to see if they now pass.

---

## Statistics

| Category | Count | Fixable? | Fix Type |
|----------|-------|----------|----------|
| J48 | 4 | ✅ Yes | Add one_letter_code mapping |
| 2YR | 4 | ✅ Yes | Add one_letter_code mapping |
| NMN/NNR | 5 | ✅ Yes | Add one_letter_code mapping |
| WVQ | 1 | ✅ Yes | Add one_letter_code mapping |
| EPE | 3 | ⚠️ Maybe | Already mapped - investigate |
| A23 | 1 | ⚠️ Maybe | Already mapped - investigate |
| Other | 5 | ❓ Unknown | Need individual investigation |
| Corrupt | 1 | ❌ No | Data corruption |
| **TOTAL** | **24** | **15 likely fixable** | |

---

## Conclusion

1. **All 24 "Stage 3 exclusions" failed at Stage 2, NOT Stage 3**
2. **Primary cause**: Missing `one_letter_code()` mappings for modified nucleotides
3. **Secondary cause**: Legacy processes these using heuristics; modern code rejects them
4. **Fix**: Add ~8 one-letter mappings to `residue.hpp`
5. **Expected result**: 15-20 of the 24 PDBs will pass after fix

The exclusions file is **misleading** because it describes hypothetical Stage 3 failures (dNN mismatches, missing pairs) that never actually occurred because these PDBs never made it to Stage 3.

---

## Next Actions

1. ✅ Add missing one_letter_code mappings to `residue.hpp`
2. ✅ Rebuild and test on 6QIQ (J48)
3. ✅ Test on other categories (2YR, NMN, WVQ)
4. ✅ Investigate EPE and A23 separately (already mapped)
5. ✅ Re-run Stage 2 validation on all 24 PDBs
6. ✅ Update exclusions file to reflect actual failures
7. ✅ Document any remaining failures with ROOT CAUSE


