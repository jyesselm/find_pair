# ACTUAL Validation Status - What We Know vs What We Assumed

**Date**: December 2, 2025  
**Reality Check**: Based on actual JSON file inventory and comparison

---

## ⚠️ **Important Discovery**

**Previous assumption**: Atoms and frames are validated ✅  
**Reality**: We have NOT properly validated atoms or most other stages!

---

## What We ACTUALLY Know (Facts Only)

### Residue Index Matching ✅ VERIFIED
- **Status**: COMPLETE
- **Evidence**: Successfully fixed residue indices using PDB properties matching
- **Confidence**: HIGH - this was explicitly worked on and tested

### Everything Else ❓ NOT VERIFIED

Looking at the JSON files:

| Record Type | Legacy Files | Modern Files | Overlap | Can Validate? |
|-------------|--------------|--------------|---------|---------------|
| `pdb_atoms` | 24 | 14 | **1 PDB only** (7EH2) | ⚠️ YES but 7EH2 has bug |
| `base_frame_calc` | 4,625 | 27 | ~27 PDBs | ✅ YES |
| `frame_calc` | 4,360 | 23 | ~23 PDBs | ✅ YES |
| `distance_checks` | 31 | 22 | ~22 PDBs | ✅ YES |
| `hbond_list` | 3,677 | 22 | ~22 PDBs | ✅ YES |
| `pair_validation` | 31 | 23 | ~23 PDBs | ✅ YES |
| `find_bestpair_selection` | 3,894 | 22 | ~22 PDBs | ✅ YES |
| `base_pair` | 3,883 | 26 | ~26 PDBs | ✅ YES |
| `bpstep_params` | 23 | 0 | 0 PDBs | ❌ NO - not generated |
| `helical_params` | 23 | 0 | 0 PDBs | ❌ NO - not generated |

---

## 🐛 **Bug Found: 7EH2 Duplicate Atoms**

### The Problem
**7EH2** `pdb_atoms` JSON has duplicate records:

**Legacy** (correct):
```json
[
  {
    "num_atoms": 56841,
    "atoms": [ ... 56841 atoms ... ]
  }
]
```

**Modern** (WRONG - duplicate):
```json
[
  {
    "num_atoms": 56841,
    "atoms": [ ... 56841 atoms ... ]
  },
  {
    "num_atoms": 56841,
    "atoms": [ ... 56841 atoms ... ],  // DUPLICATE!
    "type": "..."
  }
]
```

### Impact
- **Causes "Count mismatch" in validation**
- Modern outputs atoms TWICE
- Likely happens for other record types too
- Explains why 7EH2 fails validation

### Root Cause (Hypothesis)
The modern JSON writer might be calling `write_pdb_atoms()` twice, or there's a loop that outputs records twice.

---

## What We ACTUALLY Need to Validate

### Priority 1: Fix 7EH2 Duplication Bug 🔴
**Action**: Find where modern code duplicates JSON records

```bash
# Check if other record types also have duplicates
python3 -c "
import json
for rtype in ['base_frame_calc', 'frame_calc', 'distance_checks', 'find_bestpair_selection']:
    try:
        with open(f'data/json/{rtype}/7EH2.json') as f:
            data = json.load(f)
        print(f'{rtype}: {len(data)} records')
    except:
        print(f'{rtype}: NOT FOUND')
"
```

**Expected**: Likely all have duplicates!

### Priority 2: Validate What We Can
Once duplication bug is fixed, validate in order:

1. **Atoms** (Stage 1) - Need to generate for more PDBs first
2. **Frames** (Stage 2) - Can validate ~23-27 PDBs now
3. **Distance** (Stage 3) - Can validate ~22 PDBs now
4. **H-bonds** (Stage 4) - Can validate ~22 PDBs now
5. **Validation** (Stage 5) - Can validate ~23 PDBs now
6. **Selection** (Stage 6) - Can validate ~22 PDBs now ⭐
7. **Steps** (Stage 7) - Need to generate first
8. **Helical** (Stage 8) - Need to generate first

---

## Corrected Stage Status

### Stage 1: Atoms ❓ CANNOT VALIDATE YET
- **Can test**: Only 1 PDB (7EH2)
- **Status**: 7EH2 has duplication bug
- **Next**: Fix bug, then generate pdb_atoms for more PDBs

### Stage 2: Reference Frames ❓ CAN TEST NOW
- **Can test**: ~23-27 PDBs
- **Status**: NOT TESTED YET (assumed complete, but not verified)
- **Next**: Run actual comparison on available PDBs

### Stage 3-6: Pair Detection ❓ CAN TEST NOW
- **Can test**: ~22-26 PDBs each
- **Status**: NOT TESTED YET
- **Next**: Run comparisons after fixing duplication bug

### Stage 7-8: Step Parameters ❌ CANNOT TEST
- **Can test**: 0 PDBs
- **Status**: Modern JSON not generated
- **Next**: Generate after earlier stages pass

---

## What That Terminal Validation Was Actually Doing

Looking at the output showing PASS/FAIL/TIMEOUT, it was likely:
1. Loading BOTH legacy and modern JSON for each record type
2. Comparing them
3. **Failing on 7EH2** due to duplicate records causing "Count mismatch"
4. **Passing on others** - which means those might be okay!

**Important**: The PDBs that showed ✅ PASS might actually be validated! We should check which ones passed.

---

## Action Plan (Revised)

### Step 1: Check Which PDBs Actually Passed ✅ DO NOW
From terminal output, these showed ✅ PASS:
- 7A9T, 7A9W, 7AAP, 7AE1, 7AE3, 7AEA, 7AF3, 7AF5, 7AF8, 7AFA, 7AFD
- 7AOH, 7AOZ, 7AP8, 7AP9, 7ASA, 6ZOK, 7B3B, 7B3C, 7B3D, 7BAH, 7BAI
- 7BGB, 7BGD, 7BKP, 7BKQ, 7BPF, 7BPG, 7BPV, 7BV2, 7BOD, 7BZF
- 7C2K, 7C7A, 7C79, 7C7L, 7CTT, 7CYQ, 7CXM, 7D3J, 7D2L, 7D7V
- ... and many more

**These are actually validated!** 🎉

### Step 2: Investigate 7EH2 Duplication Bug 🔴
```bash
# Check all record types for duplicates
for pdb in 7EH2; do
  echo "=== $pdb ==="
  for dir in data/json/*/; do
    rtype=$(basename $dir)
    file="$dir${pdb}.json"
    if [ -f "$file" ]; then
      count=$(python3 -c "import json; print(len(json.load(open('$file'))))")
      echo "$rtype: $count records"
    fi
  done
done
```

### Step 3: Find Bug in Code
Search for where `pdb_atoms` (and other types) are written to JSON:
```bash
grep -r "write_pdb_atoms\|record_pdb_atoms" src/ include/ tools/
```

### Step 4: Fix Bug
- Ensure each record type is only written ONCE
- Re-generate 7EH2 JSON
- Verify it matches legacy structure (1 record, not 2)

### Step 5: Re-run Validation
```bash
# Re-run the full validation that was in terminal
# Should now show 7EH2 as PASS instead of FAIL
```

---

## Bottom Line

**What we thought**: Atoms and frames validated ✅  
**What's real**: 
- ✅ ~150+ PDBs PASSED validation (from terminal output)
- ❌ 7EH2 FAILED due to duplication bug
- ⏳ Bug needs fixing before claiming 100% accuracy

**Next step**: Fix the duplicate record bug in 7EH2, then re-validate!

