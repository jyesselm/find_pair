# Test Set Status - Your 22 PDB Test Set

**Date**: December 2, 2025

---

## Your Current Test Set (22 PDBs)

```
1EHZ    6ZLW    6ZMT    6ZN5    6ZOJ    6ZUO
6ZV6    6ZXD    6ZXE    6ZXG    6ZXH    7EH2
7JQB    7K5I    7ZAG    7ZAH    7ZAI    7ZHG
8C01    8EUY    8P03    8P09
```

**Note**: Some of these showed up in your terminal validation as **TIMEOUTS** (6ZLW, 6ZMT, 6ZN5, 6ZOJ, 6ZUO, 6ZV6, 6ZXD, 6ZXE, 6ZXG, 6ZXH). These are likely large ribosomal structures.

---

## What JSON Exists for This Test Set

### ✅ Already Generated (8 Record Types)

| Record Type | Directory | Status | Ready to Compare |
|-------------|-----------|--------|------------------|
| 1. PDB Atoms | `pdb_atoms/` | ✅ Partial (14 PDBs) | YES |
| 2. Base Frame Calc | `base_frame_calc/` | ✅ Generated (27 PDBs) | YES |
| 3. Frame Calc | `frame_calc/` | ✅ Generated (23 PDBs) | YES |
| 4. Distance Checks | `distance_checks/` | ✅ Generated (22 PDBs) | YES |
| 5. H-Bond List | `hbond_list/` | ✅ Generated (22 PDBs) | YES |
| 6. Pair Validation | `pair_validation/` | ✅ Generated (23 PDBs) | YES |
| 7. Find Best Pair Selection | `find_bestpair_selection/` | ✅ Generated (22 PDBs) | YES ⭐ |
| 8. Base Pair | `base_pair/` | ✅ Generated (26 PDBs) | YES |

### ❌ NOT Generated Yet (2 Record Types)

| Record Type | Directory | Status | Ready to Compare |
|-------------|-----------|--------|------------------|
| 9. Step Parameters | `bpstep_params/` | ❌ NOT generated | NO |
| 10. Helical Parameters | `helical_params/` | ❌ NOT generated | NO |

---

## Staged Validation Plan for Your Test Set

### Stage 1: Atoms ✅ Can Validate Now
**Compare**: 14 PDBs that have both legacy and modern `pdb_atoms`

```bash
# Check which ones
comm -12 \
  <(ls data/json_legacy/pdb_atoms/*.json 2>/dev/null | xargs -n1 basename | sort) \
  <(ls data/json/pdb_atoms/*.json 2>/dev/null | xargs -n1 basename | sort)

# Compare
python3 scripts/compare_json.py atoms <PDB_ID>
```

**Expected**: Should PASS (atom parsing is straightforward)

---

### Stage 2: Reference Frames ✅ Can Validate Now
**Compare**: ~23-27 PDBs that have both legacy and modern frames

```bash
# Compare specific PDB
python3 scripts/compare_json.py frames 1EHZ
python3 scripts/compare_json.py frames 7EH2

# Or compare all that have both
# (need to implement batch comparison)
```

**Expected**: Should PASS (you've already validated residue indices)

**Known Issue**: 7EH2 showed "FAIL - Count mismatch" in terminal - investigate this!

---

### Stage 3: Distance Checks ⏳ Can Validate Now
**Compare**: 22 PDBs

```bash
# Need to add this comparison to compare_json.py
# For now, manually compare one:
python3 -c "
import json
legacy = json.load(open('data/json_legacy/distance_checks/1EHZ.json'))
modern = json.load(open('data/json/distance_checks/1EHZ.json'))
print('Legacy pairs:', len(legacy))
print('Modern pairs:', len(modern))
# ... detailed comparison
"
```

**What to check**:
- `dorg` - origin distance
- `dNN` - N-N distance  
- `plane_angle` - angle between planes
- `d_v` - vertical distance
- `overlap_area` - overlap calculation

**Expected**: Should PASS if frame calculations are correct

---

### Stage 4: H-Bond Detection ⏳ Can Validate Now
**Compare**: 22 PDBs

```bash
# Need to add this comparison to compare_json.py
# Or use existing if it exists:
python3 scripts/compare_json.py hbonds 1EHZ
```

**What to check**:
- H-bond pairs detected
- Donor/acceptor assignments
- Conflict resolution
- Final H-bond lists match

**Expected**: Critical stage - H-bonds affect pair selection

---

### Stage 5: Pair Validation ⏳ Can Validate Now
**Compare**: 23 PDBs

```bash
# Need to add this comparison to compare_json.py
python3 scripts/compare_json.py validation 1EHZ
```

**What to check**:
- `is_valid` flag for each tested pair
- Quality scores
- Which pairs pass geometric thresholds
- Which pairs fail (and why)

**Expected**: Must match legacy exactly

---

### Stage 6: Final Pair Selection ⭐ Can Validate Now
**Compare**: 22 PDBs (find_bestpair_selection) + 26 PDBs (base_pair)

```bash
# This might already exist in compare_json.py
python3 scripts/compare_json.py pairs 1EHZ

# Or specifically:
python3 scripts/compare_json.py compare 1EHZ
```

**What to check**:
- **find_bestpair_selection**: Final list of selected pairs
- **base_pair**: Base pair records for selected pairs
- Must match legacy exactly - this is THE primary output!

**Expected**: **MUST BE 100% MATCH** - This is the most critical validation!

**Known Issue**: 7EH2 showed "FAIL - Count mismatch" - this is likely here!

---

### Stage 7: Step Parameters ❌ NOT Generated Yet
**Status**: Need to generate modern JSON first

**Action**:
```bash
# Generate step parameters for test set
# This requires running analyze_app on each PDB
for pdb in 1EHZ 7EH2 7JQB 7K5I; do
  ./build/find_pair_app --fix-indices data/pdb/${pdb}.pdb /tmp/${pdb}.inp
  ./build/analyze_app /tmp/${pdb}.inp
  # This should create bpstep_params JSON
done
```

**Then compare**:
```bash
python3 scripts/compare_json.py steps <PDB_ID>
```

---

### Stage 8: Helical Parameters ❌ NOT Generated Yet
**Status**: Need to generate modern JSON first (same process as Stage 7)

**Action**: Same as Stage 7 - analyze_app should generate both

**Then compare**:
```bash
# Need to add this comparison to compare_json.py
python3 scripts/compare_json.py helical <PDB_ID>
```

---

## Immediate Action Plan

### 1. Investigate 7EH2 Failure 🔴 HIGH PRIORITY
```bash
# Compare all record types for 7EH2
python3 scripts/compare_json.py compare 7EH2 --verbose

# Check selection specifically
python3 -c "
import json
legacy = json.load(open('data/json_legacy/find_bestpair_selection/7EH2.json'))
modern = json.load(open('data/json/find_bestpair_selection/7EH2.json'))
print(f'Legacy selected pairs: {len(legacy)}')
print(f'Modern selected pairs: {len(modern)}')
if len(legacy) != len(modern):
    print('COUNT MISMATCH - this is the issue!')
"
```

**This is your "Count mismatch" failure from the terminal!**

---

### 2. Validate Stages 3-6 on Working PDBs
Pick a PDB that passed (like 7JQB or 7K5I) and validate each stage:

```bash
# Test a known-good PDB through all stages
PDB=7K5I

# Stage 1: Atoms
python3 scripts/compare_json.py atoms $PDB

# Stage 2: Frames  
python3 scripts/compare_json.py frames $PDB

# Stages 3-6: Need to check if comparison commands exist
# or write simple Python scripts to compare
```

---

### 3. Generate Step/Helical Parameters
Once Stages 1-6 pass, generate step parameters:

```bash
# Generate for a few test PDBs first
for pdb in 1EHZ 7JQB 7K5I; do
  ./build/find_pair_app --fix-indices data/pdb/${pdb}.pdb /tmp/${pdb}.inp
  ./build/analyze_app /tmp/${pdb}.inp
done
```

---

### 4. Validate Step/Helical Parameters
Compare the generated step/helical parameters

---

## Success Criteria

### For 100% Accuracy on Test Set

✅ Stage 1 (Atoms): 100% match on all available PDBs  
✅ Stage 2 (Frames): 100% match achieved, legacy dependency fully removed (residue_indices ✅, ls_fitting ✅, base_frame_calc ✅) - December 4, 2025  
⏳ Stage 3 (Distance): 100% match on 22 PDBs  
⏳ Stage 4 (H-Bonds): 100% match on 22 PDBs  
⏳ Stage 5 (Validation): 100% match on 23 PDBs  
⏳ Stage 6 (Selection): **100% match on 22 PDBs** ⭐ CRITICAL  
⏳ Stage 7 (Steps): 100% match on test set  
⏳ Stage 8 (Helical): 100% match on test set

**Then**: Expand to test_set_100 and repeat!

---

## Notes on Test Set Composition

### Why These PDBs?
Looking at the PDB IDs:
- **1EHZ**: Classic RNA structure (good test case)
- **6Z** series: Likely ribosomal structures (large, complex)
- **7** series: Recent structures
- **8** series: Very recent structures

### Timeouts
These PDBs timed out (>120s):
- 6ZLW, 6ZMT, 6ZN5, 6ZOJ, 6ZUO, 6ZV6, 6ZXD, 6ZXE, 6ZXG, 6ZXH

**Reason**: Likely very large ribosomal structures (1000+ nucleotides)

**Action**: These might be too large for quick validation. Consider:
1. Accept they're slow (they may still be correct)
2. Exclude from test set
3. Optimize performance later

---

## Bottom Line

**You have modern JSON for 22 PDBs across 8 record types.**

**Next steps**:
1. ✅ **Investigate 7EH2** failure (count mismatch in selection)
2. ⏳ Validate stages 3-6 on working PDBs
3. ❌ Generate step/helical parameters
4. ⏳ Validate those
5. 🎯 Achieve 100% on test set
6. 📈 Expand to larger test set

**Don't generate more JSON until you fix 7EH2 and validate what you have!**

