# Validation Bugs Found and Fixed

**Date**: December 2, 2025  
**Status**: 2 bugs found, 2 bugs fixed

---

## Bug 1: RMSD Check Skipped for Standard Nucleotides ✅ FIXED

### Issue
**PDB**: 2VQF, Residue A:1331 G  
**Symptom**: Count mismatch (modern 1524, legacy 1523)

### Root Cause
Modern code was **skipping RMSD geometry check** for standard nucleotides (A, C, G, T, U).

**File**: `src/x3dna/algorithms/base_frame_calculator.cpp` (line 279)

**Before**:
```cpp
if (has_ring_atoms && needs_rmsd_check) {  // Only modified nucleotides
    auto rmsd_result = check_nt_type_by_rmsd(residue);
    if (!rmsd_result.has_value() || *rmsd_result > 0.2618) {
        return result; // RMSD check failed
    }
}
```

**Problem**: `needs_rmsd_check = false` for standard nucleotides → skipped check

### The Failure Case
**A:1331 G** in 2VQF:
- Has all ring atoms ✅
- **RMS fit = 0.268139** (bad geometry/distorted)
- **Threshold = 0.2618**  
- **Should be filtered**: 0.268139 > 0.2618
- **Legacy**: Filters it out (runs RMSD check on ALL) ✅
- **Modern** (before fix): Accepted it (skipped RMSD check) ❌

### The Fix
```cpp
if (has_ring_atoms) {  // Check ALL nucleotides, not just modified
    auto rmsd_result = check_nt_type_by_rmsd(residue);
    if (!rmsd_result.has_value() || *rmsd_result > 0.2618) {
        return result; // RMSD check failed
    }
}
```

### Validation After Fix
- 2VQF: Modern 1523 = Legacy 1523 ✅ PASS
- All subsequent PDBs continue passing

**Impact**: Critical - ensures modern applies SAME geometry filtering as legacy

---

## Bug 2: Stale Legacy JSON Files ✅ IDENTIFIED

### Issue
**PDB**: 6V9Q  
**Symptom**: Count mismatch (modern 60, legacy 14)

### Root Cause
**Stale legacy JSON** from old version of legacy code or incomplete run.

### Investigation
**Modern**: 60 nucleotides in chain K (K:1-60, with K:61 filtered)
- All have good RMS values (max 0.0181 << 0.2618)
- All have required ring atoms
- All properly identified

**Legacy JSON (old)**: Only 14 residues
- K:42, 44, 46-51, 55-60 (sparse, many gaps)
- Missing K:1-41, 43, 45, 52-54

### The Fix
**Regenerated legacy JSON** with current legacy code:
```bash
python scripts/rebuild_json.py regenerate 6V9Q --legacy-only
```

**Result**: Legacy now has 74 records (not 14)
- Includes K:1-60 (plus some duplicates)
- Modern 60 = Legacy 60 ✅ PASS

### Solution
**Before running full validation**: Regenerate all legacy JSON with current code
```bash
python scripts/rebuild_json.py regenerate --legacy-only
```

**Impact**: Ensures comparing against current legacy behavior, not old versions

---

## Validation Statistics (After Fixes)

### Testing Progress
- **Total tested**: 2,500+ PDBs
- **✅ PASS**: 2,310+ (92.4%)
- **❌ FAIL**: 1 (6V9Q - fixed but CSV not updated yet)
- **⏭️  SKIP**: 187 (7.5% - no legacy JSON)

### Bugs Found
1. ✅ **RMSD check bug**: Fixed in code
2. ✅ **Stale JSON bug**: Fixed by regeneration

### Confidence After Fixes
With both bugs fixed:
- RMSD filtering now matches legacy for ALL nucleotides
- Legacy JSON regenerated ensures accurate comparison
- Expect remaining ~1,600 PDBs to continue passing

---

## Recommendations

### Before Full-Scale Validation
1. **Regenerate ALL legacy JSON** with current legacy code
   ```bash
   python scripts/rebuild_json.py regenerate --legacy-only
   ```

2. **Clear validation CSV** to start fresh
   ```bash
   rm data/index_validation_status.csv
   ```

3. **Run comprehensive validation** with fixed code
   ```bash
   python scripts/validate_all_indices.py --batch-size 100 --threads 20 --clean
   ```

### Expected Outcome
With both fixes applied and fresh legacy JSON:
- **Target**: >95% PASS rate
- **Remaining failures**: Edge cases or genuine differences to investigate

---

## Next Steps

1. ✅ Commit RMSD fix
2. ✅ Document stale JSON issue
3. ⏳ Complete current validation run
4. ⏳ Regenerate all legacy JSON
5. ⏳ Run final comprehensive validation
6. → Proceed to Week 2 (unified comparison framework)

---

**Status**: Validation continuing with fixed code. Will report when complete!

