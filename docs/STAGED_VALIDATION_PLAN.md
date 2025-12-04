# Staged Validation Plan - Step by Step to 100% Accuracy

**Goal**: Validate modern code matches legacy in clear stages  
**Strategy**: Only generate JSON for current stage to save space  
**Current Stage**: Need to determine based on what's already validated

---

## All JSON Comparison Types (10 Total)

### Stage 1: Core Data Parsing ✅ **COMPLETE - DO NOT REGENERATE**
**Focus**: Verify we're reading PDB files correctly

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `pdb_atoms` | `pdb_atoms/` | Atom parsing from PDB | ✅ **COMPLETE** - 3602/3602 PDBs validated (100% pass rate) |

**Validation Results** (December 2, 2025):
- ✅ All 3602 fast PDBs tested and validated
- ✅ 100% success rate - all atoms match legacy exactly
- ✅ Atom indices match correctly (matched by atom_idx, not position)
- ✅ Results saved to `data/atoms_test_results.json`

**Note**: Atoms JSON generation is complete. No need to regenerate atoms JSON for future testing.

**JSON Files Needed**:
- `data/json_legacy/atoms/<PDB>.json`
- `data/json/atoms/<PDB>.json`

**Validation Command**:
```bash
python3 scripts/compare_json.py atoms <PDB_ID>
```

---

### Stage 2: Reference Frame Calculation ✅ COMPLETE
**Focus**: Verify frame calculations (rotation matrices, origins)
**Status**: 100% match achieved AND legacy dependency removed from ALL Stage 2 generation (residue_indices, ls_fitting, base_frame_calc). Verified December 4, 2025.

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `base_frame_calc` | `base_frame_calc/` | Frame calculation metadata (template matching, RMS) | ✅ DONE |
| `frame_calc` / `ref_frame` | `frame_calc/` or `ref_frame/` | Reference frames (3x3 rotation + origin) | ✅ DONE |

**JSON Files Needed**:
- `data/json_legacy/base_frame_calc/<PDB>.json`
- `data/json/base_frame_calc/<PDB>.json`
- `data/json_legacy/frame_calc/<PDB>.json` OR `ref_frame/<PDB>.json`
- `data/json/frame_calc/<PDB>.json` OR `ref_frame/<PDB>.json`

**Validation Command**:
```bash
python3 scripts/compare_json.py frames <PDB_ID>
```

---

### Stage 3: Base Pair Detection - Geometric Checks 🔄 CURRENT?
**Focus**: Verify geometric distance/angle calculations

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `distance_checks` | `distance_checks/` | dorg, dNN, plane_angle, d_v, overlap_area | ⏳ TODO |

**JSON Files Needed**:
- `data/json_legacy/distance_checks/<PDB>.json`
- `data/json/distance_checks/<PDB>.json`

**Validation Command**:
```bash
# Need to add this comparison type to compare_json.py
python3 scripts/compare_json.py distance-checks <PDB_ID>
```

**What This Tests**:
- `dorg` - Origin distance between bases
- `dNN` - N-N distance (N9 for purines, N1 for pyrimidines)
- `plane_angle` - Angle between base planes
- `d_v` - Vertical distance
- `overlap_area` - Base overlap calculation

---

### Stage 4: Base Pair Detection - Hydrogen Bonds ⏳ TODO
**Focus**: Verify H-bond detection and conflict resolution

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `hbond_list` | `hbond_list/` | H-bond detection for each pair tested | ⏳ TODO |

**JSON Files Needed**:
- `data/json_legacy/hbond_list/<PDB>.json`
- `data/json/hbond_list/<PDB>.json`

**Validation Command**:
```bash
# Need to add this comparison type to compare_json.py
python3 scripts/compare_json.py hbonds <PDB_ID>
```

**What This Tests**:
- Initial H-bond detection (distance + angle checks)
- Conflict resolution (one atom = one H-bond)
- Final H-bond assignments

---

### Stage 5: Base Pair Detection - Validation Results ⏳ TODO
**Focus**: Verify validation logic (pass/fail for each tested pair)

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `pair_validation` | `pair_validation/` | is_valid flag and quality scores for all tested pairs | ⏳ TODO |

**JSON Files Needed**:
- `data/json_legacy/pair_validation/<PDB>.json`
- `data/json/pair_validation/<PDB>.json`

**Validation Command**:
```bash
# Need to add this comparison type to compare_json.py
python3 scripts/compare_json.py validation <PDB_ID>
```

**What This Tests**:
- Which pairs pass geometric thresholds
- Which pairs fail (and why)
- Quality score calculations

---

### Stage 6: Base Pair Selection ⏳ TODO
**Focus**: Verify final selected pairs (THE MOST CRITICAL OUTPUT)

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `find_bestpair_selection` | `find_bestpair_selection/` | Final selected base pairs from greedy algorithm | ⏳ TODO |
| `base_pair` | `base_pair/` | Base pair records (only for selected pairs) | ⏳ TODO |

**JSON Files Needed**:
- `data/json_legacy/find_bestpair_selection/<PDB>.json`
- `data/json/find_bestpair_selection/<PDB>.json`
- `data/json_legacy/base_pair/<PDB>.json`
- `data/json/base_pair/<PDB>.json`

**Validation Command**:
```bash
# Likely already exists
python3 scripts/compare_json.py pairs <PDB_ID>
```

**What This Tests**:
- Greedy selection algorithm
- Mutual best pair matching
- Final pair list matches legacy exactly

**⭐ THIS IS THE PRIMARY OUTPUT - MUST BE 100% MATCH**

---

### Stage 7: Step Parameters ⏳ TODO
**Focus**: Verify step parameter calculations

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `bpstep_params` | `bpstep_params/` | Shift, Slide, Rise, Tilt, Roll, Twist | ⏳ TODO |

**JSON Files Needed**:
- `data/json_legacy/bpstep_params/<PDB>.json`
- `data/json/bpstep_params/<PDB>.json`

**Validation Command**:
```bash
python3 scripts/compare_json.py steps <PDB_ID>
```

**What This Tests**:
- 6 step parameters for consecutive base pairs
- Midpoint frame calculations
- Transformation matrices

---

### Stage 8: Helical Parameters ⏳ TODO
**Focus**: Verify helical parameter calculations

| Record Type | Directory | What It Tests | Status |
|-------------|-----------|---------------|--------|
| `helical_params` | `helical_params/` | x_displacement, y_displacement, rise, inclination, tip, twist | ⏳ TODO |

**JSON Files Needed**:
- `data/json_legacy/helical_params/<PDB>.json`
- `data/json/helical_params/<PDB>.json`

**Validation Command**:
```bash
# Need to add this comparison type to compare_json.py
python3 scripts/compare_json.py helical <PDB_ID>
```

**What This Tests**:
- Helical axis calculations
- Alternative parameter set
- Helical rise, twist, etc.

---

## Stage-by-Stage Approach

### ✅ COMPLETED STAGES

#### Stage 1: Atom Parsing
- **Status**: ✅ COMPLETE
- **What we validated**: PDB file parsing matches legacy
- **JSON generated**: `atoms/` directory

#### Stage 2: Reference Frames
- **Status**: ✅ COMPLETE
- **What we validated**: Frame calculations match legacy
- **JSON generated**: `base_frame_calc/`, `frame_calc/` or `ref_frame/`

### 🔄 NEXT STAGES TO DO (In Order)

#### Stage 3: Distance Checks
- **Status**: ⏳ TODO NEXT
- **What to validate**: Geometric measurements (dorg, dNN, angles, overlap)
- **JSON to generate**: `distance_checks/` only
- **Why this first**: Foundation for all pair detection

**Action Plan**:
1. Generate modern JSON for distance_checks only (small files)
2. Compare with legacy distance_checks
3. Fix any mismatches
4. Move to Stage 4

#### Stage 4: H-Bond Detection
- **Status**: ⏳ TODO AFTER Stage 3
- **What to validate**: H-bond detection and conflict resolution
- **JSON to generate**: `hbond_list/` only
- **Why this next**: Required for pair validation

#### Stage 5: Pair Validation
- **Status**: ⏳ TODO AFTER Stage 4
- **What to validate**: Which pairs pass/fail validation
- **JSON to generate**: `pair_validation/` only
- **Why this next**: Required for selection algorithm

#### Stage 6: Final Pair Selection ⭐
- **Status**: ⏳ TODO AFTER Stage 5
- **What to validate**: FINAL SELECTED PAIRS (most critical!)
- **JSON to generate**: `find_bestpair_selection/`, `base_pair/`
- **Why this matters**: This is THE primary output

#### Stage 7: Step Parameters
- **Status**: ⏳ TODO AFTER Stage 6
- **What to validate**: Step parameters calculation
- **JSON to generate**: `bpstep_params/` only
- **Why this next**: Depends on final pairs being correct

#### Stage 8: Helical Parameters
- **Status**: ⏳ TODO AFTER Stage 7
- **What to validate**: Helical parameters calculation
- **JSON to generate**: `helical_params/` only
- **Why this last**: Final validation

---

## How to Work Stage-by-Stage

### For Each Stage:

1. **Check What JSON Already Exists**
   ```bash
   # Check legacy
   ls data/json_legacy/<stage_directory>/
   
   # Check modern
   ls data/json/<stage_directory>/
   ```

2. **Generate ONLY This Stage's Modern JSON**
   ```bash
   # Modify generate_modern_json to output only specific record types
   # OR use appropriate flags
   ./build/generate_modern_json --output-types=<stage_types> data/pdb/<PDB>.pdb data/json/
   ```

3. **Compare This Stage Only**
   ```bash
   python3 scripts/compare_json.py <stage_command> <PDB_ID> --verbose
   ```

4. **Fix Issues**
   - Investigate mismatches
   - Update modern code
   - Rebuild
   - Regenerate this stage's JSON only
   - Re-compare

5. **Validate on Test Set**
   ```bash
   python3 scripts/compare_json.py <stage_command> --test-set 10
   ```

6. **Move to Next Stage**
   - Only after 100% match on current stage
   - Don't generate all JSON types at once

---

## Space-Saving Strategy

### Current Issue
If we generate ALL 10 JSON types for ALL PDBs:
- 1737 PDBs × 10 types × ~2 files (legacy + modern) = ~35,000 files
- Each file could be 1KB-1MB depending on PDB size
- Total: Could be several GB of JSON

### Solution: Stage-by-Stage Generation

1. **Only generate current stage's JSON**
2. **Validate and fix**
3. **Delete modern JSON after validation passes** (keep legacy as reference)
4. **Move to next stage**

### Which JSON to Keep

**Always Keep** (Reference):
- All `data/json_legacy/*` - This is our reference, never delete

**Generate as Needed** (Modern):
- Generate `data/json/<stage>/` only when working on that stage
- After stage validates 100%, can optionally delete modern JSON
- Regenerate if needed for debugging

---

## Current Status Check Needed

We need to determine:

1. **What JSON files currently exist?**
   ```bash
   # Check all directories in legacy
   ls -la data/json_legacy/
   
   # Check all directories in modern
   ls -la data/json/
   ```

2. **What validation has already passed?**
   - Stage 1 (atoms): ✅ DONE
   - Stage 2 (frames): ✅ DONE
   - Stage 3-8: Need to check what exists

3. **What should we validate next?**
   - Based on what exists
   - Following the stage order

---

## Immediate Next Steps

### Step 1: Inventory Current State
```bash
# See what JSON directories exist
ls data/json_legacy/
ls data/json/

# Count files in each
for dir in data/json_legacy/*/; do 
  echo "$(basename $dir): $(ls $dir 2>/dev/null | wc -l) files"
done
```

### Step 2: Identify Next Stage
Based on inventory, determine which stage to tackle next

### Step 3: Generate Only That Stage's JSON
Modify tools to output only the record types for that stage

### Step 4: Validate Stage
Run comparisons for only that stage

### Step 5: Repeat
Move through stages 3-8 in order

---

## Summary: The Staged Approach

```
Stage 1: Atoms          ✅ DONE → validated
Stage 2: Frames         ✅ DONE → validated
Stage 3: Distance       ⏳ NEXT → generate & validate
Stage 4: H-bonds        ⏳ TODO → then this
Stage 5: Validation     ⏳ TODO → then this
Stage 6: Selection ⭐   ⏳ TODO → then this (CRITICAL)
Stage 7: Steps          ⏳ TODO → then this
Stage 8: Helical        ⏳ TODO → finally this

100% Match Achieved! 🎯
```

**Key Principle**: One stage at a time, only generate what you need, validate before moving forward.

