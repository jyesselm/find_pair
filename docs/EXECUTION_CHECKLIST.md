# Execution Checklist - Clean Slate Validation

**Plan**: See [CLEAN_SLATE_VALIDATION_PLAN.md](CLEAN_SLATE_VALIDATION_PLAN.md)  
**Status**: Ready to execute

---

## ✅ Completed

- [x] Identified 4,123 valid PDBs in `valid_pdbs.json`
- [x] Created stage-by-stage validation framework
- [x] Baseline validation (0% - as expected, no modern JSON)
- [x] Documented complete plan
- [x] Mapped production code to each stage

---

## 🔄 Next: Implementation Steps

### Step 1: Clean Slate - Delete ALL JSON ⏳ DO NEXT

```bash
# AGGRESSIVE cleanup - delete everything, regenerate on-demand
cd /Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2

# Delete ALL JSON (can regenerate)
rm -rf data/json/
rm -rf data/json_legacy/
rm -rf data/validation_results/

# Create fresh directories
mkdir -p data/json
mkdir -p data/json_legacy  
mkdir -p data/validation_results

echo "✅ Clean slate complete!"
echo "Will generate JSON on-demand per batch/stage"
```

**Checklist**:
- [ ] Deleted `data/json/`
- [ ] Deleted `data/json_legacy/`
- [ ] Deleted `data/validation_results/`
- [ ] Created fresh directories
- [ ] Confirmed PDB files still exist: `ls data/pdb/ | wc -l`

---

### Step 2: Modify generate_modern_json for Selective Output ⏳ TODO

**Approach**: Add `--output-types` flag

**Changes needed** in `tools/generate_modern_json.cpp`:

```cpp
// 1. Add to argument parsing:
std::string output_types_str = "all";  // Default: all types
// Parse --output-types=pdb_atoms,base_frame_calc

// 2. Create enabled set:
std::set<std::string> enabled_types;
if (output_types_str == "all") {
    // Empty set = all enabled
} else {
    // Parse comma-separated list
    // enabled_types = {"pdb_atoms", "base_frame_calc"}
}

// 3. Conditional JSON writing:
if (enabled_types.empty() || enabled_types.count("pdb_atoms")) {
    writer.write_pdb_atoms(structure);
}

if (enabled_types.empty() || enabled_types.count("base_frame_calc")) {
    // Write base_frame_calc
}

// etc. for all types
```

**Test**:
```bash
# Rebuild
make release

# Test selective output
./build/generate_modern_json \
  --output-types=base_frame_calc,frame_calc \
  data/pdb/1EHZ.pdb \
  data/json/

# Verify only frame JSONs created
ls data/json/*/1EHZ.json
# Should show ONLY base_frame_calc and frame_calc
```

**Checklist**:
- [ ] Modified `generate_modern_json.cpp` with `--output-types` flag
- [ ] Rebuilt: `make release`
- [ ] Tested on 1EHZ (verified selective output works)
- [ ] Cleaned up test: `rm data/json/*/1EHZ.json`

---

### Step 3: Create Batch Validation Script ⏳ TODO

**New script**: `scripts/batch_validation_workflow.py`

**Features**:
- Generate modern JSON for batch (specific stage only)
- Validate batch
- Delete successes, keep failures
- Update result JSON incrementally
- Resume from checkpoint

**Test**:
```bash
# Test on small batch
python3 scripts/batch_validation_workflow.py \
  --stage stage2_frames \
  --batch-start 0 \
  --batch-size 10 \
  --cleanup

# Verify:
# - Modern JSON generated for 10 PDBs
# - Validation ran
# - Successes deleted
# - Failures kept
# - Results updated
```

**Checklist**:
- [ ] Created `scripts/batch_validation_workflow.py`
- [ ] Implemented batch generation
- [ ] Implemented validation
- [ ] Implemented cleanup logic
- [ ] Tested on 10 PDBs
- [ ] Verified cleanup works

---

### Step 4: Run Stage 2 Validation (Frames) ⏳ TODO

**Process**: Validate frames in batches

```bash
# Batch 1 (PDBs 0-99)
python3 scripts/batch_validation_workflow.py \
  --stage stage2_frames \
  --batch-start 0 \
  --batch-size 100 \
  --cleanup

# Check results
python3 -c "
import json
with open('data/validation_results/stage2_frames_results.json') as f:
    r = json.load(f)
print(f'Passed: {r[\"summary\"][\"passed_count\"]}')
print(f'Failed: {r[\"summary\"][\"failed_count\"]}')
if r['failed']:
    print('First failure:', r['failed'][0])
"

# If failures: investigate and fix
# Re-run batch if needed

# Continue with remaining batches...
```

**Target**: Process all ~42 batches for Stage 2

**Checklist**:
- [ ] Batch 1 complete (0-99)
- [ ] Batch 2 complete (100-199)
- [ ] ... (track progress)
- [ ] All 42 batches complete
- [x] Stage 2: 100% match rate achieved
- [x] Stage 2: Legacy dependency removed from residue_indices (December 4, 2025)
- [x] Stage 2: Legacy dependency removed from ls_fitting (December 4, 2025)
- [x] Stage 2: Legacy dependency removed from base_frame_calc (December 4, 2025)

---

### Step 5: Run Stages 3-8 ⏳ TODO

Repeat for each stage:

**Stage 3: Distance Checks**
- [ ] All batches complete
- [ ] 100% pass rate

**Stage 4: H-Bonds**
- [ ] All batches complete
- [ ] 100% pass rate

**Stage 5: Validation**
- [ ] All batches complete
- [ ] 100% pass rate

**Stage 6: Selection** ⭐
- [ ] All batches complete
- [ ] 100% pass rate
- [ ] **CRITICAL MILESTONE**

**Stage 7: Step Parameters**
- [ ] Generated modern JSON
- [ ] 100% pass rate

**Stage 8: Helical Parameters**
- [ ] Generated modern JSON
- [ ] 100% pass rate

---

## Final Validation

### After All Stages Complete

```bash
# Run full validation on all PDBs
python3 scripts/stage_by_stage_validation.py

# Should show:
# All 8 stages: 100% pass rate
```

**Success Criteria**:
```
Atoms               : 4123 passed, 0 failed, 100.00% pass rate
Reference Frames    : 4123 passed, 0 failed, 100.00% pass rate
Distance Checks     : 4123 passed, 0 failed, 100.00% pass rate
H-Bonds             : 4123 passed, 0 failed, 100.00% pass rate
Pair Validation     : 4123 passed, 0 failed, 100.00% pass rate
Pair Selection      : 4123 passed, 0 failed, 100.00% pass rate
Step Parameters     : 4123 passed, 0 failed, 100.00% pass rate
Helical Parameters  : 4123 passed, 0 failed, 100.00% pass rate
```

**Then**: 🎉 **100% ACCURACY ACHIEVED!**

---

## Current Status

**Completed**:
- ✅ Validation framework created
- ✅ Plan documented
- ✅ Ready to execute

**Next immediate action**:
- ⏳ Archive existing data (Step 1)
- ⏳ Modify generate_modern_json (Step 2)
- ⏳ Create batch script (Step 3)
- ⏳ Start Stage 2 validation (Step 4)

---

## Quick Reference Commands

### Archive Data
```bash
mkdir -p data/archive/clean_slate_$(date +%Y%m%d_%H%M%S)
mv data/json data/archive/clean_slate_$(date +%Y%m%d_%H%M%S)/json_modern
mv data/validation_results data/archive/clean_slate_$(date +%Y%m%d_%H%M%S)/validation_results
mkdir -p data/json data/validation_results
```

### Test Selective Output (After Code Changes)
```bash
./build/generate_modern_json --output-types=base_frame_calc,frame_calc data/pdb/1EHZ.pdb data/json/
```

### Run Batch Validation (After Script Created)
```bash
python3 scripts/batch_validation_workflow.py --stage stage2_frames --batch-start 0 --batch-size 100 --cleanup
```

### Check Progress
```bash
cat data/validation_results/stage2_frames_results.json | python3 -m json.tool | head -50
```

---

**Ready to execute!** Start with Step 1 when you're ready to begin.

