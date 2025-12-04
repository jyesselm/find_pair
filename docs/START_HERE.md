# START HERE - Clean Slate Staged Validation

**Date**: December 2, 2025  
**Goal**: 100% accuracy with legacy indexes across all validation types

**Current Status**:
- ✅ **Stage 1 (Atoms)**: COMPLETE - 3602/3602 PDBs validated (100% pass rate)
- 🔄 **Stage 2 (Residue Indices)**: In progress
- ⏳ **Stages 3-8**: Pending  
**Approach**: Batch-by-batch, stage-by-stage, on-demand generation with aggressive cleanup

---

## 📚 **Complete Documentation Set**

1. **[START_HERE.md](START_HERE.md)** ← You are here
2. **[CLEAN_SLATE_VALIDATION_PLAN.md](CLEAN_SLATE_VALIDATION_PLAN.md)** - Complete strategy
3. **[EXECUTION_CHECKLIST.md](EXECUTION_CHECKLIST.md)** - Step-by-step implementation
4. **[STAGED_VALIDATION_PLAN.md](STAGED_VALIDATION_PLAN.md)** - 8-stage breakdown

---

## 🎯 **What We're Doing**

### The Goal
Validate that modern C++ code produces **identical output** to legacy code for:
- **4,123 valid PDBs**
- **8 validation stages** (atoms → frames → distances → hbonds → validation → selection → steps → helical)
- **100% accuracy** on all stages

### Why Clean Slate?
- Existing JSON has count mismatches (unclear what was tested)
- Need systematic validation from scratch
- Old approach uses ~1.6GB disk space (too much!)

### New Approach: On-Demand + Aggressive Cleanup
1. **Generate on-demand**: Only create JSON for current batch/stage
2. **Validate immediately**: Compare legacy vs modern
3. **Delete aggressively**: Remove BOTH legacy AND modern after successful validation
4. **Keep failures only**: Retain JSON for failed PDBs (debugging)
5. **Track in result files**: All outcomes recorded in `data/validation_results/*.json`

**Space**: ~82MB (keep only ~5% failures) vs ~1.6GB (keep everything)

---

## 📋 **The 8 Stages**

| Stage | Name | Production Code Tested | JSON Types | Priority |
|-------|------|----------------------|------------|----------|
| 1 | Atoms | `PdbParser::parse()` | `pdb_atoms` | Medium |
| 2 | Frames | `BaseFrameCalculator` | `base_frame_calc`, `frame_calc` | **HIGH** |
| 3 | Distances | `BasePairValidator` (distances) | `distance_checks` | **HIGH** |
| 4 | H-Bonds | `HydrogenBondFinder` | `hbond_list` | **HIGH** |
| 5 | Validation | `BasePairValidator` (complete) | `pair_validation` | **HIGH** |
| 6 | Selection | `BasePairFinder` (greedy algorithm) | `find_bestpair_selection`, `base_pair` | **CRITICAL** ⭐ |
| 7 | Steps | `ParameterCalculator` (steps) | `bpstep_params` | Medium |
| 8 | Helical | `ParameterCalculator` (helical) | `helical_params` | Medium |

**Focus**: Stages 2-6 are most critical (base pair detection pipeline)

---

## ✅ **What's Already Done**

- [x] Identified 4,123 valid PDBs
- [x] Created validation framework (`scripts/stage_by_stage_validation.py`)
- [x] Baseline validation (confirmed we need to start fresh)
- [x] Documented complete plan
- [x] Mapped production code to each stage
- [x] Ready to execute

---

## 🔄 **What Needs to Be Done**

### Implementation (Code Changes)

**Priority 1**: Modify `tools/generate_modern_json.cpp`
- [ ] Add `--output-types` flag for selective JSON output
- [ ] Add `--output-stage` shorthand (e.g., `--output-stage=frames`)
- [ ] Test on 1 PDB to verify it works
- [ ] Rebuild: `make release`

**Priority 2**: Create `scripts/batch_validation_workflow.py`
- [ ] Batch processing (100 PDBs at a time)
- [ ] Generate legacy JSON for batch
- [ ] Generate modern JSON for batch
- [ ] Validate batch
- [ ] Cleanup (delete successes, keep failures)
- [ ] Update result JSON incrementally
- [ ] Test on small batch (10 PDBs)

### Execution (Run Validation)

**Stage 2: Reference Frames** (Start here - most important foundation)
- [ ] Batch 1: PDBs 1-100
- [ ] Batch 2: PDBs 101-200
- [ ] ... continue
- [ ] Batch 42: PDBs 4101-4123
- [ ] **Target: 100% pass rate**

**Stages 3-6**: Repeat batch process
- [ ] Stage 3: Distance Checks → 100%
- [ ] Stage 4: H-Bonds → 100%
- [ ] Stage 5: Validation → 100%
- [ ] Stage 6: Selection ⭐ → 100% (CRITICAL MILESTONE)

**Stages 7-8**: Step/Helical Parameters
- [ ] Stage 7: Steps → 100%
- [ ] Stage 8: Helical → 100%

---

## 🚀 **Next Immediate Actions**

### Action 1: Clean Slate ⏳ DO NOW
```bash
# Delete everything
rm -rf data/json/
rm -rf data/json_legacy/
rm -rf data/validation_results/

# Create fresh
mkdir -p data/json data/json_legacy data/validation_results
```

### Action 2: Modify generate_modern_json ⏳ DO NEXT
Add `--output-types` flag to selectively output only specific JSON types

### Action 3: Create Batch Script ⏳ DO NEXT
Create `scripts/batch_validation_workflow.py` for automated batch processing

### Action 4: Start Validation ⏳ THEN
Begin Stage 2 (frames) validation in batches

---

## 📊 **Success Criteria**

### Per Batch
- ✅ Legacy JSON generated for batch
- ✅ Modern JSON generated for batch
- ✅ Validation complete
- ✅ Successes deleted (space saved)
- ✅ Failures documented in result JSON

### Per Stage
- ✅ All batches complete (4,123 PDBs)
- ✅ 100% pass rate
- ✅ All failures investigated and fixed
- ✅ Production code verified correct

### Final
- ✅ All 8 stages at 100%
- ✅ Result JSON files document everything
- ✅ Minimal disk usage (~82MB for failures)
- ✅ **100% accuracy with legacy achieved!** 🎯

---

## 📁 **What Will Exist After Completion**

```
data/
├── pdb/                          # Input PDB files (ALWAYS KEEP)
├── json/                         # Modern JSON (only failures ~5%)
├── json_legacy/                  # Legacy JSON (only failures ~5%)
└── validation_results/           # Result JSON files (8 files, ~100KB each)
    ├── stage1_atoms_results.json
    ├── stage2_frames_results.json
    ├── stage3_distances_results.json
    ├── stage4_hbonds_results.json
    ├── stage5_validation_results.json
    ├── stage6_selection_results.json      ⭐ MOST CRITICAL
    ├── stage7_steps_results.json
    └── stage8_helical_results.json
```

Each result file contains:
```json
{
  "stage_name": "Reference Frames",
  "total_pdbs_tested": 4123,
  "passed": ["1EHZ", "100D", ...],           // List of successful PDBs
  "failed": [                                 // Details of failures
    {
      "pdb_id": "7EH2",
      "issue": "Count mismatch: legacy=88, modern=48",
      "details": {...}
    }
  ],
  "summary": {
    "passed_count": 4100,
    "failed_count": 23,
    "pass_rate": 99.44
  }
}
```

---

## 🔍 **Ensuring Production Code is Tested**

We're NOT running test mocks - we're running **actual production executables**:

### Legacy Code
```bash
cd org
./build/bin/find_pair_analyze ../data/pdb/1EHZ.pdb
# → Calls legacy C code
# → Outputs JSON to data/json_legacy/
```

### Modern Code
```bash
./build/generate_modern_json data/pdb/1EHZ.pdb data/json/
# → Calls production C++ code:
#   - PdbParser::parse() → reads PDB
#   - BaseFrameCalculator::calculate_frame() → calculates frames
#   - BasePairFinder::find_pairs() → finds pairs
#   - JsonWriter::write_*() → outputs JSON
# → Outputs JSON to data/json/
```

### Validation
```python
# Compare actual JSON outputs
legacy = json.load(open('data/json_legacy/base_frame_calc/1EHZ.json'))
modern = json.load(open('data/json/base_frame_calc/1EHZ.json'))
# → If different: production code has a bug
# → If identical: production code is correct
```

**This tests the REAL algorithm, not test code!**

---

## 📖 **Documentation Reference**

**Planning**:
- [CLEAN_SLATE_VALIDATION_PLAN.md](CLEAN_SLATE_VALIDATION_PLAN.md) - Complete strategy
- [STAGED_VALIDATION_PLAN.md](STAGED_VALIDATION_PLAN.md) - 8-stage breakdown

**Execution**:
- [EXECUTION_CHECKLIST.md](EXECUTION_CHECKLIST.md) - Step-by-step tasks
- **[START_HERE.md](START_HERE.md)** - This file (overview)

**Results** (after execution):
- `data/validation_results/stage*_results.json` - Validation outcomes

---

## ⚡ **Quick Start**

### 1. Clean Slate
```bash
rm -rf data/json/ data/json_legacy/ data/validation_results/
mkdir -p data/json data/json_legacy data/validation_results
```

### 2. Modify Code
Add `--output-types` flag to `generate_modern_json.cpp`

### 3. Create Batch Script
Implement `scripts/batch_validation_workflow.py`

### 4. Run Validation
```bash
python3 scripts/batch_validation_workflow.py --stage stage2_frames
```

### 5. Achieve 100%!
Repeat for all 8 stages

---

**Ready to begin!** See [EXECUTION_CHECKLIST.md](EXECUTION_CHECKLIST.md) for detailed steps.

