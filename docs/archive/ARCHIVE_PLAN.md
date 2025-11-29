# Documentation Archive Plan

**Date**: 2025-11-28  
**Purpose**: Reduce documentation from 117 files to ~15 essential files

---

## Current Situation

- **Total MD files**: 117
- **Essential**: ~10 files
- **Reference**: ~5 files  
- **Can archive**: ~100 files

---

## Essential Documents (Keep in root)

1. **100_PERCENT_MATCH_STATUS.md** - Main status document
2. **FIND_PAIR_VALIDATION_COMPREHENSIVE.md** - Validation results
3. **BASE_PAIR_RECORDING.md** - Base pair fix documentation
4. **TESTING_GUIDE.md** - Testing workflows
5. **FIX_INDICES_OPTION.md** - Residue indexing fix
6. **COMPARISON_COVERAGE.md** - What's being compared
7. **JSON_DATA_TYPES_AND_COMPARISONS.md** - JSON reference
8. **BUILD_INSTRUCTIONS.md** - Build guide
9. **QUICK_START.md** - Quick start guide
10. **README.md** - Documentation index

---

## Reference Documents (Keep in root)

1. **ALGORITHM_CRITICAL_GUIDE.md** - Algorithm details
2. **DEBUGGING_STRATEGIES.md** - Debugging guide
3. **CONFIGURATION_OPTIONS.md** - Configuration reference
4. **DATA_STRUCTURE.md** - Data structure reference
5. **MODERNIZATION_PLAN.md** - Modernization strategy

---

## Documents to Archive

### Superseded Documents
- FIND_PAIR_VALIDATION.md
- FIND_PAIR_VALIDATION_EXPANDED.md
- BASE_PAIR_RECORDING_FIX.md
- BASE_PAIR_RECORDING_SUCCESS.md
- BASE_PAIR_RECORDING_COMPREHENSIVE.md
- BASE_PAIR_RECORDING_EXPANDED_TEST.md
- BASE_PAIR_RECORDING_VERIFICATION_LOG.md

### Old Test Results
- FIX_INDICES_TEST_RESULTS.md
- PHASE1_TEST_RESULTS.md
- BROAD_TESTING_RESULTS.md
- TEST_RESULTS_BP_TYPE_ID_FIX.md
- FRAME_COMPARISON_RESULTS.md

### Historical Investigations
- 1T0K_VALIDATION_ANALYSIS.md
- 3G8T_analysis_updated.md
- BASE_PAIR_RECORD_ANALYSIS.md
- BASE_PAIR_FRAMES_100_PERCENT_PLAN.md
- BASE_PAIR_FRAMES_PROGRESS.md
- DORG_CALCULATION_INVESTIGATION.md
- FRAME_ORIGIN_MISMATCH.md
- HBOND_INVESTIGATION.md
- INVESTIGATION_SUMMARY.md
- ISSUES_FOUND_SUMMARY.md
- IS_NUCLEOTIDE_BUG_ANALYSIS.md
- OFF_BY_ONE_ANALYSIS.md
- OVERLAP_CALCULATION.md
- PAIR_REJECTION_ANALYSIS.md
- QUALITY_SCORE_DIFFERENCES.md
- RESIDUE_ORDERING_COMPARISON.md
- VALIDATION_DIFFERENCES.md

### Old Plans/Status
- 100_PERCENT_MATCH_PLAN.md
- ACTION_PLAN_NEXT_STEPS.md
- NEXT_STEPS.md
- NEXT_STEPS_AFTER_FIX_INDICES.md

### Completed Fixes (info in main docs)
- BP_TYPE_ID_BUG_FIX.md
- BP_TYPE_ID_REQUIREMENTS.md
- FIX_INDICES_IMPLEMENTATION_SUMMARY.md
- RESIDUE_GROUPING_FIX.md
- RESIDUE_INDEXING_COMPLETE.md
- RESIDUE_INDEXING_ISSUE.md
- RESIDUE_INDEXING_QUICK_REFERENCE.md
- RESIDUE_INDEXING_STATUS.md

### Redundant/Duplicate
- COMPARISON_CONFIG.md
- COMPARISON_WORKFLOW.md
- JSON_COMPARISON_GUIDE.md
- JSON_GENERATION_SYSTEM.md
- JSON_RECORD_ORDER.md
- JSON_SERIALIZATION.md
- KNOWN_DIFFERENCES_CATALOG.md
- LEGACY_JSON_GENERATION.md
- LEGACY_MODE_DESIGN.md
- LEGACY_STEP_PARAMETER_ANALYSIS.md
- OOP_CLASS_HIERARCHY.md
- PARAMETER_MAPPING_VERIFICATION.md
- PDB_PROPERTIES_MATCHING_APPROACH.md
- PDB_PROPERTIES_MATCHING_SUMMARY.md
- PROTOCOL_LEGACY_COMPARISON.md
- PROTOCOL_LEGACY_DETAILED_COMPARISON.md
- SIMPLIFIED_DEVELOPMENT.md
- STORAGE_OPTIMIZATION_STRATEGY.md
- TOOLS_TO_CREATE.md
- DEBUG_JSON_GUIDE.md

---

## Archive Structure

```
docs/
├── README.md (updated index)
├── Essential docs (10 files)
├── Reference docs (5 files)
└── archive/
    └── historical/
        ├── superseded/
        ├── test_results/
        ├── investigations/
        ├── old_plans/
        └── completed_fixes/
```

---

## Action Plan

1. Create archive directory structure
2. Move files to appropriate archive subdirectories
3. Update README.md to reflect new structure
4. Add note in archived files pointing to current docs

---

*This plan will reduce documentation from 117 files to ~15 essential files.*

