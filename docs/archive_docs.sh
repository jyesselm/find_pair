#!/bin/bash
# Archive unnecessary documentation files

cd "$(dirname "$0")"

# Superseded documents
mv FIND_PAIR_VALIDATION.md archive/historical/superseded/ 2>/dev/null
mv FIND_PAIR_VALIDATION_EXPANDED.md archive/historical/superseded/ 2>/dev/null
mv BASE_PAIR_RECORDING_FIX.md archive/historical/superseded/ 2>/dev/null
mv BASE_PAIR_RECORDING_SUCCESS.md archive/historical/superseded/ 2>/dev/null
mv BASE_PAIR_RECORDING_COMPREHENSIVE.md archive/historical/superseded/ 2>/dev/null
mv BASE_PAIR_RECORDING_EXPANDED_TEST.md archive/historical/superseded/ 2>/dev/null
mv BASE_PAIR_RECORDING_VERIFICATION_LOG.md archive/historical/superseded/ 2>/dev/null

# Old test results
mv FIX_INDICES_TEST_RESULTS.md archive/historical/test_results/ 2>/dev/null
mv PHASE1_TEST_RESULTS.md archive/historical/test_results/ 2>/dev/null
mv BROAD_TESTING_RESULTS.md archive/historical/test_results/ 2>/dev/null
mv TEST_RESULTS_BP_TYPE_ID_FIX.md archive/historical/test_results/ 2>/dev/null
mv FRAME_COMPARISON_RESULTS.md archive/historical/test_results/ 2>/dev/null

# Historical investigations
mv 1T0K_VALIDATION_ANALYSIS.md archive/historical/investigations/ 2>/dev/null
mv 3G8T_analysis_updated.md archive/historical/investigations/ 2>/dev/null
mv BASE_PAIR_RECORD_ANALYSIS.md archive/historical/investigations/ 2>/dev/null
mv BASE_PAIR_FRAMES_100_PERCENT_PLAN.md archive/historical/investigations/ 2>/dev/null
mv BASE_PAIR_FRAMES_PROGRESS.md archive/historical/investigations/ 2>/dev/null
mv DORG_CALCULATION_INVESTIGATION.md archive/historical/investigations/ 2>/dev/null
mv FRAME_ORIGIN_MISMATCH.md archive/historical/investigations/ 2>/dev/null
mv HBOND_INVESTIGATION.md archive/historical/investigations/ 2>/dev/null
mv INVESTIGATION_SUMMARY.md archive/historical/investigations/ 2>/dev/null
mv ISSUES_FOUND_SUMMARY.md archive/historical/investigations/ 2>/dev/null
mv IS_NUCLEOTIDE_BUG_ANALYSIS.md archive/historical/investigations/ 2>/dev/null
mv OFF_BY_ONE_ANALYSIS.md archive/historical/investigations/ 2>/dev/null
mv OVERLAP_CALCULATION.md archive/historical/investigations/ 2>/dev/null
mv PAIR_REJECTION_ANALYSIS.md archive/historical/investigations/ 2>/dev/null
mv QUALITY_SCORE_DIFFERENCES.md archive/historical/investigations/ 2>/dev/null
mv RESIDUE_ORDERING_COMPARISON.md archive/historical/investigations/ 2>/dev/null
mv VALIDATION_DIFFERENCES.md archive/historical/investigations/ 2>/dev/null

# Old plans/status
mv 100_PERCENT_MATCH_PLAN.md archive/historical/old_plans/ 2>/dev/null
mv ACTION_PLAN_NEXT_STEPS.md archive/historical/old_plans/ 2>/dev/null
mv NEXT_STEPS.md archive/historical/old_plans/ 2>/dev/null
mv NEXT_STEPS_AFTER_FIX_INDICES.md archive/historical/old_plans/ 2>/dev/null

# Completed fixes (info now in main docs)
mv BP_TYPE_ID_BUG_FIX.md archive/historical/completed_fixes/ 2>/dev/null
mv BP_TYPE_ID_REQUIREMENTS.md archive/historical/completed_fixes/ 2>/dev/null
mv FIX_INDICES_IMPLEMENTATION_SUMMARY.md archive/historical/completed_fixes/ 2>/dev/null
mv RESIDUE_GROUPING_FIX.md archive/historical/completed_fixes/ 2>/dev/null
mv RESIDUE_INDEXING_COMPLETE.md archive/historical/completed_fixes/ 2>/dev/null
mv RESIDUE_INDEXING_ISSUE.md archive/historical/completed_fixes/ 2>/dev/null
mv RESIDUE_INDEXING_QUICK_REFERENCE.md archive/historical/completed_fixes/ 2>/dev/null
mv RESIDUE_INDEXING_STATUS.md archive/historical/completed_fixes/ 2>/dev/null

# Redundant/duplicate
mv COMPARISON_CONFIG.md archive/historical/redundant/ 2>/dev/null
mv COMPARISON_WORKFLOW.md archive/historical/redundant/ 2>/dev/null
mv JSON_COMPARISON_GUIDE.md archive/historical/redundant/ 2>/dev/null
mv JSON_GENERATION_SYSTEM.md archive/historical/redundant/ 2>/dev/null
mv JSON_RECORD_ORDER.md archive/historical/redundant/ 2>/dev/null
mv JSON_SERIALIZATION.md archive/historical/redundant/ 2>/dev/null
mv KNOWN_DIFFERENCES_CATALOG.md archive/historical/redundant/ 2>/dev/null
mv LEGACY_JSON_GENERATION.md archive/historical/redundant/ 2>/dev/null
mv LEGACY_MODE_DESIGN.md archive/historical/redundant/ 2>/dev/null
mv LEGACY_STEP_PARAMETER_ANALYSIS.md archive/historical/redundant/ 2>/dev/null
mv OOP_CLASS_HIERARCHY.md archive/historical/redundant/ 2>/dev/null
mv PARAMETER_MAPPING_VERIFICATION.md archive/historical/redundant/ 2>/dev/null
mv PDB_PROPERTIES_MATCHING_APPROACH.md archive/historical/redundant/ 2>/dev/null
mv PDB_PROPERTIES_MATCHING_SUMMARY.md archive/historical/redundant/ 2>/dev/null
mv PROTOCOL_LEGACY_COMPARISON.md archive/historical/redundant/ 2>/dev/null
mv PROTOCOL_LEGACY_DETAILED_COMPARISON.md archive/historical/redundant/ 2>/dev/null
mv SIMPLIFIED_DEVELOPMENT.md archive/historical/redundant/ 2>/dev/null
mv STORAGE_OPTIMIZATION_STRATEGY.md archive/historical/redundant/ 2>/dev/null
mv TOOLS_TO_CREATE.md archive/historical/redundant/ 2>/dev/null
mv DEBUG_JSON_GUIDE.md archive/historical/redundant/ 2>/dev/null

echo "Archived unnecessary documentation files"
echo "Essential docs remain in docs/ root"
echo "Archived docs moved to docs/archive/historical/"

