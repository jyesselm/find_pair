# Git Commit Guide

**Date**: Current  
**Purpose**: Identify what files should be committed to git

## ✅ Files That SHOULD Be Committed

### Essential Root Documentation (High Priority)
These are the core project documentation files:

```bash
git add PROJECT_STATUS.md
git add README.md
git add DOCUMENTATION_INDEX.md
git add FINAL_100_PDB_COMPARISON_RESULTS.md
git add 3CME_EDGE_CASE_ANALYSIS.md
git add TIE_BREAKING_FIX_SUMMARY.md
git add TIE_BREAKING_FIX_FINAL_RESULTS.md
git add CLEANUP_PLAN_FILES.md
git add CLEANUP_COMPLETE.md
```

### Active Documentation in `docs/` (High Priority)
These are actively maintained documentation files:

```bash
git add docs/1T0K_VALIDATION_ANALYSIS.md
git add docs/BP_TYPE_ID_BUG_FIX.md
git add docs/BP_TYPE_ID_INVESTIGATION_SUMMARY.md
git add docs/BROAD_TESTING_RESULTS.md
git add docs/DIFFERENCES_SUMMARY.md
git add docs/FINAL_TESTING_SUMMARY.md
git add docs/INVESTIGATION_COMPLETE.md
git add docs/IS_NUCLEOTIDE_BUG_ANALYSIS.md
git add docs/JSON_GENERATION_SYSTEM.md
git add docs/KNOWN_DIFFERENCES_CATALOG.md
git add docs/LEGACY_JSON_GENERATION.md
git add docs/LEGACY_STEP_PARAMETER_ANALYSIS.md
git add docs/QUALITY_SCORE_DIFFERENCES.md
git add docs/STORAGE_OPTIMIZATION_STRATEGY.md
git add docs/TEST_RESULTS_BP_TYPE_ID_FIX.md
git add docs/VALIDATION_DIFFERENCES.md
```

### Archive Files (Medium Priority - Historical Reference)
These provide historical context but are not actively maintained:

```bash
git add docs/archive/
```

### Useful Scripts (Medium Priority)
These are utility scripts that may be useful:

```bash
git add scripts/CLEANUP_PLAN.md
git add scripts/analyze_find_bestpair_differences.py
git add scripts/check_legacy_json_progress.py
git add scripts/cleanup_optional_json.py
git add scripts/investigate_specific_pairs.py
git add scripts/test_multiple_pdbs.py
```

### IDE Configuration (Optional)
If you want to share IDE settings:

```bash
git add .cursorrules  # Optional - IDE configuration
```

## ❌ Files That Should NOT Be Committed

### IDE Files (Should be in .gitignore)
```bash
# Add to .gitignore if not already there:
.cursor/
```

### Temporary/Generated Files
```bash
# These should NOT be committed:
base_frame_calc/              # Temporary test directory
comparison_results_summary.json  # Can be regenerated
```

### Data Files (Already in .gitignore)
```bash
# Already ignored via .gitignore:
data/  # Contains PDB files and JSON (large, regenerable)
```

## 📝 Recommended Commit Commands

### Option 1: Commit Essential Files Only (Recommended)
```bash
# Essential root documentation
git add PROJECT_STATUS.md README.md DOCUMENTATION_INDEX.md
git add FINAL_100_PDB_COMPARISON_RESULTS.md 3CME_EDGE_CASE_ANALYSIS.md
git add TIE_BREAKING_FIX_SUMMARY.md TIE_BREAKING_FIX_FINAL_RESULTS.md
git add CLEANUP_PLAN_FILES.md CLEANUP_COMPLETE.md

# Active documentation
git add docs/*.md

# Archive files
git add docs/archive/

# Useful scripts
git add scripts/*.py scripts/CLEANUP_PLAN.md

# Commit
git commit -m "Add essential project documentation and cleanup summary

- Add project status and final comparison results
- Add tie-breaking fix documentation
- Add cleanup documentation and organization
- Add active documentation files
- Archive historical investigation files"
```

### Option 2: Commit Everything Except Temporary Files
```bash
# Add all markdown files
git add *.md
git add docs/*.md
git add docs/archive/
git add scripts/*.py scripts/*.md

# Add IDE config (optional)
git add .cursorrules

# Commit
git commit -m "Add project documentation and cleanup files"
```

## 🔧 Update .gitignore (Recommended)

Add these to `.gitignore` if not already present:

```bash
# IDE files
.cursor/

# Temporary test directories
base_frame_calc/

# Generated comparison results (can be regenerated)
comparison_results_summary.json
```

## Summary

**High Priority (Should Commit)**:
- ✅ 9 essential root documentation files
- ✅ 15 active documentation files in `docs/`
- ✅ Archive files in `docs/archive/`
- ✅ 5 useful scripts

**Optional**:
- ⚠️ `.cursorrules` (IDE config - personal preference)

**Should NOT Commit**:
- ❌ `.cursor/` (IDE files)
- ❌ `base_frame_calc/` (temporary)
- ❌ `comparison_results_summary.json` (generated)
- ❌ `data/` (already ignored, large files)

## Quick Command

To commit all essential files at once:

```bash
git add PROJECT_STATUS.md README.md DOCUMENTATION_INDEX.md \
        FINAL_100_PDB_COMPARISON_RESULTS.md 3CME_EDGE_CASE_ANALYSIS.md \
        TIE_BREAKING_FIX_SUMMARY.md TIE_BREAKING_FIX_FINAL_RESULTS.md \
        CLEANUP_PLAN_FILES.md CLEANUP_COMPLETE.md \
        docs/*.md docs/archive/ \
        scripts/*.py scripts/CLEANUP_PLAN.md && \
git commit -m "Add essential project documentation and cleanup files"
```

