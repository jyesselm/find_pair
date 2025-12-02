# Documentation Index

**Last Updated**: 2025-01-XX  
**Status**: Organized and focused on current work

---

## 🎯 Start Here - Current Work

### [MATCHING_PLAN.md](MATCHING_PLAN.md) ⭐⭐⭐ **CURRENT WORK**
**Plan to match legacy and modern results exactly**
- Strategy 1: Create minimal test cases (2 consecutive base pairs)
- Strategy 2: Break legacy code into smaller testable units
- Strategy 3: Subdivide steps into finer-grained operations
- Implementation plan with tools to create
- **Use this for**: Understanding current approach to achieve exact matching

### [TESTING_GUIDE.md](TESTING_GUIDE.md) ⭐⭐⭐ **TESTING REFERENCE**
**Complete guide for testing and validation**
- JSON comparison workflows using `compare_json.py`
- Test types and execution
- Troubleshooting guide
- Step parameter comparison
- Reference frames comparison
- **Use this for**: How to test, compare, and validate results

---

## 📚 Essential Reference Documents

### Project Overview
- **[PROJECT_SUMMARY.md](PROJECT_SUMMARY.md)** - Complete project overview and status
- **[CODE_FLOW.md](CODE_FLOW.md)** - Detailed code flow and architecture
- **[QUICK_START.md](QUICK_START.md)** - Quick start guide for new developers
- **[BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md)** - Build and compilation instructions

### Data & JSON Structure
- **[JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)** - JSON structure, field meanings, comparison methodology
- **[DATA_STRUCTURE.md](DATA_STRUCTURE.md)** - Data directory structure and organization

### Algorithm & Implementation Guides
- **[ALGORITHM_CRITICAL_GUIDE.md](ALGORITHM_CRITICAL_GUIDE.md)** - Critical algorithm details and edge cases
- **[STEP_PARAMETERS_IMPLEMENTATION.md](STEP_PARAMETERS_IMPLEMENTATION.md)** - Step parameters implementation

### Residue Indexing & Legacy Matching
- **[LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md)** - Guide to using legacy indices for accurate comparisons
- **[FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md)** - Using `--fix-indices` option for legacy matching

### Reference Frames
- **[REF_FRAMES_NEXT_STEPS.md](REF_FRAMES_NEXT_STEPS.md)** - Next steps for ref_frames comparison
- **[REF_FRAMES_DEBUGGING_APPROACH.md](REF_FRAMES_DEBUGGING_APPROACH.md)** - Systematic debugging approach
- **[REF_FRAMES_COMPARISON_FIX.md](REF_FRAMES_COMPARISON_FIX.md)** - Fix documentation for ref_frames
- **[REF_FRAMES_QUICK_START.md](REF_FRAMES_QUICK_START.md)** - Quick start for ref_frames debugging

---

## 📖 Legacy Code Reference

### [legacy/](legacy/) - Complete Legacy Code Documentation
- **[00_INDEX.md](legacy/00_INDEX.md)** - Navigation guide for all legacy docs
- **[01_ARCHITECTURE.md](legacy/01_ARCHITECTURE.md)** - Program structure and flows
- **[02_DATA_STRUCTURES.md](legacy/02_DATA_STRUCTURES.md)** - All data organization
- **[03_CORE_FUNCTIONS.md](legacy/03_CORE_FUNCTIONS.md)** - Critical function profiles
- **[04_ALGORITHMS.md](legacy/04_ALGORITHMS.md)** - Mathematical derivations and algorithms
- **[10_JSON_STRUCTURE.md](legacy/10_JSON_STRUCTURE.md)** - Legacy JSON format reference

**Use legacy docs for**: Understanding legacy code behavior, algorithm details, data structures

---

## 📦 Archived Documents

Historical, outdated, and completed documentation has been moved to `archive/`:

### Archive Structure
- `archive/status_docs/` - Old status documents and plans
- `archive/historical_fixes/` - Completed fix documentation
- `archive/old_debug/` - Old debug reports
- `archive/historical/` - Historical investigations, test results, status reports
  - `completed_fixes/` - Documentation of completed fixes
  - `investigations/` - Historical investigations
  - `investigations_2025/` - 2025 investigation documents
  - `status_reports_2025/` - 2025 status reports
  - `old_plans/` - Old action plans
  - `redundant/` - Redundant/duplicate docs
  - `superseded/` - Superseded documents

### Recently Archived
The following documents were archived to keep the main docs clean:
- `100_PERCENT_MATCH_STATUS.md` → `archive/status_docs/`
- `NEXT_STEPS.md` → `archive/status_docs/` (superseded by MATCHING_PLAN.md)
- `MODERNIZATION_PLAN.md` → `archive/status_docs/` (modernization complete)
- `FIND_PAIR_VALIDATION_COMPREHENSIVE.md` → `archive/status_docs/`
- `COMPARISON_COVERAGE.md` → `archive/status_docs/`
- `DEBUGGING_STRATEGIES.md` → `archive/status_docs/`
- `CONFIGURATION_OPTIONS.md` → `archive/status_docs/`
- `ATOM_INDEX_CONVERSION_FIX.md` → `archive/historical_fixes/`
- `BASE_PAIR_RECORDING.md` → `archive/historical_fixes/`
- `pdb_debug_reports/` → `archive/old_debug/`

---

## Quick Navigation

### By Task

**I want to...**
- **See current work plan**: [MATCHING_PLAN.md](MATCHING_PLAN.md) ⭐⭐⭐
- **Run tests/compare results**: [TESTING_GUIDE.md](TESTING_GUIDE.md) ⭐⭐⭐
- **Understand the code**: [CODE_FLOW.md](CODE_FLOW.md)
- **Get project overview**: [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md)
- **Match legacy indices**: [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md)
- **Use --fix-indices**: [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md)
- **Debug ref_frames**: [REF_FRAMES_NEXT_STEPS.md](REF_FRAMES_NEXT_STEPS.md)
- **Understand JSON structure**: [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)
- **Understand legacy code**: [legacy/00_INDEX.md](legacy/00_INDEX.md)

### By Topic

**Current Focus:**
- **Matching Plan**: [MATCHING_PLAN.md](MATCHING_PLAN.md) - Current work to achieve exact matching
- **Testing**: [TESTING_GUIDE.md](TESTING_GUIDE.md) - Testing and comparison workflows
- **Ref Frames**: [REF_FRAMES_NEXT_STEPS.md](REF_FRAMES_NEXT_STEPS.md) - Reference frames debugging

**Reference:**
- **Residue Indexing**: [LEGACY_INDICES_GUIDE.md](LEGACY_INDICES_GUIDE.md), [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md)
- **Step Parameters**: [STEP_PARAMETERS_IMPLEMENTATION.md](STEP_PARAMETERS_IMPLEMENTATION.md)
- **Algorithms**: [ALGORITHM_CRITICAL_GUIDE.md](ALGORITHM_CRITICAL_GUIDE.md)
- **JSON Structure**: [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)
- **Legacy Code**: [legacy/](legacy/) - Complete legacy documentation

---

## Document Status

| Document | Status | Priority |
|----------|--------|----------|
| MATCHING_PLAN.md | ✅ Active | ⭐⭐⭐ |
| TESTING_GUIDE.md | ✅ Active | ⭐⭐⭐ |
| REF_FRAMES_NEXT_STEPS.md | ✅ Active | ⭐⭐ |
| LEGACY_INDICES_GUIDE.md | 📚 Reference | ⭐⭐ |
| FIX_INDICES_OPTION.md | 📚 Reference | ⭐⭐ |
| CODE_FLOW.md | 📚 Reference | ⭐⭐ |
| PROJECT_SUMMARY.md | 📚 Reference | ⭐⭐ |
| JSON_DATA_TYPES_AND_COMPARISONS.md | 📚 Reference | ⭐⭐ |

**Legend:**
- ✅ Active - Currently being used/updated
- 📚 Reference - Technical reference documentation
- ⭐⭐⭐ High priority
- ⭐⭐ Medium priority

---

*This documentation is organized to focus on current work (matching legacy and modern results) while keeping essential reference material easily accessible.*
