# Documentation Index

**Last Updated**: 2025-01-XX  
**Status**: Consolidated and organized

---

## 🎯 Start Here - Essential Documents

### [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md) ⭐⭐⭐ **MAIN SUMMARY**
**Complete project overview and essential information**
- Project status and achievements
- Architecture overview
- Key components and data structures
- Usage examples
- **Use this for**: Quick overview, understanding the project, getting started

### [CODE_FLOW.md](CODE_FLOW.md) ⭐⭐⭐ **CODE ARCHITECTURE**
**Detailed code flow and architecture**
- Complete workflow from input to output
- Phase 1 (find_pair) and Phase 2 (analyze) flows
- Algorithm details and data structures
- **Use this for**: Understanding how the code works, debugging, development

### [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) ⭐ **STATUS DOCUMENT**
**Single source of truth for match status and progress**
- Current status: 97.8% perfect matches on primary output (312/319 PDBs)
- 10-PDB test set: 100% perfect matches (10/10)
- Completed fixes: Residue indexing, base_pair recording, H-bond conflict resolution, overlap calculation, bp_type_id preservation, Stage 6 implementation
- Remaining issues: 6 PDBs with selection mismatches (1TN1, 1TN2, 1TTT, 3F2T, 5V0O, 9CF3)
- **Use this for**: Current status, what's done, what's next, PDBs with differences

### [FIND_PAIR_VALIDATION_COMPREHENSIVE.md](FIND_PAIR_VALIDATION_COMPREHENSIVE.md) ⭐ **VALIDATION REPORT**
**Comprehensive validation results and status**
- Executable functionality: 100% success rate (140+ PDBs tested)
- Primary output match: 97.8% perfect match (312/319 PDBs)
- Test coverage: 10-PDB set, 30-PDB sample, 100-PDB sample, 319 PDBs comprehensive
- **Use this for**: Validation results, test statistics, known issues

### [BASE_PAIR_RECORDING.md](BASE_PAIR_RECORDING.md) ⭐ **BASE PAIR FIX**
**Complete documentation of base_pair recording fix**
- Problem: Modern recorded all validated pairs, legacy only records selected pairs
- Solution: Updated to match legacy behavior
- Results: 100% match on all tested PDBs (319/319, 25,323 pairs)
- **Use this for**: Understanding base_pair recording behavior

### [TESTING_GUIDE.md](TESTING_GUIDE.md) ⭐ **TESTING GUIDE**
**Complete guide for testing and validation**
- JSON comparison workflows
- Test types and execution
- Troubleshooting
- **Use this for**: How to test, compare, and validate

### [NEXT_STEPS.md](NEXT_STEPS.md) ⭐ **NEXT STEPS**
**What to do next after find_pair validation**
- Priority 1: Move to step parameters (analyze phase)
- Priority 2: Investigate 6 mismatched PDBs (optional)
- Priority 3: Expand testing (optional)
- **Use this for**: Clear guidance on what to work on next

### [HBOND_INVESTIGATION.md](HBOND_INVESTIGATION.md)
**H-bond investigation and fixes**
- Problem: Modified nucleotides getting incorrect H-bond types
- Solution: `get_base_type_for_hbond()` function
- Status: Complete (see [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md) for current status)
- **Use this for**: Understanding the H-bond investigation process

### [TOOLS_TO_CREATE.md](TOOLS_TO_CREATE.md)
**Recommended tools to help reach 100% match**
- Priority-ordered list of tools to create
- Time estimates and use cases
- **Use this for**: Deciding which tools to build next

---

## 📚 Reference Documents

### Data Structures & JSON
- **[JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)** - JSON structure, field meanings, comparison methodology
- **[legacy/10_JSON_STRUCTURE.md](legacy/10_JSON_STRUCTURE.md)** - Legacy JSON format reference (see [legacy/](legacy/) for all legacy documentation)
- **[JSON_SERIALIZATION.md](JSON_SERIALIZATION.md)** - Modern JSON serialization implementation

### Legacy Code Reference
- **[legacy/](legacy/)** - Complete legacy code documentation
  - [Index](legacy/00_INDEX.md) - Navigation guide for all legacy docs
  - [Architecture](legacy/01_ARCHITECTURE.md) - Program structure and flows
  - [Data Structures](legacy/02_DATA_STRUCTURES.md) - All data organization
  - [Core Functions](legacy/03_CORE_FUNCTIONS.md) - Critical function profiles
  - [Algorithms](legacy/04_ALGORITHMS.md) - Mathematical derivations
  - [Implementation Guide](legacy/08_IMPLEMENTATION_GUIDE.md) - Conversion strategies

### Algorithm Implementation
- **[ALGORITHM_CRITICAL_GUIDE.md](ALGORITHM_CRITICAL_GUIDE.md)** - Critical algorithm details and edge cases
- **[OVERLAP_CALCULATION.md](OVERLAP_CALCULATION.md)** - Overlap calculation implementation details
- **[BP_TYPE_ID_REQUIREMENTS.md](BP_TYPE_ID_REQUIREMENTS.md)** - Base pair type ID calculation requirements

### Development & Debugging
- **[BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md)** - Build and compilation instructions
- **[CONFIGURATION_OPTIONS.md](CONFIGURATION_OPTIONS.md)** - Configuration options reference
- **[DEBUGGING_STRATEGIES.md](DEBUGGING_STRATEGIES.md)** - Debugging strategies and tools
- **[QUICK_START.md](QUICK_START.md)** - Quick start guide for new developers

### Modernization
- **[MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md)** - Overall modernization strategy
- **[modernization/](modernization/)** - Stage-by-stage modernization details (10 stages)

---

## 📖 Historical/Investigation Documents

### [HBOND_INVESTIGATION.md](HBOND_INVESTIGATION.md)
**H-bond investigation process** (historical reference)
- Full investigation that led to the H-bond fix
- Detailed findings and test results
- **Status**: Complete, kept for reference
- **Use this for**: Understanding how the H-bond issue was diagnosed

---

## 📦 Archived Documents

**Note**: Many historical documents have been archived to `docs/archive/historical/` to keep the main documentation clean and focused.

### Archive Structure

- `archive/historical/superseded/` - Superseded by consolidated docs
- `archive/historical/test_results/` - Old test results
- `archive/historical/investigations/` - Historical investigations
- `archive/historical/investigations_2025/` - 2025 investigation documents
- `archive/historical/status_reports_2025/` - 2025 status reports
- `archive/historical/old_plans/` - Old action plans
- `archive/historical/completed_fixes/` - Completed fix documentation
- `archive/historical/redundant/` - Redundant/duplicate docs

**Total archived**: ~80+ files  
**Remaining in root**: ~20 essential files

### Key Archived Documents

- Superseded validation docs → `archive/historical/superseded/`
- Old test results → `archive/historical/test_results/`
- Historical investigations → `archive/historical/investigations/` and `investigations_2025/`
- Status reports → `archive/historical/status_reports_2025/`
- Old plans/status → `archive/historical/old_plans/`
- Completed fixes → `archive/historical/completed_fixes/`
- Redundant docs → `archive/historical/redundant/`

**See**: [archive/ARCHIVE_PLAN.md](archive/ARCHIVE_PLAN.md) for complete archive plan

---

## Quick Navigation

### By Task

**I want to...**
- **Get project overview**: [PROJECT_SUMMARY.md](PROJECT_SUMMARY.md) ⭐⭐⭐ **Start here**
- **Understand code flow**: [CODE_FLOW.md](CODE_FLOW.md) ⭐⭐⭐ **Architecture**
- **See current status**: [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md) ⭐ **Status**
- **Know what to do next**: [100_PERCENT_MATCH_STATUS.md#what-to-do-next](100_PERCENT_MATCH_STATUS.md#what-to-do-next) ⭐ **Next steps**
- **See validation results**: [FIND_PAIR_VALIDATION_COMPREHENSIVE.md](FIND_PAIR_VALIDATION_COMPREHENSIVE.md)
- **Understand base_pair recording**: [BASE_PAIR_RECORDING.md](BASE_PAIR_RECORDING.md)
- **Run tests**: [TESTING_GUIDE.md](TESTING_GUIDE.md)
- **Use --fix-indices**: [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md)
- **Compare JSON structures**: [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)
- **See what's compared**: [COMPARISON_COVERAGE.md](COMPARISON_COVERAGE.md)
- **Debug an issue**: [DEBUGGING_STRATEGIES.md](DEBUGGING_STRATEGIES.md)
- **Build the project**: [BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md)

### By Topic

**Topics:**
- **100% Match Goal**: [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md)
- **Validation Results**: [FIND_PAIR_VALIDATION_COMPREHENSIVE.md](FIND_PAIR_VALIDATION_COMPREHENSIVE.md)
- **Base Pair Recording**: [BASE_PAIR_RECORDING.md](BASE_PAIR_RECORDING.md)
- **Residue Indexing**: [FIX_INDICES_OPTION.md](FIX_INDICES_OPTION.md)
- **Testing**: [TESTING_GUIDE.md](TESTING_GUIDE.md)
- **H-bond Detection**: [HBOND_INVESTIGATION.md](HBOND_INVESTIGATION.md), [legacy/04_ALGORITHMS.md](legacy/04_ALGORITHMS.md#3-h-bond-detection-workflow)
- **Quality Scores**: [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md#quality-score-calculation)
- **Stage 6 (ParameterCalculator)**: [100_PERCENT_MATCH_STATUS.md](100_PERCENT_MATCH_STATUS.md#8-stage-6-parametercalculator-implementation-2025-11-26--complete)
- **JSON Structure**: [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)
- **Algorithms**: [ALGORITHM_CRITICAL_GUIDE.md](ALGORITHM_CRITICAL_GUIDE.md)
- **Modernization**: [MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md), [modernization/](modernization/)

---

## Document Status

| Document | Status | Last Updated | Priority |
|----------|--------|--------------|----------|
| PROJECT_SUMMARY.md | ✅ Active | 2025-01-XX | ⭐⭐⭐ |
| CODE_FLOW.md | ✅ Active | 2025-01-XX | ⭐⭐⭐ |
| 100_PERCENT_MATCH_STATUS.md | ✅ Active | 2025-11-28 | ⭐⭐⭐ |
| FIND_PAIR_VALIDATION_COMPREHENSIVE.md | ✅ Active | 2025-11-28 | ⭐⭐⭐ |
| BASE_PAIR_RECORDING.md | ✅ Active | 2025-11-28 | ⭐⭐⭐ |
| TESTING_GUIDE.md | ✅ Active | 2025-11-28 | ⭐⭐⭐ |
| FIX_INDICES_OPTION.md | 📚 Reference | - | ⭐⭐ |
| COMPARISON_COVERAGE.md | 📚 Reference | 2025-11-28 | ⭐⭐ |
| JSON_DATA_TYPES_AND_COMPARISONS.md | 📚 Reference | - | ⭐⭐ |
| ALGORITHM_CRITICAL_GUIDE.md | 📚 Reference | - | ⭐⭐ |
| DEBUGGING_STRATEGIES.md | 📚 Reference | - | ⭐ |
| BUILD_INSTRUCTIONS.md | 📚 Reference | - | ⭐⭐ |
| MODERNIZATION_PLAN.md | 📚 Reference | - | ⭐ |

**Legend:**
- ✅ Active - Currently being used/updated
- 📖 Reference - Historical/investigation documentation
- 📚 Reference - Technical reference documentation
- ⭐⭐⭐ High priority
- ⭐⭐ Medium priority
- ⭐ Low priority

---

## Maintenance Notes

- **Primary document**: Always keep `100_PERCENT_MATCH_STATUS.md` as the single source of truth
- **When adding new docs**: Update this README
- **When removing docs**: Mark them in "Removed Documents" section
- **Status updates**: Update "Document Status" table
