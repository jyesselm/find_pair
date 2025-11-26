# Documentation Index

**Last Updated**: 2025-11-25  
**Status**: Consolidated and organized

---

## 🎯 Primary Documents (Start Here)

### [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md) ⭐ **MAIN DOCUMENT**
**Single source of truth for reaching 100% match**
- Current status: 99.5% match (1044/1048 pairs)
- Completed fixes: H-bond type detection for modified nucleotides
- Next steps: Prioritized action plan with 5 clear steps
- **Use this for**: Current status, what's done, what's next

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

## 🗑️ Removed Documents

The following documents were removed as they were outdated or superseded:

- ~~**PAIR_MATCHING_ISSUES.md**~~ - Outdated status (contradicted current plan)
- ~~**PAIR_MATCHING_DEBUGGING.md**~~ - Superseded by HBOND_INVESTIGATION.md and 100_PERCENT_MATCH_PLAN.md
- ~~**RESIDUE_INVESTIGATION_FINDINGS.md**~~ - Information consolidated into 100_PERCENT_MATCH_PLAN.md
- ~~**RESIDUE_INVESTIGATION_TOOLS.md**~~ - Tool information moved to TOOLS_TO_CREATE.md

---

## Quick Navigation

### By Task

**I want to...**
- **See current status**: [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md#executive-summary)
- **Know what to do next**: [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md#next-steps-to-reach-100)
- **Understand the H-bond fix**: [HBOND_INVESTIGATION.md](HBOND_INVESTIGATION.md)
- **Find tools to build**: [TOOLS_TO_CREATE.md](TOOLS_TO_CREATE.md)
- **Compare JSON structures**: [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)
- **Debug an issue**: [DEBUGGING_STRATEGIES.md](DEBUGGING_STRATEGIES.md)
- **Build the project**: [BUILD_INSTRUCTIONS.md](BUILD_INSTRUCTIONS.md)

### By Topic

**Topics:**
- **100% Match Goal**: [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md)
- **H-bond Detection**: [HBOND_INVESTIGATION.md](HBOND_INVESTIGATION.md), [legacy/04_ALGORITHMS.md](legacy/04_ALGORITHMS.md#3-h-bond-detection-workflow)
- **Quality Scores**: [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md#step-2-quality-score-investigation-high-priority---1-2-days)
- **Residue Recognition**: [100_PERCENT_MATCH_PLAN.md](100_PERCENT_MATCH_PLAN.md#step-3-residue-inclusion-fix-high-priority---1-2-days)
- **JSON Structure**: [JSON_DATA_TYPES_AND_COMPARISONS.md](JSON_DATA_TYPES_AND_COMPARISONS.md)
- **Algorithms**: [ALGORITHM_CRITICAL_GUIDE.md](ALGORITHM_CRITICAL_GUIDE.md)
- **Modernization**: [MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md), [modernization/](modernization/)

---

## Document Status

| Document | Status | Last Updated | Priority |
|----------|--------|--------------|----------|
| 100_PERCENT_MATCH_PLAN.md | ✅ Active | 2025-11-25 | ⭐⭐⭐ |
| HBOND_INVESTIGATION.md | 📖 Reference | 2025-11-25 | ⭐ |
| TOOLS_TO_CREATE.md | ✅ Active | 2025-11-25 | ⭐⭐ |
| HBOND_INVESTIGATION.md | 📖 Reference | 2025-11-25 | ⭐ |
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

- **Primary document**: Always keep `100_PERCENT_MATCH_PLAN.md` as the single source of truth
- **When adding new docs**: Update this README
- **When removing docs**: Mark them in "Removed Documents" section
- **Status updates**: Update "Document Status" table
