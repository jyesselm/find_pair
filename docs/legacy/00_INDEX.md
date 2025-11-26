# Legacy Code Documentation Index

**Date**: 2025-01-XX  
**Purpose**: Comprehensive documentation of legacy X3DNA code structure, algorithms, and implementation details  
**Status**: Complete reference guide for matching legacy behavior

---

## Document Organization

This documentation is organized into focused subdocuments covering different aspects of the legacy codebase:

### Core Documents

1. **[Architecture Overview](01_ARCHITECTURE.md)**
   - Program entry points and main flows
   - File organization and structure
   - Module dependencies
   - Build system overview

2. **[Data Structures](02_DATA_STRUCTURES.md)**
   - All data structure definitions
   - Memory allocation patterns
   - Indexing conventions (1-based)
   - Array access patterns
   - Frame matrix representation

3. **[Core Functions](03_CORE_FUNCTIONS.md)**
   - Critical function profiles (base_frame, check_pair, etc.)
   - Function signatures and parameters
   - Input/output specifications
   - Dependencies and call graphs
   - Function responsibilities

4. **[Algorithms](04_ALGORITHMS.md)**
   - Mathematical derivations
   - Step-by-step pseudocode
   - Numerical examples
   - Edge cases and special handling
   - Algorithm correctness proofs

5. **[Helper Functions](05_HELPER_FUNCTIONS.md)**
   - Utility function reference
   - Mathematical operations
   - String manipulation
   - Memory management
   - File I/O helpers

6. **[Parameters & Constants](06_PARAMETERS.md)**
   - All parameter definitions
   - Default values and ranges
   - Configuration loading
   - Validation rules
   - Constants and macros

7. **[Workflows](07_WORKFLOWS.md)**
   - Complete program workflows
   - Function call sequences
   - Data flow diagrams
   - Integration patterns
   - Execution order

8. **[Implementation Guide](08_IMPLEMENTATION_GUIDE.md)**
   - How to match legacy behavior
   - Conversion guidelines (1-based ↔ 0-based)
   - Common pitfalls
   - Verification checklist
   - Testing strategy

9. **[Knowledge Base](09_KNOWLEDGE_BASE.md)**
   - Quick reference for common patterns
   - Indexing conventions
   - Residue indexing details
   - Common pitfalls and fixes
   - Recent fixes and solutions

10. **[JSON Structure](10_JSON_STRUCTURE.md)**
    - JSON file organization
    - Record types and stages
    - Base pair identification flow
    - JSON record formats
    - Indexing in JSON

11. **[Test Tools](11_TEST_TOOLS.md)**
    - Isolated component testing tools
    - Comparison workflow
    - Usage examples
    - Integration with modern code

---

## Quick Reference

### Most Critical Functions

1. **`base_frame()`** → [Core Functions](03_CORE_FUNCTIONS.md#1-base_frame)
   - Calculates base reference frames
   - Uses least-squares fitting

2. **`ls_fitting()`** → [Algorithms](04_ALGORITHMS.md#11-ls_fitting)
   - Quaternion-based fitting algorithm
   - Mathematical foundation

3. **`check_pair()`** → [Core Functions](03_CORE_FUNCTIONS.md#3-check_pair)
   - Validates base pairs
   - Core validation logic

4. **`hb_atompair()`** → [Algorithms](04_ALGORITHMS.md#5-hb_atompair)
   - H-bond conflict resolution
   - Critical for matching

5. **`find_bestpair()`** → [Core Functions](03_CORE_FUNCTIONS.md#7-find_bestpair)
   - Greedy matching algorithm
   - Base pair selection

### Key Concepts

- **1-based indexing**: All arrays use 1-based indexing → [Data Structures](02_DATA_STRUCTURES.md#indexing-conventions)
- **Frame representation**: 9-element flattened matrices → [Data Structures](02_DATA_STRUCTURES.md#base-reference-frames)
- **Residue indexing**: Sequential assignment during PDB parsing → [Data Structures](02_DATA_STRUCTURES.md#residue-indexing)
- **H-bond validation**: Three-phase process → [Algorithms](04_ALGORITHMS.md#13-h-bond-detection-workflow)

---

## Navigation by Task

### "I want to understand how..."

**...base frames are calculated**
→ Start: [Core Functions: base_frame](03_CORE_FUNCTIONS.md#1-base_frame)
→ Math: [Algorithms: ls_fitting](04_ALGORITHMS.md#11-ls_fitting)

**...base pairs are validated**
→ Start: [Core Functions: check_pair](03_CORE_FUNCTIONS.md#3-check_pair)
→ Details: [Algorithms: Validation Process](04_ALGORITHMS.md#validation-process)

**...H-bonds are detected**
→ Start: [Core Functions: get_hbond_ij](03_CORE_FUNCTIONS.md#4-get_hbond_ij)
→ Algorithm: [Algorithms: H-bond Detection](04_ALGORITHMS.md#13-h-bond-detection-workflow)
→ Conflict: [Algorithms: hb_atompair](04_ALGORITHMS.md#5-hb_atompair)

**...the greedy matching works**
→ Start: [Core Functions: find_bestpair](03_CORE_FUNCTIONS.md#7-find_bestpair)
→ Workflow: [Workflows: Base Pair Finding](07_WORKFLOWS.md#workflow-2-base-pair-finding)

**...step parameters are calculated**
→ Start: [Algorithms: bpstep_par](04_ALGORITHMS.md#12-bpstep_par)
→ Math: [Algorithms: Step Parameters](04_ALGORITHMS.md#12-bpstep_par)

**...data is organized**
→ Start: [Data Structures: Overview](02_DATA_STRUCTURES.md)
→ Memory: [Data Structures: Memory Allocation](02_DATA_STRUCTURES.md#memory-allocation-pattern)

### "I need to match..."

**...frame calculations exactly**
→ Reference: [Core Functions: base_frame](03_CORE_FUNCTIONS.md#1-base_frame)
→ Math: [Algorithms: ls_fitting](04_ALGORITHMS.md#11-ls_fitting)
→ Guide: [Implementation Guide: Frame Matching](08_IMPLEMENTATION_GUIDE.md#matching-base-frames)

**...H-bond detection exactly**
→ Reference: [Algorithms: H-bond Detection](04_ALGORITHMS.md#13-h-bond-detection-workflow)
→ Conflict: [Algorithms: hb_atompair](04_ALGORITHMS.md#5-hb_atompair)
→ Guide: [Implementation Guide: H-bond Matching](08_IMPLEMENTATION_GUIDE.md#matching-h-bonds)

**...validation logic exactly**
→ Reference: [Core Functions: check_pair](03_CORE_FUNCTIONS.md#3-check_pair)
→ Parameters: [Parameters: Validation](06_PARAMETERS.md#validation-parameters)
→ Guide: [Implementation Guide: Validation Matching](08_IMPLEMENTATION_GUIDE.md#matching-validation)

**...indexing conventions**
→ Reference: [Data Structures: Indexing](02_DATA_STRUCTURES.md#indexing-conventions)
→ Guide: [Implementation Guide: Index Conversion](08_IMPLEMENTATION_GUIDE.md#index-conversion)

---

## Document Status

| Document | Status | Last Updated |
|----------|--------|--------------|
| Architecture | ✅ Complete | 2025-01-XX |
| Data Structures | ✅ Complete | 2025-01-XX |
| Core Functions | ✅ Complete | 2025-01-XX |
| Algorithms | ✅ Complete | 2025-01-XX |
| Helper Functions | ✅ Complete | 2025-01-XX |
| Parameters | ✅ Complete | 2025-01-XX |
| Workflows | ✅ Complete | 2025-01-XX |
| Implementation Guide | ✅ Complete | 2025-01-XX |
| Knowledge Base | ✅ Complete | 2025-01-XX |
| JSON Structure | ✅ Complete | 2025-01-XX |
| Test Tools | ✅ Complete | 2025-01-XX |

---

## Related Documentation

- **[LEGACY_KNOWLEDGE.md](../LEGACY_KNOWLEDGE.md)**: Quick reference for common issues
- **[ALGORITHM_CRITICAL_GUIDE.md](../ALGORITHM_CRITICAL_GUIDE.md)**: High-level algorithm overview
- **[LEGACY_JSON_STRUCTURE.md](../LEGACY_JSON_STRUCTURE.md)**: JSON output format

---

**Next Steps**: Start with [Architecture Overview](01_ARCHITECTURE.md) for high-level understanding, then dive into specific areas as needed.

