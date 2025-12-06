# Comprehensive Stage 2 Refactoring Plan

**Date**: December 5, 2025  
**Status**: Planning Phase  
**Scope**: Refactor all Stage 2 code (frame calculation) and related infrastructure  
**Goal**: Modern OOP architecture that is clean, extensible, and maintainable

---

## Executive Summary

Stage 2 (frame calculation) is now **100% validated** (3,602/3,602 PDBs passing), but the code has become messy through iterative development. This refactoring will transform it into a clean, modern OOP architecture while preserving all functionality.

### Current Problems

1. **Mixed Responsibilities**: `BaseFrameCalculator` (933 lines) does calculation AND JSON recording
2. **Code Duplication**: Iteration logic duplicated in 3+ places
3. **Confusing Abstractions**: `LsFittingCalculator` is just a wrapper around `BaseFrameCalculator`
4. **Hard to Extend**: Adding features requires changes across multiple files
5. **Single Responsibility Violations**: Classes doing too much
6. **Poor Testability**: Can't test algorithms without JSON dependencies

### Target Architecture

```
┌─────────────────────────────────────────────────┐
│ ALGORITHMS (Pure Calculation)                   │
│  ├─ BaseFrameCalculator                         │
│  ├─ RingAtomMatcher                             │
│  ├─ StandardBaseTemplates                       │
│  └─ LeastSquaresFitter                          │
└─────────────────────────────────────────────────┘
                    │
                    │ uses
                    ▼
┌─────────────────────────────────────────────────┐
│ SERVICES (Orchestration)                        │
│  └─ FrameCalculationService                     │
│     ├─ calculate_frame(Residue&)                │
│     ├─ calculate_all_frames(Structure&)         │
│     └─ validate_template_assignment(Residue&)   │
└─────────────────────────────────────────────────┘
                    │
                    │ uses
                    ▼
┌─────────────────────────────────────────────────┐
│ RECORDERS (JSON Output)                         │
│  └─ FrameJsonRecorder                           │
│     ├─ record_ls_fitting()                      │
│     ├─ record_base_frame_calc()                 │
│     └─ record_frame_calc()                      │
└─────────────────────────────────────────────────┘
                    │
                    │ uses
                    ▼
┌─────────────────────────────────────────────────┐
│ I/O (Low-level JSON)                            │
│  └─ JsonWriter                                  │
└─────────────────────────────────────────────────┘
```

---

## Part 1: Core Refactoring (Frame Calculation)

### Problem Analysis

#### BaseFrameCalculator (933 lines) - TOO MANY RESPONSIBILITIES

**Current structure:**
```cpp
class BaseFrameCalculator {
    // GOOD: Core algorithm
    FrameCalculationResult calculate_frame(Residue&);
    
    // BAD: JSON recording (should be separate)
    void calculate_and_record_frames(Structure&, JsonWriter&);
    void calculate_and_record_frames_only(Structure&, JsonWriter&);
    
    // BAD: Iteration logic (should be in service layer)
    // ... 100+ lines of residue iteration ...
};
```

**Lines breakdown:**
- Lines 1-180: Core algorithm ✅ KEEP
- Lines 181-530: RMSD checks, type detection ✅ KEEP (but extract to helpers)
- Lines 531-650: JSON recording ❌ MOVE to FrameJsonRecorder
- Lines 651-933: Helper functions ✅ KEEP (but refactor)

#### LsFittingCalculator - UNNECESSARY WRAPPER

**Current structure:**
```cpp
class LsFittingCalculator {
    BaseFrameCalculator calculator_;
    
    // Just delegates everything!
    size_t calculate_and_record(...) {
        return calculator_.calculate_and_record_frames(...);
    }
};
```

**Verdict**: DELETE this class entirely ❌

---

## Part 2: New Architecture

### 2.1 Algorithm Layer (Pure Calculation)

#### BaseFrameCalculator (Refactored)
**File**: `include/x3dna/algorithms/base_frame_calculator.hpp`  
**Responsibility**: Calculate reference frames for nucleotides  
**Size**: ~400 lines (down from 933)

```cpp
class BaseFrameCalculator {
public:
    // Core calculation
    FrameCalculationResult calculate_frame(core::Residue& residue);
    FrameCalculationResult calculate_frame_const(const core::Residue& residue) const;
    
    // Configuration
    void set_template_path(const std::filesystem::path& path);
    void set_is_rna(bool is_rna);
    void set_legacy_mode(bool legacy_mode);
    
    // Utilities
    static bool detect_rna(const core::Structure& structure);
    
    // NO JSON recording methods! ✅
    
private:
    FrameCalculationResult calculate_frame_impl(const core::Residue& residue) const;
    
    // Extract to helper classes:
    StandardBaseTemplates templates_;
    RingAtomMatcher matcher_;  // NEW - extract ring matching logic
    ResidueTypeDetector type_detector_;  // NEW - extract RMSD type detection
};
```

**Changes:**
1. ❌ REMOVE: `calculate_and_record_frames()` → Move to `FrameJsonRecorder`
2. ❌ REMOVE: `calculate_and_record_frames_only()` → Move to `FrameJsonRecorder`
3. ✅ ADD: `RingAtomMatcher` helper class
4. ✅ ADD: `ResidueTypeDetector` helper class
5. ✅ SIMPLIFY: Focus only on frame calculation

---

#### RingAtomMatcher (NEW)
**File**: `include/x3dna/algorithms/ring_atom_matcher.hpp`  
**Responsibility**: Match experimental ring atoms to standard templates  
**Size**: ~150 lines

```cpp
struct RingMatchResult {
    std::vector<std::string> matched_atom_names;
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    bool has_purine_atoms;
    int purine_atom_count;
};

class RingAtomMatcher {
public:
    RingMatchResult match(const core::Residue& residue) const;
    
private:
    RingMatchResult match_all_atoms(const core::Residue& residue) const;
    RingMatchResult match_pyrimidine_only(const core::Residue& residue) const;
};
```

**Extracted from**: `BaseFrameCalculator::check_nt_type_by_rmsd()` (lines 52-164)

---

#### ResidueTypeDetector (NEW)
**File**: `include/x3dna/algorithms/residue_type_detector.hpp`  
**Responsibility**: Determine nucleotide type via RMSD and atom analysis  
**Size**: ~200 lines

```cpp
struct TypeDetectionResult {
    core::ResidueType detected_type;
    bool used_fallback;
    std::optional<double> rmsd;
    std::string detection_method;  // "registry", "rmsd", "atom_analysis"
};

class ResidueTypeDetector {
public:
    TypeDetectionResult detect_type(const core::Residue& residue) const;
    
private:
    bool should_use_rmsd_check(const core::Residue& residue) const;
    core::ResidueType detect_by_atoms(const core::Residue& residue) const;
    TypeDetectionResult detect_by_rmsd(const core::Residue& residue) const;
};
```

**Extracted from**: `BaseFrameCalculator::calculate_frame_impl()` (lines 187-530)

---

### 2.2 Service Layer (Orchestration)

#### FrameCalculationService (NEW)
**File**: `include/x3dna/services/frame_calculation_service.hpp`  
**Responsibility**: Orchestrate frame calculation across structure  
**Size**: ~200 lines

```cpp
class FrameCalculationService {
public:
    explicit FrameCalculationService(const std::filesystem::path& template_path);
    
    // Batch operations
    void calculate_all_frames(core::Structure& structure);
    std::vector<FrameCalculationResult> calculate_frames_const(
        const core::Structure& structure) const;
    
    // Single operations
    FrameCalculationResult calculate_frame(core::Residue& residue);
    
    // Configuration
    void set_legacy_mode(bool legacy_mode);
    void auto_detect_rna(const core::Structure& structure);
    
    // Validation
    bool validate_template_assignment(const core::Residue& residue) const;
    
private:
    BaseFrameCalculator calculator_;
    ResidueTypeDetector type_detector_;
};
```

**Purpose:**
- High-level API for frame calculation
- Handles RNA detection automatically
- Provides validation methods
- Clean interface for other components

---

### 2.3 Recording Layer (JSON Output)

#### FrameJsonRecorder (NEW)
**File**: `include/x3dna/io/frame_json_recorder.hpp`  
**Responsibility**: Record frame calculation results to JSON  
**Size**: ~250 lines

```cpp
class FrameJsonRecorder {
public:
    explicit FrameJsonRecorder(FrameCalculationService& service);
    
    // Individual record methods
    size_t record_ls_fitting(core::Structure& structure, io::JsonWriter& writer);
    size_t record_base_frame_calc(core::Structure& structure, io::JsonWriter& writer);
    size_t record_frame_calc(core::Structure& structure, io::JsonWriter& writer);
    
    // Combined recording
    size_t record_all(core::Structure& structure, io::JsonWriter& writer);
    
private:
    FrameCalculationService& service_;
    
    // Shared iteration logic (DRY principle)
    template<typename RecordFunc>
    size_t iterate_and_record(core::Structure& structure, 
                               io::JsonWriter& writer, 
                               RecordFunc record_func);
};
```

**Benefits:**
1. Eliminates duplication (iteration logic in ONE place)
2. Separates calculation from recording
3. Easy to test independently
4. Can add new record types without touching algorithms

---

## Part 3: Template Assignment Improvements

### 3.1 Current Issues

**TemplateAssignment** is good but could be better:
- Hardcoded nucleotide tables (should be data-driven)
- Manual overrides scattered in code
- Not easily extensible

### 3.2 Enhanced Template Assignment

#### TemplateAssignmentRegistry (NEW)
**File**: `include/x3dna/algorithms/template_assignment_registry.hpp`  
**Responsibility**: Centralized template assignment rules  
**Size**: ~150 lines

```cpp
struct TemplateRule {
    std::string residue_name;
    core::ResidueType assigned_type;
    std::string template_file;
    std::string reason;  // Why this override is needed
};

class TemplateAssignmentRegistry {
public:
    static TemplateAssignmentRegistry& instance();
    
    // Query
    std::optional<TemplateRule> get_rule(const std::string& residue_name) const;
    bool has_override(const std::string& residue_name) const;
    
    // Load from JSON
    void load_from_file(const std::filesystem::path& config_file);
    
    // Registration
    void register_rule(const TemplateRule& rule);
    
private:
    std::unordered_map<std::string, TemplateRule> rules_;
};
```

**Configuration file**: `resources/config/template_assignment_rules.json`

```json
{
  "rules": [
    {
      "residue_name": "A23",
      "assigned_type": "ADENINE",
      "template_file": "Atomic.a.pdb",
      "reason": "2'-deoxy-2'-fluoroadenosine - auto-detection fails due to missing purine atoms",
      "citation": "STAGE2_COMPLETE_FINAL.md - Memory ID: 11882176"
    },
    {
      "residue_name": "DI",
      "assigned_type": "INOSINE",
      "template_file": "Atomic_I.pdb",
      "reason": "2'-Deoxyinosine - was incorrectly classified as Guanine"
    }
  ]
}
```

**Benefits:**
1. ✅ Data-driven (edit JSON, not code)
2. ✅ Self-documenting (includes reasons)
3. ✅ Easy to audit and update
4. ✅ Version controllable
5. ✅ Can generate reports of all overrides

---

## Part 4: Configuration Management

### 4.1 Current Issues

Configuration is scattered:
- Template paths hardcoded
- RNA detection done multiple times
- Legacy mode flags passed around
- No central configuration

### 4.2 ConfigurationManager (Enhanced)

**File**: `include/x3dna/config/configuration_manager.hpp`

```cpp
struct FrameCalculationConfig {
    std::filesystem::path template_path = "data/templates";
    bool legacy_mode = false;
    bool auto_detect_rna = true;
    bool force_rna = false;
    double rmsd_tolerance = 0.05;
    
    // Template assignment
    std::filesystem::path template_rules_file = 
        "resources/config/template_assignment_rules.json";
    
    // Modified nucleotides
    std::filesystem::path modified_nucleotides_file = 
        "resources/config/modified_nucleotides.json";
    
    // Special residues
    std::filesystem::path special_residues_file = 
        "resources/config/special_residues.json";
};

class ConfigurationManager {
public:
    static ConfigurationManager& instance();
    
    // Load configuration
    void load_from_file(const std::filesystem::path& config_file);
    void load_defaults();
    
    // Access configuration
    const FrameCalculationConfig& frame_calculation_config() const;
    
    // Update configuration
    void set_legacy_mode(bool legacy_mode);
    void set_template_path(const std::filesystem::path& path);
    
private:
    FrameCalculationConfig frame_calc_config_;
    // ... other configs ...
};
```

**Benefits:**
1. Single source of truth for configuration
2. Easy to override from command line
3. Can load from config files
4. Environment-specific configs (dev/test/prod)

---

## Part 5: Error Handling and Validation

### 5.1 Current Issues

Error handling is inconsistent:
- Some functions throw exceptions
- Some return `std::optional`
- Some return invalid results
- Limited error context

### 5.2 Structured Error Handling

#### FrameCalculationError (NEW)
**File**: `include/x3dna/algorithms/frame_calculation_error.hpp`

```cpp
enum class FrameErrorCode {
    SUCCESS,
    INSUFFICIENT_ATOMS,
    TEMPLATE_NOT_FOUND,
    RMSD_FIT_FAILED,
    INVALID_RESIDUE_TYPE,
    NON_NUCLEOTIDE,
    BUFFER_MOLECULE,
    UNKNOWN_ERROR
};

struct FrameCalculationError {
    FrameErrorCode code;
    std::string message;
    std::string residue_name;
    char chain_id;
    int seq_num;
    
    std::string to_string() const;
    bool is_fatal() const;
};

class FrameCalculationResult {
public:
    // ... existing fields ...
    
    // Error handling
    bool is_valid = false;
    std::optional<FrameCalculationError> error;
    
    // Convenience
    bool succeeded() const { return is_valid && !error.has_value(); }
    bool failed() const { return !succeeded(); }
    std::string error_message() const;
};
```

**Benefits:**
1. Clear error codes
2. Detailed error context
3. Can log/report errors systematically
4. Easy debugging

---

## Part 6: Testing Infrastructure

### 6.1 Unit Tests (NEW)

Each class should have comprehensive unit tests:

```
tests/
├─ unit/
│  ├─ algorithms/
│  │  ├─ test_base_frame_calculator.cpp
│  │  ├─ test_ring_atom_matcher.cpp
│  │  ├─ test_residue_type_detector.cpp
│  │  └─ test_template_assignment_registry.cpp
│  ├─ services/
│  │  └─ test_frame_calculation_service.cpp
│  └─ io/
│     └─ test_frame_json_recorder.cpp
├─ integration/
│  ├─ test_stage2_complete.cpp
│  └─ test_modified_nucleotides.cpp
└─ fixtures/
   ├─ small_rna.pdb
   ├─ modified_nucleotides.pdb
   └─ edge_cases.pdb
```

### 6.2 Test Framework

Use **Google Test** or **Catch2**:

```cpp
// Example unit test
TEST(RingAtomMatcher, MatchesStandardPurine) {
    // Given
    auto residue = create_test_residue("A", {
        {" C4 ", {0, 0, 0}},
        {" N3 ", {1, 0, 0}},
        // ... 
    });
    
    RingAtomMatcher matcher;
    
    // When
    auto result = matcher.match(residue);
    
    // Then
    EXPECT_EQ(result.matched_atom_names.size(), 9);
    EXPECT_TRUE(result.has_purine_atoms);
    EXPECT_EQ(result.purine_atom_count, 3);
}
```

---

## Part 7: Documentation

### 7.1 Code Documentation

All public methods should have Doxygen comments:

```cpp
/**
 * @brief Calculate reference frame for a nucleotide residue
 * 
 * Performs least-squares fitting of experimental ring atoms to standard
 * template geometry. The resulting transformation defines the reference
 * frame for the residue.
 * 
 * @param residue The residue to calculate frame for (will be modified)
 * @return FrameCalculationResult containing frame and fit quality
 * 
 * @throws std::runtime_error if template file cannot be loaded
 * 
 * @note For modified nucleotides, template assignment is determined by
 *       ModifiedNucleotideRegistry and TemplateAssignmentRegistry
 * 
 * @see FrameCalculationResult
 * @see ModifiedNucleotideRegistry
 * @see TemplateAssignmentRegistry
 */
FrameCalculationResult calculate_frame(core::Residue& residue);
```

### 7.2 Architecture Documentation

**File**: `docs/ARCHITECTURE_STAGE2.md`

Should include:
1. High-level architecture diagram
2. Class responsibilities
3. Data flow diagrams
4. Sequence diagrams for key operations
5. Extension points
6. Common pitfalls

### 7.3 Developer Guide

**File**: `docs/DEVELOPER_GUIDE_STAGE2.md`

Should include:
1. How to add new modified nucleotides
2. How to add new template assignment rules
3. How to extend frame calculation
4. Debugging tips
5. Common issues and solutions

---

## Part 8: Migration Strategy

### Phase 1: Preparation (1 day)
1. ✅ Create comprehensive plan (this document)
2. ✅ Review plan with stakeholders
3. ✅ Set up feature branch: `refactor/stage2-modern-oop`
4. ✅ Create milestone in project tracker

### Phase 2: Extract Helper Classes (2-3 days)
1. Create `RingAtomMatcher`
   - Extract from `BaseFrameCalculator::check_nt_type_by_rmsd()`
   - Write unit tests
   - Integrate back into `BaseFrameCalculator`
   
2. Create `ResidueTypeDetector`
   - Extract type detection logic
   - Write unit tests
   - Integrate back into `BaseFrameCalculator`

3. Create `TemplateAssignmentRegistry`
   - Design JSON schema
   - Migrate hardcoded rules to JSON
   - Write unit tests
   - Update `TemplateAssignment` to use registry

**Validation**: Run Stage 2 tests - should still pass 3,602/3,602

### Phase 3: Create Service Layer (2 days)
1. Create `FrameCalculationService`
   - Implement high-level API
   - Move orchestration logic from `BaseFrameCalculator`
   - Write unit tests

2. Update `BaseFrameCalculator`
   - Remove orchestration methods
   - Keep only core algorithm
   - Simplify to ~400 lines

**Validation**: Run Stage 2 tests - should still pass 3,602/3,602

### Phase 4: Create Recording Layer (2 days)
1. Create `FrameJsonRecorder`
   - Implement iteration logic (DRY)
   - Implement individual record methods
   - Write unit tests

2. Remove recording from `BaseFrameCalculator`
   - Delete `calculate_and_record_frames()`
   - Delete `calculate_and_record_frames_only()`

3. Delete `LsFittingCalculator`
   - Remove files entirely
   - Update CMakeLists.txt

**Validation**: Run Stage 2 tests - should still pass 3,602/3,602

### Phase 5: Update Tools (1 day)
1. Update `generate_modern_json.cpp`
   - Use `FrameCalculationService`
   - Use `FrameJsonRecorder`
   - Extract helper functions
   - Reduce from 793 lines to ~200-300 lines

**Validation**: Run full Stage 2 validation - should pass 3,602/3,602

### Phase 6: Enhanced Error Handling (1 day)
1. Create `FrameCalculationError`
2. Update `FrameCalculationResult`
3. Update all classes to use structured errors

**Validation**: Run Stage 2 tests - should still pass 3,602/3,602

### Phase 7: Configuration Management (1 day)
1. Enhance `ConfigurationManager`
2. Create `template_assignment_rules.json`
3. Update classes to use configuration

**Validation**: Run Stage 2 tests - should still pass 3,602/3,602

### Phase 8: Documentation (2 days)
1. Add Doxygen comments to all public methods
2. Write `ARCHITECTURE_STAGE2.md`
3. Write `DEVELOPER_GUIDE_STAGE2.md`
4. Update existing documentation

### Phase 9: Testing (2 days)
1. Write comprehensive unit tests
2. Write integration tests
3. Achieve 90%+ code coverage

### Phase 10: Cleanup & Polish (1 day)
1. Remove dead code
2. Fix linter warnings
3. Run clang-format
4. Final review

### Phase 11: Merge & Deploy (1 day)
1. Final validation: 3,602/3,602 PDBs ✅
2. Code review
3. Merge to main
4. Tag release: `v2.0.0-stage2-refactored`
5. Update deployment documentation

**Total Estimate**: 15-18 days

---

## Part 9: Success Metrics

### Quantitative Metrics

| Metric | Before | Target | Measure |
|--------|--------|--------|---------|
| `BaseFrameCalculator` LOC | 933 | 400 | Lines of code |
| Number of classes | 2 | 6 | Class count |
| Code duplication | High | None | DRY violations |
| Test coverage | ~20% | 90%+ | Line coverage |
| Validation pass rate | 100% | 100% | Must maintain |
| Build time | Baseline | <110% | Seconds |
| `generate_modern_json.cpp` LOC | 793 | 200-300 | Lines of code |

### Qualitative Metrics

- [ ] Clear separation of concerns
- [ ] Easy to understand for new developers
- [ ] Easy to add new modified nucleotides
- [ ] Easy to add new template rules
- [ ] Well documented
- [ ] Comprehensive tests
- [ ] No tight coupling between layers
- [ ] Follows SOLID principles

---

## Part 10: Risks and Mitigation

### Risk 1: Breaking Validation
**Impact**: High  
**Likelihood**: Medium  
**Mitigation**: 
- Run validation after each phase
- Keep validation threshold at 100% (3,602/3,602)
- Automated testing in CI/CD

### Risk 2: Scope Creep
**Impact**: Medium  
**Likelihood**: High  
**Mitigation**:
- Stick to the plan
- No new features during refactoring
- Track scope changes explicitly

### Risk 3: Regression Bugs
**Impact**: High  
**Likelihood**: Medium  
**Mitigation**:
- Comprehensive unit tests
- Integration tests
- Before/after comparisons
- JSON output validation

### Risk 4: Performance Degradation
**Impact**: Medium  
**Likelihood**: Low  
**Mitigation**:
- Benchmark before/after
- Profile hot paths
- Keep algorithm core unchanged

### Risk 5: Over-Engineering
**Impact**: Medium  
**Likelihood**: Medium  
**Mitigation**:
- Follow YAGNI principle
- Keep it simple
- Focus on real extensibility needs

---

## Part 11: Post-Refactoring Benefits

### For Development
1. ✅ **Faster feature development**: Clear extension points
2. ✅ **Easier debugging**: Better error messages and logging
3. ✅ **Better testing**: Can test components independently
4. ✅ **Less cognitive load**: Each class has single responsibility

### For Maintenance
1. ✅ **Easier to understand**: Clear architecture
2. ✅ **Easier to modify**: Changes localized to single class
3. ✅ **Easier to extend**: Well-defined interfaces
4. ✅ **Better documentation**: Self-documenting code + docs

### For Stage 3+
1. ✅ **Template for future stages**: Same pattern for other algorithms
2. ✅ **Reusable components**: Services and recorders for other stages
3. ✅ **Consistent architecture**: All stages follow same pattern

---

## Part 12: Follow-up Refactorings

After Stage 2, apply same pattern to:

### Stage 3: Distance Calculations
- `DistanceCalculationService`
- `DistanceJsonRecorder`

### Stage 4: H-Bond Detection
- `HydrogenBondService`
- `HydrogenBondJsonRecorder`

### Stage 5+: Pair Validation
- `PairValidationService`
- `PairJsonRecorder`

**Goal**: Consistent architecture across ALL stages

---

## Summary

### Before Refactoring
- ❌ 933-line `BaseFrameCalculator` doing too much
- ❌ Confusing `LsFittingCalculator` wrapper
- ❌ Code duplication in 3+ places
- ❌ Mixed responsibilities (calculation + JSON)
- ❌ Hard to test
- ❌ Hard to extend
- ❌ 793-line `generate_modern_json.cpp`

### After Refactoring
- ✅ Clean separation: Algorithms → Services → Recorders → I/O
- ✅ Single Responsibility Principle throughout
- ✅ No code duplication (DRY)
- ✅ Easy to test (each layer independently)
- ✅ Easy to extend (data-driven configuration)
- ✅ Well documented (architecture + developer guide)
- ✅ Comprehensive tests (90%+ coverage)
- ✅ Modern OOP (SOLID principles)
- ✅ 100% validation maintained (3,602/3,602 PDBs)
- ✅ Template for future stages

### Key Achievement
**Production-ready, maintainable, extensible architecture that serves as the foundation for all future development.**

---

**Next Steps**: Review this plan, get approval, create feature branch, and begin Phase 1.

