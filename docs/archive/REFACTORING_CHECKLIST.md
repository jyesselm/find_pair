# Stage 2 Refactoring Implementation Checklist

**Date**: December 5, 2025  
**Feature Branch**: `refactor/stage2-modern-oop`  
**Validation Target**: 3,602/3,602 PDBs (100%) - MUST MAINTAIN

---

## Quick Reference

```
Current State              Target State
──────────────────────────────────────────────────────────
BaseFrameCalculator        → BaseFrameCalculator (pure algorithm, 400 LOC)
  (933 LOC, mixed)         → RingAtomMatcher (NEW, 150 LOC)
                          → ResidueTypeDetector (NEW, 200 LOC)
                          → FrameCalculationService (NEW, 200 LOC)
                          → FrameJsonRecorder (NEW, 250 LOC)

LsFittingCalculator       → DELETED ❌
  (wrapper, ~100 LOC)

generate_modern_json      → generate_modern_json (simplified, 200-300 LOC)
  (793 LOC, duplicated)

TemplateAssignment        → TemplateAssignment + TemplateAssignmentRegistry
  (hardcoded)               (data-driven JSON config)
```

---

## Phase 1: Preparation ✅

### Setup
- [ ] Review comprehensive plan (`docs/COMPREHENSIVE_STAGE2_REFACTORING_PLAN.md`)
- [ ] Create feature branch: `git checkout -b refactor/stage2-modern-oop`
- [ ] Ensure baseline validation passes: `pytest tests_python/ -k stage2`
  - [ ] Expected: 3,602/3,602 PDBs passing
- [ ] Create backup tag: `git tag -a v1.0-pre-refactor -m "Backup before Stage 2 refactoring"`

### Documentation
- [ ] Read existing docs:
  - [ ] `docs/REFACTOR_FRAME_CALCULATORS.md`
  - [ ] `docs/STAGE2_COMPLETE_FINAL.md`
  - [ ] `docs/CODE_FLOW.md`

**Estimated Time**: 0.5 days  
**Validation**: Baseline tests pass

---

## Phase 2: Extract Helper Classes

### 2.1: Create RingAtomMatcher

#### Files to Create
- [ ] `include/x3dna/algorithms/ring_atom_matcher.hpp`
- [ ] `src/x3dna/algorithms/ring_atom_matcher.cpp`
- [ ] `tests/unit/algorithms/test_ring_atom_matcher.cpp`

#### Implementation Steps
1. [ ] Define `RingMatchResult` struct
   ```cpp
   struct RingMatchResult {
       std::vector<std::string> matched_atom_names;
       std::vector<geometry::Vector3D> experimental_coords;
       std::vector<geometry::Vector3D> standard_coords;
       bool has_purine_atoms;
       int purine_atom_count;
   };
   ```

2. [ ] Extract logic from `BaseFrameCalculator::check_nt_type_by_rmsd()`
   - [ ] Copy ring atom matching logic (lines 52-164)
   - [ ] Refactor into clean methods
   - [ ] Add Doxygen comments

3. [ ] Write unit tests
   - [ ] Test matching standard purine (9 atoms)
   - [ ] Test matching standard pyrimidine (6 atoms)
   - [ ] Test matching incomplete residue
   - [ ] Test modified nucleotide (e.g., A23)
   - [ ] **Target**: 90%+ coverage

4. [ ] Integrate into `BaseFrameCalculator`
   - [ ] Add `RingAtomMatcher matcher_;` member
   - [ ] Replace inline logic with `matcher_.match(residue)`
   - [ ] Update `check_nt_type_by_rmsd()` to use matcher

5. [ ] **Validate**: Run Stage 2 tests
   ```bash
   cmake --build build --target generate_modern_json -j8
   pytest tests_python/ -k stage2 -v
   ```
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

6. [ ] Commit
   ```bash
   git add include/x3dna/algorithms/ring_atom_matcher.*
   git add src/x3dna/algorithms/ring_atom_matcher.cpp
   git add tests/unit/algorithms/test_ring_atom_matcher.cpp
   git commit -m "refactor: Extract RingAtomMatcher from BaseFrameCalculator

   - Separates ring atom matching into dedicated class
   - Improves testability and reusability
   - No functional changes
   - Validation: 3,602/3,602 PDBs passing"
   ```

**Estimated Time**: 1 day

---

### 2.2: Create ResidueTypeDetector

#### Files to Create
- [ ] `include/x3dna/algorithms/residue_type_detector.hpp`
- [ ] `src/x3dna/algorithms/residue_type_detector.cpp`
- [ ] `tests/unit/algorithms/test_residue_type_detector.cpp`

#### Implementation Steps
1. [ ] Define `TypeDetectionResult` struct
   ```cpp
   struct TypeDetectionResult {
       core::ResidueType detected_type;
       bool used_fallback;
       std::optional<double> rmsd;
       std::string detection_method;
   };
   ```

2. [ ] Extract type detection logic from `BaseFrameCalculator`
   - [ ] Lines 187-530 (type detection)
   - [ ] Include RMSD check logic
   - [ ] Include atom-based detection
   - [ ] Include registry check

3. [ ] Implement methods
   - [ ] `detect_type(const Residue& residue)`
   - [ ] `should_use_rmsd_check(const Residue& residue)`
   - [ ] `detect_by_atoms(const Residue& residue)`
   - [ ] `detect_by_rmsd(const Residue& residue)`

4. [ ] Write unit tests
   - [ ] Test standard nucleotides (A, C, G, T, U)
   - [ ] Test modified nucleotides (A23, DI, EPE)
   - [ ] Test non-nucleotides (buffer molecules)
   - [ ] Test edge cases (QUO, 9DG)
   - [ ] **Target**: 90%+ coverage

5. [ ] Integrate into `BaseFrameCalculator`
   - [ ] Add `ResidueTypeDetector type_detector_;` member
   - [ ] Replace inline logic with `type_detector_.detect_type(residue)`
   - [ ] Simplify `calculate_frame_impl()`

6. [ ] **Validate**: Run Stage 2 tests
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

7. [ ] Commit
   ```bash
   git commit -m "refactor: Extract ResidueTypeDetector from BaseFrameCalculator

   - Separates type detection into dedicated class
   - Handles RMSD-based and atom-based detection
   - Improves testability
   - Validation: 3,602/3,602 PDBs passing"
   ```

**Estimated Time**: 1.5 days

---

### 2.3: Create TemplateAssignmentRegistry

#### Files to Create
- [ ] `include/x3dna/algorithms/template_assignment_registry.hpp`
- [ ] `src/x3dna/algorithms/template_assignment_registry.cpp`
- [ ] `resources/config/template_assignment_rules.json`
- [ ] `tests/unit/algorithms/test_template_assignment_registry.cpp`

#### Implementation Steps
1. [ ] Define `TemplateRule` struct
   ```cpp
   struct TemplateRule {
       std::string residue_name;
       core::ResidueType assigned_type;
       std::string template_file;
       std::string reason;
   };
   ```

2. [ ] Implement singleton registry
   - [ ] `static TemplateAssignmentRegistry& instance()`
   - [ ] `get_rule(const std::string& residue_name)`
   - [ ] `has_override(const std::string& residue_name)`
   - [ ] `load_from_file(const std::filesystem::path&)`
   - [ ] `register_rule(const TemplateRule&)`

3. [ ] Create JSON configuration
   ```json
   {
     "rules": [
       {
         "residue_name": "A23",
         "assigned_type": "ADENINE",
         "template_file": "Atomic.a.pdb",
         "reason": "2'-deoxy-2'-fluoroadenosine - auto-detection fails"
       }
     ]
   }
   ```

4. [ ] Migrate existing hardcoded rules
   - [ ] From `TemplateAssignment` class
   - [ ] From memory (A23, DI, etc.)
   - [ ] Document each rule

5. [ ] Write unit tests
   - [ ] Test loading from JSON
   - [ ] Test rule lookup
   - [ ] Test missing rules
   - [ ] Test rule registration
   - [ ] **Target**: 90%+ coverage

6. [ ] Update `TemplateAssignment` to use registry
   - [ ] Check registry first before hardcoded logic
   - [ ] Fall back to auto-detection if no rule

7. [ ] **Validate**: Run Stage 2 tests
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

8. [ ] Commit
   ```bash
   git commit -m "refactor: Create data-driven TemplateAssignmentRegistry

   - Moves hardcoded template rules to JSON config
   - Makes rules easy to update and audit
   - Self-documenting with reasons for each override
   - Validation: 3,602/3,602 PDBs passing"
   ```

**Estimated Time**: 1 day

**Total Phase 2 Time**: 3.5 days

---

## Phase 3: Create Service Layer

### 3.1: Create FrameCalculationService

#### Files to Create
- [ ] `include/x3dna/services/frame_calculation_service.hpp`
- [ ] `src/x3dna/services/frame_calculation_service.cpp`
- [ ] `tests/unit/services/test_frame_calculation_service.cpp`

#### Implementation Steps
1. [ ] Create service directory structure
   ```bash
   mkdir -p include/x3dna/services
   mkdir -p src/x3dna/services
   mkdir -p tests/unit/services
   ```

2. [ ] Define service interface
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

3. [ ] Move orchestration logic from `BaseFrameCalculator`
   - [ ] Extract `calculate_all_frames()` logic
   - [ ] Add RNA auto-detection
   - [ ] Add validation methods

4. [ ] Write unit tests
   - [ ] Test batch calculation
   - [ ] Test RNA detection
   - [ ] Test legacy mode
   - [ ] Test validation
   - [ ] **Target**: 90%+ coverage

5. [ ] Update `BaseFrameCalculator`
   - [ ] Remove orchestration methods
   - [ ] Keep only core `calculate_frame()` and `calculate_frame_impl()`
   - [ ] Should be ~400 lines now

6. [ ] **Validate**: Run Stage 2 tests
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

7. [ ] Commit
   ```bash
   git commit -m "refactor: Create FrameCalculationService layer

   - Separates orchestration from algorithm
   - Provides high-level API for frame calculation
   - BaseFrameCalculator reduced to 400 LOC (from 933)
   - Validation: 3,602/3,602 PDBs passing"
   ```

**Estimated Time**: 2 days

---

## Phase 4: Create Recording Layer

### 4.1: Create FrameJsonRecorder

#### Files to Create
- [ ] `include/x3dna/io/frame_json_recorder.hpp`
- [ ] `src/x3dna/io/frame_json_recorder.cpp`
- [ ] `tests/unit/io/test_frame_json_recorder.cpp`

#### Implementation Steps
1. [ ] Define recorder interface
   ```cpp
   class FrameJsonRecorder {
   public:
       explicit FrameJsonRecorder(FrameCalculationService& service);
       
       size_t record_ls_fitting(core::Structure& structure, io::JsonWriter& writer);
       size_t record_base_frame_calc(core::Structure& structure, io::JsonWriter& writer);
       size_t record_frame_calc(core::Structure& structure, io::JsonWriter& writer);
       size_t record_all(core::Structure& structure, io::JsonWriter& writer);
       
   private:
       FrameCalculationService& service_;
       
       template<typename RecordFunc>
       size_t iterate_and_record(core::Structure& structure, 
                                  io::JsonWriter& writer, 
                                  RecordFunc record_func);
   };
   ```

2. [ ] Implement iteration logic (DRY - only ONE place)
   - [ ] Extract from `BaseFrameCalculator::calculate_and_record_frames()`
   - [ ] Extract from `LsFittingCalculator::calculate_and_record()`
   - [ ] Consolidate into `iterate_and_record()` template

3. [ ] Implement individual record methods
   - [ ] `record_ls_fitting()`: writes ls_fitting JSON
   - [ ] `record_base_frame_calc()`: writes base_frame_calc JSON
   - [ ] `record_frame_calc()`: writes frame_calc JSON
   - [ ] `record_all()`: convenience method

4. [ ] Write unit tests
   - [ ] Test ls_fitting recording
   - [ ] Test base_frame_calc recording
   - [ ] Test frame_calc recording
   - [ ] Test record_all
   - [ ] Test iteration logic
   - [ ] **Target**: 90%+ coverage

5. [ ] **Validate**: Build and test
   ```bash
   cmake --build build --target frame_json_recorder -j8
   cmake --build build --target test_frame_json_recorder -j8
   ./build/tests/unit/io/test_frame_json_recorder
   ```

6. [ ] Commit
   ```bash
   git commit -m "refactor: Create FrameJsonRecorder layer

   - Separates JSON recording from calculation
   - Consolidates iteration logic (DRY)
   - Easy to test independently
   - Validation: Unit tests passing"
   ```

**Estimated Time**: 1.5 days

---

### 4.2: Remove Recording from BaseFrameCalculator

#### Files to Modify
- [ ] `include/x3dna/algorithms/base_frame_calculator.hpp`
- [ ] `src/x3dna/algorithms/base_frame_calculator.cpp`

#### Implementation Steps
1. [ ] Remove methods from header
   - [ ] Delete `calculate_and_record_frames()` declaration
   - [ ] Delete `calculate_and_record_frames_only()` declaration

2. [ ] Remove methods from implementation
   - [ ] Delete `calculate_and_record_frames()` (lines 531-591)
   - [ ] Delete `calculate_and_record_frames_only()` (lines 593-649)
   - [ ] This should reduce file to ~400 lines ✅

3. [ ] Clean up includes
   - [ ] Remove `#include <x3dna/io/json_writer.hpp>` if not needed
   - [ ] Add any new includes for helper classes

4. [ ] **Validate**: Build
   ```bash
   cmake --build build --target base_frame_calculator -j8
   ```

5. [ ] Commit
   ```bash
   git commit -m "refactor: Remove JSON recording from BaseFrameCalculator

   - Deleted calculate_and_record_frames()
   - Deleted calculate_and_record_frames_only()
   - BaseFrameCalculator now pure algorithm (400 LOC)
   - Recording moved to FrameJsonRecorder"
   ```

**Estimated Time**: 0.5 days

---

### 4.3: Delete LsFittingCalculator

#### Files to Delete
- [ ] `include/x3dna/algorithms/ls_fitting_calculator.hpp`
- [ ] `src/x3dna/algorithms/ls_fitting_calculator.cpp`

#### Implementation Steps
1. [ ] Remove from CMakeLists.txt
   ```cmake
   # Find and remove:
   # src/x3dna/algorithms/ls_fitting_calculator.cpp
   ```

2. [ ] Delete files
   ```bash
   git rm include/x3dna/algorithms/ls_fitting_calculator.hpp
   git rm src/x3dna/algorithms/ls_fitting_calculator.cpp
   ```

3. [ ] Find all usages and replace
   ```bash
   grep -r "LsFittingCalculator" --include="*.cpp" --include="*.hpp"
   ```
   - [ ] Update `generate_modern_json.cpp`
   - [ ] Update any other files

4. [ ] **Validate**: Build
   ```bash
   cmake --build build -j8
   ```

5. [ ] Commit
   ```bash
   git commit -m "refactor: Delete LsFittingCalculator (unnecessary wrapper)

   - Was just wrapping BaseFrameCalculator
   - Functionality moved to FrameJsonRecorder
   - Simplifies architecture"
   ```

**Estimated Time**: 0.5 days

**Total Phase 4 Time**: 2.5 days

---

## Phase 5: Update Tools

### 5.1: Refactor generate_modern_json.cpp

#### Files to Modify
- [ ] `tools/generate_modern_json.cpp`

#### Implementation Steps

1. [ ] Create helper functions (top of file)
   ```cpp
   // Utility: RNA detection
   bool detect_rna_structure(const core::Structure& structure) {
       return FrameCalculationService::detect_rna(structure);
   }
   
   // Helper: Setup frame calculation
   std::tuple<FrameCalculationService, io::JsonWriter> 
   setup_frame_calculation(const std::filesystem::path& pdb_file,
                           bool legacy_mode,
                           core::Structure& structure) {
       FrameCalculationService service("data/templates");
       service.set_legacy_mode(legacy_mode);
       service.auto_detect_rna(structure);
       
       io::JsonWriter writer(pdb_file);
       writer.record_residue_indices(structure);
       
       return {std::move(service), std::move(writer)};
   }
   
   // Helper: Print stage result
   void print_stage_result(const std::string& stage_name,
                           const std::string& pdb_name,
                           size_t count) {
       std::cout << "  ✅ " << stage_name << "/" << pdb_name << ".json\n";
       std::cout << "     " << count << " records\n\n";
   }
   ```

2. [ ] Refactor Stage 3 (ls_fitting)
   ```cpp
   // Before (30 lines):
   if (stage == "ls_fitting" || stage == "all") {
       std::cout << "Stage 3: Writing ls_fitting...\n";
       JsonWriter writer(pdb_file);
       writer.record_residue_indices(structure);
       
       LsFittingCalculator calculator("data/templates");
       calculator.set_legacy_mode(legacy_mode);
       
       bool is_rna = LsFittingCalculator::detect_rna(structure);
       // ... more setup ...
       
       size_t count = calculator.calculate_and_record(structure, writer);
       writer.write_split_files(json_output_dir, true);
       // ... print result ...
   }
   
   // After (8 lines):
   if (stage == "ls_fitting" || stage == "all") {
       std::cout << "Stage 3: Writing ls_fitting...\n";
       auto [service, writer] = setup_frame_calculation(pdb_file, legacy_mode, structure);
       FrameJsonRecorder recorder(service);
       size_t count = recorder.record_ls_fitting(structure, writer);
       writer.write_split_files(json_output_dir, true);
       print_stage_result("ls_fitting", pdb_name, count);
   }
   ```

3. [ ] Refactor Stage 4 (frames)
   ```cpp
   // Similar simplification for base_frame_calc + frame_calc
   if (stage == "frames" || stage == "all") {
       std::cout << "Stage 4: Writing base_frame_calc and frame_calc...\n";
       auto [service, writer] = setup_frame_calculation(pdb_file, legacy_mode, structure);
       FrameJsonRecorder recorder(service);
       size_t count = recorder.record_base_frame_calc(structure, writer);
       count += recorder.record_frame_calc(structure, writer);
       writer.write_split_files(json_output_dir, true);
       print_stage_result("frames", pdb_name, count);
   }
   ```

4. [ ] Remove duplicate RNA detection code
   - [ ] Lines 203-209: Use `detect_rna_structure()` helper
   - [ ] Lines 236-241: Use `detect_rna_structure()` helper
   - [ ] Lines 286-297: Use `detect_rna_structure()` helper

5. [ ] Remove duplicate frame setup code
   - [ ] Use `setup_frame_calculation()` helper everywhere

6. [ ] Extract `build_paired_indices()` helper
   ```cpp
   std::set<int> build_paired_indices(const std::vector<core::BasePair>& pairs) {
       std::set<int> indices;
       for (const auto& pair : pairs) {
           // ... existing logic ...
       }
       return indices;
   }
   ```
   - [ ] Replace lines 544-572 with call to helper
   - [ ] Replace lines 574-602 with call to helper

7. [ ] Final cleanup
   - [ ] Remove commented code
   - [ ] Fix formatting
   - [ ] Add section comments

8. [ ] **Measure**: Count lines
   ```bash
   wc -l tools/generate_modern_json.cpp
   ```
   - [ ] Expected: 200-300 lines (down from 793) ✅

9. [ ] **Validate**: Run full Stage 2 validation
   ```bash
   cmake --build build --target generate_modern_json -j8
   python3 scripts/validate_frames_parallel.py
   ```
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

10. [ ] Commit
    ```bash
    git commit -m "refactor: Simplify generate_modern_json.cpp

    - Use FrameCalculationService and FrameJsonRecorder
    - Extract helper functions (DRY)
    - Reduce from 793 to ~250 lines (68% reduction)
    - No functional changes
    - Validation: 3,602/3,602 PDBs passing"
    ```

**Estimated Time**: 1 day

---

## Phase 6: Enhanced Error Handling

### 6.1: Create FrameCalculationError

#### Files to Create
- [ ] `include/x3dna/algorithms/frame_calculation_error.hpp`
- [ ] `src/x3dna/algorithms/frame_calculation_error.cpp`
- [ ] `tests/unit/algorithms/test_frame_calculation_error.cpp`

#### Implementation Steps
1. [ ] Define error codes
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
   ```

2. [ ] Define error struct
   ```cpp
   struct FrameCalculationError {
       FrameErrorCode code;
       std::string message;
       std::string residue_name;
       char chain_id;
       int seq_num;
       
       std::string to_string() const;
       bool is_fatal() const;
   };
   ```

3. [ ] Update `FrameCalculationResult`
   ```cpp
   struct FrameCalculationResult {
       // ... existing fields ...
       
       // Error handling
       std::optional<FrameCalculationError> error;
       
       // Convenience
       bool succeeded() const { return is_valid && !error.has_value(); }
       bool failed() const { return !succeeded(); }
       std::string error_message() const;
   };
   ```

4. [ ] Update all classes to use structured errors
   - [ ] `BaseFrameCalculator`
   - [ ] `ResidueTypeDetector`
   - [ ] `RingAtomMatcher`

5. [ ] Write unit tests
   - [ ] Test error creation
   - [ ] Test error messages
   - [ ] Test error codes

6. [ ] **Validate**: Run Stage 2 tests
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

7. [ ] Commit
   ```bash
   git commit -m "feat: Add structured error handling for frame calculation

    - Defined FrameErrorCode enum
    - Created FrameCalculationError struct
    - Updated FrameCalculationResult with error field
    - Better debugging and error reporting
    - Validation: 3,602/3,602 PDBs passing"
    ```

**Estimated Time**: 1 day

---

## Phase 7: Configuration Management

### 7.1: Enhance ConfigurationManager

#### Files to Modify
- [ ] `include/x3dna/config/config_manager.hpp`
- [ ] `src/x3dna/config/config_manager.cpp`

#### Files to Create
- [ ] `resources/config/frame_calculation.json`

#### Implementation Steps
1. [ ] Define `FrameCalculationConfig` struct
   ```cpp
   struct FrameCalculationConfig {
       std::filesystem::path template_path = "data/templates";
       bool legacy_mode = false;
       bool auto_detect_rna = true;
       double rmsd_tolerance = 0.05;
       std::filesystem::path template_rules_file = 
           "resources/config/template_assignment_rules.json";
   };
   ```

2. [ ] Add to `ConfigurationManager`
   ```cpp
   class ConfigurationManager {
   public:
       const FrameCalculationConfig& frame_calculation_config() const;
       void set_frame_calculation_config(const FrameCalculationConfig& config);
       
   private:
       FrameCalculationConfig frame_calc_config_;
   };
   ```

3. [ ] Create JSON config file
   ```json
   {
     "frame_calculation": {
       "template_path": "data/templates",
       "legacy_mode": false,
       "auto_detect_rna": true,
       "rmsd_tolerance": 0.05,
       "template_rules_file": "resources/config/template_assignment_rules.json"
     }
   }
   ```

4. [ ] Update classes to use configuration
   - [ ] `FrameCalculationService`
   - [ ] `BaseFrameCalculator`

5. [ ] **Validate**: Run Stage 2 tests
   - [ ] Expected: 3,602/3,602 PDBs passing ✅

6. [ ] Commit
   ```bash
   git commit -m "feat: Add centralized configuration for frame calculation

    - Created FrameCalculationConfig struct
    - Added to ConfigurationManager
    - Created JSON config file
    - Single source of truth for configuration
    - Validation: 3,602/3,602 PDBs passing"
    ```

**Estimated Time**: 1 day

---

## Phase 8: Documentation

### 8.1: Add Doxygen Comments

#### Files to Document
- [ ] `include/x3dna/algorithms/base_frame_calculator.hpp`
- [ ] `include/x3dna/algorithms/ring_atom_matcher.hpp`
- [ ] `include/x3dna/algorithms/residue_type_detector.hpp`
- [ ] `include/x3dna/algorithms/template_assignment_registry.hpp`
- [ ] `include/x3dna/services/frame_calculation_service.hpp`
- [ ] `include/x3dna/io/frame_json_recorder.hpp`

#### Documentation Standards
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
 */
FrameCalculationResult calculate_frame(core::Residue& residue);
```

**Estimated Time**: 0.5 days

---

### 8.2: Write Architecture Documentation

#### File to Create
- [ ] `docs/ARCHITECTURE_STAGE2.md`

#### Content Outline
1. [ ] Overview and objectives
2. [ ] Architecture diagram
3. [ ] Layer descriptions
   - [ ] Algorithm layer
   - [ ] Service layer
   - [ ] Recording layer
   - [ ] I/O layer
4. [ ] Class responsibilities
5. [ ] Data flow diagrams
6. [ ] Sequence diagrams
7. [ ] Extension points
8. [ ] Common pitfalls

**Estimated Time**: 0.5 days

---

### 8.3: Write Developer Guide

#### File to Create
- [ ] `docs/DEVELOPER_GUIDE_STAGE2.md`

#### Content Outline
1. [ ] Quick start
2. [ ] How to add new modified nucleotides
3. [ ] How to add new template assignment rules
4. [ ] How to extend frame calculation
5. [ ] Debugging tips
6. [ ] Common issues and solutions
7. [ ] FAQ

**Estimated Time**: 0.5 days

---

### 8.4: Update Existing Documentation

#### Files to Update
- [ ] `docs/README.md` - Add link to new docs
- [ ] `docs/START_HERE.md` - Update Stage 2 section
- [ ] `docs/CODE_FLOW.md` - Update architecture diagrams
- [ ] `docs/STAGE2_COMPLETE_FINAL.md` - Add refactoring notes

**Estimated Time**: 0.5 days

**Total Phase 8 Time**: 2 days

---

## Phase 9: Testing

### 9.1: Write Unit Tests

#### Test Files to Create
- [ ] `tests/unit/algorithms/test_base_frame_calculator.cpp`
- [ ] `tests/unit/algorithms/test_ring_atom_matcher.cpp` ✅ (Phase 2)
- [ ] `tests/unit/algorithms/test_residue_type_detector.cpp` ✅ (Phase 2)
- [ ] `tests/unit/algorithms/test_template_assignment_registry.cpp` ✅ (Phase 2)
- [ ] `tests/unit/services/test_frame_calculation_service.cpp` ✅ (Phase 3)
- [ ] `tests/unit/io/test_frame_json_recorder.cpp` ✅ (Phase 4)

#### Coverage Target
- [ ] Run coverage analysis
  ```bash
  cmake -DCMAKE_BUILD_TYPE=Coverage ..
  cmake --build build --target coverage -j8
  ```
- [ ] **Target**: 90%+ line coverage ✅

**Estimated Time**: 1 day

---

### 9.2: Write Integration Tests

#### Test Files to Create
- [ ] `tests/integration/test_stage2_complete.cpp`
- [ ] `tests/integration/test_modified_nucleotides.cpp`
- [ ] `tests/integration/test_edge_cases.cpp`

#### Test Scenarios
- [ ] Test complete pipeline (PDB → frame calculation → JSON)
- [ ] Test all modified nucleotides (A23, DI, EPE, etc.)
- [ ] Test edge cases (QUO, 9DG, CM0, etc.)
- [ ] Test buffer molecule filtering
- [ ] Test RNA vs DNA detection

**Estimated Time**: 1 day

**Total Phase 9 Time**: 2 days

---

## Phase 10: Cleanup & Polish

### 10.1: Remove Dead Code

#### Files to Check
- [ ] Search for unused functions
  ```bash
  # Use clang tooling or manual review
  find src include -name "*.cpp" -o -name "*.hpp" | xargs grep -l "TODO\|FIXME\|XXX\|DEPRECATED"
  ```
- [ ] Remove commented-out code
- [ ] Remove debug print statements

**Estimated Time**: 0.25 days

---

### 10.2: Fix Linter Warnings

#### Steps
1. [ ] Run clang-tidy
   ```bash
   cmake --build build --target clang-tidy -j8
   ```

2. [ ] Fix warnings
   - [ ] `base_frame_calculator.cpp`
   - [ ] `ring_atom_matcher.cpp`
   - [ ] `residue_type_detector.cpp`
   - [ ] `frame_calculation_service.cpp`
   - [ ] `frame_json_recorder.cpp`
   - [ ] `generate_modern_json.cpp`

3. [ ] **Target**: Zero warnings ✅

**Estimated Time**: 0.25 days

---

### 10.3: Format Code

#### Steps
1. [ ] Run clang-format on all modified files
   ```bash
   find include src -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i
   ```

2. [ ] Verify formatting
   ```bash
   git diff
   ```

3. [ ] Commit formatting separately
   ```bash
   git commit -m "style: Format code with clang-format"
   ```

**Estimated Time**: 0.25 days

---

### 10.4: Final Review

#### Checklist
- [ ] All unit tests passing
- [ ] All integration tests passing
- [ ] Coverage ≥ 90%
- [ ] Zero linter warnings
- [ ] Code formatted
- [ ] Documentation complete
- [ ] No dead code
- [ ] No TODOs in production code

**Estimated Time**: 0.25 days

**Total Phase 10 Time**: 1 day

---

## Phase 11: Merge & Deploy

### 11.1: Final Validation

#### Steps
1. [ ] Run complete Stage 2 validation
   ```bash
   python3 scripts/validate_frames_parallel.py
   ```
   - [ ] **Expected**: 3,602/3,602 PDBs passing ✅

2. [ ] Run all unit tests
   ```bash
   cmake --build build --target test -j8
   ./build/run_tests
   ```
   - [ ] **Expected**: All tests passing ✅

3. [ ] Build in release mode
   ```bash
   cmake -DCMAKE_BUILD_TYPE=Release ..
   cmake --build build -j8
   ```
   - [ ] **Expected**: Clean build ✅

4. [ ] Performance benchmark
   ```bash
   time ./build/tools/generate_modern_json --pdb-file test.pdb
   ```
   - [ ] **Expected**: ≤110% of baseline time ✅

**Estimated Time**: 0.25 days

---

### 11.2: Code Review

#### Steps
1. [ ] Self-review entire diff
   ```bash
   git diff main..refactor/stage2-modern-oop
   ```

2. [ ] Create pull request
   - [ ] Title: "Refactor Stage 2: Modern OOP Architecture"
   - [ ] Description: Link to `docs/COMPREHENSIVE_STAGE2_REFACTORING_PLAN.md`
   - [ ] Reviewers: Assign team members

3. [ ] Address review comments

4. [ ] Get approval

**Estimated Time**: 0.5 days

---

### 11.3: Merge to Main

#### Steps
1. [ ] Squash commits if needed
   ```bash
   git rebase -i main
   ```

2. [ ] Merge to main
   ```bash
   git checkout main
   git merge --no-ff refactor/stage2-modern-oop
   ```

3. [ ] Tag release
   ```bash
   git tag -a v2.0.0-stage2-refactored -m "Stage 2 refactored with modern OOP architecture

   - Separated concerns: Algorithms → Services → Recorders
   - Created 6 new classes
   - Reduced BaseFrameCalculator from 933 to 400 LOC
   - Reduced generate_modern_json from 793 to ~250 LOC
   - Eliminated code duplication
   - Added comprehensive tests (90%+ coverage)
   - 100% validation maintained (3,602/3,602 PDBs)
   "
   git push origin v2.0.0-stage2-refactored
   ```

4. [ ] Push to remote
   ```bash
   git push origin main
   ```

**Estimated Time**: 0.25 days

**Total Phase 11 Time**: 1 day

---

## Total Estimated Time

| Phase | Description | Days |
|-------|-------------|------|
| 1 | Preparation | 0.5 |
| 2 | Extract Helper Classes | 3.5 |
| 3 | Create Service Layer | 2.0 |
| 4 | Create Recording Layer | 2.5 |
| 5 | Update Tools | 1.0 |
| 6 | Enhanced Error Handling | 1.0 |
| 7 | Configuration Management | 1.0 |
| 8 | Documentation | 2.0 |
| 9 | Testing | 2.0 |
| 10 | Cleanup & Polish | 1.0 |
| 11 | Merge & Deploy | 1.0 |
| **TOTAL** | | **17.5 days** |

**With buffer**: ~20 working days (4 weeks)

---

## Success Criteria

### Must Have ✅
- [ ] All 3,602 PDBs pass validation (100%)
- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Code coverage ≥ 90%
- [ ] Zero linter warnings
- [ ] Documentation complete
- [ ] `BaseFrameCalculator` ≤ 400 LOC
- [ ] `generate_modern_json.cpp` ≤ 300 LOC
- [ ] Build time ≤ 110% of baseline

### Nice to Have ⭐
- [ ] Performance improvement (faster than baseline)
- [ ] Coverage > 95%
- [ ] Doxygen website generated
- [ ] Tutorial video/walkthrough

---

## Daily Progress Tracking

Use this section to track daily progress:

```
Day 1 (YYYY-MM-DD):
- [ ] Task 1
- [ ] Task 2
Progress: X% complete

Day 2 (YYYY-MM-DD):
- [ ] Task 1
- [ ] Task 2
Progress: X% complete

...
```

---

## Notes & Blockers

Use this section to note any issues, blockers, or important decisions:

```
YYYY-MM-DD: Note about X
YYYY-MM-DD: Blocker: Y (resolved by Z)
```

---

## Final Checklist

Before declaring Phase 11 complete:

- [ ] ✅ 3,602/3,602 PDBs pass validation
- [ ] ✅ All tests pass
- [ ] ✅ Coverage ≥ 90%
- [ ] ✅ Zero warnings
- [ ] ✅ Documentation complete
- [ ] ✅ Code formatted
- [ ] ✅ PR approved and merged
- [ ] ✅ Release tagged: `v2.0.0-stage2-refactored`
- [ ] ✅ Announced to team

🎉 **REFACTORING COMPLETE!** 🎉

