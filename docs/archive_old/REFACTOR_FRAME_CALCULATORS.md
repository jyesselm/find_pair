# Refactoring Frame Calculator Classes

## Current Problems

### 1. Mixed Responsibilities

**BaseFrameCalculator** currently does too much:
- ✅ Calculates frames (core algorithm) - **GOOD**
- ❌ Records JSON (base_frame_calc, ls_fitting, frame_calc) - **SHOULD BE SEPARATE**
- ❌ Has multiple "record" methods with duplicated logic:
  - `calculate_and_record_frames()` - records all 3 types
  - `calculate_and_record_frames_only()` - records only base_frame_calc + frame_calc
- ❌ Contains iteration logic (looping through residues, checking legacy_residue_idx)

**LsFittingCalculator** is a confusing wrapper:
- ❌ Just wraps `BaseFrameCalculator` and duplicates iteration logic
- ❌ Delegates all methods to `BaseFrameCalculator` (proxy pattern gone wrong)
- ❌ No real separation of concerns - just filters what gets recorded

### 2. Code Duplication

The iteration logic is duplicated in 3 places:
1. `BaseFrameCalculator::calculate_and_record_frames()` (lines 529-591)
2. `BaseFrameCalculator::calculate_and_record_frames_only()` (lines 593-649)
3. `LsFittingCalculator::calculate_and_record()` (lines 20-68)

All three:
- Auto-detect RNA
- Iterate through `structure.residues_in_legacy_order()`
- Skip amino acids
- Call `calculate_frame()`
- Check `legacy_residue_idx`
- Record different JSON types

### 3. Unclear Separation

- What's the difference between `BaseFrameCalculator` and `LsFittingCalculator`?
- Why does `LsFittingCalculator` exist if it just wraps `BaseFrameCalculator`?
- Why does `BaseFrameCalculator` have recording methods at all?

## Proposed Architecture

### Principle: Separation of Concerns

**Calculators** = Pure algorithms (calculate frames, no JSON)
**Recorders** = Write JSON (use calculators, don't calculate)

### New Structure

```
BaseFrameCalculator (algorithm only)
  ├── calculate_frame(Residue&) -> FrameCalculationResult
  ├── calculate_all_frames(Structure&) -> void
  └── NO JSON recording methods

FrameJsonRecorder (composition pattern)
  ├── Uses BaseFrameCalculator internally
  ├── record_base_frame_calc(Structure&, JsonWriter&)
  ├── record_ls_fitting(Structure&, JsonWriter&)
  ├── record_frame_calc(Structure&, JsonWriter&)
  └── record_all(Structure&, JsonWriter&) // convenience method
```

### Detailed Design

#### 1. BaseFrameCalculator (Algorithm Only)

**Responsibilities:**
- Calculate frames for residues
- Match ring atoms
- Perform least-squares fitting
- Return `FrameCalculationResult`

**NO JSON recording!**

```cpp
class BaseFrameCalculator {
public:
    // Core calculation methods
    FrameCalculationResult calculate_frame(Residue& residue);
    FrameCalculationResult calculate_frame_const(const Residue& residue) const;
    void calculate_all_frames(Structure& structure);
    
    // Configuration
    void set_template_path(const std::filesystem::path& path);
    void set_is_rna(bool is_rna);
    void set_legacy_mode(bool legacy_mode);
    static bool detect_rna(const Structure& structure);
    
    // NO recording methods!
};
```

#### 2. FrameJsonRecorder (Recording Only)

**Responsibilities:**
- Iterate through residues in legacy order
- Use `BaseFrameCalculator` to calculate frames
- Record JSON via `JsonWriter`
- Handle different recording scenarios

```cpp
class FrameJsonRecorder {
public:
    explicit FrameJsonRecorder(BaseFrameCalculator& calculator);
    
    // Individual record methods
    size_t record_base_frame_calc(Structure& structure, JsonWriter& writer);
    size_t record_ls_fitting(Structure& structure, JsonWriter& writer);
    size_t record_frame_calc(Structure& structure, JsonWriter& writer);
    
    // Convenience method for all three
    size_t record_all(Structure& structure, JsonWriter& writer);
    
private:
    BaseFrameCalculator& calculator_;
    
    // Shared iteration logic
    template<typename RecordFunc>
    size_t iterate_and_record(Structure& structure, JsonWriter& writer, RecordFunc record_func);
};
```

#### 3. Usage in generate_modern_json.cpp

**Stage 3: ls_fitting**
```cpp
BaseFrameCalculator calculator("data/templates");
calculator.set_legacy_mode(legacy_mode);
calculator.set_is_rna(LsFittingCalculator::detect_rna(structure));

FrameJsonRecorder recorder(calculator);
size_t count = recorder.record_ls_fitting(structure, writer);
```

**Stage 4: frames**
```cpp
BaseFrameCalculator calculator("data/templates");
calculator.set_legacy_mode(legacy_mode);
calculator.set_is_rna(BaseFrameCalculator::detect_rna(structure));

FrameJsonRecorder recorder(calculator);
size_t count = recorder.record_base_frame_calc(structure, writer);
count += recorder.record_frame_calc(structure, writer);
```

**Stage "all":**
```cpp
BaseFrameCalculator calculator("data/templates");
calculator.set_legacy_mode(legacy_mode);
calculator.set_is_rna(BaseFrameCalculator::detect_rna(structure));

FrameJsonRecorder recorder(calculator);
size_t count = recorder.record_all(structure, writer);
```

## Migration Plan

### Phase 1: Create FrameJsonRecorder
1. Create `include/x3dna/io/frame_json_recorder.hpp`
2. Create `src/x3dna/io/frame_json_recorder.cpp`
3. Move iteration logic from `BaseFrameCalculator` to `FrameJsonRecorder`
4. Implement all record methods

### Phase 2: Remove Recording from BaseFrameCalculator
1. Remove `calculate_and_record_frames()` from `BaseFrameCalculator`
2. Remove `calculate_and_record_frames_only()` from `BaseFrameCalculator`
3. Keep only calculation methods

### Phase 3: Remove LsFittingCalculator
1. Delete `include/x3dna/algorithms/ls_fitting_calculator.hpp`
2. Delete `src/x3dna/algorithms/ls_fitting_calculator.cpp`
3. Remove from `CMakeLists.txt`
4. Update `generate_modern_json.cpp` to use `FrameJsonRecorder`

### Phase 4: Update generate_modern_json.cpp
1. Replace `LsFittingCalculator` with `BaseFrameCalculator` + `FrameJsonRecorder`
2. Update Stage 3 to use `recorder.record_ls_fitting()`
3. Update Stage 4 to use `recorder.record_base_frame_calc()` + `recorder.record_frame_calc()`

## Benefits

1. **Clear Separation**: Calculators calculate, recorders record
2. **No Duplication**: Iteration logic in one place (`FrameJsonRecorder`)
3. **Flexible**: Can record any combination of JSON types
4. **Testable**: Can test calculation and recording separately
5. **Maintainable**: Changes to recording don't affect calculation logic

## File Structure After Refactoring

```
include/x3dna/algorithms/
  ├── base_frame_calculator.hpp  (algorithm only, ~150 lines)
  
include/x3dna/io/
  ├── frame_json_recorder.hpp    (new, ~80 lines)
  
src/x3dna/algorithms/
  ├── base_frame_calculator.cpp  (algorithm only, ~500 lines, removed ~120 lines)
  
src/x3dna/io/
  ├── frame_json_recorder.cpp    (new, ~150 lines)
  
tools/
  ├── generate_modern_json.cpp   (simplified, uses recorder)
```

## Example: FrameJsonRecorder Implementation

```cpp
// frame_json_recorder.hpp
class FrameJsonRecorder {
public:
    explicit FrameJsonRecorder(BaseFrameCalculator& calculator);
    
    size_t record_base_frame_calc(Structure& structure, JsonWriter& writer);
    size_t record_ls_fitting(Structure& structure, JsonWriter& writer);
    size_t record_frame_calc(Structure& structure, JsonWriter& writer);
    size_t record_all(Structure& structure, JsonWriter& writer);
    
private:
    BaseFrameCalculator& calculator_;
    
    // Helper: iterate residues and call record_func for each valid frame
    template<typename RecordFunc>
    size_t iterate_and_record(Structure& structure, JsonWriter& writer, RecordFunc record_func);
};

// frame_json_recorder.cpp
template<typename RecordFunc>
size_t FrameJsonRecorder::iterate_and_record(Structure& structure, JsonWriter& writer, RecordFunc record_func) {
    auto residues = structure.residues_in_legacy_order();
    size_t count = 0;
    
    for (const auto* residue_ptr : residues) {
        auto* residue = const_cast<Residue*>(residue_ptr);
        
        if (residue->residue_type() == ResidueType::AMINO_ACID) {
            continue;
        }
        
        FrameCalculationResult result = calculator_.calculate_frame(*residue);
        if (!result.is_valid) {
            continue;
        }
        
        int legacy_residue_idx = 0;
        if (!residue->atoms().empty()) {
            legacy_residue_idx = residue->atoms()[0].legacy_residue_idx();
        }
        
        if (legacy_residue_idx <= 0) {
            continue;
        }
        
        size_t record_idx = static_cast<size_t>(legacy_residue_idx);
        record_func(record_idx, *residue, result, writer);
        count++;
    }
    
    return count;
}

size_t FrameJsonRecorder::record_ls_fitting(Structure& structure, JsonWriter& writer) {
    return iterate_and_record(structure, writer, [](size_t idx, const Residue& res, 
                                                     const FrameCalculationResult& result, 
                                                     JsonWriter& w) {
        w.record_ls_fitting(idx, result.num_matched, result.rms_fit,
                           result.rotation_matrix, result.translation,
                           res.name(), res.chain_id(), res.seq_num(), res.insertion());
    });
}
```

## Summary

**Current State:**
- ❌ `BaseFrameCalculator` mixes calculation + recording
- ❌ `LsFittingCalculator` is a confusing wrapper
- ❌ Code duplication in 3 places
- ❌ Unclear responsibilities

**Target State:**
- ✅ `BaseFrameCalculator` = pure algorithm (calculate only)
- ✅ `FrameJsonRecorder` = recording only (uses calculator)
- ✅ No duplication (iteration logic in one place)
- ✅ Clear separation of concerns
- ✅ Flexible and maintainable

---

# Extended Refactoring: All Algorithm Classes

## Additional Problems Found

### 1. BasePairFinder - Mixed Responsibilities

**Current Issues:**
- Has `find_pairs_with_recording()` method that mixes finding pairs with JSON recording
- Takes `JsonWriter*` as optional parameter (nullptr to skip recording)
- Recording logic embedded in algorithm code
- Makes algorithm depend on `JsonWriter` (tight coupling)

**Location:** `src/x3dna/algorithms/base_pair_finder.cpp` (lines 32-1060+)

**Problem Pattern:**
```cpp
std::vector<BasePair> find_pairs_with_recording(Structure& structure, JsonWriter* writer) const {
    // ... algorithm logic ...
    if (writer) {
        writer->record_pair_validation(...);
        writer->record_distance_checks(...);
        writer->record_hbond_list(...);
        // ... more recording ...
    }
    // ... more algorithm logic ...
}
```

### 2. generate_modern_json.cpp - Massive Code Duplication

**Current Issues:**
- **793 lines** - way too long (should be ~200-300 lines)
- RNA detection duplicated in **3+ places**:
  - Lines 286-297 (manual loop)
  - Lines 203-209 (using LsFittingCalculator::detect_rna)
  - Lines 236-241 (using BaseFrameCalculator::detect_rna)
- Frame calculation/recording logic duplicated in **4+ places**:
  - Lines 373-531 (old "all" stage logic)
  - Lines 189-218 (Stage 3: ls_fitting)
  - Lines 220-251 (Stage 4: frames)
  - Lines 477-515 (duplicate frame recording)
- Duplicate code for building `paired_legacy_indices` (lines 544-572 and 574-602)
- Too many responsibilities:
  - Parsing command-line arguments
  - Loading PDB files
  - Calculating frames
  - Recording JSON
  - Finding pairs
  - Calculating parameters
  - Writing files
  - Index validation

### 3. ParameterCalculator - Good Example! ✅

**Status:** Clean separation
- Only calculates parameters
- No JSON recording
- Pure algorithm class
- **This is the model to follow!**

### 4. JsonWriter - Too Many Responsibilities?

**Current State:**
- Has 17+ `record_*` methods
- Each method knows how to format specific JSON types
- This is actually OK - it's a writer, it should know how to write

**Potential Issue:**
- Algorithms shouldn't call `JsonWriter` directly
- Should use recorder classes as intermediaries

## Comprehensive Refactoring Plan

### Principle: Clear Separation of Concerns

```
┌─────────────────────────────────────────────────────────┐
│ Algorithms (Pure Calculation)                         │
│  ├── BaseFrameCalculator                                │
│  ├── BasePairFinder                                     │
│  ├── ParameterCalculator ✅ (already clean)            │
│  └── HydrogenBondFinder                                 │
└─────────────────────────────────────────────────────────┘
                        │
                        │ uses
                        ▼
┌─────────────────────────────────────────────────────────┐
│ Recorders (JSON Writing)                                │
│  ├── FrameJsonRecorder                                  │
│  ├── PairJsonRecorder                                   │
│  ├── ParameterJsonRecorder                               │
│  └── HydrogenBondJsonRecorder                           │
└─────────────────────────────────────────────────────────┘
                        │
                        │ uses
                        ▼
┌─────────────────────────────────────────────────────────┐
│ JsonWriter (Low-level JSON Writing)                    │
│  └── record_*(...) methods                             │
└─────────────────────────────────────────────────────────┘
```

### Detailed Refactoring

#### 1. BasePairFinder - Remove Recording

**Current:**
```cpp
class BasePairFinder {
    std::vector<BasePair> find_pairs_with_recording(Structure&, JsonWriter*);
};
```

**Target:**
```cpp
class BasePairFinder {
    // Pure algorithm - no JsonWriter dependency
    std::vector<BasePair> find_pairs(Structure&);
    
    // Return validation results for recording
    struct PairFindingResult {
        std::vector<BasePair> pairs;
        std::vector<ValidationResult> validations;
        std::vector<DistanceCheck> distance_checks;
        std::vector<HydrogenBond> hbonds;
        // ... other data needed for JSON ...
    };
    
    PairFindingResult find_pairs_with_details(Structure&);
};
```

**New Recorder:**
```cpp
class PairJsonRecorder {
public:
    explicit PairJsonRecorder(BasePairFinder& finder);
    
    // Record all pair-related JSON
    size_t record_all(Structure& structure, JsonWriter& writer);
    
    // Individual record methods
    size_t record_pair_validation(Structure& structure, JsonWriter& writer);
    size_t record_distance_checks(Structure& structure, JsonWriter& writer);
    size_t record_hbond_list(Structure& structure, JsonWriter& writer);
    size_t record_base_pairs(const std::vector<BasePair>& pairs, JsonWriter& writer);
    
private:
    BasePairFinder& finder_;
};
```

#### 2. ParameterCalculator - Add Recorder

**Current:** ✅ Already clean (only calculates)

**Add:**
```cpp
class ParameterJsonRecorder {
public:
    explicit ParameterJsonRecorder(ParameterCalculator& calculator);
    
    size_t record_bpstep_params(const std::vector<BasePair>& pairs, JsonWriter& writer);
    size_t record_helical_params(const std::vector<BasePair>& pairs, JsonWriter& writer);
    
private:
    ParameterCalculator& calculator_;
};
```

#### 3. generate_modern_json.cpp - Massive Simplification

**Current:** 793 lines with massive duplication

**Target:** ~200-300 lines, clean and simple

**Strategy:**
1. Extract RNA detection to utility function
2. Extract frame setup to helper function
3. Use recorder classes instead of direct algorithm calls
4. Remove all duplicate code

**Before (Stage 3 - 30 lines):**
```cpp
// Stage 3: LS Fitting
if (stage == "ls_fitting" || stage == "all") {
    std::cout << "Stage 3: Writing ls_fitting...\n";
    JsonWriter writer(pdb_file);
    writer.record_residue_indices(structure);
    
    LsFittingCalculator calculator("data/templates");
    calculator.set_legacy_mode(legacy_mode);
    
    bool is_rna = LsFittingCalculator::detect_rna(structure);
    if (is_rna) {
        std::cout << "Detected RNA structure (O2' atoms found)\n";
    } else {
        std::cout << "Detected DNA structure (no O2' atoms)\n";
    }
    
    size_t records_count = calculator.calculate_and_record(structure, writer);
    writer.write_split_files(json_output_dir, true);
    std::cout << "  ✅ ls_fitting/" << pdb_name << ".json\n";
    std::cout << "     " << records_count << " ls_fitting records calculated\n\n";
}
```

**After (Stage 3 - 10 lines):**
```cpp
// Stage 3: LS Fitting
if (stage == "ls_fitting" || stage == "all") {
    auto [calculator, writer] = setup_frame_calculation(pdb_file, legacy_mode, structure);
    FrameJsonRecorder recorder(calculator);
    size_t count = recorder.record_ls_fitting(structure, writer);
    writer.write_split_files(json_output_dir, true);
    print_stage_result("ls_fitting", pdb_name, count);
}
```

**Helper Functions:**
```cpp
// Utility: RNA detection (used everywhere)
bool detect_rna_structure(const Structure& structure);

// Helper: Setup frame calculation (used in stages 3, 4, and "all")
std::pair<BaseFrameCalculator, JsonWriter> 
setup_frame_calculation(const std::filesystem::path& pdb_file, 
                       bool legacy_mode, 
                       Structure& structure);

// Helper: Print stage result
void print_stage_result(const std::string& stage_name, 
                       const std::string& pdb_name, 
                       size_t count);
```

#### 4. Remove Duplicate Code Patterns

**Pattern 1: RNA Detection**
- **Current:** Duplicated in 3+ places
- **Fix:** Single utility function `detect_rna_structure()`

**Pattern 2: Frame Setup**
- **Current:** Duplicated setup code (create calculator, set legacy mode, detect RNA)
- **Fix:** `setup_frame_calculation()` helper function

**Pattern 3: Iteration Logic**
- **Current:** Duplicated in BaseFrameCalculator methods and LsFittingCalculator
- **Fix:** Move to `FrameJsonRecorder::iterate_and_record()`

**Pattern 4: paired_legacy_indices Building**
- **Current:** Duplicated in lines 544-572 and 574-602
- **Fix:** Extract to helper function `build_paired_indices()`

## Complete Refactoring Checklist

### Phase 1: Frame Calculators (Already Documented)
- [ ] Create `FrameJsonRecorder`
- [ ] Remove recording from `BaseFrameCalculator`
- [ ] Remove `LsFittingCalculator`
- [ ] Update `generate_modern_json.cpp` stages 3 & 4

### Phase 2: Pair Finding
- [ ] Refactor `BasePairFinder::find_pairs_with_recording()`
  - Remove `JsonWriter*` parameter
  - Return `PairFindingResult` with all data
- [ ] Create `PairJsonRecorder` class
- [ ] Update `generate_modern_json.cpp` to use recorder

### Phase 3: Parameters
- [ ] Create `ParameterJsonRecorder` class
- [ ] Update `generate_modern_json.cpp` to use recorder

### Phase 4: generate_modern_json.cpp Cleanup
- [ ] Extract RNA detection to utility function
- [ ] Extract frame setup to helper function
- [ ] Extract `paired_legacy_indices` building to helper
- [ ] Remove all duplicate code
- [ ] Simplify stage implementations
- [ ] Target: Reduce from 793 lines to ~200-300 lines

### Phase 5: Hydrogen Bonds (Future)
- [ ] Review `HydrogenBondFinder` for similar issues
- [ ] Create `HydrogenBondJsonRecorder` if needed

## File Structure After Complete Refactoring

```
include/x3dna/algorithms/
  ├── base_frame_calculator.hpp      (algorithm only, ~150 lines)
  ├── base_pair_finder.hpp            (algorithm only, remove recording)
  ├── parameter_calculator.hpp        (already clean ✅)
  
include/x3dna/io/
  ├── frame_json_recorder.hpp        (new, ~80 lines)
  ├── pair_json_recorder.hpp          (new, ~100 lines)
  ├── parameter_json_recorder.hpp    (new, ~60 lines)
  ├── json_writer.hpp                 (unchanged, low-level writer)
  
src/x3dna/algorithms/
  ├── base_frame_calculator.cpp      (algorithm only, ~500 lines)
  ├── base_pair_finder.cpp            (algorithm only, remove ~200 lines)
  
src/x3dna/io/
  ├── frame_json_recorder.cpp        (new, ~150 lines)
  ├── pair_json_recorder.cpp          (new, ~200 lines)
  ├── parameter_json_recorder.cpp    (new, ~100 lines)
  
tools/
  ├── generate_modern_json.cpp        (simplified, ~200-300 lines, down from 793)
  └── generate_modern_json_utils.cpp  (new, helper functions, ~100 lines)
```

## Benefits of Complete Refactoring

1. **Clear Architecture:**
   - Algorithms = pure calculation
   - Recorders = JSON writing
   - Tools = orchestration

2. **No Duplication:**
   - RNA detection: 1 function
   - Frame setup: 1 function
   - Iteration logic: in recorders
   - paired_legacy_indices: 1 function

3. **Maintainable:**
   - `generate_modern_json.cpp` reduced from 793 to ~200-300 lines
   - Each class has single responsibility
   - Easy to test algorithms independently

4. **Flexible:**
   - Can use algorithms without JSON
   - Can record JSON in different formats
   - Can combine algorithms in different ways

5. **Testable:**
   - Test algorithms without JSON
   - Test recorders with mock algorithms
   - Test tools with mock recorders

## Summary: Complete Refactoring

**Current State:**
- ❌ `BaseFrameCalculator` mixes calculation + recording
- ❌ `LsFittingCalculator` is confusing wrapper
- ❌ `BasePairFinder` has `find_pairs_with_recording()`
- ❌ `generate_modern_json.cpp` is 793 lines with massive duplication
- ❌ RNA detection duplicated 3+ times
- ❌ Frame setup duplicated 4+ times
- ❌ `paired_legacy_indices` building duplicated 2+ times

**Target State:**
- ✅ `BaseFrameCalculator` = pure algorithm
- ✅ `BasePairFinder` = pure algorithm (no JsonWriter dependency)
- ✅ `ParameterCalculator` = pure algorithm (already ✅)
- ✅ `FrameJsonRecorder` = records frame JSON
- ✅ `PairJsonRecorder` = records pair JSON
- ✅ `ParameterJsonRecorder` = records parameter JSON
- ✅ `generate_modern_json.cpp` = ~200-300 lines, clean orchestration
- ✅ No duplication (utilities for common patterns)
- ✅ Clear separation of concerns throughout

