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

