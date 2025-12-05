# Residue Architecture Refactoring Plan

## Current Problems

1. ❌ Circular dependencies between `residue.hpp` and `modified_nucleotide_registry.hpp`
2. ❌ `ResidueType` enum defined in multiple places causing redefinition errors
3. ❌ `one_letter_code()` is a method with giant if-statements instead of stored value
4. ❌ No clear separation of concerns
5. ❌ Modified nucleotide info scattered across multiple files
6. ❌ Residue computes properties on-demand instead of storing them

## Target Architecture

### 1. **ResidueType Enum** (Stand-alone)
**File**: `include/x3dna/core/residue_type.hpp`
- Pure enum definition
- No dependencies
- Can be included anywhere

### 2. **ModifiedNucleotideRegistry** (Data Provider)
**File**: `include/x3dna/core/modified_nucleotide_registry.hpp`
- Loads from `resources/config/modified_nucleotides.json`
- Provides lookup methods
- Dependencies: ResidueType enum only

### 3. **ResidueFactory** (Creator Pattern)
**File**: `include/x3dna/core/residue_factory.hpp`
- Creates Residue objects with all properties initialized
- Uses ModifiedNucleotideRegistry for lookups
- Centralizes all residue creation logic

### 4. **Residue Class** (Data Holder)
**File**: `include/x3dna/core/residue.hpp`
- **Stores** properties (not computes):
  - `std::string three_letter_code_` (the name)
  - `char one_letter_code_` (stored value)
  - `ResidueType type_`
  - `bool is_purine_`
  - `bool is_pyrimidine_`
- Simple getters (no complex logic)
- Atoms and other data

### 5. **RingAtomMatcher** (Uses Residue Properties)
- Reads residue properties
- Determines matching strategy based on stored flags

## Refactoring Steps

### Step 1: Extract ResidueType Enum ✅
```cpp
// include/x3dna/core/residue_type.hpp
#pragma once
namespace x3dna { namespace core {
enum class ResidueType {
    UNKNOWN = -2,
    ADENINE = 1,
    CYTOSINE = 2,
    // ... etc
};
}} // namespace x3dna::core
```

### Step 2: Fix ModifiedNucleotideRegistry ✅
- Include only `residue_type.hpp`
- Load from JSON
- No circular dependencies

### Step 3: Create ResidueFactory 🔨
```cpp
class ResidueFactory {
public:
    static Residue create_from_pdb_data(
        const std::string& name,
        int seq_num,
        char chain_id,
        char insertion_code,
        const std::vector<Atom>& atoms
    );
    
private:
    static char determine_one_letter_code(const std::string& name);
    static ResidueType determine_type(const std::string& name, char one_letter);
    static bool is_purine(const std::string& name, ResidueType type);
};
```

### Step 4: Update Residue Class 🔨
**Remove:**
- `one_letter_code()` method with if-statements
- `residue_type()` method with if-statements
- `is_nucleotide()` complex logic

**Add:**
- `char one_letter_code_` member
- `ResidueType type_` member
- `bool is_purine_` member
- Simple getters

**Constructor** - takes all properties:
```cpp
Residue(std::string name, char one_letter_code, ResidueType type, 
        bool is_purine, int seq_num, char chain_id, ...)
```

### Step 5: Update PDB Parser 🔨
- Use `ResidueFactory::create_from_pdb_data()` instead of direct Residue construction

### Step 6: Update All Residue Users 🔨
- Change from `residue.one_letter_code()` computation to `residue.one_letter_code()` getter
- No logic changes, just using stored values

### Step 7: Clean Up Template Assignment 🔨
- Simplify to use ResidueFactory
- Remove old MODIFIED_PURINES/PYRIMIDINES maps

### Step 8: Test Everything 🧪
```bash
# Rebuild
make clean-all && make release

# Test on known working PDBs
./build/generate_modern_json data/pdb/6QIQ.pdb /tmp/test_refactor/

# Compare with legacy
python3 scripts/compare_json.py compare 6QIQ

# Run stage 2 validation on subset
python3 scripts/validate_stage2.py --test-set 100
```

## Implementation Order

1. ✅ Create `include/x3dna/core/residue_type.hpp`
2. ✅ Update `modified_nucleotide_registry.hpp` to use it
3. 🔨 Create `include/x3dna/core/residue_factory.hpp`
4. 🔨 Create `src/x3dna/core/residue_factory.cpp`
5. 🔨 Update `Residue` class to store properties
6. 🔨 Update `PdbParser` to use `ResidueFactory`
7. 🔨 Update all code using residue properties
8. 🔨 Clean up old code
9. 🧪 Build and test
10. 🧪 Validate on Stage 2 PDBs

## Benefits

1. ✅ No circular dependencies
2. ✅ Single source of truth for modified nucleotides (JSON file)
3. ✅ Clean separation of concerns
4. ✅ Easy to add new modified nucleotides (just update JSON)
5. ✅ Residue is simple data holder
6. ✅ Factory pattern centralizes creation logic
7. ✅ Better performance (compute once, store, not compute every access)

## Files to Create

- `include/x3dna/core/residue_type.hpp` (new)
- `include/x3dna/core/residue_factory.hpp` (new)
- `src/x3dna/core/residue_factory.cpp` (new)

## Files to Modify

- `include/x3dna/core/modified_nucleotide_registry.hpp`
- `include/x3dna/core/residue.hpp`
- `src/x3dna/io/pdb_parser.cpp`
- `src/x3dna/algorithms/template_assignment.cpp`
- `CMakeLists.txt`

## Success Criteria

1. ✅ Build succeeds with no errors
2. ✅ No circular dependencies
3. ✅ All Stage 3 exclusions still pass (23/23)
4. ✅ Stage 2 validation maintains or improves success rate
5. ✅ Code is cleaner and more maintainable

---

**Status**: Ready to implement
**Estimated Time**: 1-2 hours
**Priority**: HIGH - fixes architectural issues

