# Legacy Index Strategy

**Date**: December 2, 2025  
**Principle**: All objects store both modern (0-based) and legacy (1-based) indices

---

## Why We Need Both Index Systems

### Modern Indices (0-based)
- **Purpose**: Internal C++ array/vector access
- **Range**: [0, n-1]
- **Usage**: `vector[modern_idx]`, loops, algorithms
- **Standard**: C++ convention

### Legacy Indices (1-based)
- **Purpose**: JSON output, comparison with legacy
- **Range**: [1, n]
- **Usage**: JSON serialization, validation
- **Standard**: Legacy code convention (Fortran-style)

---

## Index Storage Strategy

### Every Object Stores Both

#### Atom
```cpp
class Atom {
private:
    size_t atom_idx_;         // 0-based, modern
    int legacy_atom_idx_;     // 1-based, for JSON
};

// JSON output:
j["atom_idx"] = legacy_atom_idx_;  // Use 1-based
```

#### Residue
```cpp
class Residue {
private:
    size_t residue_idx_;      // 0-based, modern
    int legacy_residue_idx_;  // 1-based, for JSON
};

// JSON output:
j["residue_idx"] = legacy_residue_idx_;  // Use 1-based
```

#### BasePair
```cpp
class BasePair {
private:
    // Modern indices (0-based)
    size_t residue_idx1_;
    size_t residue_idx2_;
    size_t basepair_idx_;     // Modern basepair index
    
    // Legacy indices (1-based, for JSON)
    int legacy_basepair_idx_;    // Which base pair # (1, 2, 3...)
    int legacy_residue_idx1_;    // Residue 1's legacy index
    int legacy_residue_idx2_;    // Residue 2's legacy index
};

// JSON output:
j["basepair_idx"] = legacy_basepair_idx_;  // Use 1-based
j["base_i"] = legacy_residue_idx1_;        // Use 1-based
j["base_j"] = legacy_residue_idx2_;        // Use 1-based
```

---

## Assignment Strategy

### During PDB Parsing

```cpp
// PdbParser::parse_file()
int atom_counter = 1;        // 1-based for legacy
int residue_counter = 1;     // 1-based for legacy

for each atom in PDB file order:
    atom.set_atom_idx(modern_idx);           // 0-based
    atom.set_legacy_atom_idx(atom_counter++);  // 1-based
    
for each residue in PDB file order:
    res.set_residue_idx(modern_idx);              // 0-based
    res.set_legacy_residue_idx(residue_counter++);  // 1-based
```

**Key**: Process in PDB file order → legacy indices match legacy code

### During Base Pair Creation

```cpp
// BasePairFinder::find_pairs()
int basepair_counter = 1;  // 1-based for legacy

for each selected pair:
    BasePair pair(idx1, idx2);  // Modern indices 0-based
    
    // Get residues' legacy indices
    int leg_idx1 = structure.residue(idx1).legacy_residue_idx();
    int leg_idx2 = structure.residue(idx2).legacy_residue_idx();
    
    // Assign legacy indices to pair
    pair.set_legacy_basepair_idx(basepair_counter++);  // 1-based
    pair.set_legacy_residue_indices(leg_idx1, leg_idx2);  // Already 1-based
    
    pairs.push_back(pair);
```

---

## JSON Serialization (No Conversion)

### Atom JSON
```cpp
nlohmann::json Atom::to_json() const {
    nlohmann::json j;
    j["atom_idx"] = legacy_atom_idx_;  // Already 1-based, no +1
    j["atom_name"] = name_;
    // ...
    return j;
}
```

### Residue Frame JSON
```cpp
nlohmann::json Residue::frame_to_json() const {
    nlohmann::json j;
    j["residue_idx"] = legacy_residue_idx_;  // Already 1-based, no +1
    j["base_type"] = base_type_;
    // ...
    return j;
}
```

### BasePair JSON
```cpp
nlohmann::json BasePair::to_json() const {
    nlohmann::json j;
    j["basepair_idx"] = legacy_basepair_idx_;  // Already 1-based, no +1
    j["base_i"] = legacy_residue_idx1_;        // Already 1-based, no +1
    j["base_j"] = legacy_residue_idx2_;        // Already 1-based, no +1
    // ...
    return j;
}
```

**Key Rule**: **NO arithmetic during serialization**. Just write the stored value.

---

## Why This Matters

### Correct (Current Approach)
```cpp
// Store legacy index during creation
atom.set_legacy_atom_idx(counter);  // 1-based

// Write to JSON (no conversion)
j["atom_idx"] = atom.legacy_atom_idx();  // Just write it
```

### Wrong (Don't Do This)
```cpp
// Store modern index only
atom.set_atom_idx(i);  // 0-based

// Convert during JSON writing
j["atom_idx"] = atom.atom_idx() + 1;  // Math during serialization ❌
```

**Why wrong**:
- Conversion scattered across code
- Error-prone (forget +1 somewhere)
- Unclear which is which
- Hard to debug

**Why correct**:
- Clear separation: modern for C++, legacy for JSON
- No conversion/math
- Explicit what's being used
- Easy to debug

---

## Validation Strategy

### During Generation
```cpp
// Just generate with legacy indices
atom.set_legacy_atom_idx(counter++);
// Write to JSON
j["atom_idx"] = atom.legacy_atom_idx();
```

### During Comparison
```python
# Compare the written values
assert legacy_json["atom_idx"] == modern_json["atom_idx"]
# Both should be 1-based, no conversion
```

---

## Implementation Checklist

### Verify All Classes Store Legacy Indices

- [x] **Atom**: Has `legacy_atom_idx_` ✅
- [x] **Residue**: Has `legacy_residue_idx_` ✅
- [ ] **BasePair**: Add `legacy_basepair_idx_`, `legacy_residue_idx1_`, `legacy_residue_idx2_`
- [ ] **Chain**: Add `legacy_chain_idx_` (if needed)

### Verify Assignment

- [ ] PdbParser assigns legacy indices in PDB file order
- [ ] BasePairFinder assigns legacy basepair indices
- [ ] Counters start at 1 (not 0)

### Verify JSON Writing

- [ ] All `to_json()` methods use legacy indices
- [ ] No `+ 1` or `- 1` during serialization
- [ ] Direct write: `j["idx"] = legacy_idx_`

---

## Example: Complete BasePair

### Class Definition
```cpp
// include/x3dna/models/base_pair.hpp
class BasePair {
public:
    // Getters for modern (0-based)
    size_t residue_idx1() const { return residue_idx1_; }
    size_t residue_idx2() const { return residue_idx2_; }
    size_t basepair_idx() const { return basepair_idx_; }
    
    // Getters for legacy (1-based)
    int legacy_basepair_idx() const { return legacy_basepair_idx_; }
    int legacy_residue_idx1() const { return legacy_residue_idx1_; }
    int legacy_residue_idx2() const { return legacy_residue_idx2_; }
    
    // Setters
    void set_legacy_basepair_idx(int idx) { legacy_basepair_idx_ = idx; }
    void set_legacy_residue_indices(int idx1, int idx2) {
        legacy_residue_idx1_ = idx1;
        legacy_residue_idx2_ = idx2;
    }
    
    // JSON serialization
    nlohmann::json to_json() const;
    
private:
    // Modern indices (0-based, for C++ arrays)
    size_t basepair_idx_;
    size_t residue_idx1_;
    size_t residue_idx2_;
    
    // Legacy indices (1-based, for JSON)
    int legacy_basepair_idx_;
    int legacy_residue_idx1_;
    int legacy_residue_idx2_;
    
    // Other data...
};
```

### JSON Serialization
```cpp
// src/x3dna/models/base_pair.cpp
nlohmann::json BasePair::to_json() const {
    nlohmann::json j;
    
    // All indices are legacy (1-based), no conversion
    j["basepair_idx"] = legacy_basepair_idx_;
    j["base_i"] = legacy_residue_idx1_;
    j["base_j"] = legacy_residue_idx2_;
    j["bp_type"] = bp_type_;
    
    // Frames
    j["orien_i"] = frame1_.rotation_matrix_to_json();
    j["orien_j"] = frame2_.rotation_matrix_to_json();
    j["org_i"] = frame1_.origin_to_json();
    j["org_j"] = frame2_.origin_to_json();
    
    return j;
}
```

### Creation in BasePairFinder
```cpp
// src/x3dna/algorithms/base_pair_finder.cpp
std::vector<BasePair> BasePairFinder::find_pairs(Structure& structure) {
    std::vector<BasePair> pairs;
    int basepair_counter = 1;  // 1-based counter for legacy
    
    // ... greedy selection algorithm ...
    
    for (auto [idx1, idx2] : selected_pairs) {
        BasePair pair(idx1, idx2);  // Modern indices
        
        // Assign legacy indices
        pair.set_legacy_basepair_idx(basepair_counter++);
        
        auto& res1 = structure.residue(idx1);
        auto& res2 = structure.residue(idx2);
        pair.set_legacy_residue_indices(
            res1.legacy_residue_idx(),  // Already 1-based
            res2.legacy_residue_idx()   // Already 1-based
        );
        
        pairs.push_back(pair);
    }
    
    return pairs;
}
```

---

## Summary: Universal Legacy Index Pattern

**Every indexed object needs**:
1. Modern index (0-based) - for C++ internal use
2. Legacy index (1-based) - for JSON output

**Objects that need this**:
- ✅ Atom: `atom_idx_` (0-based) + `legacy_atom_idx_` (1-based)
- ✅ Residue: `residue_idx_` (0-based) + `legacy_residue_idx_` (1-based)
- ⏳ BasePair: Need to add `legacy_basepair_idx_`, `legacy_residue_idx1_`, `legacy_residue_idx2_`
- ⏳ Chain: `chain_idx_` (0-based) + `legacy_chain_idx_` (1-based) if tracked

**JSON serialization rule**: Always use legacy indices, never convert

