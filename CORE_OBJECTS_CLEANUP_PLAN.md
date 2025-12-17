# Core Objects Cleanup Plan

**Branch**: `refactor/clean-code-principles`
**Scope**: Core data structures (Atom, Residue, Chain, Structure, BasePair, ReferenceFrame)
**Goal**: Clean container-like behavior, simple constructors, centralized ring atom/modified residue handling

## Completion Status

| Phase | Status | Commits |
|-------|--------|---------|
| Phase 1: Ring Atom Registry | ✅ Complete | `0f99d4f` |
| Phase 2: Consolidate is_purine() | ✅ Complete | `def80f8` |
| Phase 3: Residue Builder | ✅ Complete | `c999294` |
| Phase 4: BasePair Frames Required | ✅ Complete | `669c5b2` |
| Phase 5: Container Iteration | ✅ Complete | `2425e04` |
| Phase 6: JSON Serializers | ✅ Complete | `2425e04` |
| Additional: Named Constants | ✅ Complete | `c262bb8` |
| Additional: Flatten Nesting | ✅ Complete | `b4c63bf`, `2fd9eb2`, `31e3ae4`, `22c9cf7` |

**Validation**: test-set-100 = 99% (1 pre-existing failure in 8RUJ)

---

## Current State Analysis

### Issues Identified

1. **Scattered Ring Atom Logic**
   - `Atom::is_ring_atom()` hardcodes 9 atom names (N1, C2, N3, C4, C5, C6, N7, C8, N9)
   - `RingAtomMatcher::get_ring_atom_names()` has similar logic
   - `Residue::ring_atoms()` delegates to `Atom::is_ring_atom()`
   - No single source of truth for ring atom definitions

2. **Scattered Purine/Pyrimidine Classification**
   - `ModifiedNucleotideRegistry::is_purine()` - for modified nucleotides
   - `RingAtomMatcher::is_purine()` - for residue types
   - `ResidueTypeDetector::is_purine()` - duplicate of above
   - `Residue::is_purine_` - stored value, set by factory
   - `ResidueFactory::determine_is_purine()` - actual logic

3. **Residue Construction Complexity**
   - Two constructors: simple (4 params) and full (8 params)
   - Full constructor has 8 parameters - violates 3-4 max guideline
   - Properties (one_letter_code, type, is_purine) set at construction via factory
   - But from_json methods bypass factory

4. **BasePair Uses `std::optional` for Required Data**
   - `frame1_`, `frame2_` are `std::optional<ReferenceFrame>` but always required
   - Precondition documented: "residues must have frames calculated"
   - Every usage does `.value()` or checks `.has_value()` then `.value()`
   - Creates false impression that frames are optional

5. **Container Behavior Inconsistencies**
   - `Structure` returns `std::vector<const Residue*>` - raw pointers
   - `Chain::nucleotides()` returns by value (copy)
   - `Residue::ring_atoms()` returns by value (copy)
   - No consistent iteration interface

6. **Mutable Non-const Access**
   - `Structure::chains()` returns mutable reference
   - `Chain::residues()` returns mutable reference
   - `Residue::atoms()` returns mutable reference
   - Breaks encapsulation - external code can modify internals

---

## Phase 1: Centralize Ring Atom Definitions

**Goal**: Single source of truth for what constitutes a ring atom.

### Step 1.1: Create `core/ring_atom_registry.hpp`

```cpp
// include/x3dna/core/ring_atom_registry.hpp
namespace x3dna::core {

/**
 * @brief Single source of truth for ring atom definitions
 *
 * Ring atoms are the base ring atoms used for least-squares fitting.
 * Purines have 9 ring atoms, pyrimidines have 6.
 */
class RingAtomRegistry {
public:
    /**
     * @brief Check if an atom name is a ring atom
     * @param atom_name Atom name (will be trimmed)
     * @return true if ring atom
     */
    [[nodiscard]] static bool is_ring_atom(std::string_view atom_name);

    /**
     * @brief Get ring atom names for purine bases (A, G, I)
     * @return 9 atom names: N1, C2, N3, C4, C5, C6, N7, C8, N9
     */
    [[nodiscard]] static const std::vector<std::string>& purine_atoms();

    /**
     * @brief Get ring atom names for pyrimidine bases (C, U, T, P)
     * @return 6 atom names: N1, C2, N3, C4, C5, C6
     */
    [[nodiscard]] static const std::vector<std::string>& pyrimidine_atoms();

    /**
     * @brief Get ring atom names for a residue type
     * @param type ResidueType
     * @return Ring atom names for that type
     */
    [[nodiscard]] static const std::vector<std::string>& atoms_for_type(ResidueType type);

private:
    static const std::vector<std::string> PURINE_RING_ATOMS;
    static const std::vector<std::string> PYRIMIDINE_RING_ATOMS;
};

} // namespace x3dna::core
```

### Step 1.2: Update Atom to Use Registry

```cpp
// In atom.hpp
[[nodiscard]] bool is_ring_atom() const {
    return RingAtomRegistry::is_ring_atom(name_);
}
```

### Step 1.3: Update RingAtomMatcher to Use Registry

```cpp
// In ring_atom_matcher.cpp
std::vector<std::string> RingAtomMatcher::get_ring_atom_names(ResidueType type) {
    return RingAtomRegistry::atoms_for_type(type);
}
```

### Files Changed
- NEW: `include/x3dna/core/ring_atom_registry.hpp`
- NEW: `src/x3dna/core/ring_atom_registry.cpp`
- UPDATE: `include/x3dna/core/atom.hpp` (delegate to registry)
- UPDATE: `src/x3dna/algorithms/ring_atom_matcher.cpp` (delegate to registry)

### Unit Test
- `tests/unit/core/ring_atom_registry_test.cpp`

---

## Phase 2: Consolidate Purine/Pyrimidine Classification

**Goal**: Single source of truth for purine vs pyrimidine classification.

### Step 2.1: Extend `ModifiedNucleotideRegistry`

Add method for standard nucleotides too:

```cpp
// In modified_nucleotide_registry.hpp
/**
 * @brief Check if a residue type is a purine
 * @param type ResidueType enum value
 * @return true for ADENINE, GUANINE, INOSINE; false otherwise
 */
[[nodiscard]] static bool is_purine(ResidueType type);
```

### Step 2.2: Remove Duplicate `is_purine()` Methods

Remove from:
- `RingAtomMatcher::is_purine()` - replace with `ModifiedNucleotideRegistry::is_purine()`
- `ResidueTypeDetector::is_purine()` - replace with `ModifiedNucleotideRegistry::is_purine()`

### Files Changed
- UPDATE: `include/x3dna/core/modified_nucleotide_registry.hpp`
- UPDATE: `src/x3dna/core/modified_nucleotide_registry.cpp`
- UPDATE: `src/x3dna/algorithms/ring_atom_matcher.cpp`
- UPDATE: `src/x3dna/algorithms/residue_type_detector.cpp`

---

## Phase 3: Clean Up Residue Construction

**Goal**: Simple constructors, clear builder pattern.

### Step 3.1: Add Residue::Builder

```cpp
class Residue {
public:
    class Builder;

    /**
     * @brief Create a builder for fluent construction
     * @param name Residue name (required)
     * @param seq_num Sequence number (required)
     * @param chain_id Chain ID (required)
     */
    [[nodiscard]] static Builder create(std::string name, int seq_num, char chain_id);

    // Keep simple constructor for compatibility
    Residue(const std::string& name, int seq_num, char chain_id, char insertion = ' ');

    // Remove the 8-parameter constructor
    // Residue(name, one_letter_code, type, is_purine, seq_num, chain_id, insertion, atoms)
    // -> Replace with Builder pattern

private:
    friend class Builder;
    friend class ResidueFactory;  // Factory can still set private members

    // ...members...
};

class Residue::Builder {
public:
    Builder(std::string name, int seq_num, char chain_id);

    Builder& insertion(char ins);
    Builder& one_letter_code(char code);
    Builder& type(ResidueType type);
    Builder& is_purine(bool purine);
    Builder& atoms(std::vector<Atom> atoms);
    Builder& add_atom(Atom atom);

    [[nodiscard]] Residue build() const;

private:
    Residue residue_;
};
```

### Step 3.2: Update ResidueFactory to Use Builder

```cpp
[[nodiscard]] static Residue create(...) {
    return Residue::create(name, seq_num, chain_id)
        .insertion(insertion_code)
        .one_letter_code(determine_one_letter_code(name))
        .type(determine_type(name, code))
        .is_purine(determine_is_purine(name, type))
        .atoms(atoms)
        .build();
}
```

### Step 3.3: Fix `from_json` to Use Factory

```cpp
[[nodiscard]] static Residue from_json_legacy(const nlohmann::json& j) {
    // Parse atoms first
    std::vector<Atom> atoms;
    if (j.contains("atoms")) {
        for (const auto& atom_json : j["atoms"]) {
            atoms.push_back(Atom::from_json_legacy(atom_json));
        }
    }

    // Use factory to create with proper properties
    return ResidueFactory::create(name, seq_num, chain_id, ' ', atoms);
}
```

### Files Changed
- UPDATE: `include/x3dna/core/residue.hpp`
- UPDATE: `src/x3dna/core/residue_factory.cpp`

### Unit Test
- UPDATE: `tests/unit/core/residue_test.cpp`

---

## Phase 4: Make BasePair Frames Required

**Goal**: Frames are always required - make the API honest about this.

### Step 4.1: Change BasePair Frame Storage

```cpp
class BasePair {
public:
    // Constructor REQUIRES frames
    BasePair(size_t idx1, size_t idx2,
             const ReferenceFrame& frame1, const ReferenceFrame& frame2,
             BasePairType type = BasePairType::UNKNOWN);

    // Frames are non-optional
    [[nodiscard]] const ReferenceFrame& frame1() const { return frame1_; }
    [[nodiscard]] const ReferenceFrame& frame2() const { return frame2_; }

    // Computed values no longer need has_value() checks
    [[nodiscard]] double origin_distance() const {
        return frame1_.origin().distance_to(frame2_.origin());
    }

private:
    ReferenceFrame frame1_;  // NOT std::optional
    ReferenceFrame frame2_;  // NOT std::optional
    // ...
};
```

### Step 4.2: Update BasePair Construction Sites

All places that create BasePair must now provide frames:

```cpp
// In base_pair_finder.cpp or similar
// BEFORE:
BasePair bp(idx1, idx2, type);
bp.set_frame1(frame1);
bp.set_frame2(frame2);

// AFTER:
BasePair bp(idx1, idx2, frame1, frame2, type);
```

### Step 4.3: Extract hydrogen_bond to Separate File

```cpp
// include/x3dna/core/hydrogen_bond.hpp
namespace x3dna::core {

/**
 * @brief Represents a hydrogen bond in a base pair
 */
struct HydrogenBond {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;  // '-' for standard, ' ' for non-standard
    std::optional<size_t> hbond_idx;

    [[nodiscard]] bool is_standard() const { return type == '-'; }
};

} // namespace x3dna::core
```

### Files Changed
- NEW: `include/x3dna/core/hydrogen_bond.hpp`
- UPDATE: `include/x3dna/core/base_pair.hpp`
- UPDATE: `src/x3dna/algorithms/base_pair_finder.cpp`
- UPDATE: `src/x3dna/algorithms/base_pair_validator.cpp`

### Unit Test
- UPDATE: `tests/unit/core/base_pair_test.cpp`

---

## Phase 5: Clean Container Interfaces

**Goal**: Consistent, safe container access patterns.

### Step 5.1: Remove Mutable Reference Accessors

```cpp
// BEFORE (allows external modification)
std::vector<Residue>& residues() { return residues_; }

// AFTER (read-only access, internal modification via methods)
[[nodiscard]] const std::vector<Residue>& residues() const { return residues_; }
```

Apply to:
- `Structure::chains()`
- `Chain::residues()`
- `Residue::atoms()`

### Step 5.2: Add Iteration Support

```cpp
// In Chain
[[nodiscard]] auto begin() const { return residues_.begin(); }
[[nodiscard]] auto end() const { return residues_.end(); }
[[nodiscard]] size_t size() const { return residues_.size(); }
[[nodiscard]] bool empty() const { return residues_.empty(); }
[[nodiscard]] const Residue& operator[](size_t idx) const { return residues_[idx]; }
```

### Step 5.3: Use Spans Instead of Raw Pointers

```cpp
// BEFORE
[[nodiscard]] std::vector<const Residue*> all_residues() const;

// AFTER (C++20 span or just return const ref to internal vector)
// For filtering, use algorithm + ranges
[[nodiscard]] std::vector<std::reference_wrapper<const Residue>> nucleotides() const;
```

### Files Changed
- UPDATE: `include/x3dna/core/structure.hpp`
- UPDATE: `include/x3dna/core/chain.hpp`
- UPDATE: `include/x3dna/core/residue.hpp`
- UPDATE: Various call sites that used mutable access

---

## Phase 6: Simplify JSON Serialization

**Goal**: Move JSON logic out of core objects.

### Step 6.1: Create Dedicated Serializers

```cpp
// include/x3dna/io/residue_serializer.hpp
namespace x3dna::io {

class ResidueSerializer {
public:
    [[nodiscard]] static nlohmann::json to_legacy_json(const core::Residue& r);
    [[nodiscard]] static nlohmann::json to_modern_json(const core::Residue& r);
    [[nodiscard]] static core::Residue from_legacy_json(const nlohmann::json& j);
    [[nodiscard]] static core::Residue from_modern_json(const nlohmann::json& j);
};

} // namespace x3dna::io
```

### Step 6.2: Keep Simple `to_json()` in Classes

For convenience, keep simple delegation:

```cpp
// In Residue
[[nodiscard]] nlohmann::json to_json() const {
    return io::ResidueSerializer::to_modern_json(*this);
}
```

### Files Changed
- NEW: `include/x3dna/io/residue_serializer.hpp`
- NEW: `include/x3dna/io/atom_serializer.hpp`
- NEW: `include/x3dna/io/base_pair_serializer.hpp`
- NEW: `src/x3dna/io/residue_serializer.cpp`
- NEW: `src/x3dna/io/atom_serializer.cpp`
- NEW: `src/x3dna/io/base_pair_serializer.cpp`
- UPDATE: Core classes to delegate to serializers

---

## Implementation Order

**Recommended sequence** (lowest to highest risk):

1. **Phase 1** (ring atoms) - New file, then update references
2. **Phase 2** (purine/pyrimidine) - Consolidate duplicates
3. **Phase 5.1** (remove mutable refs) - May break code, fix call sites
4. **Phase 3** (Residue builder) - Construction pattern change
5. **Phase 4** (BasePair frames) - API change, update call sites
6. **Phase 5.2-5.3** (iteration/spans) - Container improvements
7. **Phase 6** (JSON serializers) - Extract to separate files

---

## Validation After Each Phase

```bash
make release && make test && fp2-validate validate core --test-set 10
```

---

## Summary of Changes

| Phase | New Files | Files Updated | Key Change |
|-------|-----------|---------------|------------|
| 1 | 2 | 2 | Centralize ring atom definitions |
| 2 | 0 | 4 | Consolidate is_purine() |
| 3 | 0 | 2 | Residue::Builder pattern |
| 4 | 1 | 4 | BasePair frames required, HydrogenBond extracted |
| 5 | 0 | 3+ | Container interface cleanup |
| 6 | 6 | 3 | JSON serializers extracted |
| **Total** | **~9** | **~18** | |

---

## Principles Applied

- **Functions do one thing**: Ring atom check is one function
- **Low complexity**: No deep nesting in new code
- **Good function names**: `is_ring_atom()`, `purine_atoms()`, `atoms_for_type()`
- **3-4 max parameters**: Builder pattern for Residue
- **Max 3 levels indentation**: Early returns, extracted helpers
- **Single Responsibility**: RingAtomRegistry only knows ring atoms
- **Composition over inheritance**: Builder composes Residue
- **const everywhere**: Immutable after construction where possible
- **Strong Encapsulation**: Private data, no mutable refs returned
- **[[nodiscard]]**: On all computed values
- **Early returns**: For validation checks
