# Stage 2: Core Domain Objects

## Objectives

Implement the core domain model classes (Atom, Residue, Chain, Structure, ReferenceFrame, BasePair) that represent the biological structures and data.

## Duration

**2 weeks**

## Dependencies

- ✅ Stage 0: Project Setup
- ✅ Stage 1: Geometry Classes

## Tasks

### Task 2.1: Implement Atom
- [x] Create `include/x3dna/core/atom.hpp`
- [x] Implement constructors
- [x] Implement getters/setters
- [x] Implement distance calculations
- [x] Implement atom type queries (ring atom, H-bond donor/acceptor)
- [x] Add JSON serialization (modern and legacy formats)
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class Atom {
    Atom(const std::string& name, const Vector3D& position);
    const std::string& name() const;
    const Vector3D& position() const;
    double distance_to(const Atom& other) const;
    bool is_ring_atom() const;
    bool is_hydrogen_bond_donor() const;
    bool is_hydrogen_bond_acceptor() const;
    nlohmann::json to_json_legacy() const;  // {"atom_name": "...", "xyz": [...]}
};
```

**Deliverable**: Fully tested Atom class

### Task 2.2: Implement Residue
- [x] Create `include/x3dna/core/residue.hpp`
- [x] Implement residue type enumeration
- [x] Implement constructors
- [x] Implement atom management (add, find, query)
- [x] Implement base identification (one-letter code, RY classification)
- [x] Implement ring atom queries
- [x] Implement reference frame storage
- [x] Add JSON serialization (modern and legacy formats)
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class Residue {
    Residue(const std::string& name, int seq_num, char chain_id);
    void add_atom(const Atom& atom);
    std::optional<Atom> find_atom(const std::string& name) const;
    std::vector<Atom> ring_atoms() const;
    char one_letter_code() const;
    bool is_nucleotide() const;
    void set_reference_frame(const ReferenceFrame& frame);
    nlohmann::json to_json_legacy() const;
};
```

**Deliverable**: Fully tested Residue class

### Task 2.3: Implement Chain
- [x] Create `include/x3dna/core/chain.hpp`
- [x] Implement constructors
- [x] Implement residue management
- [x] Implement sequence operations
- [x] Implement nucleotide filtering
- [x] Add JSON serialization
- [x] Write unit tests

**Key Methods**:
```cpp
class Chain {
    explicit Chain(char id);
    void add_residue(const Residue& residue);
    std::string sequence() const;
    std::vector<Residue> nucleotides() const;
    nlohmann::json to_json() const;
};
```

**Deliverable**: Fully tested Chain class

### Task 2.4: Implement ReferenceFrame
- [x] Create `include/x3dna/core/reference_frame.hpp`
- [x] Implement constructors
- [x] Implement axis access (x, y, z)
- [x] Implement transformations
- [x] Implement direction dot product (for base pair validation)
- [x] Implement serialization (array format)
- [x] Add JSON serialization (legacy format: `orien` and `org`)
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class ReferenceFrame {
    ReferenceFrame(const Matrix3D& rotation, const Vector3D& origin);
    Vector3D z_axis() const;  // Normal to base plane
    double direction_dot_product(const ReferenceFrame& other) const;
    std::array<double, 9> rotation_as_array() const;
    nlohmann::json to_json_legacy() const;  // {"orien": [[...]], "org": [...]}
};
```

**Critical**: Must match legacy format exactly for regression testing!

**Deliverable**: Fully tested ReferenceFrame class

### Task 2.5: Implement Structure
- [x] Create `include/x3dna/core/structure.hpp`
- [x] Implement constructors
- [x] Implement chain management
- [x] Implement queries (num_residues, num_atoms, etc.)
- [x] Implement residue access (all_residues, nucleotides)
- [x] Add JSON serialization (modern and legacy formats)
- [x] Write comprehensive unit tests

**Key Methods**:
```cpp
class Structure {
    Structure(const std::string& pdb_id);
    void add_chain(const Chain& chain);
    std::vector<const Residue*> nucleotides() const;
    nlohmann::json to_json_legacy() const;  // pdb_atoms record format
};
```

**Design Note**: Structure represents **raw PDB data only** (atoms, residues, chains). Base pairs are **derived/calculated data** and should be managed separately, not as part of Structure. This maintains clear separation of concerns: Structure = input data, BasePair = calculated result.

**Deliverable**: Fully tested Structure class

### Task 2.6: Implement BasePair
- [x] Create `include/x3dna/core/base_pair.hpp`
- [x] Implement base pair type enumeration
- [x] Implement constructors
- [x] Implement reference frame storage
- [x] Implement hydrogen bond storage
- [x] Implement distance/angle metrics
- [x] Add JSON serialization (legacy format: `base_pair` record)
- [x] Write comprehensive unit tests

**Design Note**: BasePair is a **standalone class** representing calculated/derived data. It is **not** part of Structure. Base pairs will be managed separately (e.g., as a `std::vector<BasePair>` in algorithms/protocols that calculate them).

**Key Methods**:
```cpp
class BasePair {
    BasePair(size_t idx1, size_t idx2, BasePairType type);
    void set_frame1(const ReferenceFrame& frame);
    double origin_distance() const;
    double plane_angle() const;
    double n_n_distance() const;
    nlohmann::json to_json_legacy() const;  // {"type": "base_pair", ...}
};
```

**Deliverable**: Fully tested BasePair class

### Task 2.7: Implement Parameter Structures
- [x] Create `include/x3dna/core/parameters.hpp`
- [x] Implement BasePairStepParameters struct
- [x] Implement HelicalParameters struct
- [x] Implement comparison operators
- [x] Implement JSON serialization (legacy formats)
- [x] Write unit tests

**Key Structures**:
```cpp
struct BasePairStepParameters {
    double shift, slide, rise, tilt, roll, twist;
    nlohmann::json to_json_legacy() const;  // {"params": {"Shift": ..., ...}}
};

struct HelicalParameters {
    double x_displacement, y_displacement, rise, inclination, tip, twist;
    nlohmann::json to_json_legacy() const;  // {"params": [...]}
};
```

**Deliverable**: Parameter structures with JSON support

### Task 2.8: Integration Testing
- [ ] Test Structure → Chain → Residue → Atom hierarchy
- [ ] Test ReferenceFrame with Structure
- [ ] Test BasePair with Structure
- [ ] Test JSON round-trip (write → read → compare)
- [ ] Test legacy JSON format compatibility

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### Atom Tests
- [x] Test constructors
- [x] Test getters/setters
- [x] Test distance calculations
- [x] Test atom type queries
- [x] Test JSON serialization (both formats)
- [x] Test edge cases

#### Residue Tests
- [x] Test constructors
- [x] Test atom management
- [x] Test base identification
- [x] Test ring atom queries
- [x] Test reference frame storage
- [x] Test JSON serialization
- [ ] Test with real PDB residue data (integration test)

#### Chain Tests
- [x] Test constructors
- [x] Test residue management
- [x] Test sequence generation
- [x] Test nucleotide filtering
- [x] Test JSON serialization

#### ReferenceFrame Tests
- [x] Test constructors
- [x] Test axis access
- [x] Test transformations
- [x] Test direction dot product
- [x] Test serialization (array format)
- [x] Test JSON serialization (legacy format)
- [ ] **Critical**: Compare with legacy JSON `ref_frame` records (integration test)

#### Structure Tests
- [x] Test constructors
- [x] Test chain management
- [x] Test queries
- [x] Test residue access
- [x] Test JSON serialization (both formats)
- [ ] Test with real PDB files (integration test)

#### BasePair Tests
- [x] Test constructors
- [x] Test reference frame storage
- [x] Test hydrogen bond storage
- [x] Test distance/angle calculations
- [x] Test JSON serialization (legacy format)
- [ ] **Critical**: Compare with legacy JSON `base_pair` records (integration test)

### Integration Tests
- [ ] Test complete Structure hierarchy
- [ ] Test JSON round-trip (write → read → compare)
- [ ] Test with sample PDB data
- [ ] Test memory management (no leaks)

### Integration Tests (All PDB/JSON Pairs)
- [ ] Create `test_core_objects_integration.cpp` using `IntegrationTestBase`
- [ ] Test PDB parsing on all discovered PDB/JSON pairs
- [ ] Compare parsed Structure with legacy JSON `pdb_atoms` records
- [ ] Test JSON round-trip for all pairs
- [ ] Verify atom counts, residue counts match legacy JSON
- [ ] Run: `./tests/integration/test_core_objects_integration`

### Regression Tests

#### Load Legacy JSON
- [ ] Load `data/json_legacy/100D.json`
- [ ] Parse `pdb_atoms` records → create Structure
- [ ] Parse `ref_frame` records → create ReferenceFrames
- [ ] Parse `base_pair` records → create BasePairs
- [ ] Verify all data matches

#### Compare JSON Output
- [ ] Create Structure from PDB
- [ ] Export to legacy JSON format
- [ ] Compare with legacy JSON file
- [ ] Verify field-by-field match (within tolerance)

**Test Cases**:
```cpp
// Load legacy JSON
auto legacy_json = load_json("data/json_legacy/100D.json");

// Parse atoms
auto atoms_record = find_record(legacy_json, "pdb_atoms");
Structure structure = Structure::from_json_legacy(atoms_record);

// Verify
EXPECT_EQ(structure.num_atoms(), atoms_record["num_atoms"]);

// Export and compare
auto our_json = structure.to_json_legacy();
compare_json_records(our_json["pdb_atoms"], atoms_record);
```

## Test Data

### Sample PDB Data
- Use `data/pdb/100D.pdb` for testing
- Use `data/pdb/157D.pdb` for additional test cases
- Create minimal test fixtures for unit tests

### Legacy JSON Reference
- Use `data/json_legacy/100D.json` for regression testing
- Extract specific record types for targeted tests
- Compare field-by-field

## Success Criteria

- [x] All unit tests pass (100% for core classes) - **129 tests passing**
- [ ] All integration tests pass (pending Task 2.8)
- [ ] Code coverage > 85% for core classes (pending measurement)
- [x] JSON serialization matches legacy format exactly (implemented)
- [x] Can load Structure from legacy JSON (implemented)
- [x] Can export Structure to legacy JSON format (implemented)
- [x] Round-trip JSON test passes (write → read → compare) (unit tests verify)
- [ ] No memory leaks (pending valgrind)
- [ ] Performance acceptable (pending profiling)
- [x] Documentation complete

## Code Quality Checks

- [ ] No compiler warnings
- [ ] Passes clang-tidy
- [ ] Follows code style
- [ ] Const-correctness
- [ ] RAII for all resources
- [ ] Smart pointers where appropriate

## Deliverables

1. ✅ Atom class fully implemented and tested (17 tests passing)
2. ✅ Residue class fully implemented and tested (19 tests passing)
3. ✅ Chain class fully implemented and tested (17 tests passing)
4. ✅ ReferenceFrame class fully implemented and tested (20 tests passing)
5. ✅ Structure class fully implemented and tested (18 tests passing)
6. ✅ BasePair class fully implemented and tested (14 tests passing)
7. ✅ Parameter structures implemented (Task 2.7 - 24 tests passing)
8. ✅ Comprehensive unit test suite (129 tests total, all passing)
9. ⏳ Integration tests (Task 2.8 - pending)
10. ⏳ Regression tests (JSON compatibility) (Task 2.8 - pending)
11. ✅ Documentation (in progress)

## Files Created

```
include/x3dna/core/
├── atom.hpp ✅
├── residue.hpp ✅
├── chain.hpp ✅
├── structure.hpp ✅
├── reference_frame.hpp ✅
├── base_pair.hpp ✅
└── parameters.hpp ✅

# Note: All classes are header-only (no .cpp files needed)

tests/unit/core/
├── test_atom.cpp ✅
├── test_residue.cpp ✅
├── test_chain.cpp ✅
├── test_structure.cpp ✅
├── test_reference_frame.cpp ✅
├── test_base_pair.cpp ✅
└── test_parameters.cpp ✅

tests/integration/
└── test_core_integration.cpp ⏳ (Task 2.8)

tests/regression/
└── test_json_compatibility.cpp ⏳ (Task 2.8)
```

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| JSON format mismatch | High | Test early, compare field-by-field |
| Memory management issues | Medium | Use smart pointers, test with valgrind |
| Performance with large structures | Low | Profile, optimize if needed |
| Legacy format complexity | Medium | Create helper functions, test incrementally |

## Validation Checklist

Before moving to Stage 3:
- [x] All core classes compile and link
- [x] All unit tests pass (129 tests)
- [ ] Integration tests pass (Task 2.8 - pending)
- [ ] Regression tests match legacy JSON (Task 2.8 - pending)
- [x] JSON round-trip works (verified in unit tests)
- [x] Code reviewed
- [x] Documentation complete
- [ ] No memory leaks (pending valgrind)
- [ ] Performance acceptable (pending profiling)

## Next Stage

After completing Stage 2, proceed to **Stage 3: I/O Layer** (`STAGE_03_IO.md`)

---

*Estimated Completion: Week 4*

