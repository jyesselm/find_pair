# Stage 4: Algorithms Part 1 - Frame Calculation

## Objectives

Implement the base frame calculation algorithm, which is critical for all subsequent parameter calculations. This includes loading standard base templates, matching ring atoms, and performing least-squares fitting.

## Duration

**2 weeks**

## Dependencies

- ✅ Stage 0: Project Setup
- ✅ Stage 1: Geometry Classes (especially LeastSquaresFitter)
- ✅ Stage 2: Core Domain Objects
- ✅ Stage 3: I/O Layer

## Tasks

### Task 4.1: Implement Standard Base Template Loader
- [ ] Create `include/x3dna/algorithms/StandardBaseTemplates.hpp`
- [ ] Implement loading of standard base PDB files (Atomic_A.pdb, etc.)
- [ ] Cache loaded templates
- [ ] Handle template path configuration
- [ ] Write unit tests

**Key Methods**:
```cpp
class StandardBaseTemplates {
    Structure load_template(ResidueType type);
    void set_template_path(const std::filesystem::path& path);
    std::filesystem::path get_template_path(ResidueType type) const;
};
```

**Deliverable**: Template loader

### Task 4.2: Implement Ring Atom Matcher
- [ ] Create `include/x3dna/algorithms/RingAtomMatcher.hpp`
- [ ] Implement matching experimental atoms to standard template atoms
- [ ] Handle different base types (purine vs pyrimidine)
- [ ] Handle missing atoms gracefully
- [ ] Return matched atom pairs
- [ ] Write unit tests

**Key Methods**:
```cpp
class RingAtomMatcher {
    struct MatchedAtoms {
        std::vector<Atom> experimental;
        std::vector<Atom> standard;
        std::vector<std::string> atom_names;
    };
    
    MatchedAtoms match(const Residue& residue, const Structure& standard_template);
};
```

**Deliverable**: Ring atom matcher

### Task 4.3: Implement BaseFrameCalculator
- [ ] Create `include/x3dna/algorithms/BaseFrameCalculator.hpp`
- [ ] Implement frame calculation for a residue
- [ ] Use LeastSquaresFitter for alignment
- [ ] Calculate RMS fit quality
- [ ] Store frame in residue
- [ ] Support batch calculation for all residues
- [ ] Add JSON recording (matches legacy `base_frame_calc`, `ls_fitting`, `frame_calc`)
- [ ] Write comprehensive unit tests

**Key Methods**:
```cpp
class BaseFrameCalculator {
    struct FrameCalculationResult {
        ReferenceFrame frame;
        double rms_fit;
        std::vector<std::string> matched_atoms;
        std::string template_file;
    };
    
    ReferenceFrame calculate_frame(Residue& residue);
    FrameCalculationResult calculate_with_metrics(Residue& residue);
    void calculate_all_frames(Structure& structure);
    
    void set_template_path(const std::filesystem::path& path);
};
```

**Critical**: Must match original `base_frame()` algorithm exactly!

**Deliverable**: Fully tested BaseFrameCalculator

### Task 4.4: Integration with Structure
- [ ] Integrate BaseFrameCalculator with Structure
- [ ] Ensure frames are stored correctly
- [ ] Test with real PDB files
- [ ] Compare frames with legacy JSON

**Deliverable**: Integrated frame calculation

### Task 4.5: Performance Optimization (if needed)
- [ ] Profile frame calculation
- [ ] Optimize hot paths
- [ ] Cache template loading
- [ ] Verify performance is acceptable

**Deliverable**: Optimized implementation

## Testing Plan

### Unit Tests

#### StandardBaseTemplates Tests
- [ ] Test loading each base type (A, C, G, T, U)
- [ ] Test template path configuration
- [ ] Test caching
- [ ] Test error handling (missing files)

#### RingAtomMatcher Tests
- [ ] Test matching for each base type
- [ ] Test with complete atom sets
- [ ] Test with missing atoms
- [ ] Test atom name matching
- [ ] Verify matched atoms are correct

#### BaseFrameCalculator Tests
- [ ] Test frame calculation for each base type
- [ ] Test RMS fit calculation
- [ ] Test matched atoms list
- [ ] Test batch calculation
- [ ] Test error handling

### Integration Tests
- [ ] Test complete frame calculation pipeline
- [ ] Test with real PDB files
- [ ] Test with Structure containing multiple residues
- [ ] Verify frames are stored in residues

### Integration Tests (All PDB/JSON Pairs)
- [ ] Create `test_frame_calculation_integration.cpp` using `IntegrationTestBase`
- [ ] Test frame calculation on all discovered PDB/JSON pairs
- [ ] Compare `base_frame_calc` records with legacy JSON
- [ ] Compare `ls_fitting` records (rotation matrices, translations, RMS)
- [ ] Compare `frame_calc` records (complete frame data)
- [ ] Verify rotation matrices match within 0.001
- [ ] Verify origins match within 0.001
- [ ] Verify RMS values match within 0.001
- [ ] Run: `./tests/integration/test_frame_calculation_integration`

### Regression Tests

#### Compare with Legacy JSON
- [ ] Load `data/json_legacy/100D.json`
- [ ] Extract `base_frame_calc` records
- [ ] Extract `ls_fitting` records
- [ ] Extract `frame_calc` records
- [ ] Calculate frames for same residues
- [ ] Compare rotation matrices (within 0.001 tolerance)
- [ ] Compare origins (within 0.001 tolerance)
- [ ] Compare RMS values (within 0.001 tolerance)

**Test Cases**:
```cpp
TEST(FrameCalculationRegression, MatchLegacy) {
    // Load legacy JSON
    auto legacy = load_json("data/json_legacy/100D.json");
    auto frame_calc_records = find_records(legacy, "frame_calc");
    
    // Load PDB and calculate frames
    PdbParser parser;
    Structure structure = parser.parse_file("data/pdb/100D.pdb");
    
    BaseFrameCalculator calculator;
    calculator.calculate_all_frames(structure);
    
    // Compare
    for (size_t i = 0; i < frame_calc_records.size(); ++i) {
        auto legacy_record = frame_calc_records[i];
        size_t residue_idx = legacy_record["residue_idx"];
        
        auto residue = structure.get_residue_by_index(residue_idx);
        auto frame = residue->reference_frame();
        
        // Compare rotation matrix
        auto legacy_orien = legacy_record["orien"];  // From ls_fitting record
        auto our_orien = frame.rotation().to_json_legacy()["orien"];
        compare_matrices(our_orien, legacy_orien, 0.001);
        
        // Compare origin
        auto legacy_org = legacy_record["org"];
        auto our_org = frame.origin().to_json_legacy()["org"];
        compare_vectors(our_org, legacy_org, 0.001);
        
        // Compare RMS
        double legacy_rms = legacy_record["rms_fit"];
        // Get RMS from calculation result
        // ...
    }
}
```

## Success Criteria

- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Code coverage > 85% for algorithms
- [ ] **Critical**: Frame calculations match legacy JSON within 0.001 tolerance
- [ ] Rotation matrices match legacy format
- [ ] Origins match legacy format
- [ ] RMS values match legacy format
- [ ] Performance acceptable
- [ ] Documentation complete

## Code Quality Checks

- [ ] No compiler warnings
- [ ] Passes clang-tidy
- [ ] Follows code style
- [ ] Proper error handling
- [ ] Const-correctness
- [ ] Memory safety (no leaks)

## Deliverables

1. ✅ StandardBaseTemplates implemented
2. ✅ RingAtomMatcher implemented
3. ✅ BaseFrameCalculator fully implemented and tested
4. ✅ Integration with Structure
5. ✅ Comprehensive unit test suite
6. ✅ Integration tests
7. ✅ Regression tests (compare with legacy JSON)
8. ✅ Documentation

## Files Created

```
include/x3dna/algorithms/
├── BaseFrameCalculator.hpp
├── StandardBaseTemplates.hpp
└── RingAtomMatcher.hpp

src/x3dna/algorithms/
├── BaseFrameCalculator.cpp
├── StandardBaseTemplates.cpp
└── RingAtomMatcher.cpp

tests/unit/algorithms/
├── test_base_frame_calculator.cpp
├── test_standard_base_templates.cpp
└── test_ring_atom_matcher.cpp

tests/integration/
└── test_frame_calculation.cpp

tests/regression/
└── test_frame_calculation_regression.cpp
```

## Risks & Mitigation

| Risk | Impact | High | Mitigation |
|------|--------|------|------------|
| Frame calculation doesn't match original | **CRITICAL** | Compare with legacy JSON early, debug step-by-step |
| Least-squares fitting issues | High | Test LeastSquaresFitter thoroughly in Stage 1 |
| Template file path issues | Medium | Use ConfigManager, test path resolution |
| Performance issues | Low | Profile, optimize if needed |

## Validation Checklist

Before moving to Stage 5:
- [ ] All algorithm classes compile and link
- [ ] All unit tests pass
- [ ] Integration tests pass
- [ ] **CRITICAL**: Regression tests match legacy JSON frames
- [ ] Rotation matrices match within 0.001
- [ ] Origins match within 0.001
- [ ] RMS values match within 0.001
- [ ] Code reviewed
- [ ] Documentation complete
- [ ] No memory leaks
- [ ] Performance acceptable

## Next Stage

After completing Stage 4, proceed to **Stage 5: Algorithms Part 2** (`STAGE_05_ALGORITHMS_2.md`)

---

*Estimated Completion: Week 8*

