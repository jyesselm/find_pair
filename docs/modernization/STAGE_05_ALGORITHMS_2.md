# Stage 5: Algorithms Part 2 - Base Pair Finding

## Objectives

Implement base pair detection algorithms, including hydrogen bond validation, base pair finding strategies, and helix detection.

## Duration

**2 weeks**

## Dependencies

- ✅ Stage 0-4: All previous stages
- ✅ Stage 4: BaseFrameCalculator (frames must be calculated first)

## Tasks

### Task 5.1: Implement HydrogenBondValidator
- [x] Create `include/x3dna/algorithms/HydrogenBondValidator.hpp` (integrated into BasePairValidator)
- [x] Implement H-bond detection between two residues
- [x] Implement distance/angle validation
- [x] Identify donor/acceptor pairs
- [x] Count H-bonds
- [x] Return detailed H-bond information
- [ ] Add JSON recording (matches legacy `hbond_list` format) - pending
- [ ] Write comprehensive unit tests - pending

**Key Methods**:
```cpp
class HydrogenBondValidator {
    struct HydrogenBond {
        std::string donor_atom;
        std::string acceptor_atom;
        double distance;
        double angle;
        bool is_standard;
        char type;  // '-' or ' '
        int linkage_type;
    };
    
    std::vector<HydrogenBond> find_hydrogen_bonds(
        const Residue& res1,
        const Residue& res2
    ) const;
    
    bool has_sufficient_hbonds(
        const Residue& res1,
        const Residue& res2,
        int min_count = 1
    ) const;
    
    void set_max_distance(double distance);
    void set_min_angle(double angle);
};
```

**Deliverable**: Fully tested HydrogenBondValidator

### Task 5.2: Implement BasePairValidator
- [x] Create `include/x3dna/algorithms/BasePairValidator.hpp`
- [x] Implement `check_pair()` equivalent
- [x] Check distance constraints (dorg, dNN)
- [x] Check plane angle
- [x] Check vertical distance (d_v)
- [x] Validate hydrogen bonds
- [x] Calculate direction vectors (dir_x, dir_y, dir_z)
- [x] Return validation result with details
- [ ] Add JSON recording (matches legacy `pair_validation` format) - pending
- [ ] Write comprehensive unit tests - pending

**Key Methods**:
```cpp
class BasePairValidator {
    struct ValidationResult {
        bool is_valid;
        BasePairType type;
        double dir_z;
        double dorg;
        double plane_angle;
        double dNN;
        std::vector<HydrogenBond> hbonds;
        // ... validation checks
    };
    
    ValidationResult validate(
        const Residue& res1,
        const Residue& res2
    ) const;
};
```

**Critical**: Must match original `check_pair()` algorithm exactly!

**Deliverable**: Fully tested BasePairValidator

### Task 5.3: Implement BasePairFinder
- [x] Create `include/x3dna/algorithms/BasePairFinder.hpp`
- [x] Implement best pair strategy (greedy mutual best match)
- [x] Implement all pairs strategy (exhaustive search)
- [ ] Implement distance-based strategy - pending
- [x] Find best partner for a residue
- [x] Find all valid pairs
- [x] Return vector of BasePair objects
- [x] Integration with generate_modern_json.cpp
- [ ] Write comprehensive unit tests - pending

**Key Methods**:
```cpp
class BasePairFinder {
    enum class Strategy {
        BEST_PAIR,
        ALL_PAIRS,
        DISTANCE_BASED
    };
    
    std::vector<BasePair> find_pairs(Structure& structure);
    void set_strategy(Strategy strategy);
    
private:
    std::vector<BasePair> find_best_pairs(Structure& structure);
    std::vector<BasePair> find_all_pairs(const Structure& structure);
    BasePair find_best_partner(const Residue& residue, const Structure& structure) const;
};
```

**Deliverable**: Fully tested BasePairFinder

### Task 5.4: Implement HelixDetector
- [ ] Create `include/x3dna/algorithms/HelixDetector.hpp`
- [ ] Implement helix detection algorithm
- [ ] Implement base pair context analysis
- [ ] Implement 5'→3' reordering
- [ ] Identify helical regions
- [ ] Handle circular structures
- [ ] Write comprehensive unit tests

**Key Methods**:
```cpp
class HelixDetector {
    struct Helix {
        std::vector<size_t> base_pair_indices;
        size_t start_index;
        size_t end_index;
        bool is_circular;
    };
    
    std::vector<Helix> detect_helices(const Structure& structure);
    void reorder_base_pairs(Structure& structure);
    void ensure_five_to_three_ordering(std::vector<BasePair>& pairs);
};
```

**Deliverable**: Fully tested HelixDetector

### Task 5.5: Integration Testing
- [x] Test complete base pair finding pipeline (basic testing done)
- [x] Test with real PDB files (integrated in generate_modern_json)
- [ ] Compare results with legacy JSON - pending (need to check for base_pair records or .inp files)
- [ ] Verify helix detection - pending (requires HelixDetector)

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### HydrogenBondValidator Tests
- [ ] Test H-bond detection for known pairs (A-T, G-C)
- [ ] Test distance validation
- [ ] Test angle validation
- [ ] Test donor/acceptor identification
- [ ] Test H-bond counting
- [ ] Test edge cases (no H-bonds, multiple H-bonds)

#### BasePairValidator Tests
- [ ] Test validation for valid pairs
- [ ] Test validation for invalid pairs
- [ ] Test distance checks (dorg, dNN)
- [ ] Test plane angle check
- [ ] Test H-bond validation
- [ ] Test direction vector calculation
- [ ] **Critical**: Compare with legacy `pair_validation` records

#### BasePairFinder Tests
- [ ] Test best pair strategy
- [ ] Test all pairs strategy
- [ ] Test with known structures
- [ ] Test mutual best match logic
- [ ] Test edge cases (no pairs, single pair)

#### HelixDetector Tests
- [ ] Test helix detection
- [ ] Test reordering
- [ ] Test circular structure handling
- [ ] Test context analysis

### Integration Tests
- [ ] Test complete base pair finding workflow
- [ ] Test with real PDB files
- [ ] Verify pairs are stored in Structure
- [ ] Verify helix information

### Integration Tests (All PDB/JSON Pairs)
- [ ] Create `test_base_pair_integration.cpp` using `IntegrationTestBase`
- [ ] Test base pair finding on all discovered PDB/JSON pairs
- [ ] Compare `base_pair` records with legacy JSON
- [ ] Compare `pair_validation` records (validation results)
- [ ] Compare `hbond_list` records (hydrogen bond information)
- [ ] Compare `distance_checks` records
- [ ] Verify base pair indices match
- [ ] Verify frames match within 0.001
- [ ] Verify H-bond counts and distances match
- [ ] Run: `./tests/integration/test_base_pair_integration`

### Regression Tests

#### Compare with Legacy JSON
- [ ] Load `data/json_legacy/100D.json`
- [ ] Extract `base_pair` records
- [ ] Extract `pair_validation` records
- [ ] Extract `hbond_list` records
- [ ] Find base pairs with our code
- [ ] Compare base pair indices
- [ ] Compare validation results
- [ ] Compare H-bond information

**Test Cases**:
```cpp
TEST(BasePairRegression, MatchLegacy) {
    // Load legacy JSON
    auto legacy = load_json("data/json_legacy/100D.json");
    auto legacy_pairs = extract_base_pairs(legacy);
    
    // Find pairs with our code
    Structure structure = load_and_process("data/pdb/100D.pdb");
    BasePairFinder finder;
    auto our_pairs = finder.find_pairs(structure);
    
    // Compare
    EXPECT_EQ(our_pairs.size(), legacy_pairs.size());
    
    for (size_t i = 0; i < our_pairs.size(); ++i) {
        EXPECT_EQ(our_pairs[i].residue_index1(), legacy_pairs[i]["base_i"]);
        EXPECT_EQ(our_pairs[i].residue_index2(), legacy_pairs[i]["base_j"]);
        
        // Compare frames
        compare_frames(our_pairs[i].frame1(), legacy_pairs[i]["orien_i"], 0.001);
    }
}
```

## Success Criteria

- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Code coverage > 85%
- [ ] **Critical**: Base pairs match legacy JSON
- [ ] Validation results match legacy JSON
- [ ] H-bond information matches legacy JSON
- [ ] Performance acceptable
- [ ] Documentation complete

## Deliverables

1. ⏳ HydrogenBondValidator implemented and tested (integrated into BasePairValidator)
2. ✅ BasePairValidator implemented (matches legacy check_pair)
3. ✅ BasePairFinder implemented (matches legacy find_bestpair)
4. ⏳ HelixDetector implemented and tested (pending)
5. ⏳ Comprehensive unit test suite (pending)
6. ⏳ Integration tests (pending)
7. ⏳ Regression tests (pending)
8. ✅ Documentation (in progress)

## Files Created

```
include/x3dna/algorithms/
├── BasePairValidator.hpp ✅ (H-bond validation integrated)
├── BasePairFinder.hpp ✅
└── HelixDetector.hpp ⏳ (pending)

src/x3dna/algorithms/
├── BasePairValidator.cpp ✅
├── BasePairFinder.cpp ✅
└── HelixDetector.cpp ⏳ (pending)

tests/unit/algorithms/
├── test_base_pair_validator.cpp ⏳ (pending)
├── test_base_pair_finder.cpp ⏳ (pending)
└── test_helix_detector.cpp ⏳ (pending)

tests/integration/
└── test_base_pair_integration.cpp ⏳ (pending)
```

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Validation doesn't match original | High | Compare with legacy JSON early, debug step-by-step |
| H-bond detection issues | High | Test with known pairs, compare with legacy |
| Performance with large structures | Medium | Profile, optimize if needed |

## Next Stage

After completing Stage 5, proceed to **Stage 6: Algorithms Part 3** (`STAGE_06_ALGORITHMS_3.md`)

---

*Estimated Completion: Week 10*

