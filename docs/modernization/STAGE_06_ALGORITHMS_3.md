# Stage 6: Algorithms Part 3 - Parameter Calculation

## Objectives

Implement parameter calculation algorithms for base pair step parameters and helical parameters. This is the final output of the analysis pipeline.

## Duration

**2 weeks**

## Dependencies

- ✅ Stage 0-5: All previous stages
- ✅ Stage 5: Base pairs must be found first

## Tasks

### Task 6.1: Implement ParameterCalculator
- [ ] Create `include/x3dna/algorithms/ParameterCalculator.hpp`
- [ ] Implement `bpstep_par()` equivalent (6 step parameters)
- [ ] Implement `helical_par()` equivalent (6 helical parameters)
- [ ] Calculate midstep frames
- [ ] Calculate helical midstep frames
- [ ] Support batch calculation for all steps
- [ ] Add JSON recording (matches legacy `bpstep_params`, `helical_params`)
- [ ] Write comprehensive unit tests

**Key Methods**:
```cpp
class ParameterCalculator {
    BasePairStepParameters calculate_step_parameters(
        const BasePair& pair1,
        const BasePair& pair2
    );
    
    HelicalParameters calculate_helical_parameters(
        const BasePair& pair1,
        const BasePair& pair2
    );
    
    ReferenceFrame calculate_midstep_frame(
        const ReferenceFrame& frame1,
        const ReferenceFrame& frame2
    );
    
    std::vector<BasePairStepParameters> calculate_all_step_parameters(
        const std::vector<BasePair>& pairs
    );
    
private:
    // Core algorithm (matches original bpstep_par)
    void bpstep_par_impl(
        const Matrix3D& r1, const Vector3D& o1,
        const Matrix3D& r2, const Vector3D& o2,
        BasePairStepParameters& params,
        ReferenceFrame& midstep_frame
    );
};
```

**CRITICAL**: Must match original `bpstep_par()` algorithm exactly!

**Deliverable**: Fully tested ParameterCalculator

### Task 6.2: Verify Algorithm Correctness
- [ ] Compare with original implementation step-by-step
- [ ] Debug any discrepancies
- [ ] Verify all 6 parameters match
- [ ] Verify midstep frames match

**Deliverable**: Verified correct implementation

### Task 6.3: Integration Testing
- [ ] Test complete parameter calculation pipeline
- [ ] Test with real PDB files
- [ ] Compare results with legacy JSON
- [ ] Verify all step parameters

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### ParameterCalculator Tests
- [ ] Test step parameter calculation
- [ ] Test helical parameter calculation
- [ ] Test midstep frame calculation
- [ ] Test with known frame pairs
- [ ] Test batch calculation
- [ ] **Critical**: Compare with legacy `bpstep_params` records
- [ ] **Critical**: Compare with legacy `helical_params` records

**Test Cases**:
```cpp
TEST(ParameterCalculation, StepParameters) {
    // Create test base pairs with known frames
    BasePair pair1 = create_test_pair1();
    BasePair pair2 = create_test_pair2();
    
    ParameterCalculator calculator;
    auto params = calculator.calculate_step_parameters(pair1, pair2);
    
    // Verify all 6 parameters are calculated
    EXPECT_NE(params.shift, 0.0);
    EXPECT_NE(params.slide, 0.0);
    // ... etc
}
```

### Integration Tests
- [ ] Test complete parameter calculation workflow
- [ ] Test with real PDB files
- [ ] Verify parameters are stored correctly
- [ ] Test with multiple base pairs

### Integration Tests (All PDB/JSON Pairs)
- [ ] Create `test_parameter_integration.cpp` using `IntegrationTestBase`
- [ ] Test parameter calculation on all discovered PDB/JSON pairs
- [ ] Compare `bpstep_params` records with legacy JSON
- [ ] Compare `helical_params` records with legacy JSON
- [ ] Verify all 6 step parameters match within 0.001 (Shift, Slide, Rise, Tilt, Roll, Twist)
- [ ] Verify all 6 helical parameters match within 0.001
- [ ] Verify midstep frames match within 0.001
- [ ] **CRITICAL**: This is the final output - must match exactly!
- [ ] Run: `./tests/integration/test_parameter_integration`

### Regression Tests

#### Compare with Legacy JSON
- [ ] Load `data/json_legacy/100D.json`
- [ ] Extract `bpstep_params` records
- [ ] Extract `helical_params` records
- [ ] Calculate parameters with our code
- [ ] Compare all 6 step parameters (within 0.001)
- [ ] Compare all 6 helical parameters (within 0.001)
- [ ] Compare midstep frames (within 0.001)

**Test Cases**:
```cpp
TEST(ParameterRegression, MatchLegacy) {
    // Load legacy JSON
    auto legacy = load_json("data/json_legacy/100D.json");
    auto legacy_params = extract_bpstep_params(legacy);
    
    // Calculate with our code
    Structure structure = load_and_process("data/pdb/100D.pdb");
    ParameterCalculator calculator;
    auto our_params = calculator.calculate_all_step_parameters(
        structure.base_pairs()
    );
    
    // Compare
    EXPECT_EQ(our_params.size(), legacy_params.size());
    
    for (size_t i = 0; i < our_params.size(); ++i) {
        EXPECT_NEAR(our_params[i].shift, legacy_params[i]["params"]["Shift"], 0.001);
        EXPECT_NEAR(our_params[i].slide, legacy_params[i]["params"]["Slide"], 0.001);
        EXPECT_NEAR(our_params[i].rise, legacy_params[i]["params"]["Rise"], 0.001);
        EXPECT_NEAR(our_params[i].tilt, legacy_params[i]["params"]["Tilt"], 0.001);
        EXPECT_NEAR(our_params[i].roll, legacy_params[i]["params"]["Roll"], 0.001);
        EXPECT_NEAR(our_params[i].twist, legacy_params[i]["params"]["Twist"], 0.001);
        
        // Compare midstep frames
        compare_frames(
            calculator.calculate_midstep_frame(pair1.frame(), pair2.frame()),
            legacy_params[i]["mst_orien"],
            0.001
        );
    }
}
```

## Success Criteria

- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Code coverage > 85%
- [ ] **CRITICAL**: All 6 step parameters match legacy JSON within 0.001
- [ ] **CRITICAL**: All 6 helical parameters match legacy JSON within 0.001
- [ ] Midstep frames match legacy JSON
- [ ] Performance acceptable
- [ ] Documentation complete

## Deliverables

1. ✅ ParameterCalculator fully implemented and tested
2. ✅ Step parameter calculation verified
3. ✅ Helical parameter calculation verified
4. ✅ Comprehensive unit test suite
5. ✅ Integration tests
6. ✅ Regression tests (compare with legacy JSON)
7. ✅ Documentation

## Files Created

```
include/x3dna/algorithms/
└── ParameterCalculator.hpp

src/x3dna/algorithms/
└── ParameterCalculator.cpp

tests/unit/algorithms/
└── test_parameter_calculator.cpp

tests/integration/
└── test_parameter_calculation.cpp

tests/regression/
└── test_parameter_regression.cpp
```

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| Parameters don't match original | **CRITICAL** | Compare step-by-step, debug carefully |
| Algorithm implementation errors | High | Test with known values, verify math |
| Frame extraction issues | High | Verify frames are correct from Stage 4 |

## Validation Checklist

Before moving to Stage 7:
- [ ] All algorithm classes compile and link
- [ ] All unit tests pass
- [ ] Integration tests pass
- [ ] **CRITICAL**: Regression tests match legacy JSON parameters
- [ ] All 6 step parameters match within 0.001
- [ ] All 6 helical parameters match within 0.001
- [ ] Midstep frames match
- [ ] Code reviewed
- [ ] Documentation complete
- [ ] No memory leaks
- [ ] Performance acceptable

## Next Stage

After completing Stage 6, proceed to **Stage 7: Protocols** (`STAGE_07_PROTOCOLS.md`)

---

*Estimated Completion: Week 12*

