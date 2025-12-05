# Overall Testing Strategy

## Testing Philosophy

1. **Test-Driven Development**: Write tests alongside implementation
2. **Incremental Testing**: Test each stage before moving to the next
3. **Regression Testing**: Compare with legacy JSON files at every stage
4. **Comprehensive Coverage**: Aim for >80% code coverage
5. **Real Data Testing**: Use actual PDB and JSON files from the project

## Testing Levels

### 1. Unit Tests

**Purpose**: Test individual classes and functions in isolation

**Framework**: Google Test

**Location**: `tests/unit/`

**Coverage Target**: >90% for core classes, >80% overall

**When to Write**: During implementation of each class

**Example Structure**:
```cpp
TEST(AtomTest, Construction) {
    Atom atom(" C1'", 1.0, 2.0, 3.0);
    EXPECT_EQ(atom.name(), " C1'");
    EXPECT_DOUBLE_EQ(atom.x(), 1.0);
}
```

### 2. Integration Tests

**Purpose**: Test interactions between multiple components

**Framework**: Google Test

**Location**: `tests/integration/`

**Coverage Target**: All major workflows

**When to Write**: After completing a stage

**Example Structure**:
```cpp
TEST(IntegrationTest, PdbToStructureToJson) {
    PdbParser parser;
    Structure structure = parser.parse_file("data/pdb/100D.pdb");
    
    JsonWriter writer;
    writer.record_pdb_atoms(structure);
    auto json = writer.get_json();
    
    // Verify structure preserved
    EXPECT_GT(json["calculations"].size(), 0);
}
```

### 3. Regression Tests

**Purpose**: Ensure output matches legacy JSON files exactly

**Framework**: Google Test + custom comparison utilities

**Location**: `tests/regression/`

**Coverage Target**: All calculation types

**When to Write**: After implementing algorithms

**Example Structure**:
```cpp
TEST(RegressionTest, FrameCalculation) {
    // Load legacy JSON
    auto legacy = load_json("data/json_legacy/100D.json");
    auto legacy_frames = extract_frames(legacy);
    
    // Calculate with our code
    Structure structure = load_and_process("data/pdb/100D.pdb");
    auto our_frames = extract_frames(structure);
    
    // Compare
    for (size_t i = 0; i < legacy_frames.size(); ++i) {
        compare_frames(our_frames[i], legacy_frames[i], 0.001);
    }
}
```

## Test Data

### PDB Files
- `data/pdb/100D.pdb` - Primary test file
- `data/pdb/157D.pdb` - Secondary test file
- Additional files as needed

### Legacy JSON Files
- `data/json_legacy/100D.json` - Primary reference
- `data/json_legacy/100D_globals.json` - Global variables
- `data/json_legacy/157D.json` - Secondary reference
- `data/json_legacy/157D_globals.json` - Global variables

### Test Fixtures
- Minimal test data for unit tests
- Synthetic data for edge cases
- Real data for integration/regression tests

## Testing Tools

### Test Framework
- **Google Test**: Primary testing framework
- **Google Mock**: For mocking (if needed)

### Code Coverage
- **gcov/lcov**: Code coverage analysis
- **Target**: >80% overall, >90% for critical algorithms

### Memory Checking
- **valgrind**: Memory leak detection
- **AddressSanitizer**: Runtime error detection
- **Target**: Zero memory leaks, zero errors

### Static Analysis
- **clang-tidy**: Static code analysis
- **clang-format**: Code formatting
- **Target**: Zero warnings, consistent style

### Performance Testing
- **Google Benchmark**: Performance benchmarks
- **Target**: Comparable or better than original

## Regression Testing Strategy

### Comparison Tolerances

| Data Type | Tolerance | Notes |
|-----------|-----------|-------|
| Rotation matrices | 0.001 | Element-wise comparison |
| Coordinates | 0.001 Å | For origins, positions |
| Angles | 0.001° | For tilt, roll, twist |
| Distances | 0.001 Å | For dorg, dNN, etc. |
| RMS values | 0.001 | For fit quality |

### Comparison Functions

```cpp
// Compare rotation matrices
bool compare_matrices(const Matrix3D& m1, const Matrix3D& m2, double tolerance) {
    for (size_t i = 0; i < 9; ++i) {
        if (std::abs(m1.as_array()[i] - m2.as_array()[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

// Compare vectors
bool compare_vectors(const Vector3D& v1, const Vector3D& v2, double tolerance) {
    return (v1 - v2).length() < tolerance;
}

// Compare JSON records
void compare_json_records(const nlohmann::json& j1, const nlohmann::json& j2, double tolerance) {
    // Field-by-field comparison
    // ...
}
```

### Test Categories

#### 1. Frame Calculation Tests
- Compare `base_frame_calc` records
- Compare `ls_fitting` records
- Compare `frame_calc` records
- Verify rotation matrices match
- Verify origins match
- Verify RMS values match

#### 2. Base Pair Tests
- Compare `base_pair` records
- Compare `pair_validation` records
- Verify base indices match
- Verify frames match
- Verify validation results match

#### 3. Parameter Tests
- Compare `bpstep_params` records
- Compare `helical_params` records
- Verify all 6 parameters match
- Verify midstep frames match

#### 4. Hydrogen Bond Tests
- Compare `hbond_list` records
- Verify H-bond counts match
- Verify distances match
- Verify atom pairs match

## Test Execution

### Running Tests

```bash
# Build tests
cd build
cmake ..
make

# Run all tests
ctest

# Run specific test
./tests/unit/test_atom

# Run with coverage
make coverage
```

### Continuous Integration

- Run tests on every commit
- Run regression tests nightly
- Generate coverage reports
- Check for memory leaks

## Test Maintenance

### Adding New Tests
1. Create test file in appropriate directory
2. Follow naming convention: `test_<class_name>.cpp`
3. Use descriptive test names
4. Add to CMakeLists.txt
5. Document test purpose

### Updating Tests
- Update when API changes
- Update when requirements change
- Keep tests in sync with implementation

### Test Data Management
- Version control test data
- Document test data sources
- Keep test data minimal but representative

## Success Metrics

### Coverage Metrics
- Overall code coverage: >80%
- Critical algorithm coverage: >90%
- Core class coverage: >90%

### Quality Metrics
- Zero memory leaks (valgrind clean)
- Zero undefined behavior (AddressSanitizer clean)
- Zero compiler warnings
- All tests pass

### Regression Metrics
- 100% of legacy JSON records match
- All numerical values within tolerance
- All record types validated

## Test Documentation

### Test Plan Documents
- Each stage has detailed testing plan
- Document test cases
- Document test data
- Document expected results

### Test Results
- Track test results over time
- Document failures and fixes
- Maintain test history

## Best Practices

1. **Write Tests First**: When possible, write tests before implementation
2. **Test Edge Cases**: Don't just test happy path
3. **Use Real Data**: Test with actual PDB/JSON files
4. **Compare with Legacy**: Always compare with legacy JSON
5. **Document Tests**: Explain what each test validates
6. **Keep Tests Fast**: Unit tests should run quickly
7. **Isolate Tests**: Tests should not depend on each other
8. **Clean Up**: Tests should clean up after themselves

## Troubleshooting

### Tests Fail
1. Check test data is correct
2. Check tolerance values
3. Check comparison logic
4. Debug step-by-step
5. Compare intermediate values

### Performance Issues
1. Profile slow tests
2. Optimize hot paths
3. Use test fixtures efficiently
4. Consider test parallelization

### Memory Issues
1. Run valgrind
2. Check for leaks
3. Verify RAII usage
4. Check smart pointer usage

---

*This testing strategy ensures the modernized codebase matches the original implementation exactly while maintaining high code quality.*

