# Integration Testing Framework

## Overview

Comprehensive integration testing framework that automatically discovers all PDB files with corresponding JSON files and runs comparison tests at each stage of development.

## Philosophy

- **Run on every PDB/JSON pair**: Automatically discover and test all available pairs
- **Stage-by-stage validation**: Test at each stage where comparison is possible
- **Comprehensive comparison**: Compare all calculation types with legacy JSON
- **Automated discovery**: No manual test case maintenance needed

## Test Discovery

### PDB/JSON Pair Discovery

The framework automatically discovers PDB files that have corresponding JSON files:

```
data/pdb/100D.pdb  →  data/json_legacy/100D.json  ✅ Test
data/pdb/157D.pdb  →  data/json_legacy/157D.json  ✅ Test
data/pdb/1EHZ.pdb  →  data/json_legacy/1EHZ.json  ❌ Skip (no JSON)
```

### Discovery Utility

```cpp
// tests/integration/TestDataDiscovery.hpp
class TestDataDiscovery {
public:
    struct PdbJsonPair {
        std::filesystem::path pdb_file;
        std::filesystem::path json_file;
        std::filesystem::path globals_file;
        std::string pdb_name;
    };
    
    static std::vector<PdbJsonPair> discover_pairs(
        const std::filesystem::path& pdb_dir,
        const std::filesystem::path& json_dir
    );
    
    static bool has_json(const std::filesystem::path& pdb_file,
                        const std::filesystem::path& json_dir);
};
```

## Integration Test Structure

### Base Integration Test Class

```cpp
// tests/integration/IntegrationTestBase.hpp
class IntegrationTestBase : public ::testing::Test {
protected:
    void SetUp() override {
        // Discover all PDB/JSON pairs
        pairs_ = TestDataDiscovery::discover_pairs(
            "data/pdb",
            "data/json_legacy"
        );
        
        ASSERT_GT(pairs_.size(), 0) 
            << "No PDB/JSON pairs found for testing";
    }
    
    std::vector<TestDataDiscovery::PdbJsonPair> pairs_;
    
    // Helper methods
    nlohmann::json load_legacy_json(const std::filesystem::path& json_file);
    Structure load_and_process_pdb(const std::filesystem::path& pdb_file);
    void compare_calculations(const nlohmann::json& legacy,
                             const nlohmann::json& ours,
                             double tolerance = 0.001);
};
```

## Stage-Specific Integration Tests

### Stage 2: Core Objects Integration

```cpp
// tests/integration/test_core_objects_integration.cpp
class CoreObjectsIntegrationTest : public IntegrationTestBase {};

TEST_F(CoreObjectsIntegrationTest, ParsePDBMatchesLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        auto legacy_atoms = find_record(legacy_json, "pdb_atoms");
        
        // Parse PDB with our code
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        // Compare
        EXPECT_EQ(structure.num_atoms(), legacy_atoms["num_atoms"]);
        
        // Compare atom data
        compare_atoms(structure, legacy_atoms, 0.001);
    }
}

TEST_F(CoreObjectsIntegrationTest, JSONRoundTrip) {
    for (const auto& pair : pairs_) {
        // Parse PDB
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        // Export to JSON
        auto our_json = structure.to_json_legacy();
        
        // Load back
        Structure restored = Structure::from_json_legacy(our_json);
        
        // Verify round-trip
        EXPECT_EQ(structure.num_atoms(), restored.num_atoms());
        EXPECT_EQ(structure.num_residues(), restored.num_residues());
    }
}
```

### Stage 3: I/O Integration

```cpp
// tests/integration/test_io_integration.cpp
class IOIntegrationTest : public IntegrationTestBase {};

TEST_F(IOIntegrationTest, PDBToJSONMatchesLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        
        // Parse PDB and export to JSON
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        JsonWriter writer;
        writer.init(pair.pdb_file);
        writer.record_pdb_atoms(structure);
        auto our_json = writer.get_json();
        
        // Compare
        auto legacy_atoms = find_record(legacy_json, "pdb_atoms");
        auto our_atoms = find_record(our_json, "pdb_atoms");
        
        compare_json_records(our_atoms, legacy_atoms, 0.001);
    }
}
```

### Stage 4: Frame Calculation Integration

```cpp
// tests/integration/test_frame_calculation_integration.cpp
class FrameCalculationIntegrationTest : public IntegrationTestBase {};

TEST_F(FrameCalculationIntegrationTest, FramesMatchLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        auto legacy_frames = extract_frames(legacy_json);
        
        // Calculate frames with our code
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        BaseFrameCalculator calculator;
        calculator.calculate_all_frames(structure);
        
        // Compare each frame
        for (const auto& legacy_frame : legacy_frames) {
            size_t residue_idx = legacy_frame["residue_idx"];
            auto residue = structure.get_residue_by_index(residue_idx);
            
            if (residue && residue->reference_frame()) {
                auto our_frame = residue->reference_frame().value();
                
                // Compare rotation matrix
                compare_matrices(
                    our_frame.rotation(),
                    legacy_frame["orien"],
                    0.001
                );
                
                // Compare origin
                compare_vectors(
                    our_frame.origin(),
                    legacy_frame["org"],
                    0.001
                );
            }
        }
    }
}
```

### Stage 5: Base Pair Finding Integration

```cpp
// tests/integration/test_base_pair_integration.cpp
class BasePairIntegrationTest : public IntegrationTestBase {};

TEST_F(BasePairIntegrationTest, BasePairsMatchLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        auto legacy_pairs = extract_base_pairs(legacy_json);
        
        // Find base pairs with our code
        Structure structure = load_and_process_pdb(pair.pdb_file);
        BasePairFinder finder;
        auto our_pairs = finder.find_pairs(structure);
        
        // Compare
        EXPECT_EQ(our_pairs.size(), legacy_pairs.size());
        
        for (size_t i = 0; i < our_pairs.size(); ++i) {
            EXPECT_EQ(our_pairs[i].residue_index1(), legacy_pairs[i]["base_i"]);
            EXPECT_EQ(our_pairs[i].residue_index2(), legacy_pairs[i]["base_j"]);
            
            // Compare frames
            compare_frames(our_pairs[i].frame1(), legacy_pairs[i]["orien_i"], 0.001);
            compare_frames(our_pairs[i].frame2(), legacy_pairs[i]["orien_j"], 0.001);
        }
    }
}
```

### Stage 6: Parameter Calculation Integration

```cpp
// tests/integration/test_parameter_integration.cpp
class ParameterIntegrationTest : public IntegrationTestBase {};

TEST_F(ParameterIntegrationTest, StepParametersMatchLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        auto legacy_params = extract_bpstep_params(legacy_json);
        
        // Calculate parameters with our code
        Structure structure = load_and_process_pdb(pair.pdb_file);
        ParameterCalculator calculator;
        auto our_params = calculator.calculate_all_step_parameters(
            structure.base_pairs()
        );
        
        // Compare all parameters
        EXPECT_EQ(our_params.size(), legacy_params.size());
        
        for (size_t i = 0; i < our_params.size(); ++i) {
            const auto& our = our_params[i];
            const auto& legacy = legacy_params[i]["params"];
            
            EXPECT_NEAR(our.shift, legacy["Shift"], 0.001);
            EXPECT_NEAR(our.slide, legacy["Slide"], 0.001);
            EXPECT_NEAR(our.rise, legacy["Rise"], 0.001);
            EXPECT_NEAR(our.tilt, legacy["Tilt"], 0.001);
            EXPECT_NEAR(our.roll, legacy["Roll"], 0.001);
            EXPECT_NEAR(our.twist, legacy["Twist"], 0.001);
        }
    }
}
```

### Stage 7: Protocol Integration

```cpp
// tests/integration/test_protocol_integration.cpp
class ProtocolIntegrationTest : public IntegrationTestBase {};

TEST_F(ProtocolIntegrationTest, FindPairProtocolMatchesLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        
        // Run find_pair protocol
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        FindPairProtocol protocol;
        protocol.execute(structure);
        
        // Export to JSON
        JsonWriter writer;
        writer.init(pair.pdb_file);
        // ... record all data ...
        auto our_json = writer.get_json();
        
        // Compare all calculation types
        compare_all_calculations(legacy_json["calculations"], 
                                our_json["calculations"], 
                                0.001);
    }
}

TEST_F(ProtocolIntegrationTest, AnalyzeProtocolMatchesLegacy) {
    for (const auto& pair : pairs_) {
        // Load legacy JSON
        auto legacy_json = load_legacy_json(pair.json_file);
        
        // Run analyze protocol
        Structure structure = load_structure_with_pairs(pair.pdb_file);
        
        AnalyzeProtocol protocol;
        protocol.execute(structure);
        
        // Compare parameters
        compare_parameters(protocol.step_parameters(), 
                          extract_bpstep_params(legacy_json),
                          0.001);
    }
}
```

## Comparison Utilities

### JSON Comparison Functions

```cpp
// tests/integration/JsonComparison.hpp
namespace x3dna::test {

// Compare two JSON records
bool compare_json_record(const nlohmann::json& record1,
                         const nlohmann::json& record2,
                         double tolerance = 0.001);

// Compare arrays of calculations
void compare_calculations(const std::vector<nlohmann::json>& legacy,
                         const std::vector<nlohmann::json>& ours,
                         double tolerance = 0.001);

// Compare specific record types
void compare_frames(const ReferenceFrame& frame,
                   const nlohmann::json& legacy_frame,
                   double tolerance = 0.001);

void compare_matrices(const Matrix3D& matrix,
                     const nlohmann::json& legacy_matrix,
                     double tolerance = 0.001);

void compare_vectors(const Vector3D& vector,
                    const nlohmann::json& legacy_vector,
                    double tolerance = 0.001);

void compare_parameters(const BasePairStepParameters& params,
                       const nlohmann::json& legacy_params,
                       double tolerance = 0.001);

// Extract records by type
std::vector<nlohmann::json> extract_records_by_type(
    const nlohmann::json& json,
    const std::string& type
);

} // namespace x3dna::test
```

## Test Execution

### Running All Integration Tests

```bash
# Run all integration tests
cd build
ctest -R integration

# Run specific integration test
./tests/integration/test_frame_calculation_integration

# Run with verbose output
ctest -R integration -VV
```

### Running on Specific PDB Files

```cpp
// Filter tests to specific PDB files
TEST_F(FrameCalculationIntegrationTest, FramesMatchLegacy_100D) {
    auto pairs = TestDataDiscovery::discover_pairs("data/pdb", "data/json_legacy");
    auto it = std::find_if(pairs.begin(), pairs.end(),
        [](const auto& p) { return p.pdb_name == "100D"; });
    
    if (it != pairs.end()) {
        // Test only 100D
        test_single_pair(*it);
    }
}
```

## Test Data Management

### Current Test Data

- `data/pdb/100D.pdb` → `data/json_legacy/100D.json` ✅
- `data/pdb/157D.pdb` → `data/json_legacy/157D.json` ✅

### Adding More Test Data

As more JSON files are generated, they will be automatically discovered:

1. Generate JSON for a PDB file using original code
2. Place JSON in `data/json_legacy/`
3. Integration tests will automatically discover and test it

### Test Data Validation

```cpp
// Validate test data before running tests
TEST_F(IntegrationTestBase, ValidateTestData) {
    for (const auto& pair : pairs_) {
        // Verify PDB file exists
        ASSERT_TRUE(std::filesystem::exists(pair.pdb_file))
            << "PDB file missing: " << pair.pdb_file;
        
        // Verify JSON file exists
        ASSERT_TRUE(std::filesystem::exists(pair.json_file))
            << "JSON file missing: " << pair.json_file;
        
        // Verify JSON is valid
        ASSERT_NO_THROW({
            std::ifstream file(pair.json_file);
            nlohmann::json json;
            file >> json;
        }) << "Invalid JSON: " << pair.json_file;
    }
}
```

## Performance Considerations

### Parallel Execution

```cpp
// Run tests in parallel for multiple PDB files
TEST_F(IntegrationTestBase, ParallelExecution) {
    std::vector<std::future<void>> futures;
    
    for (const auto& pair : pairs_) {
        futures.push_back(std::async(std::launch::async, [&pair]() {
            test_single_pair(pair);
        }));
    }
    
    // Wait for all tests
    for (auto& future : futures) {
        future.wait();
    }
}
```

### Test Timeout

```cpp
// Set timeout for individual tests
TEST_F(IntegrationTestBase, FramesMatchLegacy_WithTimeout) {
    for (const auto& pair : pairs_) {
        // Set timeout (e.g., 30 seconds per PDB)
        auto start = std::chrono::steady_clock::now();
        
        test_single_pair(pair);
        
        auto duration = std::chrono::steady_clock::now() - start;
        EXPECT_LT(duration.count(), 30'000'000'000LL)  // 30 seconds
            << "Test took too long for " << pair.pdb_name;
    }
}
```

## Reporting

### Test Results Summary

```cpp
// Generate test summary
class IntegrationTestReporter {
public:
    struct TestResult {
        std::string pdb_name;
        bool passed;
        std::vector<std::string> failures;
        double execution_time;
    };
    
    void report_results(const std::vector<TestResult>& results);
    void generate_html_report(const std::vector<TestResult>& results);
};
```

### Failure Analysis

```cpp
// Detailed failure reporting
void report_failure(const std::string& pdb_name,
                   const std::string& test_name,
                   const std::string& failure_details) {
    std::cerr << "FAILURE: " << pdb_name << " - " << test_name << std::endl;
    std::cerr << "Details: " << failure_details << std::endl;
    
    // Log to file
    std::ofstream log("test_failures.log", std::ios::app);
    log << pdb_name << "," << test_name << "," << failure_details << std::endl;
}
```

## Integration with CI/CD

### Continuous Integration

```yaml
# .github/workflows/integration_tests.yml
name: Integration Tests

on: [push, pull_request]

jobs:
  integration-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Build
        run: |
          mkdir build && cd build
          cmake ..
          make
      - name: Run Integration Tests
        run: |
          cd build
          ctest -R integration -VV
      - name: Upload Test Results
        uses: actions/upload-artifact@v2
        with:
          name: test-results
          path: build/test_results/
```

## Best Practices

1. **Run Early**: Start integration tests as soon as I/O is implemented (Stage 3)
2. **Run Often**: Run integration tests after each stage
3. **Fix Immediately**: Don't let failures accumulate
4. **Document Failures**: Keep detailed failure logs
5. **Add More Data**: Generate JSON for more PDB files as needed

## Success Criteria

- [ ] All integration tests pass for all PDB/JSON pairs
- [ ] Tests run automatically on discovery
- [ ] All calculation types compared
- [ ] All numerical values match within tolerance
- [ ] Performance acceptable
- [ ] Tests run in CI/CD

---

*This framework ensures comprehensive validation against legacy JSON files at every stage.*

