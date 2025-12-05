# Stage 3: I/O Layer

## Objectives

Implement the I/O layer for reading PDB files, writing PDB files, and handling JSON input/output (both modern and legacy formats).

## Duration

**2 weeks**

## Dependencies

- ✅ Stage 0: Project Setup
- ✅ Stage 1: Geometry Classes
- ✅ Stage 2: Core Domain Objects

## Tasks

### Task 3.1: Implement PdbParser
- [ ] Create `include/x3dna/io/PdbParser.hpp`
- [ ] Implement PDB file format parsing
- [ ] Parse ATOM records
- [ ] Parse HETATM records (optional, configurable)
- [ ] Handle chain identification
- [ ] Handle residue numbering
- [ ] Handle alternate conformations
- [ ] Error handling and reporting
- [ ] Create Structure objects from PDB
- [ ] Write comprehensive unit tests

**Key Methods**:
```cpp
class PdbParser {
    Structure parse_file(const std::filesystem::path& path);
    Structure parse_stream(std::istream& stream);
    Structure parse_string(const std::string& content);
    
    void set_include_hetatm(bool value);
    void set_include_waters(bool value);
    
    class ParseError : public std::runtime_error {
        size_t line_number() const;
    };
};
```

**Deliverable**: Fully tested PdbParser

### Task 3.2: Implement PdbWriter
- [ ] Create `include/x3dna/io/PdbWriter.hpp`
- [ ] Implement Structure → PDB format conversion
- [ ] Write ATOM records
- [ ] Write HETATM records (if needed)
- [ ] Format coordinates correctly
- [ ] Handle chain/residue numbering
- [ ] Write to file/stream
- [ ] Write unit tests

**Key Methods**:
```cpp
class PdbWriter {
    void write_file(const Structure& structure, const std::filesystem::path& path);
    void write_stream(const Structure& structure, std::ostream& stream);
    std::string to_string(const Structure& structure);
};
```

**Deliverable**: Fully tested PdbWriter

### Task 3.3: Implement JsonReader
- [ ] Create `include/x3dna/io/JsonReader.hpp`
- [ ] Implement reading modern JSON format
- [ ] Implement reading legacy JSON format
- [ ] Parse Structure from JSON
- [ ] Parse individual records
- [ ] Error handling
- [ ] Write unit tests

**Key Methods**:
```cpp
class JsonReader {
    Structure read_structure(const std::filesystem::path& path);
    Structure read_structure(const nlohmann::json& json);
    Structure read_structure_legacy(const std::filesystem::path& path);
    
    // Individual record parsers
    std::vector<BasePair> read_base_pairs(const nlohmann::json& json);
    std::vector<ReferenceFrame> read_ref_frames(const nlohmann::json& json);
};
```

**Deliverable**: Fully tested JsonReader

### Task 3.4: Implement JsonWriter
- [ ] Create `include/x3dna/io/JsonWriter.hpp`
- [ ] Implement writing modern JSON format
- [ ] Implement writing legacy JSON format
- [ ] Match original json_writer.c functionality
- [ ] Support all calculation record types
- [ ] Pretty printing options
- [ ] Write unit tests

**Key Methods**:
```cpp
class JsonWriter {
    void init(const std::filesystem::path& pdb_file);
    void close();
    
    // Record functions (matching original API)
    void record_base_frame_calc(size_t residue_idx, char base_type, ...);
    void record_ls_fitting(size_t residue_idx, ...);
    void record_base_pair(const BasePair& pair);
    void record_bpstep_params(size_t bp_idx1, size_t bp_idx2, ...);
    void record_helical_params(size_t bp_idx1, size_t bp_idx2, ...);
    void record_pair_validation(...);
    void record_hbond_list(const BasePair& pair);
    void record_pdb_atoms(const Structure& structure);
    void record_all_ref_frames(const Structure& structure);
    // ... etc
    
    nlohmann::json get_json() const;
    void write_to_file(const std::filesystem::path& path) const;
};
```

**Critical**: Must match legacy JSON format exactly!

**Deliverable**: Fully tested JsonWriter

### Task 3.5: Implement Input File Parser (for analyze phase)
- [ ] Create `include/x3dna/io/InputFileParser.hpp`
- [ ] Parse .inp files (X3DNA input format)
- [ ] Extract PDB file path
- [ ] Extract base pair data
- [ ] Extract configuration options
- [ ] Create Structure with base pairs
- [ ] Write unit tests

**Key Methods**:
```cpp
class InputFileParser {
    struct InputData {
        std::filesystem::path pdb_file;
        std::vector<BasePair> base_pairs;
        // ... other data
    };
    
    InputData parse(const std::filesystem::path& input_file);
};
```

**Deliverable**: Input file parser

### Task 3.6: Integration Testing
- [ ] Test PDB → Structure → PDB round-trip
- [ ] Test PDB → Structure → JSON → Structure round-trip
- [ ] Test legacy JSON → Structure → legacy JSON round-trip
- [ ] Test with real PDB files
- [ ] Test with legacy JSON files

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### PdbParser Tests
- [ ] Test parsing simple PDB file
- [ ] Test parsing ATOM records
- [ ] Test parsing HETATM records (with/without flag)
- [ ] Test chain identification
- [ ] Test residue numbering
- [ ] Test coordinate parsing
- [ ] Test error handling (malformed files)
- [ ] Test with real PDB files (`data/pdb/100D.pdb`, etc.)
- [ ] Compare parsed Structure with expected structure

**Test Cases**:
```cpp
// Simple PDB parsing
PdbParser parser;
Structure structure = parser.parse_file("data/pdb/100D.pdb");
EXPECT_GT(structure.num_atoms(), 0);
EXPECT_GT(structure.num_residues(), 0);

// Verify specific atoms
auto chain = structure.get_chain('A');
auto residue = chain[0];
auto atom = residue.find_atom(" C1'");
EXPECT_TRUE(atom.has_value());
```

#### PdbWriter Tests
- [ ] Test writing Structure to PDB format
- [ ] Test coordinate formatting
- [ ] Test chain/residue numbering
- [ ] Test round-trip (parse → write → parse → compare)
- [ ] Compare output with original PDB (if available)

#### JsonReader Tests
- [ ] Test reading modern JSON format
- [ ] Test reading legacy JSON format
- [ ] Test parsing Structure from JSON
- [ ] Test parsing individual records
- [ ] Test error handling (invalid JSON)
- [ ] **Critical**: Load `data/json_legacy/100D.json` and verify

#### JsonWriter Tests
- [ ] Test writing modern JSON format
- [ ] Test writing legacy JSON format
- [ ] Test all record types
- [ ] Test pretty printing
- [ ] **Critical**: Compare output with legacy JSON files

**Test Cases**:
```cpp
// Write legacy JSON
Structure structure = load_structure();
JsonWriter writer;
writer.init("data/pdb/100D.pdb");
writer.record_pdb_atoms(structure);
auto json = writer.get_json();

// Compare with legacy
auto legacy_json = load_json("data/json_legacy/100D.json");
auto legacy_atoms = find_record(legacy_json, "pdb_atoms");
compare_json_records(json["calculations"][0], legacy_atoms);
```

### Integration Tests

#### Round-Trip Tests
- [ ] PDB → Structure → PDB → Structure (compare)
- [ ] PDB → Structure → JSON → Structure (compare)
- [ ] Legacy JSON → Structure → Legacy JSON (compare)
- [ ] Verify data integrity through round-trips

#### Integration Tests (All PDB/JSON Pairs)
- [ ] Create `test_io_integration.cpp` using `IntegrationTestBase`
- [ ] Test PDB → JSON conversion on all discovered pairs
- [ ] Compare `pdb_atoms` records with legacy JSON
- [ ] Test JSON reading for all pairs
- [ ] Verify all atom data matches (names, coordinates, etc.)
- [ ] Test round-trip: PDB → JSON → Structure → JSON (compare)
- [ ] Run: `./tests/integration/test_io_integration`

### Regression Tests

#### Compare with Legacy JSON
- [ ] Load legacy JSON file
- [ ] Parse all record types
- [ ] Recreate Structure
- [ ] Export to legacy JSON
- [ ] Compare field-by-field
- [ ] Verify numerical values match (within tolerance)

**Test Script**:
```cpp
TEST(JsonRegression, MatchLegacyFormat) {
    // Load legacy
    auto legacy = load_json("data/json_legacy/100D.json");
    
    // Parse
    Structure structure = JsonReader::read_structure_legacy(legacy);
    
    // Export
    JsonWriter writer;
    writer.init("data/pdb/100D.pdb");
    // ... record all data ...
    auto our_json = writer.get_json();
    
    // Compare
    compare_calculations(legacy["calculations"], our_json["calculations"]);
}
```

## Test Data

### PDB Files
- `data/pdb/100D.pdb`
- `data/pdb/157D.pdb`
- Create minimal test PDB files for unit tests

### Legacy JSON Files
- `data/json_legacy/100D.json`
- `data/json_legacy/100D_globals.json`
- `data/json_legacy/157D.json`
- `data/json_legacy/157D_globals.json`

## Success Criteria

- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Code coverage > 80% for I/O classes
- [ ] PDB parser correctly parses test files
- [ ] JSON writer matches legacy format exactly
- [ ] Round-trip tests pass (no data loss)
- [ ] Regression tests match legacy JSON within tolerance
- [ ] Error handling works correctly
- [ ] Performance acceptable (parsing speed)
- [ ] Documentation complete

## Code Quality Checks

- [ ] No compiler warnings
- [ ] Passes clang-tidy
- [ ] Follows code style
- [ ] Proper error handling
- [ ] RAII for file handles
- [ ] Exception safety

## Deliverables

1. ✅ PdbParser fully implemented and tested
2. ✅ PdbWriter fully implemented and tested
3. ✅ JsonReader fully implemented and tested
4. ✅ JsonWriter fully implemented and tested
5. ✅ InputFileParser implemented
6. ✅ Comprehensive unit test suite
7. ✅ Integration tests
8. ✅ Regression tests (JSON compatibility)
9. ✅ Documentation

## Files Created

```
include/x3dna/io/
├── PdbParser.hpp
├── PdbWriter.hpp
├── JsonReader.hpp
├── JsonWriter.hpp
└── InputFileParser.hpp

src/x3dna/io/
├── PdbParser.cpp
├── PdbWriter.cpp
├── JsonReader.cpp
├── JsonWriter.cpp
└── InputFileParser.cpp

tests/unit/io/
├── test_pdb_parser.cpp
├── test_pdb_writer.cpp
├── test_json_reader.cpp
├── test_json_writer.cpp
└── test_input_file_parser.cpp

tests/integration/
└── test_io_integration.cpp

tests/regression/
└── test_json_format.cpp
```

## Risks & Mitigation

| Risk | Impact | Mitigation |
|------|--------|------------|
| PDB format variations | High | Test with multiple PDB files, handle edge cases |
| JSON format mismatch | High | Compare field-by-field, test early |
| Performance issues | Medium | Profile, optimize parsing if needed |
| Memory issues with large files | Medium | Stream parsing where possible, test with large files |

## Validation Checklist

Before moving to Stage 4:
- [ ] All I/O classes compile and link
- [ ] All unit tests pass
- [ ] Integration tests pass
- [ ] Regression tests match legacy JSON
- [ ] Round-trip tests pass
- [ ] Can parse real PDB files
- [ ] Can read/write legacy JSON
- [ ] Code reviewed
- [ ] Documentation complete
- [ ] No memory leaks
- [ ] Performance acceptable

## Next Stage

After completing Stage 3, proceed to **Stage 4: Algorithms Part 1** (`STAGE_04_ALGORITHMS_1.md`)

---

*Estimated Completion: Week 6*

