# Stage 8: Application Executables

## Objectives

Implement the application executables (find_pair_app and analyze_app) that provide command-line interfaces to the protocols.

## Duration

**1 week**

## Dependencies

- ✅ Stage 0-7: All previous stages
- ✅ Protocols must be implemented

## Tasks

### Task 8.1: Implement Command-Line Parser
- [ ] Create `include/x3dna/apps/CommandLineParser.hpp`
- [ ] Parse find_pair options
- [ ] Parse analyze options
- [ ] **Parse `--legacy-mode` flag** (for exact legacy compatibility)
- [ ] Handle input/output files
- [ ] Validate arguments
- [ ] Write unit tests

**Key Methods**:
```cpp
class CommandLineParser {
    struct FindPairOptions {
        std::filesystem::path pdb_file;
        bool single_strand = false;
        bool find_all_pairs = false;
        bool divide_helices = false;
        bool legacy_mode = false;  // Enable legacy compatibility mode
        // ... other options
    };
    
    FindPairOptions parse_find_pair(int argc, char* argv[]);
    AnalyzeOptions parse_analyze(int argc, char* argv[]);
};
```

**Legacy Mode Flag**:
- `--legacy-mode`: Enable exact compatibility with legacy code
- Breaks some OOP principles for exact matching
- Essential for regression testing and comparison
- See `docs/LEGACY_MODE_DESIGN.md` for details

**Deliverable**: Command-line parser

### Task 8.2: Implement find_pair_app
- [ ] Create `apps/find_pair_app.cpp`
- [ ] Parse command-line arguments (including `--legacy-mode`)
- [ ] Load configuration
- [ ] **Set legacy mode in ConfigManager if flag is present**
- [ ] Parse PDB file
- [ ] Execute FindPairProtocol
- [ ] Write output files
- [ ] Handle errors gracefully
- [ ] Write integration tests

**Key Features**:
- Command-line interface matching original
- **`--legacy-mode` flag for exact legacy compatibility**
- Output file generation
- Error reporting
- Progress reporting (optional)

**Example Usage**:
```bash
# Modern mode (default)
./find_pair_app data/pdb/1H4S.pdb

# Legacy mode (exact match)
./find_pair_app --legacy-mode data/pdb/1H4S.pdb
```

**Deliverable**: find_pair_app executable

### Task 8.3: Implement analyze_app
- [ ] Create `apps/analyze_app.cpp`
- [ ] Parse command-line arguments
- [ ] Load input file (.inp format)
- [ ] Load PDB file
- [ ] Execute AnalyzeProtocol
- [ ] Write output files
- [ ] Handle errors gracefully
- [ ] Write integration tests

**Key Features**:
- Command-line interface matching original
- Input file parsing
- Output file generation
- Error reporting

**Deliverable**: analyze_app executable

### Task 8.4: Output Formatting
- [ ] Implement .inp file writer (for find_pair output)
- [ ] Implement parameter file writers
- [ ] Match original output formats
- [ ] Write tests

**Deliverable**: Output formatters

### Task 8.5: Integration Testing
- [ ] Test find_pair_app with real PDB files
- [ ] Test analyze_app with real input files
- [ ] Compare outputs with original executables
- [ ] Test error handling

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### CommandLineParser Tests
- [ ] Test option parsing
- [ ] Test file path handling
- [ ] Test validation
- [ ] Test error handling

#### Output Formatter Tests
- [ ] Test .inp file format
- [ ] Test parameter file formats
- [ ] Compare with original formats

### Integration Tests
- [ ] Test find_pair_app end-to-end
- [ ] Test analyze_app end-to-end
- [ ] Test with real PDB files
- [ ] Compare outputs with original

### Regression Tests
- [ ] Run find_pair_app on test PDB
- [ ] Compare .inp file with original
- [ ] Run analyze_app on .inp file
- [ ] Compare parameter files with original

## Success Criteria

- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Applications work correctly
- [ ] Outputs match original format
- [ ] Error handling works
- [ ] Documentation complete

## Deliverables

1. ✅ CommandLineParser implemented
2. ✅ find_pair_app executable
3. ✅ analyze_app executable
4. ✅ Output formatters
5. ✅ Comprehensive test suite
6. ✅ Documentation

## Files Created

```
apps/
├── find_pair_app.cpp
├── analyze_app.cpp
└── CommandLineParser.hpp

tests/integration/
└── test_applications.cpp
```

## Next Stage

After completing Stage 8, proceed to **Stage 9: Testing & Validation** (`STAGE_09_TESTING.md`)

---

*Estimated Completion: Week 14*

