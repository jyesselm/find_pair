# Stage 7: Protocol Implementations

## Objectives

Implement the protocol layer that orchestrates the algorithms to perform complete workflows (find_pair and analyze).

## Duration

**1 week**

## Dependencies

- ✅ Stage 0-6: All previous stages
- ✅ All algorithms must be implemented

## Tasks

### Task 7.1: Implement ProtocolBase
- [ ] Create `include/x3dna/protocols/ProtocolBase.hpp`
- [ ] Define protocol interface
- [ ] Implement configuration management
- [ ] Implement error handling
- [ ] Write unit tests

**Key Methods**:
```cpp
class ProtocolBase {
public:
    virtual ~ProtocolBase() = default;
    virtual void execute(Structure& structure) = 0;
    
    void set_config_manager(ConfigManager& config);
    ConfigManager& config() const;
    
protected:
    ConfigManager* config_ = nullptr;
};
```

**Deliverable**: ProtocolBase class

### Task 7.2: Implement FindPairProtocol
- [ ] Create `include/x3dna/protocols/FindPairProtocol.hpp`
- [ ] Implement find_pair workflow
- [ ] Orchestrate frame calculation
- [ ] Orchestrate base pair finding
- [ ] Orchestrate helix detection
- [ ] Handle options (single strand, all pairs, etc.)
- [ ] **Support `--legacy-mode` flag** (for exact legacy compatibility)
- [ ] Add JSON recording
- [ ] Write comprehensive unit tests

**Key Methods**:
```cpp
class FindPairProtocol : public ProtocolBase {
public:
    void execute(Structure& structure) override;
    
    void set_single_strand_mode(bool value);
    void set_find_all_pairs(bool value);
    void set_divide_helices(bool value);
    void set_legacy_mode(bool value);  // Enable legacy compatibility mode
    
    const std::vector<BasePair>& base_pairs() const;
    const std::vector<Helix>& helices() const;
    
private:
    void calculate_frames(Structure& structure);
    void find_pairs(Structure& structure);
    void detect_helices(Structure& structure);
    void reorder_pairs(Structure& structure);
    
    std::unique_ptr<BaseFrameCalculator> frame_calculator_;
    std::unique_ptr<BasePairFinder> pair_finder_;
    std::unique_ptr<HelixDetector> helix_detector_;
    bool legacy_mode_ = false;  // Legacy compatibility mode
};
```

**Legacy Mode Support**:
- When `legacy_mode_` is true, use legacy iteration order
- Use 1-based indexing where legacy does
- Match legacy algorithm behavior exactly
- Output JSON in exact legacy format
- See `docs/LEGACY_MODE_DESIGN.md` for details

**Deliverable**: Fully tested FindPairProtocol

### Task 7.3: Implement AnalyzeProtocol
- [ ] Create `include/x3dna/protocols/AnalyzeProtocol.hpp`
- [ ] Implement analyze workflow
- [ ] Recalculate frames
- [ ] Calculate parameters
- [ ] Handle options (torsions, simple params, etc.)
- [ ] Add JSON recording
- [ ] Write comprehensive unit tests

**Key Methods**:
```cpp
class AnalyzeProtocol : public ProtocolBase {
public:
    void execute(Structure& structure) override;
    
    void set_calculate_torsions(bool value);
    void set_simple_parameters(bool value);
    void set_circular_structure(bool value);
    
    const std::vector<BasePairStepParameters>& step_parameters() const;
    const std::vector<HelicalParameters>& helical_parameters() const;
    
private:
    void recalculate_frames(Structure& structure);
    void calculate_parameters(Structure& structure);
    
    std::unique_ptr<BaseFrameCalculator> frame_calculator_;
    std::unique_ptr<ParameterCalculator> param_calculator_;
};
```

**Deliverable**: Fully tested AnalyzeProtocol

### Task 7.4: Integration Testing
- [ ] Test complete find_pair workflow
- [ ] Test complete analyze workflow
- [ ] Test with real PDB files
- [ ] Compare outputs with legacy JSON

### Integration Tests (All PDB/JSON Pairs)
- [ ] Create `test_protocol_integration.cpp` using `IntegrationTestBase`
- [ ] Test FindPairProtocol on all discovered PDB/JSON pairs
- [ ] Test AnalyzeProtocol on all discovered PDB/JSON pairs
- [ ] Compare ALL calculation types with legacy JSON
- [ ] Verify complete end-to-end workflow matches legacy
- [ ] Test with different protocol options
- [ ] Verify JSON output structure matches legacy format
- [ ] Run: `./tests/integration/test_protocol_integration`

**Deliverable**: Integration tests passing

## Testing Plan

### Unit Tests

#### ProtocolBase Tests
- [ ] Test configuration management
- [ ] Test error handling
- [ ] Test interface

#### FindPairProtocol Tests
- [ ] Test complete workflow
- [ ] Test with different options
- [ ] Test frame calculation integration
- [ ] Test base pair finding integration
- [ ] Test helix detection integration
- [ ] Test JSON recording

#### AnalyzeProtocol Tests
- [ ] Test complete workflow
- [ ] Test with different options
- [ ] Test frame recalculation
- [ ] Test parameter calculation
- [ ] Test JSON recording

### Integration Tests
- [ ] Test find_pair → analyze pipeline
- [ ] Test with real PDB files
- [ ] Test end-to-end workflow
- [ ] Compare with original executables

### Regression Tests
- [ ] Run find_pair protocol
- [ ] Export to JSON
- [ ] Compare with legacy JSON
- [ ] Run analyze protocol
- [ ] Export to JSON
- [ ] Compare with legacy JSON

## Success Criteria

- [ ] All unit tests pass
- [ ] All integration tests pass
- [ ] Protocols produce correct results
- [ ] JSON output matches legacy format
- [ ] Performance acceptable
- [ ] Documentation complete

## Deliverables

1. ✅ ProtocolBase implemented
2. ✅ FindPairProtocol fully implemented and tested
3. ✅ AnalyzeProtocol fully implemented and tested
4. ✅ Comprehensive unit test suite
5. ✅ Integration tests
6. ✅ Regression tests
7. ✅ Documentation

## Files Created

```
include/x3dna/protocols/
├── ProtocolBase.hpp
├── FindPairProtocol.hpp
└── AnalyzeProtocol.hpp

src/x3dna/protocols/
├── ProtocolBase.cpp
├── FindPairProtocol.cpp
└── AnalyzeProtocol.cpp

tests/unit/protocols/
├── test_protocol_base.cpp
├── test_find_pair_protocol.cpp
└── test_analyze_protocol.cpp

tests/integration/
└── test_protocols.cpp
```

## Next Stage

After completing Stage 7, proceed to **Stage 8: Applications** (`STAGE_08_APPLICATIONS.md`)

---

*Estimated Completion: Week 13*

