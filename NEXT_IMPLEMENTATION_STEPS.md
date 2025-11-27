# Next Implementation Steps - Protocols & Applications

**Date**: Current  
**Status**: Ready to implement - all dependencies complete

## Summary

**Current State**: 
- ✅ Core algorithms: 100% complete
- ✅ 100% match rate with legacy code
- ❌ Protocols: 0% (not implemented)
- ❌ Applications: 0% (not implemented)
- ❌ Helix Detection: 0% (not implemented)

**Next Priority**: Implement Protocols (Stage 7) with `--legacy-mode` support

## Implementation Plan

### Phase 1: ProtocolBase (1-2 days)

**File**: `include/x3dna/protocols/ProtocolBase.hpp`

**Simple base class**:
```cpp
class ProtocolBase {
public:
    virtual ~ProtocolBase() = default;
    virtual void execute(Structure& structure) = 0;
    void set_config_manager(ConfigManager& config);
protected:
    ConfigManager* config_ = nullptr;
};
```

**Steps**:
1. Create header file
2. Implement base class
3. Add to CMakeLists.txt
4. Write basic tests

### Phase 2: FindPairProtocol (3-4 days)

**File**: `include/x3dna/protocols/FindPairProtocol.hpp`

**Orchestrates**:
- Frame calculation (using `BaseFrameCalculator`)
- Base pair finding (using `BasePairFinder`)
- Helix detection (using `HelixDetector` - when available)
- JSON recording (using `JsonWriter`)

**Key Implementation**:
```cpp
class FindPairProtocol : public ProtocolBase {
public:
    void execute(Structure& structure) override;
    void set_legacy_mode(bool value);
    
private:
    void calculate_frames(Structure& structure);
    void find_pairs(Structure& structure);
    void detect_helices(Structure& structure);
    
    bool legacy_mode_ = false;
    std::unique_ptr<BaseFrameCalculator> frame_calculator_;
    std::unique_ptr<BasePairFinder> pair_finder_;
};
```

**Legacy Mode Support**:
```cpp
void FindPairProtocol::find_pairs(Structure& structure) {
    std::vector<const Residue*> residues;
    if (legacy_mode_) {
        // Use legacy iteration order
        residues = get_residues_in_legacy_order(structure);
    } else {
        // Use modern order
        residues = structure.all_residues();
    }
    
    // ... pair finding logic
}
```

### Phase 3: AnalyzeProtocol (2-3 days)

**File**: `include/x3dna/protocols/AnalyzeProtocol.hpp`

**Orchestrates**:
- Frame recalculation
- Parameter calculation (using `ParameterCalculator`)
- JSON recording

### Phase 4: Applications (3-4 days)

**Files**: 
- `apps/find_pair_app.cpp`
- `apps/analyze_app.cpp`
- `include/x3dna/apps/CommandLineParser.hpp`

**Key Features**:
- Parse `--legacy-mode` flag
- Execute protocols
- Write output files

## Legacy Mode Integration

### ConfigManager Update

Add to `ConfigManager` (when implementing):
```cpp
bool legacy_mode() const { return legacy_mode_; }
void set_legacy_mode(bool value) { legacy_mode_ = value; }
```

### Protocol Integration

```cpp
void FindPairProtocol::execute(Structure& structure) {
    if (config_ && config_->legacy_mode()) {
        legacy_mode_ = true;
    }
    // ... rest of execution
}
```

### Command-Line Integration

```cpp
// In find_pair_app.cpp
bool legacy_mode = false;
for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "--legacy-mode") {
        legacy_mode = true;
    }
}

auto& config = ConfigManager::instance();
config.set_legacy_mode(legacy_mode);
```

## Testing Strategy

### Unit Tests

```cpp
TEST(FindPairProtocol, ModernMode) {
    ConfigManager::instance().set_legacy_mode(false);
    FindPairProtocol protocol;
    Structure structure = load_pdb("data/pdb/1H4S.pdb");
    protocol.execute(structure);
    // Verify modern behavior
}

TEST(FindPairProtocol, LegacyMode) {
    ConfigManager::instance().set_legacy_mode(true);
    FindPairProtocol protocol;
    Structure structure = load_pdb("data/pdb/1H4S.pdb");
    protocol.execute(structure);
    
    // Compare with legacy JSON - should match exactly
    auto legacy_json = load_json("data/json_legacy/find_bestpair_selection/1H4S.json");
    // ... verify exact match
}
```

## Files to Create

### Priority 1: Protocols

```
include/x3dna/protocols/
├── ProtocolBase.hpp          [NEW]
└── FindPairProtocol.hpp       [NEW]

src/x3dna/protocols/
├── ProtocolBase.cpp            [NEW]
└── FindPairProtocol.cpp        [NEW]
```

### Priority 2: Helix Detection

```
include/x3dna/algorithms/
└── HelixDetector.hpp           [NEW]

src/x3dna/algorithms/
└── HelixDetector.cpp           [NEW]
```

### Priority 3: Applications

```
include/x3dna/apps/
└── CommandLineParser.hpp       [NEW]

apps/
├── find_pair_app.cpp           [NEW]
└── analyze_app.cpp             [NEW]
```

## Quick Start

### To Begin Implementation

1. **Create ProtocolBase**:
   ```bash
   # Create directory
   mkdir -p include/x3dna/protocols src/x3dna/protocols
   
   # Create ProtocolBase.hpp
   # Simple abstract base class
   ```

2. **Create FindPairProtocol**:
   ```bash
   # Use existing algorithms:
   # - BaseFrameCalculator (already implemented)
   # - BasePairFinder (already implemented)
   # - JsonWriter (already implemented)
   ```

3. **Add Legacy Mode**:
   ```cpp
   // Check config.legacy_mode()
   // Use get_residues_in_legacy_order() when enabled
   ```

4. **Test**:
   ```bash
   # Test modern mode
   # Test legacy mode (compare with legacy JSON)
   ```

## Resources

- **Design**: `docs/modernization/STAGE_07_PROTOCOLS.md`
- **Legacy Mode**: `docs/LEGACY_MODE_DESIGN.md`
- **Status**: `MODERNIZATION_STATUS.md`
- **Roadmap**: `IMPLEMENTATION_ROADMAP.md`

## Estimated Timeline

- **ProtocolBase**: 1-2 days
- **FindPairProtocol**: 3-4 days
- **AnalyzeProtocol**: 2-3 days
- **Helix Detection**: 3-5 days
- **Applications**: 3-4 days

**Total**: 2-3 weeks

## Success Criteria

- [ ] Protocols implemented and tested
- [ ] `--legacy-mode` flag works correctly
- [ ] Legacy mode matches legacy JSON exactly
- [ ] Modern mode works correctly
- [ ] Applications provide command-line interface
- [ ] All tests pass

