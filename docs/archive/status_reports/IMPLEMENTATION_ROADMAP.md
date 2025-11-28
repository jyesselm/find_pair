# Implementation Roadmap

**Date**: Current  
**Status**: Core algorithms complete, protocols and applications needed

## Current Status

### ✅ Completed (Stages 0-6)
- **100% match rate** with legacy code achieved
- All core algorithms implemented and tested
- All data structures in place
- I/O layer complete

### ❌ Missing (Stages 7-8)
- **Protocols** - High-level workflow orchestration
- **Applications** - Command-line executables
- **Helix Detection** - Helix detection and organization

## Next Implementation Steps

### Step 1: Implement Protocols (Stage 7) - **START HERE**

**Priority**: HIGH  
**Estimated Time**: 1 week  
**Dependencies**: All algorithms (✅ Complete)

#### Task 1.1: ProtocolBase
**File**: `include/x3dna/protocols/ProtocolBase.hpp`

```cpp
namespace x3dna::protocols {

class ProtocolBase {
public:
    virtual ~ProtocolBase() = default;
    virtual void execute(Structure& structure) = 0;
    
    void set_config_manager(ConfigManager& config);
    ConfigManager& config() const;
    
protected:
    ConfigManager* config_ = nullptr;
};

} // namespace x3dna::protocols
```

**Implementation Steps**:
1. Create header file
2. Implement base class
3. Add configuration management
4. Add error handling
5. Write unit tests

#### Task 1.2: FindPairProtocol
**File**: `include/x3dna/protocols/FindPairProtocol.hpp`

**Key Features**:
- Orchestrate frame calculation
- Orchestrate base pair finding
- Orchestrate helix detection (when available)
- Support `--legacy-mode` flag
- JSON recording

**Implementation Steps**:
1. Create header file
2. Implement `execute()` method
3. Add private helper methods:
   - `calculate_frames()`
   - `find_pairs()`
   - `detect_helices()`
   - `reorder_pairs()`
4. Add `--legacy-mode` support
5. Integrate with existing algorithms
6. Add JSON recording
7. Write comprehensive tests

**Legacy Mode Support**:
- Use `get_residues_in_legacy_order()` when legacy mode enabled
- Use 1-based indexing where legacy does
- Match legacy iteration order exactly
- See `docs/LEGACY_MODE_DESIGN.md` for details

#### Task 1.3: AnalyzeProtocol
**File**: `include/x3dna/protocols/AnalyzeProtocol.hpp`

**Key Features**:
- Recalculate frames
- Calculate step parameters
- Calculate helical parameters
- Support `--legacy-mode` flag
- JSON recording

**Implementation Steps**:
1. Create header file
2. Implement `execute()` method
3. Add private helper methods:
   - `recalculate_frames()`
   - `calculate_parameters()`
4. Add `--legacy-mode` support
5. Integrate with existing algorithms
6. Add JSON recording
7. Write comprehensive tests

### Step 2: Implement Helix Detection (Stage 5) - **MEDIUM PRIORITY**

**Priority**: MEDIUM  
**Estimated Time**: 3-5 days  
**Dependencies**: Base pair finding (✅ Complete)

**File**: `include/x3dna/algorithms/HelixDetector.hpp`

**Key Features**:
- Detect helices from base pairs
- Organize base pairs into helices
- Reorder pairs (5' to 3')
- Detect circular structures

**Implementation Steps**:
1. Create header file
2. Implement helix detection algorithm
3. Implement pair reordering
4. Integrate with FindPairProtocol
5. Write tests

### Step 3: Implement Applications (Stage 8) - **MEDIUM PRIORITY**

**Priority**: MEDIUM  
**Estimated Time**: 1 week  
**Dependencies**: Protocols (Step 1)

#### Task 3.1: CommandLineParser
**File**: `include/x3dna/apps/CommandLineParser.hpp`

**Key Features**:
- Parse `--legacy-mode` flag
- Parse find_pair options
- Parse analyze options
- Validate arguments

#### Task 3.2: find_pair_app
**File**: `apps/find_pair_app.cpp`

**Key Features**:
- Parse command-line arguments
- Load configuration
- Set legacy mode if flag present
- Execute FindPairProtocol
- Write output files

#### Task 3.3: analyze_app
**File**: `apps/analyze_app.cpp`

**Key Features**:
- Parse command-line arguments
- Load .inp file
- Execute AnalyzeProtocol
- Write parameter files

## Implementation Order

### Recommended Sequence

1. **ProtocolBase** (1-2 days)
   - Simple base class
   - Sets foundation

2. **FindPairProtocol** (3-4 days)
   - Most important protocol
   - Uses existing algorithms
   - High value

3. **Helix Detection** (3-5 days)
   - Needed for complete find_pair
   - Can be done in parallel with AnalyzeProtocol

4. **AnalyzeProtocol** (2-3 days)
   - Simpler than FindPairProtocol
   - Uses existing ParameterCalculator

5. **Applications** (3-4 days)
   - Command-line interfaces
   - Depends on protocols

**Total Estimated Time**: 2-3 weeks

## Legacy Mode Implementation

### Integration Points

1. **ConfigManager**: Add `legacy_mode` flag
2. **Protocols**: Check legacy mode and adjust behavior
3. **Algorithms**: Add legacy-mode versions where needed
4. **JSON Output**: Support legacy format in legacy mode

### Key Considerations

- **Separate Functions**: Use separate legacy vs modern functions for clarity
- **Testing**: Test both modes independently
- **Documentation**: Clearly document what breaks OOP in legacy mode
- **Default**: Modern mode is default, legacy mode is opt-in

## Testing Strategy

### Protocol Tests

```cpp
TEST(FindPairProtocol, ModernMode) {
    ConfigManager::instance().set_legacy_mode(false);
    FindPairProtocol protocol;
    protocol.execute(structure);
    // Test modern behavior
}

TEST(FindPairProtocol, LegacyMode) {
    ConfigManager::instance().set_legacy_mode(true);
    FindPairProtocol protocol;
    protocol.execute(structure);
    // Compare with legacy JSON - should match exactly
}
```

### Integration Tests

- Test complete workflows
- Test with real PDB files
- Compare outputs with legacy JSON
- Test both modern and legacy modes

## Success Criteria

### Protocols (Stage 7)
- [ ] ProtocolBase implemented
- [ ] FindPairProtocol implemented
- [ ] AnalyzeProtocol implemented
- [ ] `--legacy-mode` support added
- [ ] All unit tests pass
- [ ] Integration tests pass
- [ ] Legacy mode matches legacy JSON exactly

### Applications (Stage 8)
- [ ] CommandLineParser implemented
- [ ] find_pair_app executable works
- [ ] analyze_app executable works
- [ ] `--legacy-mode` flag parsed correctly
- [ ] Output files match original format

## Files to Create

### Protocols
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
```

### Helix Detection
```
include/x3dna/algorithms/
└── HelixDetector.hpp

src/x3dna/algorithms/
└── HelixDetector.cpp

tests/unit/algorithms/
└── test_helix_detector.cpp
```

### Applications
```
include/x3dna/apps/
└── CommandLineParser.hpp

apps/
├── find_pair_app.cpp
└── analyze_app.cpp

tests/integration/
└── test_applications.cpp
```

## Quick Start Guide

### To Implement Protocols

1. **Start with ProtocolBase**:
   ```cpp
   // Create include/x3dna/protocols/ProtocolBase.hpp
   // Simple abstract base class
   ```

2. **Implement FindPairProtocol**:
   ```cpp
   // Use existing:
   // - BaseFrameCalculator (calculate frames)
   // - BasePairFinder (find pairs)
   // - JsonWriter (record JSON)
   // - get_residues_in_legacy_order() (for legacy mode)
   ```

3. **Add Legacy Mode Support**:
   ```cpp
   if (config.legacy_mode()) {
       residues = get_residues_in_legacy_order(structure);
   } else {
       residues = structure.all_residues();
   }
   ```

4. **Test**:
   ```cpp
   // Test modern mode
   // Test legacy mode (compare with legacy JSON)
   ```

## Resources

- **Design Documents**:
  - `docs/modernization/STAGE_07_PROTOCOLS.md` - Protocol implementation guide
  - `docs/LEGACY_MODE_DESIGN.md` - Legacy mode design
  - `docs/MODERNIZATION_PLAN.md` - Overall plan

- **Existing Code**:
  - `include/x3dna/algorithms/` - All algorithms ready to use
  - `include/x3dna/core/` - All data structures ready
  - `include/x3dna/io/` - I/O layer ready

- **Reference**:
  - `org/src/find_pair.c` - Legacy find_pair implementation
  - `org/src/analyze.c` - Legacy analyze implementation

## Summary

**Ready to Implement**:
- ✅ All dependencies complete (algorithms, data structures, I/O)
- ✅ Design documents complete
- ✅ Legacy mode design documented
- ✅ Testing strategy defined

**Next Action**: Start implementing `ProtocolBase` and `FindPairProtocol`

**Estimated Completion**: 2-3 weeks for protocols and applications

