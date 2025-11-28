# Implementation Guide - Protocols & Applications

**Quick Reference**: What to implement next and how

## 🎯 Current Status

✅ **Complete**: Stages 0-6 (all core algorithms)  
❌ **Missing**: Stages 7-8 (protocols and applications)  
⚠️ **Partial**: Helix detection

## 📋 What to Implement

### 1. Protocols (Stage 7) - **START HERE**

**Why**: Enables complete workflows using existing algorithms

**Files to Create**:
```
include/x3dna/protocols/
├── ProtocolBase.hpp
└── FindPairProtocol.hpp

src/x3dna/protocols/
├── ProtocolBase.cpp
└── FindPairProtocol.cpp
```

**Key Features**:
- Orchestrate frame calculation → pair finding → helix detection
- Support `--legacy-mode` flag for exact legacy compatibility
- JSON recording

**See**: `docs/modernization/STAGE_07_PROTOCOLS.md`

### 2. Helix Detection (Stage 5) - **MEDIUM PRIORITY**

**Why**: Needed for complete find_pair workflow

**Files to Create**:
```
include/x3dna/algorithms/HelixDetector.hpp
src/x3dna/algorithms/HelixDetector.cpp
```

**See**: `docs/modernization/STAGE_05_ALGORITHMS_2.md`

### 3. Applications (Stage 8) - **AFTER PROTOCOLS**

**Why**: Provides command-line interface

**Files to Create**:
```
apps/find_pair_app.cpp
apps/analyze_app.cpp
include/x3dna/apps/CommandLineParser.hpp
```

**Key Features**:
- Parse `--legacy-mode` flag
- Execute protocols
- Write output files

**See**: `docs/modernization/STAGE_08_APPLICATIONS.md`

## 🔧 Legacy Mode Support

**Requirement**: Add `--legacy-mode` flag that breaks some OOP for exact compatibility

**Design**: See `docs/LEGACY_MODE_DESIGN.md`

**Key Points**:
- Default: Modern mode (clean OOP)
- Legacy mode: Exact match with legacy (breaks some OOP)
- Opt-in: Only when flag is set
- Essential: For regression testing and comparison

**Implementation**:
```cpp
// In ConfigManager
bool legacy_mode() const { return legacy_mode_; }
void set_legacy_mode(bool value) { legacy_mode_ = value; }

// In Protocols
if (config_->legacy_mode()) {
    residues = get_residues_in_legacy_order(structure);
} else {
    residues = structure.all_residues();
}
```

## 📚 Documentation Created

1. **`MODERNIZATION_STATUS.md`** - Current implementation status
2. **`IMPLEMENTATION_ROADMAP.md`** - Detailed roadmap
3. **`NEXT_IMPLEMENTATION_STEPS.md`** - Step-by-step guide
4. **`docs/LEGACY_MODE_DESIGN.md`** - Legacy mode design
5. **`docs/modernization/LEGACY_MODE_REQUIREMENT.md`** - Requirement summary

## 🚀 Quick Start

### To Implement Protocols

1. **Read the design**:
   - `docs/modernization/STAGE_07_PROTOCOLS.md`
   - `docs/LEGACY_MODE_DESIGN.md`

2. **Start with ProtocolBase**:
   - Simple abstract base class
   - 1-2 days

3. **Implement FindPairProtocol**:
   - Use existing `BaseFrameCalculator`
   - Use existing `BasePairFinder`
   - Use existing `JsonWriter`
   - Add legacy mode support
   - 3-4 days

4. **Test**:
   - Test modern mode
   - Test legacy mode (compare with legacy JSON)

## 📊 Summary

**What's Done**:
- ✅ All algorithms (100% match rate)
- ✅ All data structures
- ✅ I/O layer
- ✅ Documentation organized

**What's Next**:
- ❌ Protocols (Stage 7)
- ❌ Applications (Stage 8)
- ❌ Helix Detection

**Estimated Time**: 2-3 weeks

**Priority**: Start with Protocols (Stage 7)

