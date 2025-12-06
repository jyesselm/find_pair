# Stage 0 & Stage 1 Refactoring Plan

**Scope**: PDB Parsing, Core Objects, Utilities, Configuration  
**Status**: Analysis & Planning  
**Priority**: Medium (Stage 2 is higher priority, but these should follow)

---

## Executive Summary

While **Stage 2 (frame calculation) is messy and needs immediate refactoring**, the earlier stages (PDB parsing, core objects, utilities) are actually in **better shape** but could still benefit from modern OOP improvements.

### Current State Assessment

| Component | LOC | Status | Refactoring Need |
|-----------|-----|--------|------------------|
| **PdbParser** | 766 | ✅ Good | LOW - Minor cleanup |
| **ResidueFactory** | 126 | ✅ Good | LOW - Already clean |
| **ModifiedNucleotideRegistry** | ~200 | ✅ Good | LOW - Already clean |
| **Structure** | 345 | ⚠️ OK | MEDIUM - JSON mixing |
| **JsonWriter** | 825 | ⚠️ OK | MEDIUM - Too many responsibilities |
| **ConfigurationManager** | ~150 | ⚠️ OK | MEDIUM - Needs expansion |

**Overall**: Stage 0/1 is **70% good**, Stage 2 is **30% good**  
**Recommendation**: Refactor Stage 2 first, then apply lessons learned to Stage 0/1

---

## Part 1: What's Already Good ✅

### 1.1: PdbParser - Clean and Well-Structured

**File**: `src/x3dna/io/pdb_parser.cpp` (766 LOC)

**Why it's good:**
```cpp
class PdbParser {
public:
    // Clean public API
    Structure parse_file(const std::filesystem::path& path);
    Structure parse_stream(std::istream& stream);
    Structure parse_string(const std::string& content);
    
    // Configuration
    void set_include_hetatm(bool value);
    void set_include_waters(bool value);
    
private:
    // Well-organized private helpers
    std::string parse_atom_name(const std::string& line);
    std::string parse_residue_name(const std::string& line);
    char parse_chain_id(const std::string& line);
    // ... etc
};
```

**Strengths:**
- ✅ Single Responsibility: Only parses PDB files
- ✅ Good separation: Public API vs private helpers
- ✅ Custom exception: `ParseError` with line numbers
- ✅ Well tested: 100% validation (3,602/3,602 PDBs)
- ✅ No JSON dependencies (pure parsing)

**Minor improvements needed:**
- Extract column definitions to constants
- Extract helper functions to separate utility class
- Add more comprehensive error messages

**Priority**: LOW (it works well)

---

### 1.2: ResidueFactory - Modern Factory Pattern ✅

**File**: `src/x3dna/core/residue_factory.cpp` (126 LOC)

**Why it's good:**
```cpp
class ResidueFactory {
public:
    static Residue create(
        const std::string& name,
        int sequence_number,
        char chain_id,
        char insertion_code,
        const std::vector<Atom>& atoms
    );
    
private:
    static char determine_one_letter_code(const std::string& name);
    static ResidueType determine_type(const std::string& name, char one_letter_code);
    static bool determine_is_purine(const std::string& name, ResidueType type);
};
```

**Strengths:**
- ✅ Static factory pattern (perfect for this use case)
- ✅ Uses ModifiedNucleotideRegistry (data-driven!)
- ✅ Clean separation of concerns
- ✅ Simple and focused

**No changes needed!** This is a model for good design.

**Priority**: NONE (perfect as-is)

---

### 1.3: ModifiedNucleotideRegistry - Data-Driven ✅

**Files**: `src/x3dna/core/modified_nucleotide_registry.cpp`

**Why it's good:**
```cpp
class ModifiedNucleotideRegistry {
public:
    static char get_one_letter_code(const std::string& residue_name);
    static std::optional<ResidueType> get_base_type(const std::string& residue_name);
    static bool is_purine(const std::string& residue_name);
    static bool contains(const std::string& residue_name);
    
    static void load_from_file(const std::filesystem::path& config_file);
    
private:
    // Singleton pattern
    static ModifiedNucleotideRegistry& instance();
    
    // Data loaded from JSON
    std::unordered_map<std::string, ModifiedNucleotideInfo> registry_;
};
```

**Strengths:**
- ✅ Singleton pattern (appropriate here)
- ✅ Data-driven: Loads from `modified_nucleotides.json`
- ✅ Clean static API
- ✅ Already has 85 nucleotides configured
- ✅ Easy to extend (edit JSON)

**Minor improvements:**
- Add validation on load (check for duplicates, invalid types)
- Add ability to reload configuration
- Better error messages

**Priority**: LOW (works great)

---

## Part 2: What Needs Improvement ⚠️

### 2.1: Structure - Mixed Responsibilities

**File**: `include/x3dna/core/structure.hpp` (345 LOC)

**Current problems:**
```cpp
class Structure {
public:
    // Good: Core domain model
    const std::vector<Chain>& chains() const;
    size_t num_residues() const;
    size_t num_atoms() const;
    
    // Bad: JSON writing (should be separate!)
    void write_atoms_json(const std::filesystem::path& output_dir) const;
    nlohmann::json to_json() const;
    void to_legacy_json(const std::filesystem::path& output_file) const;
    
    // Bad: Legacy order tracking (should be separate!)
    std::vector<const Residue*> residues_in_legacy_order() const;
    void set_legacy_indices(...);
};
```

**Problem**: Violates Single Responsibility Principle
- Domain model (chains/residues/atoms) ✅
- JSON serialization ❌
- Legacy ordering ❌

**Solution**: Extract to separate classes

```cpp
// AFTER: Clean separation

// 1. Pure domain model
class Structure {
public:
    const std::vector<Chain>& chains() const;
    void add_chain(const Chain& chain);
    size_t num_residues() const;
    size_t num_atoms() const;
    // NO JSON, NO LEGACY ORDER
};

// 2. JSON serialization (NEW)
class StructureJsonWriter {
public:
    explicit StructureJsonWriter(const Structure& structure);
    
    void write_atoms_json(const std::filesystem::path& output_dir) const;
    nlohmann::json to_json() const;
    void to_legacy_json(const std::filesystem::path& output_file) const;
    
private:
    const Structure& structure_;
};

// 3. Legacy order tracking (NEW or enhance existing)
class StructureLegacyOrder {
public:
    explicit StructureLegacyOrder(const Structure& structure);
    
    std::vector<const Residue*> residues_in_legacy_order() const;
    const Residue* get_residue_by_legacy_idx(int legacy_idx) const;
    int get_legacy_idx(const Residue& residue) const;
    
private:
    const Structure& structure_;
    // Cache legacy order
    mutable std::vector<const Residue*> cached_order_;
};
```

**Benefits:**
- Structure becomes pure domain model
- Easy to test each component
- Clear responsibilities
- Can serialize to different formats (JSON, XML, etc.)

**Priority**: MEDIUM

---

### 2.2: JsonWriter - Too Many Responsibilities

**File**: `src/x3dna/io/json_writer.cpp` (825 LOC!)

**Current problems:**
```cpp
class JsonWriter {
public:
    // Records for EVERY stage (17+ methods!)
    void record_residue_indices(...);
    void record_ls_fitting(...);
    void record_base_frame_calc(...);
    void record_frame_calc(...);
    void record_distance_checks(...);
    void record_hbond_list(...);
    void record_pair_validation(...);
    void record_base_pair(...);
    void record_bestpair_selection(...);
    void record_bpstep_params(...);
    void record_helical_params(...);
    // ... 6 more methods ...
    
    // File writing
    void write_split_files(...);
    void write_combined_file(...);
};
```

**Problem**: God object doing too much
- Knows about ALL stages (residues, frames, pairs, params)
- 825 lines - too long!
- Hard to test
- Violates Single Responsibility

**Solution**: Split by stage + composition

```cpp
// AFTER: Layered approach

// 1. Low-level JSON utilities
class JsonFormatter {
public:
    static nlohmann::json format_vector3d(const geometry::Vector3D& vec);
    static nlohmann::json format_matrix3d(const geometry::Matrix3D& mat);
    static nlohmann::json format_residue_info(const Residue& res);
    // ... formatting helpers
};

// 2. JSON accumulator (generic)
class JsonAccumulator {
public:
    void add_record(const std::string& type, const nlohmann::json& record);
    void write_split_files(const std::filesystem::path& output_dir) const;
    void write_combined_file(const std::filesystem::path& output_file) const;
    
private:
    std::map<std::string, std::vector<nlohmann::json>> records_;
};

// 3. Stage-specific recorders (created in Phase 4 of Stage 2 refactoring)
// These already planned in COMPREHENSIVE_STAGE2_REFACTORING_PLAN.md:
// - FrameJsonRecorder
// - PairJsonRecorder
// - ParameterJsonRecorder
// - etc.
```

**This refactoring aligns with Stage 2 refactoring!**

**Priority**: MEDIUM (will be addressed as part of Stage 2 refactoring)

---

### 2.3: ConfigurationManager - Needs Expansion

**File**: `src/x3dna/config/config_manager.cpp` (~150 LOC)

**Current state:**
```cpp
class ConfigManager {
public:
    static ConfigManager& instance();
    
    // Limited configuration
    const ParameterThresholds& thresholds() const;
    void set_x3dna_home(const std::filesystem::path& path);
    bool include_hetatm() const;
    bool include_waters() const;
    bool legacy_mode() const;
    
    void load_from_file(const std::filesystem::path& config_path);
};
```

**What's missing:**
- No frame calculation config
- No template assignment config
- No validation tolerances
- No logging configuration
- No performance tuning

**Solution**: Expand configuration (as planned in Stage 2 refactoring)

```cpp
// AFTER: Comprehensive configuration

struct FrameCalculationConfig {
    std::filesystem::path template_path = "data/templates";
    bool legacy_mode = false;
    bool auto_detect_rna = true;
    double rmsd_tolerance = 0.05;
    std::filesystem::path template_rules_file = 
        "resources/config/template_assignment_rules.json";
};

struct ParsingConfig {
    bool include_hetatm = false;
    bool include_waters = false;
    bool filter_by_occupancy = false;
    double min_occupancy = 0.5;
};

struct ValidationConfig {
    double coordinate_tolerance = 0.001;
    double angle_tolerance = 0.1;
    double distance_tolerance = 0.01;
    bool strict_mode = false;
};

struct LoggingConfig {
    std::string log_level = "INFO";  // DEBUG, INFO, WARNING, ERROR
    std::filesystem::path log_file = "";  // Empty = stdout
    bool log_to_file = false;
};

class ConfigurationManager {
public:
    static ConfigurationManager& instance();
    
    // Getters for all config sections
    const FrameCalculationConfig& frame_calculation_config() const;
    const ParsingConfig& parsing_config() const;
    const ValidationConfig& validation_config() const;
    const LoggingConfig& logging_config() const;
    const ParameterThresholds& thresholds() const;
    
    // Load from file
    void load_from_file(const std::filesystem::path& config_path);
    void load_defaults();
    
    // Update sections
    void set_frame_calculation_config(const FrameCalculationConfig& config);
    void set_parsing_config(const ParsingConfig& config);
    // ... etc
    
private:
    FrameCalculationConfig frame_calc_config_;
    ParsingConfig parsing_config_;
    ValidationConfig validation_config_;
    LoggingConfig logging_config_;
    ParameterThresholds thresholds_;
};
```

**Configuration file**: `resources/config/x3dna.json`

```json
{
  "frame_calculation": {
    "template_path": "data/templates",
    "legacy_mode": false,
    "auto_detect_rna": true,
    "rmsd_tolerance": 0.05,
    "template_rules_file": "resources/config/template_assignment_rules.json"
  },
  "parsing": {
    "include_hetatm": false,
    "include_waters": false,
    "filter_by_occupancy": false,
    "min_occupancy": 0.5
  },
  "validation": {
    "coordinate_tolerance": 0.001,
    "angle_tolerance": 0.1,
    "distance_tolerance": 0.01,
    "strict_mode": false
  },
  "logging": {
    "log_level": "INFO",
    "log_file": "",
    "log_to_file": false
  },
  "thresholds": {
    "max_distance": 10.0,
    "min_overlap": 0.5
  }
}
```

**Priority**: MEDIUM (will be addressed in Stage 2 refactoring)

---

## Part 3: Refactoring Strategy

### Approach: Bottom-Up Refactoring

**Phase 1**: Refactor Stage 2 (PRIORITY - see main refactoring plan)
- Extract helpers from BaseFrameCalculator
- Create service layer
- Create recorder layer
- Establish patterns

**Phase 2**: Apply lessons to Stage 0/1 (This document)
- Extract JSON from Structure
- Refactor JsonWriter using recorder pattern
- Expand ConfigurationManager
- Add utilities

**Why this order?**
1. Stage 2 is messier (more urgent)
2. Stage 2 refactoring establishes patterns
3. Stage 0/1 is 70% good already
4. Can apply Stage 2 lessons to Stage 0/1

---

## Part 4: Detailed Refactoring Plan

### 4.1: Extract JSON from Structure

#### Files to Create
```
include/x3dna/io/
└─ structure_json_writer.hpp (NEW)

src/x3dna/io/
└─ structure_json_writer.cpp (NEW)
```

#### Implementation
```cpp
class StructureJsonWriter {
public:
    explicit StructureJsonWriter(const core::Structure& structure)
        : structure_(structure) {}
    
    // JSON output
    void write_atoms_json(const std::filesystem::path& output_dir) const;
    nlohmann::json to_json() const;
    void to_legacy_json(const std::filesystem::path& output_file) const;
    
private:
    const core::Structure& structure_;
    
    nlohmann::json atom_to_json(const core::Atom& atom) const;
    nlohmann::json residue_to_json(const core::Residue& residue) const;
};
```

#### Migration
1. Create `StructureJsonWriter`
2. Move JSON methods from `Structure`
3. Update callers
4. Remove JSON methods from `Structure`

**Estimated time**: 1 day

---

### 4.2: Refactor JsonWriter

This will be addressed as part of **Stage 2 refactoring Phase 4** (creating recorder layer).

The plan in `COMPREHENSIVE_STAGE2_REFACTORING_PLAN.md` already covers this:
- Create `FrameJsonRecorder`
- Create `PairJsonRecorder`
- Create `ParameterJsonRecorder`
- Extract `JsonAccumulator` for low-level operations

**Estimated time**: 2.5 days (already planned in Stage 2)

---

### 4.3: Expand ConfigurationManager

This will be addressed as part of **Stage 2 refactoring Phase 7** (configuration management).

Already planned:
- Add `FrameCalculationConfig`
- Create `resources/config/x3dna.json`
- Expand ConfigurationManager with new sections

**Additional work for Stage 0/1:**
- Add `ParsingConfig`
- Add `ValidationConfig`
- Add `LoggingConfig`

**Estimated time**: 1 day (additional, on top of Stage 2's 1 day)

---

### 4.4: Add Utilities

#### JsonFormatter (NEW)

**File**: `include/x3dna/io/json_formatter.hpp`

```cpp
class JsonFormatter {
public:
    // Vector/Matrix formatting
    static nlohmann::json format_vector3d(const geometry::Vector3D& vec);
    static nlohmann::json format_matrix3d(const geometry::Matrix3D& mat);
    
    // Residue formatting
    static nlohmann::json format_residue_info(const core::Residue& res);
    static nlohmann::json format_residue_minimal(const core::Residue& res);
    
    // Atom formatting
    static nlohmann::json format_atom(const core::Atom& atom);
    static nlohmann::json format_atom_minimal(const core::Atom& atom);
    
    // Type formatting
    static std::string format_residue_type(core::ResidueType type);
    static std::string format_record_type(core::RecordType type);
};
```

**Purpose**: Centralize JSON formatting logic (DRY)

**Estimated time**: 0.5 days

---

#### PdbColumnDefinitions (NEW)

**File**: `include/x3dna/io/pdb_column_definitions.hpp`

```cpp
namespace pdb_format {

// PDB column ranges (1-based as per PDB spec)
constexpr int ATOM_NAME_START = 13;
constexpr int ATOM_NAME_END = 16;
constexpr int RESIDUE_NAME_START = 18;
constexpr int RESIDUE_NAME_END = 20;
constexpr int CHAIN_ID_COL = 22;
constexpr int RESIDUE_SEQ_START = 23;
constexpr int RESIDUE_SEQ_END = 26;
constexpr int X_START = 31;
constexpr int X_END = 38;
constexpr int Y_START = 39;
constexpr int Y_END = 46;
constexpr int Z_START = 47;
constexpr int Z_END = 54;
// ... etc

/**
 * @brief Extract substring using PDB column range (1-based)
 */
inline std::string extract_column(const std::string& line, int start, int end) {
    if (line.length() < static_cast<size_t>(end)) {
        return "";
    }
    // Convert to 0-based indexing
    return line.substr(start - 1, end - start + 1);
}

} // namespace pdb_format
```

**Purpose**: Document PDB format, avoid magic numbers

**Use in PdbParser:**
```cpp
// BEFORE:
std::string name = line.substr(12, 4);

// AFTER:
std::string name = pdb_format::extract_column(line, 
    pdb_format::ATOM_NAME_START, 
    pdb_format::ATOM_NAME_END);
```

**Estimated time**: 0.5 days

---

## Part 5: Priority Summary

### Immediate (Do with Stage 2)
1. ✅ Refactor Stage 2 (see main plan) - 17-20 days
2. ⏳ Extract JsonWriter recorders (part of Stage 2) - included above
3. ⏳ Expand ConfigurationManager (part of Stage 2) - included above

### Short-term (After Stage 2)
4. Extract JSON from Structure - 1 day
5. Add JsonFormatter utility - 0.5 days
6. Add PdbColumnDefinitions - 0.5 days
7. Expand ConfigurationManager further - 1 day

**Total additional time after Stage 2**: 3 days

### Long-term (Nice to have)
- Add comprehensive logging system
- Add performance profiling hooks
- Add memory pool for atom allocation
- Add parallel PDB parsing

---

## Part 6: Validation Strategy

Each refactoring must maintain:
- ✅ **Stage 1 validation**: 3,602/3,602 PDBs (atoms)
- ✅ **Stage 2 validation**: 3,602/3,602 PDBs (frames)

### After Each Change
```bash
# Build
cmake --build build -j8

# Run Stage 1 validation
pytest tests_python/ -k "test_atoms"

# Run Stage 2 validation
pytest tests_python/ -k "test_stage2"

# Expected: 100% pass rate
```

---

## Part 7: Architecture After All Refactoring

```
┌─────────────────────────────────────────────────────────┐
│ DOMAIN MODEL (Pure)                                     │
│  ├─ Atom                                                │
│  ├─ Residue                                             │
│  ├─ Chain                                               │
│  └─ Structure ← PURE (no JSON, no legacy order)         │
└─────────────────────────────────────────────────────────┘
                    │
┌───────────────────┼───────────────────────────────────────┐
│                   │                                       │
│ FACTORIES         │ REGISTRIES                            │
│  └─ ResidueFactory│  ├─ ModifiedNucleotideRegistry        │
│     (already ✅)  │  └─ TemplateAssignmentRegistry (NEW)  │
│                   │                                       │
└───────────────────┴───────────────────────────────────────┘
                    │
┌─────────────────────────────────────────────────────────┐
│ I/O LAYER                                               │
│  ├─ PdbParser (already ✅)                              │
│  ├─ StructureJsonWriter (NEW)                           │
│  ├─ JsonFormatter (NEW)                                 │
│  ├─ JsonAccumulator (NEW)                               │
│  └─ PdbColumnDefinitions (NEW)                          │
└─────────────────────────────────────────────────────────┘
                    │
┌─────────────────────────────────────────────────────────┐
│ ALGORITHMS                                              │
│  ├─ BaseFrameCalculator (refactored)                    │
│  ├─ RingAtomMatcher (NEW)                               │
│  ├─ ResidueTypeDetector (NEW)                           │
│  └─ ... (other algorithms)                              │
└─────────────────────────────────────────────────────────┘
                    │
┌─────────────────────────────────────────────────────────┐
│ SERVICES                                                │
│  ├─ FrameCalculationService (NEW)                       │
│  ├─ StructureLegacyOrder (NEW or enhanced)              │
│  └─ ... (other services)                                │
└─────────────────────────────────────────────────────────┘
                    │
┌─────────────────────────────────────────────────────────┐
│ RECORDERS                                               │
│  ├─ FrameJsonRecorder (NEW)                             │
│  ├─ PairJsonRecorder (NEW)                              │
│  ├─ ParameterJsonRecorder (NEW)                         │
│  └─ ... (other recorders)                               │
└─────────────────────────────────────────────────────────┘
                    │
┌─────────────────────────────────────────────────────────┐
│ CONFIGURATION                                           │
│  └─ ConfigurationManager (expanded)                     │
│     ├─ FrameCalculationConfig                           │
│     ├─ ParsingConfig                                    │
│     ├─ ValidationConfig                                 │
│     └─ LoggingConfig                                    │
└─────────────────────────────────────────────────────────┘
```

**Clean layers, clear responsibilities, easy to extend!**

---

## Part 8: Comparison with Stage 2

| Aspect | Stage 0/1 | Stage 2 |
|--------|-----------|---------|
| **Current quality** | 70% good | 30% good |
| **Urgency** | Medium | HIGH |
| **Refactoring effort** | 3 days | 17-20 days |
| **Validation** | ✅ 100% | ✅ 100% |
| **Main issues** | Minor mixing | Major mixing |
| **Priority** | 2nd | 1st |

**Recommendation**: Do Stage 2 first, then Stage 0/1

---

## Summary

### Stage 0/1 Assessment

**What's Good ✅:**
- PdbParser (766 LOC) - clean, well-structured
- ResidueFactory (126 LOC) - perfect factory pattern
- ModifiedNucleotideRegistry - data-driven, extensible
- 100% validation passing

**What Needs Work ⚠️:**
- Structure - mixed with JSON (need to extract)
- JsonWriter - too many responsibilities (825 LOC)
- ConfigurationManager - needs expansion
- Missing utilities (JsonFormatter, PdbColumnDefinitions)

**Refactoring Effort:**
- 3 additional days after Stage 2 refactoring
- Follows same patterns as Stage 2
- Maintains 100% validation

**Priority:**
1. Stage 2 first (17-20 days) - more urgent
2. Stage 0/1 second (3 days) - apply lessons learned

**Result:** Consistent, modern OOP architecture across ALL stages

---

**Next Steps:**
1. Proceed with Stage 2 refactoring (see main plan)
2. Use this document when refactoring Stage 0/1
3. Apply same patterns and principles throughout

