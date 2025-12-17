# Comprehensive C++ Refactoring Plan

**Branch**: `refactor/clean-code-principles`
**Scope**: ALL modern code (57 files, ~20,000 lines)
**Validation baseline**: 100% pass rate on core stages (100-PDB set), 99.1% on fast set

---

## Progress Summary

**Last Updated**: December 2024

| Phase | Status | Key Accomplishments |
|-------|--------|---------------------|
| 0 | ⬜ Not Started | Portability foundation |
| 1 | ✅ Mostly Done | Named constants, [[nodiscard]] added to most headers |
| 2 | ⬜ Not Started | Builder patterns for core data structures |
| 3 | 🔶 Partial | OverlapCalculator extracted, helix_organizer helpers added |
| 4 | ✅ Mostly Done | GEMMI parser refactoring, ResidueKey struct |
| 5 | ⬜ Not Started | Protocol configuration consolidation |
| 6 | ✅ Done | Early returns, named booleans applied |
| 7 | ⬜ Not Started | Final cleanup |

### Completed Commits (Recent)

| Commit | Description |
|--------|-------------|
| `d4fb211` | Add update_direction_count helper in helix_organizer.cpp |
| `2e6b306` | Refactor helix_organizer.cpp with FrameAlignment helper |
| `a157efb` | Fix GEMMI parser compatibility issues |
| `29c76f8` | Extract OverlapCalculator from base_pair_validator.cpp |
| `42be9d4` | Refactor parsers: add ResidueKey struct, remove deprecated methods |
| `fb2e91c` | Remove ResidueFactory and RingAtomRegistry classes |
| `b5a9d2d` | Consolidate modified nucleotide lookups to ModifiedNucleotideRegistry |
| `8c2c434` | Refactor PDB/CIF parsers to use GEMMI library |
| `c262bc8` | Add named constants for magic numbers |
| `b4c63bf` | Flatten nesting with early returns and named booleans |

### File Size Improvements

| File | Before | After | Change |
|------|--------|-------|--------|
| `base_pair_validator.cpp` | 911 lines | 496 lines | -46% |
| `pdb_parser.cpp` | 751 lines | 340 lines | -55% |

---

## Key Goals

1. **Portability**: Easy to drop into another project's `external/` directory
2. **Clean Code**: Follow all principles below
3. **Testability**: Unit tests for all extracted classes
4. **Encapsulation**: Private data by default

## Guiding Principles

- Functions should do one thing
- Limit function parameters to 3-4 max
- Max 3 levels of indentation
- Single Responsibility: one class, one reason to change
- Prefer composition over inheritance
- Use `const` everywhere you can
- Use `enum class` not `enum`
- Use `[[nodiscard]]` for return values that matter
- **Strong Encapsulation** - data is PRIVATE by default
- Extract complex conditionals into named booleans or functions
- Reduce getters/setters but keep data private (use builder patterns, factory methods)
- **Avoid `std::optional` for required members** - If data is always expected to exist after construction, don't make it optional. Optional types should only be used for truly optional data.

---

## Codebase Overview

| Category | Files | Lines | Key Issues |
|----------|-------|-------|------------|
| Core Data | 6 headers, 4 src | 3,200 | 70+ getters/setters, duplicate JSON |
| Geometry | 3 headers, 2 src | 800 | Missing `[[nodiscard]]` |
| Algorithms | 17 headers, 12 src | 7,500 | Monster functions, deep nesting |
| I/O | 6 headers, 5 src | 3,200 | Magic numbers, complex parsing |
| Protocols | 3 headers, 2 src | 900 | Setter proliferation |
| Config | 2 headers, 1 src | 300 | No validation |
| Apps/Tools | 5 files | 600 | Long main(), hardcoded paths |
| **TOTAL** | **57 files** | **~20,000** | |

---

## Phase 0: Portability Foundation (Critical for Embedding)

Make the library easy to embed in `external/x3dna/` of any project.

### Step 0.1: Fix Resource Path Resolution
**Problem**: Resources (templates, config) require specific file paths or environment variables
**Current**: Hardcoded paths + X3DNA_SOURCE_DIR macro + getenv() fallbacks

**Solution**: Configurable resource locator
```cpp
// include/x3dna/config/resource_locator.hpp
namespace x3dna::config {

class ResourceLocator {
public:
    // Set once at library initialization
    static void set_resources_path(const std::filesystem::path& path);

    // Get paths for specific resource types
    [[nodiscard]] static std::filesystem::path templates_dir();
    [[nodiscard]] static std::filesystem::path config_dir();
    [[nodiscard]] static std::filesystem::path template_file(const std::string& name);

    // Check if resources are available
    [[nodiscard]] static bool is_initialized();

private:
    static std::filesystem::path resources_path_;
};

} // namespace x3dna::config
```

**Files to update**:
- `src/x3dna/algorithms/standard_base_templates.cpp`
- `src/x3dna/core/modified_nucleotide_registry.cpp`
- `src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.cpp`
- `src/x3dna/config/config_manager.cpp`

**Unit test**: `tests/unit/config/resource_locator_test.cpp`

### Step 0.2: Remove X3DNA_SOURCE_DIR Compile-Time Macro
**Problem**: `X3DNA_SOURCE_DIR="${CMAKE_SOURCE_DIR}"` bakes build path into binary
**File**: `CMakeLists.txt` line 82

**Solution**: Remove macro, use runtime ResourceLocator instead
```cmake
# REMOVE this:
target_compile_definitions(x3dna PRIVATE X3DNA_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
```

### Step 0.3: Add CMake Package Export
**Problem**: Can't use `find_package(x3dna)` or proper `add_subdirectory()`
**Current**: Install targets commented out

**Create**:
```
cmake/
├── x3dnaConfig.cmake.in
├── x3dnaConfigVersion.cmake.in
└── x3dnaTargets.cmake (generated)
```

**CMakeLists.txt additions**:
```cmake
# Export library with namespace
install(TARGETS x3dna
    EXPORT x3dnaTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

install(EXPORT x3dnaTargets
    FILE x3dnaTargets.cmake
    NAMESPACE x3dna::
    DESTINATION lib/cmake/x3dna
)

# Support add_subdirectory usage
add_library(x3dna::x3dna ALIAS x3dna)
```

**Unit test**: Test CMake integration in `tests/cmake/`

### Step 0.4: Document External Usage
**Create**: `docs/EMBEDDING.md`
```markdown
# Using x3dna in Your Project

## Option 1: add_subdirectory (Recommended)
```cmake
add_subdirectory(external/x3dna)
target_link_libraries(myapp PRIVATE x3dna::x3dna)
```

```cpp
// In your code - initialize once at startup
x3dna::config::ResourceLocator::set_resources_path("external/x3dna/resources");
```

## Option 2: find_package
```cmake
find_package(x3dna REQUIRED)
target_link_libraries(myapp PRIVATE x3dna::x3dna)
```
```

### Step 0.5: Create Library Initialization API
**Problem**: Users don't know what setup is required
**Solution**: Single initialization function

```cpp
// include/x3dna/x3dna.hpp (add to main header)
namespace x3dna {

struct InitOptions {
    std::filesystem::path resources_path;
    bool enable_logging = false;
    // Future options...
};

// Call once before using library
void initialize(const InitOptions& options);

// Check if initialized
[[nodiscard]] bool is_initialized();

} // namespace x3dna
```

**Unit test**: `tests/unit/initialization_test.cpp`

---

## Phase 1: Global Modernization (Low Risk)

Apply across ALL files. Each step is one commit.

### Step 1.1: Add `[[nodiscard]]` Everywhere ✅ MOSTLY DONE
**Files**: All 50 headers
**Scope**: ~60 methods that return values that shouldn't be ignored
**Status**: Most geometry and algorithm headers already have [[nodiscard]] on key methods.
**Key targets**:
- `geometry/vector3d.hpp`: ✅ All methods annotated
- `geometry/matrix3d.hpp`: ✅ All methods annotated
- `core/residue.hpp`: ✅ Most methods annotated
- `algorithms/*.hpp`: ✅ Most methods annotated (15+ headers have [[nodiscard]])

### Step 1.2: Const Correctness Pass
**Files**: All source files
**Scope**: Parameters, member functions, local variables
**Status**: Ongoing - applied during other refactoring
**Focus areas**:
- Non-const references returned (e.g., `residue.hpp:77 atoms()`)
- Mutable members in const-marked classes

### Step 1.3: Extract Magic Numbers to Named Constants ✅ DONE
**Files**:
- `base_pair_validator.cpp`: `BOND_DISTANCE=2.0`, `XBIG=1e18`, `GAMUT=5e8`
- `pdb_parser.cpp`: Column indices (12, 16, 21, 26, 55-60, 77-78)
- `helix_organizer.cpp`: Threshold values
**Created**: `include/x3dna/core/constants.hpp`
**Commit**: `c262bc8`

### Step 1.4: Replace Debug getenv() with Proper Logging
**Files**: `helix_organizer.cpp` (5+ locations)
**Status**: Uses ConfigManager for debug settings (centralized)

---

## Phase 2: Core Data Structures (Medium Risk)

**Principle**: Keep data PRIVATE. Use Builder pattern for construction, reduce setters to minimum.

### Step 2.1: Refactor `atom.hpp` - Builder Pattern + Immutable Design
**Current**: 32+ getter/setter pairs (373 lines)
**Target**: Immutable after construction, Builder for complex creation

```cpp
// BEFORE (70+ lines of getters/setters)
class Atom {
    std::string name() const { return name_; }
    void set_name(const std::string& n) { name_ = n; }
    // ... 30 more pairs
};

// AFTER - Immutable with Builder
class Atom {
public:
    class Builder;  // Forward declare

    // Construction via builder only
    static Builder create(std::string name, geometry::Vector3D position);

    // Const accessors (keep minimal set needed by external code)
    [[nodiscard]] const std::string& name() const { return name_; }
    [[nodiscard]] const geometry::Vector3D& position() const { return position_; }
    [[nodiscard]] char chain_id() const { return chain_id_; }
    [[nodiscard]] int residue_seq() const { return residue_seq_; }

    // Computed/derived values
    [[nodiscard]] double distance_to(const Atom& other) const;
    [[nodiscard]] bool is_ring_atom() const;
    [[nodiscard]] bool is_backbone_atom() const;
    [[nodiscard]] std::string full_id() const;  // "A.GLY.42.CA"

private:
    friend class Builder;
    Atom() = default;  // Only Builder can construct

    std::string name_;
    geometry::Vector3D position_;
    std::string residue_name_;
    char chain_id_ = ' ';
    int residue_seq_ = 0;
    char alt_loc_ = ' ';
    char insertion_ = ' ';
    double occupancy_ = 1.0;
    double b_factor_ = 0.0;
    std::string element_;
    size_t legacy_atom_idx_ = 0;
    size_t legacy_residue_idx_ = 0;
};

class Atom::Builder {
public:
    Builder(std::string name, geometry::Vector3D position);

    Builder& residue_name(std::string name);
    Builder& chain_id(char id);
    Builder& residue_seq(int seq);
    Builder& alt_loc(char loc);
    Builder& insertion(char ins);
    Builder& occupancy(double occ);
    Builder& b_factor(double bf);
    Builder& element(std::string elem);
    Builder& legacy_indices(size_t atom_idx, size_t residue_idx);

    [[nodiscard]] Atom build() const;

private:
    Atom atom_;
};

// Usage:
auto atom = Atom::create("CA", pos)
    .residue_name("ALA")
    .chain_id('A')
    .residue_seq(42)
    .build();
```

**Unit test**: `tests/unit/core/atom_test.cpp`

### Step 2.2: Refactor `residue.hpp` - Same Pattern
**Current**: 18+ getter/setter pairs (356 lines)
**Target**: Builder pattern, immutable after construction

```cpp
class Residue {
public:
    class Builder;

    // Minimal accessors for what's actually needed
    [[nodiscard]] const std::string& name() const;
    [[nodiscard]] char chain_id() const;
    [[nodiscard]] int seq_num() const;
    [[nodiscard]] char one_letter_code() const;

    // Atom access (read-only)
    [[nodiscard]] const std::vector<Atom>& atoms() const;
    [[nodiscard]] std::optional<std::reference_wrapper<const Atom>> find_atom(const std::string& name) const;

    // Computed values
    [[nodiscard]] bool is_nucleotide() const;
    [[nodiscard]] std::string ry_classification() const;
    [[nodiscard]] std::vector<std::reference_wrapper<const Atom>> ring_atoms() const;

    // Frame (optional, set by frame calculator)
    [[nodiscard]] std::optional<ReferenceFrame> reference_frame() const;
    void set_reference_frame(ReferenceFrame frame);  // Only setter needed

private:
    friend class Builder;
    // ... private data
};
```

**Unit test**: `tests/unit/core/residue_test.cpp`

### Step 2.3: Refactor `base_pair.hpp` - Split Responsibilities + Encapsulation
**Current**: 469 lines (LARGEST HEADER), 20+ getters, tangled JSON
**Issues**:
- `hydrogen_bond` struct embedded
- Frame logic mixed in (`get_step_frame()`)
- Dual JSON formats duplicated
- **Frames are `std::optional` but always required** - Code calls `.value()` assuming frames exist, but API suggests they're optional. This is confusing and error-prone.

**Extract to new files**:
1. `core/hydrogen_bond.hpp` - Immutable H-bond class
2. `core/base_pair.hpp` - Core class with Builder
3. `io/base_pair_serializer.hpp` - JSON serialization (both formats)

**Key change: Make frames REQUIRED**:
- Current: `std::optional<ReferenceFrame> frame1_, frame2_`
- Problem: Every call site does `frame1_.value()` or `if (frame1_.has_value())` followed by `.value()`
- Precondition (base_pair_finder.hpp:59): "residues must have frames calculated"
- Solution: Pass frames to constructor, store as `ReferenceFrame` not `std::optional<ReferenceFrame>`

```cpp
// core/hydrogen_bond.hpp
class HydrogenBond {
public:
    HydrogenBond(std::string donor_atom, std::string acceptor_atom,
                 double distance, char type);

    [[nodiscard]] const std::string& donor_atom() const;
    [[nodiscard]] const std::string& acceptor_atom() const;
    [[nodiscard]] double distance() const;
    [[nodiscard]] char type() const;  // '-' for standard, '*' otherwise

private:
    std::string donor_atom_;
    std::string acceptor_atom_;
    double distance_;
    char type_;
};

// core/base_pair.hpp
class BasePair {
public:
    // Constructor requires frames - they're not optional!
    BasePair(size_t idx1, size_t idx2,
             const ReferenceFrame& frame1, const ReferenceFrame& frame2,
             BasePairType type = BasePairType::UNKNOWN);

    // Core identity (immutable)
    [[nodiscard]] size_t residue_idx1() const;
    [[nodiscard]] size_t residue_idx2() const;
    [[nodiscard]] BasePairType type() const;

    // Frames - REQUIRED, not optional
    [[nodiscard]] const ReferenceFrame& frame1() const { return frame1_; }
    [[nodiscard]] const ReferenceFrame& frame2() const { return frame2_; }

    // Validation metrics (computed from frames)
    [[nodiscard]] double origin_distance() const;
    [[nodiscard]] double plane_angle() const;
    [[nodiscard]] double direction_dot_product() const;

    // H-bonds (can be set after construction)
    [[nodiscard]] const std::vector<HydrogenBond>& hydrogen_bonds() const;
    void set_hydrogen_bonds(std::vector<HydrogenBond> hbonds);

private:
    size_t residue_idx1_;
    size_t residue_idx2_;
    ReferenceFrame frame1_;  // NOT optional - always required
    ReferenceFrame frame2_;  // NOT optional - always required
    BasePairType type_;
    std::vector<HydrogenBond> hbonds_;
};
```

**Unit test**: `tests/unit/core/base_pair_test.cpp`, `tests/unit/core/hydrogen_bond_test.cpp`

### Step 2.4: Consolidate JSON Conversion Logic
**Problem**: Legacy vs modern JSON conversion duplicated in:
- `atom.hpp` (2 methods)
- `residue.hpp` (2 methods)
- `base_pair.hpp` (2 methods)
- `structure.hpp` (2 methods)
- `parameters.hpp` (2 methods)

**Solution**: Create `io/json_serializer.hpp` with template pattern:
```cpp
template<typename T>
nlohmann::json to_json(const T& obj, JsonFormat format);

template<typename T>
T from_json(const nlohmann::json& j, JsonFormat format);
```

---

## Phase 3: Algorithm Classes - Break Up Monster Functions

### Step 3.1: Split `base_pair_finder.cpp` (821 lines)
**Problem**: `find_best_pairs()` is 300+ lines with 5-6 nesting levels
**Status**: Not started - file already has good internal structure with helper methods

**Extract to separate classes**:
```
algorithms/pairing/
├── pair_candidate_generator.hpp/cpp (~150 lines)
│   └── Phase 1: Generate all valid candidates
├── mutual_best_matcher.hpp/cpp (~150 lines)
│   └── Greedy mutual selection algorithm
├── quality_scorer.hpp/cpp (~100 lines) [exists partially]
│   └── Score calculation and adjustment
└── base_pair_finder.hpp/cpp (~200 lines)
    └── Orchestration via composition
```

**Unit tests**:
- `tests/unit/algorithms/pairing/pair_candidate_generator_test.cpp`
- `tests/unit/algorithms/pairing/mutual_best_matcher_test.cpp`
- `tests/unit/algorithms/pairing/quality_scorer_test.cpp`

### Step 3.2: Split `base_pair_validator.cpp` ✅ PARTIAL (911→496 lines)
**Problem**: Ring atom logic duplicated, polygon code embedded, deep nesting
**Status**: OverlapCalculator extracted (commit `29c76f8`)

**Extracted**:
```
algorithms/validation/
├── overlap_calculator.hpp/cpp ✅ DONE (~250 lines)
│   └── All polygon intersection code
```

**Remaining**:
```
├── distance_validator.hpp/cpp (~100 lines)
├── angle_validator.hpp/cpp (~50 lines)
├── ring_atom_finder.hpp/cpp (~150 lines)
│   └── Deduplicated ring atom matching
└── base_pair_validator.hpp/cpp (~200 lines)
    └── Composition of above validators
```

**Unit tests**:
- `tests/unit/algorithms/validation/distance_validator_test.cpp`
- `tests/unit/algorithms/validation/angle_validator_test.cpp`
- `tests/unit/algorithms/validation/overlap_calculator_test.cpp`
- `tests/unit/algorithms/validation/ring_atom_finder_test.cpp`

### Step 3.3: Split `helix_organizer.cpp` 🔶 PARTIAL (1,277 lines)
**Problem**: God class with 5+ responsibilities
**Status**: Helper structs/functions added inline:
- `FrameAlignment` struct with `is_aligned()` and `angle_sum()` (commit `2e6b306`)
- `compute_frame_alignment()` helper function
- `update_direction_count()` helper function (commit `d4fb211`)
- Early returns and named booleans applied (commits `22c9cf7`, `31e3ae4`)

**Remaining extraction**:
```
algorithms/helix/
├── backbone_linkage_checker.hpp/cpp (~80 lines)
│   └── O3'-P linkage detection
├── pair_geometry_helper.hpp/cpp (~100 lines)
│   └── Origin/z-axis calculations
├── helix_context_calculator.hpp/cpp (~200 lines)
│   └── Neighbor detection, endpoint finding
├── strand_direction_checker.hpp/cpp (~300 lines)
│   └── All wc_bporien, check_* functions
├── five_to_three_orderer.hpp/cpp (~200 lines)
│   └── Main ordering algorithm
└── helix_organizer.hpp/cpp (~150 lines)
    └── Facade composing all above
```

**Unit tests**:
- `tests/unit/algorithms/helix/backbone_linkage_checker_test.cpp`
- `tests/unit/algorithms/helix/pair_geometry_helper_test.cpp`
- `tests/unit/algorithms/helix/helix_context_calculator_test.cpp`
- `tests/unit/algorithms/helix/strand_direction_checker_test.cpp`
- `tests/unit/algorithms/helix/five_to_three_orderer_test.cpp`

### Step 3.4: Reduce Parameter Counts
**Problem**: 15+ functions with 5+ parameters

**Solution**: Group into parameter structs:
```cpp
// BEFORE
bool wc_bporien(const BasePair& pair_m, const BasePair& pair_n,
                bool swap_m, bool swap_n, const BackboneData& backbone);

// AFTER
struct PairComparisonContext {
    const BasePair& pair_m;
    const BasePair& pair_n;
    bool swap_m;
    bool swap_n;
    const BackboneData& backbone;
};

[[nodiscard]] bool wc_bporien(const PairComparisonContext& ctx);
```

---

## Phase 4: I/O Classes

### Step 4.1: Refactor `pdb_parser.cpp` ✅ DONE (751→340 lines)
**Problems**:
- Magic column indices scattered
- Tuple types as parameters (should be structs)
- Deep nesting in atom processing

**Completed**:
- Refactored to use GEMMI library for parsing (commit `8c2c434`)
- Added `ResidueKey` struct to replace tuples (commit `42be9d4`)
- Removed deprecated methods and ResidueFactory (commit `fb2e91c`)
- Fixed GEMMI parser compatibility (commit `a157efb`)

**ResidueKey struct** (now in `io/pdb_parser.hpp`):
```cpp
struct ResidueKey {
    std::string name;
    char chain_id;
    int seq_num;
    char insertion;

    bool operator==(const ResidueKey& other) const;
    bool operator<(const ResidueKey& other) const;
};
```

### Step 4.2: Refactor `json_writer.cpp` (891 lines)
**Problems**:
- `write_split_files()` is 100+ lines
- Deep nesting in type-to-directory mapping
- Repeated JSON construction patterns

**Extract**:
```
json_writer.cpp (891 lines)
├── io/json_directory_manager.cpp (~100 lines)
│   └── Directory creation, file path management
├── io/json_record_builder.cpp (~200 lines)
│   └── Repeated JSON construction patterns
└── json_writer.cpp (~400 lines)
    └── Core write logic
```

### Step 4.3: Extract Format Helpers from `input_file_writer.cpp` (507 lines)
**Problem**: Repeated `std::setw`, `std::fixed`, `std::setprecision`

**Create**: `io/format_helpers.hpp`
```cpp
namespace format {
    std::string fixed_width(double val, int width, int precision);
    std::string padded_int(int val, int width);
}
```

---

## Phase 5: Protocol and App Classes

### Step 5.1: Consolidate Protocol Configuration
**Problem**: 10+ setter methods in each protocol class

**Solution**: Builder pattern or config struct
```cpp
// BEFORE (in find_pair_protocol.hpp)
void set_legacy_mode(bool);
void set_split_output(bool);
void set_output_dir(const std::string&);
// ... 10 more setters

// AFTER
struct FindPairConfig {
    bool legacy_mode = false;
    bool split_output = false;
    std::string output_dir = "output";
    // ... all options with defaults

    // Validation
    [[nodiscard]] bool is_valid() const;
};

class FindPairProtocol {
public:
    explicit FindPairProtocol(FindPairConfig config);
};
```

### Step 5.2: Refactor App Main Functions
**Problem**: `find_pair_app.cpp` main() is 127 lines with deep nesting

**Extract**:
```cpp
// apps/find_pair_app.cpp
int main(int argc, char* argv[]) {
    auto config = parse_arguments(argc, argv);
    if (!config) return 1;

    return run_find_pair(*config);
}

// Extracted functions
std::optional<FindPairConfig> parse_arguments(int argc, char* argv[]);
int run_find_pair(const FindPairConfig& config);
```

### Step 5.3: Add Config Validation to ConfigManager
**Problem**: No validation in setters - invalid state possible

**Add**:
```cpp
class ConfigManager {
public:
    void set_threshold(const std::string& name, double value) {
        validate_threshold(name, value);  // NEW
        thresholds_[name] = value;
    }

private:
    void validate_threshold(const std::string& name, double value);
};
```

---

## Phase 6: Flatten Deep Nesting ✅ DONE

### Step 6.1: Apply Early Return Pattern ✅ DONE
**Target**: All functions with 4+ nesting levels
**Status**: Applied to key algorithm files

**Commits**:
- `b4c63bf` - Flatten nesting with early returns and named booleans (multiple files)
- `2fd9eb2` - Refactor find_best_partner() with early returns and named booleans
- `31e3ae4` - Refactor check_strand2() with named booleans and early continue
- `22c9cf7` - Refactor check_direction() with named booleans and early return

**Example transformation**:
```cpp
// BEFORE (4+ levels)
void process() {
    if (condition1) {
        if (condition2) {
            for (auto& item : items) {
                if (condition3) {
                    // actual work
                }
            }
        }
    }
}

// AFTER (max 2 levels)
void process() {
    if (!condition1) return;
    if (!condition2) return;

    for (auto& item : items) {
        process_item(item);
    }
}

void process_item(const Item& item) {
    if (!condition3) return;
    // actual work
}
```

### Step 6.2: Extract Complex Conditionals ✅ DONE
**Target**: All multi-clause conditionals
**Status**: Named booleans added throughout algorithm files

```cpp
// BEFORE
if (result.distance_check && result.d_v_check &&
    result.plane_angle_check && result.dNN_check && result.overlap_check) {

// AFTER
const bool all_geometric_checks_pass =
    result.distance_check &&
    result.d_v_check &&
    result.plane_angle_check &&
    result.dNN_check &&
    result.overlap_check;

if (all_geometric_checks_pass) {
```

---

## Phase 7: Final Cleanup

### Step 7.1: Remove Dead Code
- Unused includes
- Commented-out code
- Unreachable branches

### Step 7.2: Consistent Naming
- Verify all names follow conventions
- Fix any inconsistencies found

### Step 7.3: Documentation Pass
- Update CLAUDE.md with new architecture
- Add architecture diagram
- Document extracted class responsibilities

### Step 7.4: Final Validation
```bash
fp2-validate validate core --test-set fast
fp2-validate validate steps --test-set fast
make test
```

---

## Commit Strategy

Each step = one commit. Format:
```
refactor(<area>): <what changed>

- Specific change 1
- Specific change 2
- Validation: <test command and result>

🤖 Generated with Claude Code
```

---

## Execution Order

**Recommended sequence** (safest to riskiest):

1. **Phase 0** (portability) - Critical foundation for embedding
2. **Phase 1** (all steps) - Zero risk, establishes patterns
3. **Phase 6** (nesting) - Logic-preserving transformations
4. **Phase 2.1-2.2** (atom, residue) - Core data classes with Builder
5. **Phase 4.1** (pdb_parser) - Self-contained I/O
6. **Phase 3.2** (validator) - Extract validators
7. **Phase 3.3** (helix_organizer) - Largest extraction
8. **Phase 3.1** (finder) - Complex state
9. **Phase 2.3-2.4** (base_pair, JSON) - Cross-cutting
10. **Phase 4.2-4.3** (json_writer, format) - Output formatting
11. **Phase 5** (protocols, apps) - Top-level changes
12. **Phase 7** (cleanup) - Final polish

---

## Unit Test Structure

All extracted classes require unit tests. New structure:

```
tests/
├── unit/                              # NEW - Unit tests for extracted classes
│   ├── core/
│   │   ├── atom_test.cpp
│   │   ├── residue_test.cpp
│   │   ├── base_pair_test.cpp
│   │   └── hydrogen_bond_test.cpp
│   ├── config/
│   │   ├── resource_locator_test.cpp
│   │   └── config_manager_test.cpp
│   ├── algorithms/
│   │   ├── pairing/
│   │   │   ├── pair_candidate_generator_test.cpp
│   │   │   ├── mutual_best_matcher_test.cpp
│   │   │   └── quality_scorer_test.cpp
│   │   ├── validation/
│   │   │   ├── distance_validator_test.cpp
│   │   │   ├── angle_validator_test.cpp
│   │   │   ├── overlap_calculator_test.cpp
│   │   │   └── ring_atom_finder_test.cpp
│   │   └── helix/
│   │       ├── backbone_linkage_checker_test.cpp
│   │       ├── pair_geometry_helper_test.cpp
│   │       ├── helix_context_calculator_test.cpp
│   │       ├── strand_direction_checker_test.cpp
│   │       └── five_to_three_orderer_test.cpp
│   ├── io/
│   │   ├── pdb_column_parser_test.cpp
│   │   └── json_serializer_test.cpp
│   └── initialization_test.cpp
├── integration/                       # EXISTING - End-to-end tests
│   └── ... (existing tests)
└── cmake/                             # NEW - CMake integration tests
    └── test_add_subdirectory/
        └── CMakeLists.txt
```

**Test framework**: Google Test (already configured)

**Coverage target**: 80%+ for all extracted classes

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Breaking legacy match | Run `fp2-validate validate core --test-set 100` after every step |
| Merge conflicts | Small, focused commits; rebase frequently |
| Performance regression | Profile before/after Phase 3 extractions |
| Incomplete extraction | Each extracted class must be independently testable |

---

## Estimated Scope

| Phase | Files Changed | New Files | New Tests | Lines Changed |
|-------|---------------|-----------|-----------|---------------|
| Phase 0 | 6 | 3 | 3 | ~400 |
| Phase 1 | 50 | 1 | 0 | ~500 |
| Phase 2 | 6 | 4 | 4 | ~1,800 |
| Phase 3 | 3 | 15 | 12 | ~3,500 |
| Phase 4 | 3 | 4 | 2 | ~1,200 |
| Phase 5 | 5 | 0 | 0 | ~400 |
| Phase 6 | 10 | 0 | 0 | ~600 |
| Phase 7 | 20 | 0 | 0 | ~200 |
| **TOTAL** | **~60** | **~27** | **~21** | **~8,600** |

---

## Final Directory Structure

After refactoring:
```
include/x3dna/
├── x3dna.hpp                    # Main entry point with initialize()
├── core/
│   ├── atom.hpp                 # Builder pattern
│   ├── residue.hpp              # Builder pattern
│   ├── base_pair.hpp            # Builder pattern
│   ├── hydrogen_bond.hpp        # NEW - extracted
│   ├── constants.hpp            # NEW - magic numbers
│   └── ...
├── config/
│   ├── config_manager.hpp
│   └── resource_locator.hpp     # NEW - portable resource loading
├── algorithms/
│   ├── pairing/                 # NEW subdir
│   │   ├── pair_candidate_generator.hpp
│   │   ├── mutual_best_matcher.hpp
│   │   ├── quality_scorer.hpp
│   │   └── base_pair_finder.hpp
│   ├── validation/              # NEW subdir
│   │   ├── distance_validator.hpp
│   │   ├── angle_validator.hpp
│   │   ├── overlap_calculator.hpp
│   │   ├── ring_atom_finder.hpp
│   │   └── base_pair_validator.hpp
│   ├── helix/                   # NEW subdir
│   │   ├── backbone_linkage_checker.hpp
│   │   ├── pair_geometry_helper.hpp
│   │   ├── helix_context_calculator.hpp
│   │   ├── strand_direction_checker.hpp
│   │   ├── five_to_three_orderer.hpp
│   │   └── helix_organizer.hpp
│   └── ...
├── io/
│   ├── pdb_parser.hpp
│   ├── pdb_column_parser.hpp    # NEW
│   ├── json_serializer.hpp      # NEW - unified JSON
│   └── ...
└── ...
```

---

## Next Action

**Recommended next steps** (in priority order):

1. **Phase 0, Step 0.1**: Create `ResourceLocator` for portable resource path resolution
   - Critical for embedding the library in other projects
   - Removes hardcoded paths and X3DNA_SOURCE_DIR macro

2. **Phase 3.2**: Continue extracting validators from `base_pair_validator.cpp`
   - Already extracted: OverlapCalculator
   - Remaining: distance_validator, angle_validator, ring_atom_finder

3. **Phase 3.3**: Continue extracting from `helix_organizer.cpp` (1,277 lines)
   - Already added: FrameAlignment, compute_frame_alignment, update_direction_count
   - Could extract: backbone_linkage_checker, pair_geometry_helper

4. **Phase 2**: Builder patterns for core data structures
   - Atom, Residue, BasePair classes
   - Requires careful API changes

```bash
# After each step:
make release && make test && fp2-validate validate core --test-set 100
```

**Current validation status**:
- 100% pass rate on 100-PDB test set
- 99.1% pass rate on fast set (3569/3602, 33 known failures in complex structures)
