# Naming Conventions

## File Naming

All files follow **snake_case** naming convention:

### Header Files
- ✅ `test_helpers.hpp`
- ✅ `test_fixtures.hpp`
- ✅ `test_data_discovery.hpp`
- ✅ `integration_test_base.hpp`
- ✅ `json_comparison.hpp`
- ✅ `vector3d.hpp`
- ✅ `matrix3d.hpp`
- ✅ `pdb_parser.hpp`

### Source Files
- ✅ `test_basic.cpp`
- ✅ `test_main.cpp`
- ✅ `vector3d.cpp`
- ✅ `pdb_parser.cpp`

### Test Files
- ✅ `test_atom.cpp`
- ✅ `test_residue.cpp`
- ✅ `test_pdb_parser.cpp`
- ✅ `test_frame_calculation_integration.cpp`

## Class Naming

Classes follow **PascalCase** (standard C++ convention):

- ✅ `Vector3D`
- ✅ `Matrix3D`
- ✅ `PdbParser`
- ✅ `BaseFrameCalculator`
- ✅ `test_data_discovery` (namespace-level utility class)

## Variable Naming

All variables follow **snake_case**:

```cpp
// Good
std::string pdb_name;
std::filesystem::path json_file;
std::vector<pdb_json_pair> pairs;
double rms_fit;
size_t residue_idx;

// Bad
std::string pdbName;      // camelCase
std::string PdbName;      // PascalCase
double RMSFit;            // PascalCase
```

## Function Naming

All functions follow **snake_case**:

```cpp
// Good
void calculate_frame();
std::vector<BasePair> find_pairs();
bool has_json();
double distance_to();

// Bad
void calculateFrame();    // camelCase
void CalculateFrame();    // PascalCase
```

## Constant Naming

Constants follow **UPPER_SNAKE_CASE**:

```cpp
// Good
constexpr double DEFAULT_TOLERANCE = 0.001;
constexpr int MAX_ATOMS = 10000;
constexpr const char* JSON_EXTENSION = ".json";

// Bad
constexpr double defaultTolerance = 0.001;  // camelCase
constexpr double DefaultTolerance = 0.001;   // PascalCase
```

## Enum Naming

Enum types use **PascalCase**, enum values use **UPPER_SNAKE_CASE**:

```cpp
// Good
enum class ResidueType {
    ADENINE,
    CYTOSINE,
    GUANINE
};

enum class BasePairType {
    WATSON_CRICK,
    WOBBLE,
    HOOGSTEEN
};

// Bad
enum class residueType {      // camelCase
    adenine,                  // camelCase
    CYTOSINE                  // Mixed
};
```

## Namespace Naming

Namespaces use **snake_case**:

```cpp
// Good
namespace x3dna::core {}
namespace x3dna::io {}
namespace x3dna::test {}

// Bad
namespace x3dna::Core {}      // PascalCase
namespace x3dna::IO {}        // UPPER_CASE
```

## Struct Naming

Structs follow **snake_case**:

```cpp
// Good
struct pdb_json_pair {
    std::filesystem::path pdb_file;
    std::filesystem::path json_file;
};

struct base_pair_step_parameters {
    double shift;
    double slide;
    double rise;
};

// Bad
struct PdbJsonPair {};        // PascalCase
struct basePairStepParams {}; // camelCase
```

## Code Style Rules

### Avoid Ternary Operators

**Do NOT use ternary operators** (`? :`). Use explicit if-else statements instead for better readability and maintainability.

```cpp
// Bad - Ternary operator
size_t max_pairs = test_all_pdbs ? pairs_.size() : std::min(static_cast<size_t>(10), pairs_.size());
int value = condition ? 5 : 10;
std::string name = is_valid ? "valid" : "invalid";

// Good - Explicit if-else
size_t max_pairs;
if (test_all_pdbs) {
    max_pairs = pairs_.size();
} else {
    max_pairs = std::min(static_cast<size_t>(10), pairs_.size());
}

int value;
if (condition) {
    value = 5;
} else {
    value = 10;
}

std::string name;
if (is_valid) {
    name = "valid";
} else {
    name = "invalid";
}
```

**Rationale**: Explicit if-else statements are:
- More readable, especially for complex conditions
- Easier to debug (can set breakpoints on branches)
- Easier to modify (can add logging or additional logic)
- More consistent with the codebase style

### Limit Function Length

**Keep functions concise** - aim for functions **no longer than 100 lines**. If a function exceeds this, refactor by:
- Extracting logical sections into helper functions
- Breaking complex operations into smaller, focused functions
- Separating validation, processing, and output logic
- Using composition to combine smaller functions

**Rationale**: Shorter functions are:
- Easier to read and understand
- Easier to test and debug
- More reusable
- Less prone to bugs
- Easier to maintain and modify

### Limit Indentation Levels

**Avoid deep nesting** - limit indentation to **4 levels maximum**. If code requires more than 4 levels of indentation, refactor by:
- Extracting functions
- Using early returns
- Inverting conditions
- Breaking complex logic into smaller functions

```cpp
// Bad - 5 levels of indentation
void process_data() {
    if (condition1) {
        if (condition2) {
            for (auto& item : items) {
                if (condition3) {
                    if (condition4) {
                        // Too deeply nested!
                        do_something();
                    }
                }
            }
        }
    }
}

// Good - Extract function to reduce nesting
void process_data() {
    if (!condition1 || !condition2) {
        return;  // Early return
    }
    for (auto& item : items) {
        process_item(item);
    }
}

void process_item(const Item& item) {
    if (!condition3 || !condition4) {
        return;  // Early return
    }
    do_something();
}
```

**Guidelines**:
- Use early returns to reduce nesting
- Extract complex nested logic into separate functions
- Invert conditions when possible (`if (!condition) return;`)
- Consider using helper functions for deeply nested loops
- Maximum recommended: 4 levels of indentation

**Rationale**: Deep nesting:
- Makes code harder to read and understand
- Increases cognitive load
- Makes debugging more difficult
- Increases risk of bugs
- Makes refactoring harder

### Factory Methods for Non-Trivial Construction

**Provide factory methods** for all classes that don't have trivial construction. A class has non-trivial construction if:
- It requires more than 2 parameters in its constructor
- Construction involves complex initialization logic
- Construction requires validation or error handling
- Multiple constructors exist with overlapping functionality

**Use static factory methods** instead of (or in addition to) complex constructors:

```cpp
// Bad - Constructor with many parameters (error-prone)
Atom atom(" C1'", Vector3D(1.0, 2.0, 3.0), "  C", 'A', 42, 'A');
Residue residue("  C", 42, 'A');
ReferenceFrame frame(matrix, origin);
BasePair bp(0, 1, BasePairType::WATSON_CRICK);

// Good - Factory methods with descriptive names
Atom atom = Atom::create(" C1'", Vector3D(1.0, 2.0, 3.0), "  C", 'A', 42, 'A');
Residue residue = Residue::create("  C", 42, 'A');
ReferenceFrame frame = ReferenceFrame::create(matrix, origin);
BasePair bp = BasePair::create(0, 1, BasePairType::WATSON_CRICK);

// Even better - Named parameter factories for complex cases
Atom atom = Atom::create_with_metadata({
    .name = " C1'",
    .position = Vector3D(1.0, 2.0, 3.0),
    .residue_name = "  C",
    .chain_id = 'A',
    .residue_seq = 42,
    .record_type = 'A'
});
```

**Factory Method Naming Conventions**:
- Use `create()` for the most common construction pattern
- Use `create_from_<source>()` for construction from specific sources (e.g., `create_from_json()`)
- Use descriptive names for alternative construction patterns (e.g., `create_with_metadata()`)
- Keep constructors available but prefer factories for complex cases

**Examples from the codebase**:
- ✅ `ReferenceFrame::from_json()` - Factory for JSON deserialization
- ✅ `Atom::from_json()` - Factory for JSON deserialization
- ✅ `Structure::from_json()` - Factory for JSON deserialization
- ✅ `Residue::from_json()` - Factory for JSON deserialization

**When to Use Factories**:
- ✅ Classes with 3+ constructor parameters
- ✅ Classes with complex initialization logic
- ✅ Classes that need validation during construction
- ✅ Classes with multiple construction patterns
- ✅ Classes used frequently throughout the codebase

**When Constructors are Fine**:
- ✅ Trivial construction (0-2 simple parameters)
- ✅ Default constructors
- ✅ Copy/move constructors (auto-generated is fine)

**Rationale**: Factory methods provide:
- Better readability (descriptive method names)
- Easier refactoring (centralize construction logic)
- Validation and error handling at construction time
- More flexible construction patterns
- Better encapsulation of construction complexity

## Summary

| Element | Convention | Example |
|---------|-----------|---------|
| Files | snake_case | `test_helpers.hpp` |
| Classes | PascalCase | `Vector3D` |
| Structs | snake_case | `pdb_json_pair` |
| Functions | snake_case | `calculate_frame()` |
| Variables | snake_case | `pdb_name` |
| Constants | UPPER_SNAKE_CASE | `DEFAULT_TOLERANCE` |
| Enums (type) | PascalCase | `ResidueType` |
| Enums (values) | UPPER_SNAKE_CASE | `ADENINE` |
| Namespaces | snake_case | `x3dna::core` |
| Operators | Avoid ternary (`? :`) | Use explicit if-else |
| Factories | `create()`, `create_from_<source>()` | `Atom::create()`, `Structure::from_json()` |

---

*This convention ensures consistency throughout the codebase.*

