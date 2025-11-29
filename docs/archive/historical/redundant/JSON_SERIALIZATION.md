# JSON Serialization Design

## Overview

All core structure classes in the modernized X3DNA library support JSON serialization using the `nlohmann/json` library. This enables:

1. **Regression Testing**: Compare outputs with `data/json_legacy/*.json` files
2. **Data Persistence**: Save/load structures to/from JSON
3. **Debugging**: Export structures at any processing stage
4. **Interoperability**: Exchange data with other tools

## JSON Format Support

### Dual Format Strategy

Each class supports two JSON formats:

1. **Modern Format** (`to_json()` / `from_json()`)
   - Clean, structured JSON optimized for the new C++ API
   - Uses consistent naming conventions
   - Optimized for readability and maintainability

2. **Legacy Format** (`to_json_legacy()` / `from_json_legacy()`)
   - Exact match with `data/json_legacy/*.json` files
   - Preserves all original field names and structure
   - Enables direct comparison for regression testing

## Core Classes with JSON Support

### ✅ Atom
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // {"atom_name": "...", "xyz": [x, y, z]}
static Atom from_json(const nlohmann::json& j);
static Atom from_json_legacy(const nlohmann::json& j);
```

### ✅ Residue
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;
static Residue from_json(const nlohmann::json& j);
static Residue from_json_legacy(const nlohmann::json& j);
```

### ✅ Chain
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;
static Chain from_json(const nlohmann::json& j);
static Chain from_json_legacy(const nlohmann::json& j);
```

### ✅ Structure
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // Full structure with calculations array
static Structure from_json(const nlohmann::json& j);
static Structure from_json_legacy(const nlohmann::json& j);
```

### ✅ ReferenceFrame
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // {"orien": [[...]], "org": [x, y, z]}
static ReferenceFrame from_json(const nlohmann::json& j);
static ReferenceFrame from_json_legacy(const nlohmann::json& j);
```

### ✅ BasePair
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // {"type": "base_pair", "base_i": 1, ...}
static BasePair from_json(const nlohmann::json& j);
static BasePair from_json_legacy(const nlohmann::json& j);
```

### ✅ BasePairStepParameters
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // {"params": {"Shift": ..., "Slide": ..., ...}}
static BasePairStepParameters from_json(const nlohmann::json& j);
static BasePairStepParameters from_json_legacy(const nlohmann::json& j);
```

### ✅ HelicalParameters
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // {"params": [x, y, z, ...]}
static HelicalParameters from_json(const nlohmann::json& j);
static HelicalParameters from_json_legacy(const nlohmann::json& j);
```

### ✅ HydrogenBond
```cpp
nlohmann::json to_json() const;
nlohmann::json to_json_legacy() const;  // {"hbond_idx": 1, "donor_atom": "...", ...}
static HydrogenBond from_json(const nlohmann::json& j);
static HydrogenBond from_json_legacy(const nlohmann::json& j);
```

## Legacy JSON Structure Mapping

### Calculation Record Types

The legacy JSON format uses a `calculations` array containing various record types:

| Record Type | Class | Method |
|------------|-------|--------|
| `base_frame_calc` | BaseFrameCalculator::FrameCalculationResult | `to_json_legacy()` |
| `ls_fitting` | LeastSquaresFitter::FitResult | `to_json_legacy()` |
| `frame_calc` | BaseFrameCalculator::FrameCalculationResult | `to_json_legacy()` |
| `base_pair` | BasePair | `to_json_legacy()` |
| `pair_validation` | BasePairFinder::ValidationResult | `to_json_legacy()` |
| `hbond_list` | BasePair | `to_json_legacy()` (includes hbonds) |
| `bpstep_params` | BasePairStepParameters | `to_json_legacy()` |
| `helical_params` | HelicalParameters | `to_json_legacy()` |
| `pdb_atoms` | Structure | `to_json_legacy()` (atoms section) |
| `ref_frame` | ReferenceFrame | `to_json_legacy()` |
| `base_pairs` | Structure | `to_json_legacy()` (base pairs section) |
| `bp_sequence` | Structure | `to_json_legacy()` (sequence section) |
| `ry_classification` | Structure | `to_json_legacy()` (RY section) |
| `ring_atoms` | Residue | `to_json_legacy()` (ring atoms) |
| `distance_checks` | BasePairFinder::ValidationResult | `to_json_legacy()` |

## Usage Examples

### Reading Legacy JSON

```cpp
#include <x3dna/core/Structure.hpp>
#include <x3dna/io/JsonReader.hpp>
#include <fstream>
#include <nlohmann/json.hpp>

// Load legacy JSON file
std::ifstream file("data/json_legacy/100D.json");
nlohmann::json legacy_json;
file >> legacy_json;

// Parse structure from legacy format
Structure structure = Structure::from_json_legacy(legacy_json);

// Access data
auto base_pairs = structure.base_pairs();
for (const auto& pair : base_pairs) {
    std::cout << "Pair: " << pair.residue_index1() 
              << " - " << pair.residue_index2() << std::endl;
}
```

### Writing Legacy JSON

```cpp
// Create structure and process
Structure structure = parser.parse_file("data/pdb/100D.pdb");
FindPairProtocol protocol;
protocol.execute(structure);

// Export in legacy format
nlohmann::json output = structure.to_json_legacy();

// Write to file
std::ofstream out("output.json");
out << output.dump(2);  // Pretty print with 2-space indent
```

### Regression Testing

```cpp
#include <gtest/gtest.h>

TEST(JsonRegression, CompareWithLegacy) {
    // Load legacy reference
    std::ifstream ref_file("data/json_legacy/100D.json");
    nlohmann::json legacy_json;
    ref_file >> legacy_json;
    
    // Process structure
    PdbParser parser;
    Structure structure = parser.parse_file("data/pdb/100D.pdb");
    FindPairProtocol protocol;
    protocol.execute(structure);
    
    // Generate our output
    nlohmann::json our_json = structure.to_json_legacy();
    
    // Compare specific sections
    auto legacy_calcs = legacy_json["calculations"];
    auto our_calcs = our_json["calculations"];
    
    // Compare base pairs
    for (size_t i = 0; i < our_calcs.size(); ++i) {
        if (our_calcs[i]["type"] == "base_pair") {
            // Find corresponding legacy record
            auto legacy_pair = find_by_type(legacy_calcs, "base_pair", i);
            
            // Compare values
            EXPECT_NEAR(
                our_calcs[i]["org_i"][0],
                legacy_pair["org_i"][0],
                0.001
            );
        }
    }
}
```

### Using JsonWriter Utility

```cpp
#include <x3dna/io/JsonWriter.hpp>

JsonWriter writer;
writer.init("data/pdb/100D.pdb");

// Record calculations as they happen
for (const auto& residue : structure.nucleotides()) {
    auto frame = frame_calculator.calculate_with_metrics(*residue);
    writer.record_base_frame_calc(
        residue->index(),
        residue->one_letter_code(),
        frame.template_file,
        frame.rms_fit,
        frame.matched_atoms
    );
}

// Record base pairs
for (const auto& pair : structure.base_pairs()) {
    writer.record_base_pair(pair);
}

// Get final JSON
auto json_output = writer.get_json();
writer.write_to_file("output.json");
```

## Implementation Notes

### Matrix Format

Legacy format uses nested arrays for 3x3 matrices:
```json
"orien": [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]]
```

Modern format can use flattened array:
```json
"rotation": [r11, r12, r13, r21, r22, r23, r31, r32, r33]
```

### Numeric Precision

Legacy format uses 6 decimal places:
```cpp
fprintf(fp, "%.6f", value);
```

JSON serialization should match:
```cpp
json["value"] = std::round(value * 1000000.0) / 1000000.0;
// Or use nlohmann::json's default precision (which is sufficient)
```

### Field Name Mapping

| Legacy Field | Modern Field | Notes |
|-------------|--------------|-------|
| `atom_name` | `name` | Atom class |
| `xyz` | `position` | Vector3D |
| `orien` | `rotation` | Matrix3D |
| `org` | `origin` | Vector3D |
| `base_i`, `base_j` | `residue_index1`, `residue_index2` | BasePair |
| `bp_type` | `type` | BasePairType enum |

## Testing Strategy

### Unit Tests

Test each class's JSON serialization independently:

```cpp
TEST(AtomJson, Serialization) {
    Atom atom(" C1'", 1.0, 2.0, 3.0);
    auto json = atom.to_json_legacy();
    
    EXPECT_EQ(json["atom_name"], " C1'");
    EXPECT_DOUBLE_EQ(json["xyz"][0], 1.0);
    
    Atom restored = Atom::from_json_legacy(json);
    EXPECT_EQ(restored.name(), atom.name());
}
```

### Integration Tests

Test full structure serialization:

```cpp
TEST(StructureJson, FullSerialization) {
    Structure structure = load_test_structure();
    
    auto json = structure.to_json_legacy();
    Structure restored = Structure::from_json_legacy(json);
    
    EXPECT_EQ(restored.num_residues(), structure.num_residues());
    EXPECT_EQ(restored.base_pairs().size(), structure.base_pairs().size());
}
```

### Regression Tests

Compare with legacy JSON files:

```cpp
TEST(JsonRegression, MatchLegacyFormat) {
    // Load legacy JSON
    auto legacy = load_json("data/json_legacy/100D.json");
    
    // Generate our JSON
    Structure structure = process_structure("data/pdb/100D.pdb");
    auto our_json = structure.to_json_legacy();
    
    // Compare key sections
    compare_calculations(legacy["calculations"], our_json["calculations"]);
}
```

## Benefits

1. ✅ **Complete Compatibility**: Can read/write legacy JSON format exactly
2. ✅ **Regression Testing**: Direct comparison with `data/json_legacy/*.json`
3. ✅ **Data Persistence**: Save/load structures easily
4. ✅ **Debugging**: Export at any processing stage
5. ✅ **Interoperability**: Exchange data with other tools
6. ✅ **Documentation**: JSON serves as data format documentation

## Dependencies

- **nlohmann/json**: Header-only JSON library
  - Add to CMakeLists.txt: `find_package(nlohmann_json REQUIRED)`
  - Or use FetchContent for header-only version

---

*All structure classes are fully JSON-serializable with support for both modern and legacy formats.*

