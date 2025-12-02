# Refactor Plan: generate_modern_json.cpp

**Date**: December 2, 2025  
**Goal**: Simplify tool to only call production code, let objects write their own JSON  
**Philosophy**: Tools orchestrate, objects do the work

---

## Current Problems

### Problem 1: Too Much Custom Logic (500+ lines)
`generate_modern_json.cpp` currently has:
- Custom residue iteration loops (lines 211-367)
- Custom pair finding logic (lines 186-209, 380-438)
- Custom index validation (lines 481-618)
- Duplicate code for --only-paired mode

**Should be**: Call protocols, let them do the work

### Problem 2: Flags That Shouldn't Exist
- `--legacy`: Should be default behavior (always match legacy)
- `--fix-indices`: Should happen automatically
- `--only-paired`: Not needed if we match legacy by default

**Should be**: Tool has no flags, just runs production code

### Problem 3: JsonWriter Does Too Much
JsonWriter is passed around and called from many places:
- Tool calls it directly (line 121, 124, etc.)
- Protocol calls it
- Objects don't write themselves

**Should be**: Each object writes its own JSON

---

## Target Architecture

### Simple Tool (~100 lines with proper output)

```cpp
int main(int argc, char* argv[]) {
    // 1. Parse arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <output_dir>\n";
        return 1;
    }
    
    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path output_dir = argv[2];
    
    // Validate input
    if (!std::filesystem::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }
    
    std::filesystem::create_directories(output_dir);
    std::string pdb_id = pdb_file.stem().string();
    
    std::cout << "Processing: " << pdb_id << "\n";
    std::cout << "Input: " << pdb_file << "\n";
    std::cout << "Output: " << output_dir << "\n\n";
    
    try {
        // 2. Parse PDB
        std::cout << "Step 1: Parsing PDB...\n";
        PdbParser parser;
        Structure structure = parser.parse_file(pdb_file);
        std::cout << "  ✅ Parsed " << structure.num_atoms() << " atoms, "
                  << structure.num_residues() << " residues\n\n";
        
        // 3. Write atoms JSON
        std::cout << "Step 2: Writing atoms...\n";
        structure.write_atoms_json(output_dir);
        std::cout << "  ✅ Wrote pdb_atoms/" << pdb_id << ".json\n\n";
        
        // 4. Run find_pair protocol (calculates frames, finds pairs, writes JSON)
        std::cout << "Step 3: Running find_pair protocol...\n";
        FindPairProtocol protocol(output_dir);
        auto result = protocol.execute(structure);
        std::cout << "  ✅ Calculated " << result.num_frames << " frames\n";
        std::cout << "  ✅ Found " << result.num_pairs << " base pairs\n";
        std::cout << "  ✅ Wrote frame JSON files\n";
        std::cout << "  ✅ Wrote pair JSON files\n\n";
        
        // 5. Optional: Run analyze protocol for step parameters
        // std::cout << "Step 4: Calculating step parameters...\n";
        // AnalyzeProtocol analyze(output_dir);
        // auto step_result = analyze.execute(structure);
        // std::cout << "  ✅ Calculated " << step_result.num_steps << " step parameters\n\n";
        
        std::cout << "✅ Success! Generated all JSON for " << pdb_id << "\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
        return 1;
    }
}
```

**Result**: ~100 lines with clear output at each step, but all logic is in objects/protocols.

---

## Refactoring Steps

### Step 1: Move JSON Writing to Objects

#### 1a. Structure writes its own atoms
```cpp
// In src/x3dna/core/structure.cpp
void Structure::write_atoms_json(const std::filesystem::path& output_dir) const {
    nlohmann::json record;
    record["num_atoms"] = num_atoms();
    record["atoms"] = nlohmann::json::array();
    
    for (const auto& chain : chains_) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                record["atoms"].push_back(atom.to_json());  // Atom writes itself
            }
        }
    }
    
    // Write file
    std::filesystem::path file = output_dir / "pdb_atoms" / (pdb_id_ + ".json");
    std::filesystem::create_directories(file.parent_path());
    std::ofstream out(file);
    out << record.dump(2);
}
```

#### 1b. Atom writes itself
```cpp
// In src/x3dna/models/atom.cpp
nlohmann::json Atom::to_json() const {
    nlohmann::json j;
    j["atom_idx"] = legacy_atom_idx_;  // Use legacy index
    j["atom_name"] = name_;
    j["residue_name"] = residue_name_;
    j["chain_id"] = std::string(1, chain_id_);
    j["residue_seq"] = residue_seq_;
    j["xyz"] = {position_.x(), position_.y(), position_.z()};
    j["record_type"] = std::string(1, record_type_);
    if (line_number_ > 0) {
        j["line_number"] = line_number_;
    }
    return j;
}
```

#### 1c. Residue writes its frame
```cpp
// In src/x3dna/models/residue.cpp
nlohmann::json Residue::frame_to_json() const {
    if (!reference_frame_.has_value()) {
        return nlohmann::json();  // No frame
    }
    
    nlohmann::json j;
    j["residue_idx"] = legacy_residue_idx_;  // Use legacy index
    // ... write frame data
    return j;
}
```

#### 1d. BasePair writes itself
```cpp
// In src/x3dna/models/base_pair.cpp
nlohmann::json BasePair::to_json() const {
    nlohmann::json j;
    j["base_i"] = residue_idx1_ + 1;  // Convert to 1-based
    j["base_j"] = residue_idx2_ + 1;
    // ... write pair data
    return j;
}
```

### Step 2: Protocols Write Their Results

#### FindPairProtocol writes all find_pair outputs
```cpp
// In src/x3dna/protocols/find_pair_protocol.cpp
void FindPairProtocol::execute(Structure& structure) {
    // 1. Calculate frames
    BaseFrameCalculator calculator;
    calculator.calculate_all_frames(structure);
    
    // 2. Write frame JSON
    calculator.write_json(structure, output_dir_);
    
    // 3. Find pairs
    BasePairFinder finder;
    auto pairs = finder.find_pairs(structure);
    
    // 4. Write pair JSON
    finder.write_json(pairs, structure, output_dir_);
    
    // 5. Write selection JSON
    write_selection_json(pairs, output_dir_);
}
```

#### AnalyzeProtocol writes step parameters
```cpp
// In src/x3dna/protocols/analyze_protocol.cpp
void AnalyzeProtocol::execute(Structure& structure) {
    // 1. Calculate parameters
    ParameterCalculator calc;
    auto params = calc.calculate_step_parameters(base_pairs_, structure);
    
    // 2. Write parameter JSON
    calc.write_json(params, output_dir_);
}
```

### Step 3: Remove JsonWriter Class (Eventually)

**Current**: Centralized JsonWriter that everything calls  
**Target**: Each object writes its own JSON

**Migration**:
1. Add `to_json()` methods to all classes
2. Add `write_json()` methods to algorithms/protocols
3. Deprecate JsonWriter
4. Remove JsonWriter once all migrated

---

## Why Remove the Flags?

### --legacy Flag
**Current use**: Exclude C4 from matching  
**Why remove**: We should ALWAYS match legacy behavior  
**Solution**: Make it default, no flag needed

### --fix-indices Flag
**Current use**: Load legacy indices from JSON  
**Why remove**: Indices should always match  
**Solution**: 
- Parse PDB in same order as legacy
- Assign indices in same order
- No need to "fix" if we do it right from the start

### --only-paired Flag
**Current use**: Only record frames for paired residues  
**Why remove**: Legacy behavior should be default  
**Solution**: 
- Check what legacy actually does
- If it only records paired residues, make that default
- If it records all, remove this flag

---

## Implementation Plan

### Phase 1: Add to_json() Methods ✅ DO FIRST

**Files to modify**:
- `src/x3dna/models/atom.cpp` - Add `Atom::to_json()`
- `src/x3dna/models/residue.cpp` - Add `Residue::to_json()`
- `src/x3dna/models/base_pair.cpp` - Add `BasePair::to_json()`
- `src/x3dna/models/chain.cpp` - Add `Chain::to_json()`
- `src/x3dna/core/structure.cpp` - Add `Structure::write_atoms_json()`

**Checklist**:
- [ ] Add `to_json()` to Atom
- [ ] Add `to_json()` to Residue (for frames)
- [ ] Add `to_json()` to BasePair
- [ ] Add `write_atoms_json()` to Structure
- [ ] Test: verify JSON output matches current format

### Phase 2: Protocols Write JSON ✅ DO NEXT

**Files to modify**:
- `src/x3dna/protocols/find_pair_protocol.cpp` - Add output_dir member, write JSON
- `src/x3dna/protocols/analyze_protocol.cpp` - Add output_dir member, write JSON

**Changes**:
```cpp
// Constructor takes output directory
FindPairProtocol::FindPairProtocol(const std::filesystem::path& output_dir)
    : output_dir_(output_dir) {}

// Execute writes JSON
void FindPairProtocol::execute(Structure& structure) {
    // Do calculations
    calculate_frames(structure);
    auto pairs = find_pairs(structure);
    
    // Write JSON (objects write themselves)
    write_all_json(structure, pairs);
}
```

**Checklist**:
- [ ] Add output_dir to FindPairProtocol constructor
- [ ] Add write_all_json() to FindPairProtocol
- [ ] Move JSON writing from tools to protocol
- [ ] Test: verify same JSON is generated

### Phase 3: Simplify generate_modern_json ✅ DO NEXT

**Target code** (30 lines):
```cpp
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <output_dir>\n";
        return 1;
    }
    
    // Parse arguments
    std::filesystem::path pdb_file = argv[1];
    std::filesystem::path output_dir = argv[2];
    std::filesystem::create_directories(output_dir);
    
    // Parse PDB
    PdbParser parser;
    Structure structure = parser.parse_file(pdb_file);
    
    // Write atoms
    structure.write_atoms_json(output_dir);
    
    // Run find_pair (calculates frames, finds pairs, writes JSON)
    FindPairProtocol protocol(output_dir);
    protocol.execute(structure);
    
    std::cout << "Success: " << structure.pdb_id() << "\n";
    return 0;
}
```

**Checklist**:
- [ ] Remove all custom loops
- [ ] Remove flag parsing (--legacy, --fix-indices, --only-paired)
- [ ] Remove index validation code
- [ ] Keep only: parse PDB, run protocol, done
- [ ] Test: verify same output

### Phase 4: Remove JsonWriter ⏳ FUTURE

Once all objects write themselves:
- [ ] Deprecate JsonWriter class
- [ ] Remove from codebase
- [ ] Update all references

---

## Object Responsibilities

### What Each Object Writes

**Atom** (`to_json()`):
- atom_idx, atom_name, residue_name, chain_id, residue_seq
- xyz, record_type, line_number

**Residue** (`frame_to_json()`):
- base_frame_calc: residue info + RMS + matched atoms
- frame_calc: rotation matrix + origin

**BasePair** (`to_json()`):
- base_i, base_j, bp_type
- orien_i, orien_j, org_i, org_j
- H-bonds, validation results

**Structure** (`write_atoms_json()`):
- Collects all atoms
- Writes pdb_atoms/PDB_ID.json

**BaseFrameCalculator** (`write_json()`):
- Collects frame calculations
- Writes base_frame_calc/PDB_ID.json
- Writes frame_calc/PDB_ID.json

**BasePairFinder** (`write_json()`):
- Collects validation, distance checks, H-bonds
- Writes pair_validation/PDB_ID.json
- Writes distance_checks/PDB_ID.json
- Writes hbond_list/PDB_ID.json
- Writes find_bestpair_selection/PDB_ID.json
- Writes base_pair/PDB_ID.json

**ParameterCalculator** (`write_json()`):
- Writes bpstep_params/PDB_ID.json
- Writes helical_params/PDB_ID.json

---

## Benefits of Refactoring

### Simplicity
- Tools: 30 lines (just orchestration)
- Objects: Self-contained (data + serialization)
- Protocols: Workflow only

### Maintainability
- Each class owns its JSON format
- Changes are local (not scattered)
- Easy to add new JSON types

### Testability
- Test object serialization directly
- Mock-free testing (objects write real JSON)
- Clear what's being tested

### Alignment with OOP
- Objects responsible for their data
- No "god class" (JsonWriter)
- Single responsibility principle

---

## Migration Strategy

### Week 1: Add to_json() Methods
- Implement `to_json()` for all classes
- Keep existing JsonWriter (parallel implementation)
- Test both produce same output

### Week 2: Protocols Write JSON
- Add `write_json()` to protocols
- Protocols use objects' `to_json()` methods
- Test protocols produce same output

### Week 3: Simplify Tools
- Update generate_modern_json to use protocols
- Remove custom logic
- Remove unnecessary flags
- Test produces same output

### Week 4: Deprecate JsonWriter
- Remove JsonWriter class
- Clean up references
- Final testing

---

## Specific Refactoring Tasks

### Task 1: Remove --legacy Flag

**Current behavior**: Excludes C4 from matching  
**Analysis**: Should always match legacy  
**Action**:
```cpp
// In BaseFrameCalculator, make this default:
set_legacy_mode(true);  // Always

// Remove flag from tools
// Remove flag from protocol constructors
```

### Task 2: Remove --fix-indices Flag

**Current behavior**: Loads legacy indices from JSON  
**Analysis**: Why is this needed?  
**Root cause**: Modern might assign indices differently  
**Better solution**:
```cpp
// In PdbParser::parse_file()
// Assign legacy indices in same order as legacy
// Process atoms in PDB file order
// Assign residue indices based on first atom appearance
// No need to "fix" later
```

**Action**:
- Make PdbParser assign correct indices from the start
- Remove fix_residue_indices_from_json() calls
- Remove flag

### Task 3: Remove --only-paired Flag

**Current behavior**: Only record frames for paired residues  
**Analysis**: Check what legacy does  
**Action**:
- Read legacy code to see if it records all or only paired
- Make modern match by default
- Remove flag

### Task 4: Move Iteration Logic to Protocol

**Current**: Tool loops through residues (lines 211-367)  
**Should be**: Protocol handles iteration

```cpp
// In FindPairProtocol::execute()
void FindPairProtocol::execute(Structure& structure) {
    // Calculate frames for all nucleotides
    calculator_.calculate_all_frames(structure);
    
    // Write frame JSON
    write_frame_json(structure);
    
    // Find pairs
    auto pairs = finder_.find_pairs(structure);
    
    // Write pair JSON
    write_pair_json(pairs, structure);
}
```

**Tool just calls**:
```cpp
protocol.execute(structure);  // Done!
```

### Task 5: Move Index Validation Out

**Current**: Tool does index validation (lines 481-618)  
**Should be**: Not needed if indices are correct from start

**Options**:
- A) Remove entirely (if indices are always correct)
- B) Move to test suite (unit test that indices match)
- C) Make it a separate validation tool (not in generation path)

**Recommendation**: Remove from generate_modern_json, create separate `validate_indices` tool

---

## File Structure After Refactoring

### Tools (Minimal)
```
tools/
├── generate_modern_json.cpp        # 30 lines: parse PDB, run protocol
├── validate_indices.cpp            # NEW: separate validation tool
└── check_residue_indices.cpp       # Existing
```

### Core Objects (Self-Serializing)
```
src/x3dna/models/
├── atom.cpp                  # + to_json()
├── residue.cpp              # + to_json(), frame_to_json()
├── base_pair.cpp            # + to_json()
└── chain.cpp                # + to_json()
```

### Algorithms (Write Results)
```
src/x3dna/algorithms/
├── base_frame_calculator.cpp    # + write_json()
├── base_pair_finder.cpp         # + write_json()
└── parameter_calculator.cpp     # + write_json()
```

### Protocols (Orchestrate)
```
src/x3dna/protocols/
├── find_pair_protocol.cpp       # + constructor(output_dir), write_all_json()
└── analyze_protocol.cpp         # + constructor(output_dir), write_all_json()
```

---

## Implementation Checklist

### Preparation
- [ ] Document current generate_modern_json behavior
- [ ] Identify all JSON types it outputs
- [ ] Map each to responsible class
- [ ] Create test to verify output doesn't change

### Phase 1: Objects (Parallel Implementation)
- [ ] Atom::to_json()
- [ ] Residue::frame_to_json()
- [ ] BasePair::to_json()
- [ ] Structure::write_atoms_json()
- [ ] Test: Objects produce correct JSON

### Phase 2: Algorithms (Parallel Implementation)
- [ ] BaseFrameCalculator::write_json()
- [ ] BasePairFinder::write_json()
- [ ] ParameterCalculator::write_json()
- [ ] Test: Algorithms produce correct JSON

### Phase 3: Protocols (Parallel Implementation)
- [ ] Add output_dir to FindPairProtocol
- [ ] FindPairProtocol::write_all_json()
- [ ] Add output_dir to AnalyzeProtocol
- [ ] AnalyzeProtocol::write_all_json()
- [ ] Test: Protocols produce correct JSON

### Phase 4: Simplify Tool (Cut Over)
- [ ] Update generate_modern_json to use protocols
- [ ] Remove custom loops
- [ ] Remove flags
- [ ] Remove index validation
- [ ] Test: Output matches exactly
- [ ] Delete old JsonWriter references

### Phase 5: Cleanup
- [ ] Remove JsonWriter class
- [ ] Update documentation
- [ ] Final testing

---

## Testing Strategy

### Regression Tests
For each phase:
```bash
# Before changes
./build/generate_modern_json data/pdb/1EHZ.pdb /tmp/before

# After changes
./build/generate_modern_json data/pdb/1EHZ.pdb /tmp/after

# Compare
diff -r /tmp/before /tmp/after
# Should be identical!
```

### Unit Tests
```cpp
TEST(AtomTest, ToJson) {
    Atom atom("N", 1.0, 2.0, 3.0);
    atom.set_legacy_atom_idx(42);
    auto json = atom.to_json();
    EXPECT_EQ(json["atom_idx"], 42);
    EXPECT_EQ(json["atom_name"], "N");
}
```

---

## Quick Start: Minimal Working Example

### Before (Current - 626 lines)
```cpp
// generate_modern_json.cpp - HUGE FILE
// - 500+ lines of custom logic
// - Flag parsing (--legacy, --fix-indices, --only-paired)
// - Custom residue iteration loops
// - Custom pair finding logic
// - Index validation code
// - Duplicate --only-paired logic
// - etc.
```

### After (Target - ~100 lines with output)
```cpp
// generate_modern_json.cpp - CLEAN FILE
int main(int argc, char* argv[]) {
    // Parse args, validate input
    auto [pdb_file, output_dir] = parse_args(argc, argv);
    
    // Parse PDB
    std::cout << "Parsing PDB...\n";
    Structure structure = PdbParser().parse_file(pdb_file);
    std::cout << "  ✅ " << structure.num_atoms() << " atoms, " 
              << structure.num_residues() << " residues\n\n";
    
    // Write atoms
    std::cout << "Writing atoms...\n";
    structure.write_atoms_json(output_dir);
    std::cout << "  ✅ Wrote pdb_atoms\n\n";
    
    // Run find_pair
    std::cout << "Running find_pair...\n";
    auto result = FindPairProtocol(output_dir).execute(structure);
    std::cout << "  ✅ " << result.num_frames << " frames\n";
    std::cout << "  ✅ " << result.num_pairs << " pairs\n\n";
    
    std::cout << "✅ Success!\n";
    return 0;
}
```

**That's the goal!** ~100 lines with clear output, all logic in objects/protocols.

---

## Timeline

### Quick Win (This Week)
- Implement Atom::to_json()
- Implement Structure::write_atoms_json()
- Update tool to use it
- Verify pdb_atoms JSON matches

### Full Refactor (2-3 Weeks)
- Week 1: All to_json() methods
- Week 2: Protocols write JSON
- Week 3: Simplify tools, remove JsonWriter

### Priority
**Medium** - Not blocking validation, but will make code much cleaner

---

## Related Documents

- [CLEAN_SLATE_VALIDATION_PLAN.md](CLEAN_SLATE_VALIDATION_PLAN.md) - Validation strategy
- [CODE_FLOW.md](CODE_FLOW.md) - Current architecture
- [docs/modernization/STAGE_03_IO.md](docs/modernization/STAGE_03_IO.md) - I/O design

---

**Next Step**: Start with Phase 1 - add `to_json()` methods to objects. This can be done in parallel with validation work.

