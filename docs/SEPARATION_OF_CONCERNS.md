# Separation of Concerns: Generation vs Comparison

**Date**: December 2, 2025  
**Principle**: Generation tools should ONLY generate. Comparison tools should ONLY compare.

---

## The Problem

Currently `generate_modern_json.cpp` does BOTH:
1. **Generation**: Creates modern JSON ✅ Correct
2. **Comparison**: Reads legacy JSON, validates indices, compares ❌ Wrong

**This is wrong because**:
- Generation depends on legacy files existing (circular dependency)
- Mixes two separate concerns (generation + validation)
- Can't generate modern JSON independently
- Comparison logic scattered across generation tool

---

## Correct Architecture

### Generation Tools (No Legacy Dependencies)

**Purpose**: Generate JSON from PDB files only

```
Input:  PDB file
        ↓
    [Parse PDB]
        ↓
    [Run Algorithms]
        ↓
   [Write JSON]
        ↓
Output: Modern JSON files
```

**Files**:
- `tools/generate_modern_json.cpp` - Generate modern JSON
- `org/build/bin/find_pair_analyze` - Generate legacy JSON

**Dependencies**: 
- ✅ PDB files
- ✅ Templates
- ❌ NO legacy JSON files
- ❌ NO comparison logic

---

### Comparison Tools (Read-Only)

**Purpose**: Compare legacy vs modern JSON

```
Input:  Legacy JSON + Modern JSON
        ↓
   [Load Both]
        ↓
   [Compare Values]
        ↓
   [Check Tolerances]
        ↓
Output: Validation results
```

**Files**:
- `scripts/compare_json.py` - Python comparison
- `scripts/stage_by_stage_validation.py` - Batch validation
- `scripts/batch_validation_workflow.py` - Workflow automation

**Dependencies**:
- ✅ Legacy JSON (read-only)
- ✅ Modern JSON (read-only)
- ❌ NO generation
- ❌ NO modification

---

## What to Remove from generate_modern_json

### 1. Legacy File Reading

**Remove** (lines 107-113):
```cpp
// Try to find legacy JSON file for PDB line caching (optional)
std::filesystem::path legacy_json_file;
std::filesystem::path legacy_dir = json_output_dir.parent_path() / "json_legacy";
std::filesystem::path legacy_file = legacy_dir / "pdb_atoms" / (pdb_name + ".json");
if (std::filesystem::exists(legacy_file)) {
    legacy_json_file = legacy_file;
}
```

**Why**: Generation shouldn't read legacy. Get PDB lines from original PDB file if needed.

### 2. Index Fixing from Legacy

**Remove** (lines 86-105):
```cpp
if (fix_indices) {
    std::filesystem::path legacy_base_frame = json_output_dir.parent_path() /
                                              "json_legacy" / "base_frame_calc" /
                                              (pdb_name + ".json");
    if (std::filesystem::exists(legacy_base_frame)) {
        int fixed = fix_residue_indices_from_json(structure, legacy_base_frame.string());
        // ...
    }
}
```

**Why**: Modern should generate indices independently. Comparison checks if they match.

### 3. Index Validation Against Legacy

**Remove** (lines 481-618):
```cpp
// ========== Index Validation with ResidueTracker ==========
std::cout << "\n=== Residue Index Validation ===\n";

// Load legacy indices from JSON
std::filesystem::path legacy_json_path = json_output_dir.parent_path() / "json_legacy" /
                                         "base_frame_calc" / (pdb_name + ".json");

if (std::filesystem::exists(legacy_json_path)) {
    // ... loads legacy JSON ...
    // ... validates indices match ...
    // ... exits if mismatch ...
}
```

**Why**: This is comparison, not generation. Should be in separate validation tool.

### 4. --fix-indices Flag

**Remove**:
```cpp
bool fix_indices = false;
if (flag == "--fix-indices") {
    fix_indices = true;
}
```

**Why**: Generation shouldn't "fix" anything. It generates fresh.

---

## Separate Tools, Separate Concerns

### generate_modern_json (Pure Generation)

```cpp
// Input: PDB file
// Output: Modern JSON
// Dependencies: PDB file, templates only
// No legacy file reading!

int main(int argc, char* argv[]) {
    auto [pdb_file, output_dir, stage] = parse_args(argc, argv);
    
    // Parse PDB (from file, not from legacy JSON)
    Structure structure = PdbParser().parse_file(pdb_file);
    
    // Generate JSON based on stage
    if (stage == "atoms" || stage == "all") {
        structure.write_atoms_json(output_dir);
    }
    
    if (is_findpair_stage(stage)) {
        FindPairProtocol(output_dir, stage).execute(structure);
    }
    
    // NO legacy file reading anywhere!
    // NO index validation!
    // NO comparison!
    
    return 0;
}
```

### validate_indices (Separate Tool)

```cpp
// Input: Legacy JSON + Modern JSON
// Output: Validation report
// Pure comparison, no generation

int main(int argc, char* argv[]) {
    auto [legacy_json, modern_json] = parse_args(argc, argv);
    
    // Load both
    auto legacy_data = load_json(legacy_json);
    auto modern_data = load_json(modern_json);
    
    // Compare indices
    auto result = compare_indices(legacy_data, modern_data);
    
    // Report
    std::cout << result.summary();
    
    return result.success ? 0 : 1;
}
```

### compare_json.py (Existing - Good Example)

Already does this correctly:
```python
# Loads legacy JSON (read-only)
legacy = json.load(open('data/json_legacy/frames/1EHZ.json'))

# Loads modern JSON (read-only)
modern = json.load(open('data/json/frames/1EHZ.json'))

# Compares (no generation)
matches = compare(legacy, modern)

# Reports
print(f"Match: {matches}")
```

---

## Benefits of Separation

### Independent Generation
```bash
# Generate modern JSON (no legacy needed!)
./build/generate_modern_json data/pdb/1EHZ.pdb data/json --stage=frames

# Generate legacy JSON (separate process)
cd org && ./build/bin/find_pair_analyze ../data/pdb/1EHZ.pdb
```

### Independent Comparison
```bash
# Compare after both generated
python3 scripts/compare_json.py frames 1EHZ

# Or batch comparison
python3 scripts/stage_by_stage_validation.py --stage stage2_frames
```

### No Circular Dependencies
- Modern generation doesn't need legacy to exist
- Can generate modern JSON for new PDBs without legacy
- Can regenerate modern without touching legacy
- Clean separation of concerns

---

## Refactored Tool Structure

### generate_modern_json.cpp (Generation ONLY)

**Inputs**:
- ✅ PDB file (data/pdb/1EHZ.pdb)
- ✅ Templates (resources/templates/)
- ✅ --stage flag (which JSON to output)

**Outputs**:
- ✅ Modern JSON files

**NOT Allowed**:
- ❌ Reading data/json_legacy/
- ❌ Comparing with legacy
- ❌ Validating indices
- ❌ Any reference to legacy files

**Dependencies**: PDB file + templates ONLY

### Legacy Generation (Separate)

```bash
# Legacy has its own tool
cd org && ./build/bin/find_pair_analyze ../data/pdb/1EHZ.pdb
# → Outputs to data/json_legacy/
```

### Comparison (Separate)

```bash
# Python scripts compare
python3 scripts/compare_json.py frames 1EHZ
# → Loads data/json_legacy/frame_calc/1EHZ.json
# → Loads data/json/frame_calc/1EHZ.json
# → Compares
# → Reports
```

---

## Updated Refactoring Checklist

### Remove Legacy File Reading
- [ ] Remove legacy_json_file parameter from JsonWriter constructor
- [ ] Remove PDB line caching from legacy JSON (get from PDB file instead)
- [ ] Remove --fix-indices flag and all index fixing code
- [ ] Remove index validation code (move to separate tool)
- [ ] Remove all references to data/json_legacy/ from generation code
- [ ] Verify tool runs without any legacy files present

### Remove Comparison Logic
- [ ] Remove index validation (lines 481-618)
- [ ] Move to separate `validate_indices` tool if needed
- [ ] Or just use existing comparison scripts

### Simplify to Pure Generation
- [ ] Parse PDB
- [ ] Run protocol with stage
- [ ] Write JSON
- [ ] Done
- [ ] No legacy file access anywhere

---

## Testing Independence

### Test 1: Generate Without Legacy
```bash
# Delete all legacy JSON
rm -rf data/json_legacy/

# Modern generation should still work!
./build/generate_modern_json data/pdb/1EHZ.pdb data/json --stage=frames
# Should succeed with NO errors about missing legacy files
```

### Test 2: Compare After Generation
```bash
# Generate legacy
cd org && ./build/bin/find_pair_analyze ../data/pdb/1EHZ.pdb

# Generate modern (independent!)
./build/generate_modern_json data/pdb/1EHZ.pdb data/json --stage=frames

# Compare (separate step!)
python3 scripts/compare_json.py frames 1EHZ
```

---

## Summary

**Current**: generate_modern_json reads legacy files, compares, validates (WRONG)  
**Target**: generate_modern_json only generates modern JSON (CORRECT)

**Separation**:
- **Generation**: `generate_modern_json` (no legacy dependencies)
- **Comparison**: `compare_json.py` (reads both, compares)
- **Validation**: `stage_by_stage_validation.py` (orchestrates both)

**Benefit**: Clean architecture, no circular dependencies, tools do one thing well

