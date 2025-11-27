# Legacy Mode Design

**Purpose**: Enable direct comparison with legacy code by allowing mode that breaks some OOP principles for exact compatibility

## Overview

The `--legacy-mode` flag enables a compatibility mode that prioritizes exact matching with legacy code over clean OOP design. This is essential for:
- Regression testing
- Direct comparison with legacy output
- Debugging discrepancies
- Ensuring 100% compatibility

## Design Principles

### Modern Mode (Default)
- Clean OOP design
- 0-based indexing internally
- Modern C++ best practices
- Optimized data structures

### Legacy Mode (`--legacy-mode`)
- Exact match with legacy behavior
- 1-based indexing where needed
- Legacy iteration order
- Legacy data structure organization
- May break encapsulation for compatibility

## Configuration

### ConfigManager Integration

```cpp
namespace x3dna::config {

class ConfigManager {
public:
    // Legacy mode flag
    bool legacy_mode() const { return legacy_mode_; }
    void set_legacy_mode(bool value) { legacy_mode_ = value; }
    
    // Legacy mode affects:
    // - Indexing (1-based vs 0-based)
    // - Iteration order
    // - Output format
    // - Algorithm behavior
    
private:
    bool legacy_mode_ = false;
};

} // namespace x3dna::config
```

### Command-Line Interface

```cpp
// find_pair_app
int main(int argc, char* argv[]) {
    bool legacy_mode = false;
    
    // Parse --legacy-mode flag
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--legacy-mode") {
            legacy_mode = true;
        }
    }
    
    auto& config = ConfigManager::instance();
    config.set_legacy_mode(legacy_mode);
    
    // ... rest of application
}
```

## What Legacy Mode Affects

### 1. Indexing

**Modern Mode** (Default):
- Internal: 0-based indexing
- JSON output: Convert to 1-based for compatibility
- Algorithms: Use 0-based consistently

**Legacy Mode**:
- Internal: Use 1-based indexing where legacy does
- JSON output: 1-based (matches legacy exactly)
- Algorithms: May use 1-based indexing to match legacy iteration

**Implementation**:
```cpp
class BasePairFinder {
    std::vector<BasePair> find_pairs(Structure& structure) {
        auto& config = ConfigManager::instance();
        
        if (config.legacy_mode()) {
            // Use 1-based iteration to match legacy exactly
            return find_pairs_legacy_mode(structure);
        } else {
            // Use 0-based iteration (modern)
            return find_pairs_modern_mode(structure);
        }
    }
    
private:
    std::vector<BasePair> find_pairs_legacy_mode(Structure& structure) {
        // Iterate using legacy indices (1-based)
        auto residues = get_residues_in_legacy_order(structure);
        for (size_t i = 0; i < residues.size(); i++) {
            int legacy_idx = i + 1;  // Convert to 1-based
            // Use legacy_idx for comparison/iteration
        }
    }
    
    std::vector<BasePair> find_pairs_modern_mode(Structure& structure) {
        // Use 0-based indexing (clean OOP)
        for (size_t i = 0; i < structure.residues().size(); i++) {
            // Use i directly (0-based)
        }
    }
};
```

### 2. Iteration Order

**Modern Mode**:
- May optimize iteration order
- Use efficient data structures
- May skip unnecessary iterations

**Legacy Mode**:
- Exact legacy iteration order
- Process residues in PDB file order
- Match legacy sequence exactly

**Implementation**:
```cpp
std::vector<BasePair> find_pairs(Structure& structure) {
    auto& config = ConfigManager::instance();
    
    std::vector<const Residue*> residues;
    if (config.legacy_mode()) {
        // Use legacy order (PDB file order)
        residues = get_residues_in_legacy_order(structure);
    } else {
        // Use optimized order (modern)
        residues = structure.all_residues();
    }
    
    // ... pair finding logic
}
```

### 3. Data Structure Organization

**Modern Mode**:
- Clean OOP encapsulation
- Private members
- Access through methods

**Legacy Mode**:
- May expose internal data for compatibility
- May use legacy data structure layout
- May break encapsulation to match legacy exactly

**Example**:
```cpp
class Residue {
private:
    std::vector<Atom> atoms_;  // Modern: private
    
public:
    // Modern mode: Access through methods
    const std::vector<Atom>& atoms() const { return atoms_; }
    
    // Legacy mode: May need direct access for compatibility
    std::vector<Atom>& atoms_legacy() {
        auto& config = ConfigManager::instance();
        if (config.legacy_mode()) {
            // Break encapsulation for legacy compatibility
            return atoms_;
        }
        throw std::runtime_error("Direct access only in legacy mode");
    }
};
```

### 4. Algorithm Behavior

**Modern Mode**:
- Optimized algorithms
- May use different data structures
- May skip unnecessary steps

**Legacy Mode**:
- Exact legacy algorithm behavior
- Match legacy step-by-step
- Use same data structures as legacy

**Example**:
```cpp
class BasePairFinder {
    BasePair find_best_partner(const Residue& residue, const Structure& structure) {
        auto& config = ConfigManager::instance();
        
        if (config.legacy_mode()) {
            // Use exact legacy algorithm
            return find_best_partner_legacy(residue, structure);
        } else {
            // Use optimized modern algorithm
            return find_best_partner_modern(residue, structure);
        }
    }
    
private:
    BasePair find_best_partner_legacy(const Residue& residue, const Structure& structure) {
        // Match legacy iteration exactly
        // Use 1-based indexing
        // Match legacy tie-breaking exactly
        // Use legacy data structures
    }
    
    BasePair find_best_partner_modern(const Residue& residue, const Structure& structure) {
        // Optimized algorithm
        // Use 0-based indexing
        // Modern data structures
    }
};
```

### 5. Output Format

**Modern Mode**:
- Clean JSON structure
- May include additional metadata
- Optimized format

**Legacy Mode**:
- Exact legacy JSON format
- Match legacy field names exactly
- Match legacy field order exactly
- Match legacy data types exactly

**Implementation**:
```cpp
class JsonWriter {
    nlohmann::json to_json(const Structure& structure) {
        auto& config = ConfigManager::instance();
        
        if (config.legacy_mode()) {
            return to_json_legacy(structure);
        } else {
            return to_json_modern(structure);
        }
    }
    
private:
    nlohmann::json to_json_legacy(const Structure& structure) {
        // Match legacy JSON format exactly
        // Use 1-based indices
        // Match legacy field names
        // Match legacy field order
    }
    
    nlohmann::json to_json_modern(const Structure& structure) {
        // Modern JSON format
        // May include additional fields
        // Optimized structure
    }
};
```

## Implementation Strategy

### Phase 1: Add Legacy Mode Flag

1. Add `legacy_mode` to `ConfigManager`
2. Add `--legacy-mode` command-line flag
3. Pass flag through protocols

### Phase 2: Indexing Support

1. Ensure `get_residues_in_legacy_order()` works correctly
2. Add legacy index conversion utilities
3. Update algorithms to use legacy indices when in legacy mode

### Phase 3: Algorithm Compatibility

1. Add legacy-mode versions of key algorithms
2. Ensure iteration order matches legacy
3. Ensure tie-breaking matches legacy

### Phase 4: Output Compatibility

1. Ensure JSON output matches legacy exactly in legacy mode
2. Match field names and order
3. Match data types

## Code Organization

### Option 1: Separate Legacy Functions

```cpp
class BasePairFinder {
public:
    std::vector<BasePair> find_pairs(Structure& structure) {
        auto& config = ConfigManager::instance();
        if (config.legacy_mode()) {
            return find_pairs_legacy(structure);
        } else {
            return find_pairs_modern(structure);
        }
    }
    
private:
    std::vector<BasePair> find_pairs_legacy(Structure& structure);
    std::vector<BasePair> find_pairs_modern(Structure& structure);
};
```

**Pros**:
- Clear separation
- Easy to maintain both versions
- Can optimize modern version independently

**Cons**:
- Code duplication
- Must maintain two versions

### Option 2: Conditional Logic

```cpp
class BasePairFinder {
    std::vector<BasePair> find_pairs(Structure& structure) {
        auto& config = ConfigManager::instance();
        bool legacy_mode = config.legacy_mode();
        
        // Use conditional logic throughout
        std::vector<const Residue*> residues;
        if (legacy_mode) {
            residues = get_residues_in_legacy_order(structure);
        } else {
            residues = structure.all_residues();
        }
        
        // ... rest of algorithm with conditionals
    }
};
```

**Pros**:
- Single code path
- Less duplication
- Easier to keep in sync

**Cons**:
- More complex code
- Conditionals throughout
- Harder to optimize modern path

### Option 3: Strategy Pattern

```cpp
class PairFindingStrategy {
public:
    virtual ~PairFindingStrategy() = default;
    virtual std::vector<BasePair> find_pairs(Structure& structure) = 0;
};

class LegacyPairFindingStrategy : public PairFindingStrategy {
    std::vector<BasePair> find_pairs(Structure& structure) override;
};

class ModernPairFindingStrategy : public PairFindingStrategy {
    std::vector<BasePair> find_pairs(Structure& structure) override;
};

class BasePairFinder {
    std::unique_ptr<PairFindingStrategy> strategy_;
    
public:
    BasePairFinder(bool legacy_mode) {
        if (legacy_mode) {
            strategy_ = std::make_unique<LegacyPairFindingStrategy>();
        } else {
            strategy_ = std::make_unique<ModernPairFindingStrategy>();
        }
    }
    
    std::vector<BasePair> find_pairs(Structure& structure) {
        return strategy_->find_pairs(structure);
    }
};
```

**Pros**:
- Clean separation
- Easy to test
- Follows OOP principles

**Cons**:
- More classes
- More complex setup

**Recommendation**: Use **Option 1** (Separate Functions) for simplicity and clarity.

## Testing Strategy

### Legacy Mode Tests

```cpp
TEST(LegacyMode, FindPairsMatchesLegacy) {
    ConfigManager::instance().set_legacy_mode(true);
    
    Structure structure = load_pdb("data/pdb/1H4S.pdb");
    BasePairFinder finder;
    auto pairs = finder.find_pairs(structure);
    
    // Compare with legacy JSON
    auto legacy_json = load_json("data/json_legacy/find_bestpair_selection/1H4S.json");
    
    // Should match exactly
    ASSERT_EQ(pairs.size(), legacy_json.size());
    for (size_t i = 0; i < pairs.size(); i++) {
        EXPECT_EQ(pairs[i].residue_index1() + 1, legacy_json[i]["base_i"]);
        EXPECT_EQ(pairs[i].residue_index2() + 1, legacy_json[i]["base_j"]);
    }
}
```

### Modern Mode Tests

```cpp
TEST(ModernMode, FindPairsWorks) {
    ConfigManager::instance().set_legacy_mode(false);
    
    Structure structure = load_pdb("data/pdb/1H4S.pdb");
    BasePairFinder finder;
    auto pairs = finder.find_pairs(structure);
    
    // Should work correctly (may not match legacy exactly)
    ASSERT_GT(pairs.size(), 0);
    // ... test modern behavior
}
```

## Usage Examples

### Command Line

```bash
# Legacy mode - exact match with legacy
./find_pair_app --legacy-mode data/pdb/1H4S.pdb

# Modern mode - optimized (default)
./find_pair_app data/pdb/1H4S.pdb
```

### Programmatic

```cpp
auto& config = ConfigManager::instance();
config.set_legacy_mode(true);  // Enable legacy mode

FindPairProtocol protocol;
protocol.execute(structure);

// Output will match legacy exactly
```

## Benefits

1. **Regression Testing**: Easy to compare with legacy
2. **Debugging**: Can isolate differences
3. **Compatibility**: Ensures 100% match when needed
4. **Flexibility**: Modern mode can optimize independently

## Trade-offs

### What Breaks OOP

1. **Encapsulation**: May expose internal data in legacy mode
2. **Abstraction**: May use concrete types instead of interfaces
3. **Single Responsibility**: May mix legacy and modern concerns
4. **Open/Closed**: May need to modify classes for legacy mode

### Acceptable Because

1. **Testing**: Essential for regression testing
2. **Compatibility**: Required for exact matching
3. **Isolated**: Only affects behavior when flag is set
4. **Documented**: Clear that it breaks OOP for compatibility

## Implementation Checklist

- [ ] Add `legacy_mode` flag to `ConfigManager`
- [ ] Add `--legacy-mode` command-line flag
- [ ] Update `BasePairFinder` to support legacy mode
- [ ] Update `BasePairValidator` to support legacy mode
- [ ] Update `JsonWriter` to support legacy mode
- [ ] Update iteration order logic
- [ ] Update indexing logic
- [ ] Add legacy mode tests
- [ ] Document legacy mode behavior
- [ ] Update protocols to respect legacy mode

## Conclusion

The `--legacy-mode` flag is essential for:
- Ensuring 100% compatibility with legacy code
- Regression testing
- Debugging discrepancies
- Maintaining exact match capability

While it may break some OOP principles, this is acceptable because:
- It's opt-in (default is modern mode)
- It's isolated (only affects behavior when enabled)
- It's documented (clear purpose and trade-offs)
- It's essential (needed for compatibility verification)

**Recommendation**: Implement as separate functions (Option 1) for clarity and maintainability.

