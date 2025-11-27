# Legacy Mode Requirement

**Status**: Design Document  
**Purpose**: Document the `--legacy-mode` flag requirement for exact legacy compatibility

## Requirement

Add a `--legacy-mode` flag that allows direct comparison to legacy code, even if it breaks some OOP principles. Everything else should be modern.

## Rationale

1. **Regression Testing**: Essential for comparing modern code with legacy
2. **Debugging**: Helps isolate differences when they occur
3. **Compatibility**: Ensures 100% match capability when needed
4. **Flexibility**: Modern mode can optimize independently

## Design

See `docs/LEGACY_MODE_DESIGN.md` for complete design document.

## Key Points

### What Legacy Mode Does

1. **Indexing**: Uses 1-based indexing where legacy does
2. **Iteration Order**: Matches legacy iteration order exactly
3. **Algorithm Behavior**: Matches legacy step-by-step
4. **Output Format**: Matches legacy JSON format exactly
5. **Data Structures**: May break encapsulation for compatibility

### What Breaks OOP

- **Encapsulation**: May expose internal data
- **Abstraction**: May use concrete types
- **Single Responsibility**: May mix concerns
- **Open/Closed**: May need modifications

### Why It's Acceptable

- **Opt-in**: Default is modern mode
- **Isolated**: Only affects behavior when enabled
- **Documented**: Clear purpose and trade-offs
- **Essential**: Needed for compatibility verification

## Implementation

### ConfigManager

```cpp
class ConfigManager {
    bool legacy_mode() const { return legacy_mode_; }
    void set_legacy_mode(bool value) { legacy_mode_ = value; }
    
private:
    bool legacy_mode_ = false;
};
```

### Command Line

```bash
# Modern mode (default)
./find_pair_app data/pdb/1H4S.pdb

# Legacy mode
./find_pair_app --legacy-mode data/pdb/1H4S.pdb
```

### Protocols

```cpp
class FindPairProtocol {
    void set_legacy_mode(bool value) { legacy_mode_ = value; }
    
private:
    bool legacy_mode_ = false;
};
```

## Integration Points

1. **Stage 7 (Protocols)**: Add legacy mode support
2. **Stage 8 (Applications)**: Add `--legacy-mode` flag parsing
3. **Algorithms**: Add legacy-mode versions where needed
4. **JSON Output**: Support legacy format in legacy mode

## Testing

```cpp
TEST(LegacyMode, ExactMatch) {
    ConfigManager::instance().set_legacy_mode(true);
    // ... test that output matches legacy exactly
}
```

## Status

- ✅ Design documented (`docs/LEGACY_MODE_DESIGN.md`)
- ✅ Modernization plan updated
- ⏳ Implementation pending (Stage 7)

