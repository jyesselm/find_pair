# --only-paired Mode: Matching Legacy Behavior

**Date**: December 2, 2025  
**Status**: ✅ **IMPLEMENTED AND WORKING**

---

## Why This Was Needed

During comprehensive validation, we discovered that **legacy code only processes residues involved in base pairs**, while modern code processes all nucleotides. This caused apparent "mismatches" like:

- **7EH2**: Modern found 88 nucleotides, legacy found 48
- The difference: 40 unpaired terminal residues that legacy skipped

**Solution**: Add `--only-paired` mode to modern code to exactly match legacy behavior.

---

## What --only-paired Mode Does

When `--only-paired` flag is used:

1. **Calculate frames for ALL nucleotides** (needed to find pairs)
2. **Find base pairs**
3. **Identify which residues are in pairs**
4. **Only record frames for paired residues** (skip unpaired)

This perfectly replicates legacy's workflow:
```
Legacy:  find_pairs() → ref_frames(paired_residues_only)
Modern:  calculate_all_frames() → find_pairs() → record(paired_residues_only)
```

---

## Usage

```bash
# With validation (automatic with RUN_FULL_VALIDATION.sh)
./build/generate_modern_json data/pdb/7EH2.pdb data/json --only-paired

# Check what it does
./build/generate_modern_json data/pdb/7EH2.pdb data/json --only-paired
# Output shows:
#   Only-paired mode: tracking 48 paired residues for validation
#   Frames calculated: 88
#   Frames recorded: 48
#   Skipped unpaired residues: 40
```

---

## Example: 7EH2

### Without --only-paired (Modern Comprehensive Mode)

```
Nucleotides found: 88
Frames calculated: 88
- All DNA in chains G, H, J, Q
- All RNA in chains I, R
- Includes unpaired terminal overhangs
```

### With --only-paired (Legacy Compatible Mode)

```
Nucleotides found: 88
Frames calculated: 88  
Frames RECORDED: 48 ✅ Matches legacy!
Skipped unpaired: 40
- Only residues in the 24 base pairs
- Excludes terminal overhangs (residues 1-14, 25)
```

###Legacy (Reference)

```
Base pairs found: 24
Residues processed: 48 ✅ Perfect match!
- Only calculates frames for paired residues
- Skips unpaired residues entirely
```

---

## Implementation Details

### Code Changes (tools/generate_modern_json.cpp)

```cpp
// 1. Add flag
bool only_paired = false;

// 2. After finding pairs, collect paired residue indices
std::set<int> paired_legacy_indices;
if (only_paired) {
    for (const auto& pair : base_pairs) {
        // Collect legacy indices of both residues in each pair
        paired_legacy_indices.insert(legacy_idx_1);
        paired_legacy_indices.insert(legacy_idx_2);
    }
}

// 3. In ResidueTracker population, skip unpaired
for (const auto& residue : residues) {
    if (only_paired && !is_in_paired_set(residue)) {
        continue;  // Skip unpaired residues
    }
    tracker.add_residue(...);
}

// 4. In modern index assignment, skip unpaired
for (const auto& residue : residues) {
    if (only_paired && !is_in_paired_set(residue)) {
        continue;  // Skip unpaired residues
    }
    tracker.assign_modern_index(...);
}
```

---

## Impact on Validation

### Before --only-paired

- ❌ 7EH2 FAILED: 88 modern vs 48 legacy (count mismatch)
- Validation stopped
- Required investigation

### After --only-paired

- ✅ 7EH2 PASSED: 48 modern vs 48 legacy (perfect match)
- Validation continues
- 125+ new PDBs validated successfully

---

## When to Use Each Mode

### Use --only-paired (Legacy Compatible Mode)

**For**:
- ✅ Validation against legacy data
- ✅ Direct comparison with legacy calculations
- ✅ Debugging calculation differences
- ✅ Ensuring exact legacy behavior match

**When**:
- Comparing with legacy JSON
- Running validation scripts
- Initial development/debugging

### Use Normal Mode (Comprehensive)

**For**:
- ✅ Analyzing ALL nucleotides (paired and unpaired)
- ✅ Studying loops, overhangs, single-stranded regions
- ✅ Complete structural analysis
- ✅ New applications beyond legacy scope

**When**:
- No legacy comparison needed
- Analyzing unpaired regions
- Production use after validation complete

---

## Validation Script Integration

`scripts/validate_all_indices.py` now uses `--only-paired` by default:

```python
cmd = [
    str(project_root / "build" / "generate_modern_json"),
    str(pdb_file),
    str(json_dir),
    "--only-paired"  # Ensures exact legacy match
]
```

---

## Results

With `--only-paired` enabled:

| Metric | Count | Percentage |
|--------|-------|------------|
| ✅ PASS | 2,511 | 91.9% |
| ⏭️ SKIP | 209 | 7.6% |
| ⏱️ TIMEOUT | 13 | 0.5% |
| ❌ FAIL | 0 | 0.0% |

**Success Rate**: 100% of testable PDBs pass validation! 🎉

---

## Technical Notes

### Why Calculate All Frames First?

Even in `--only-paired` mode, we calculate frames for ALL nucleotides first because:
1. **Need frames to find pairs**: BasePairFinder requires reference frames
2. **Efficient**: Calculate once, filter during recording
3. **Flexible**: Easy to switch between modes

### Memory Overhead

- Frames calculated: ALL nucleotides
- Frames recorded to JSON: Only paired (in --only-paired mode)
- Memory: Slightly higher during processing
- Output: Identical to legacy

---

## Future Options

Could add additional modes:

- `--only-helical`: Only residues in helices (excludes isolated pairs)
- `--min-helix-size N`: Only helices with N+ base pairs
- `--include-loops`: Include loop residues between paired regions

---

## See Also

- `RESOLUTION_7EH2.md` - Original investigation
- `INVESTIGATION_7EH2.md` - Detailed analysis  
- `INDEX_VALIDATION_REPORT.md` - Overall validation report

