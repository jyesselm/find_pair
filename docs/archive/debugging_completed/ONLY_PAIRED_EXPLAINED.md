# --only-paired Mode: Detailed Explanation

**Date**: December 2, 2025  
**Purpose**: Understand what --only-paired does and its impacts

---

## What --only-paired Does (Step by Step)

### Normal Mode (WITHOUT --only-paired)

```
1. Parse PDB file → Get all atoms
2. Identify nucleotides → Find all residues with ring atoms
3. Calculate reference frames → ALL nucleotides (88 in 7EH2)
4. Find base pairs → Using the frames
5. Record everything to JSON:
   - base_frame_calc: 88 entries
   - base_pair: 24 pairs
6. Validate indices → Check all 88 against legacy
```

### --only-paired Mode

```
1. Parse PDB file → Get all atoms (same)
2. Identify nucleotides → Find all residues with ring atoms (same)
3. Calculate reference frames → ALL nucleotides (same as normal!)
4. Find base pairs → Using the frames (same)
5. **FILTER**: Identify which residues are in pairs
6. Record ONLY paired residues to JSON:
   - base_frame_calc: 48 entries (only paired)
   - base_pair: 24 pairs (same)
7. Validate indices → Check only 48 against legacy
```

### Key Difference

**Frames are CALCULATED for all nucleotides in both modes!**

The difference is **what gets RECORDED** to JSON:
- Normal mode: Records all 88
- --only-paired: Records only 48 (those in pairs)

---

## Code-Level Explanation

### Where It Applies (tools/generate_modern_json.cpp)

```cpp
// Step 1: Calculate frames for ALL nucleotides (lines 183-232)
if (only_paired) {
    calculator.calculate_all_frames(structure);  // ALL nucleotides!
    
    // Find pairs
    auto temp_pairs = temp_finder.find_pairs(structure);
    
    // Build set of paired residue indices
    paired_legacy_indices = {6995, 6996, ..., 7022};  // 48 indices
} else {
    calculator.calculate_all_frames(structure);  // Same!
}

// Step 2: Record frames to JSON (lines 256-392)
for (auto* residue : residues_in_order) {
    // ...
    
    // **THIS IS WHERE THE FILTER APPLIES**
    if (only_paired && !is_in_paired_set(residue)) {
        skipped_unpaired++;
        continue;  // Don't record to JSON
    }
    
    // Record to JSON (base_frame_calc, ls_fitting, etc.)
    writer.record_base_frame_calc(...);
}

// Step 3: Validation (lines 449-560)
for (auto& residue : structure) {
    // **FILTER APPLIES HERE TOO**
    if (only_paired && !is_in_paired_set(residue)) {
        continue;  // Don't add to ResidueTracker
    }
    
    tracker.add_residue(...);
}
```

---

## What It Affects

### ✅ AFFECTED: JSON Output

**Files that change**:
- `data/json/base_frame_calc/<PDB>.json` - Only paired residues
- `data/json/frame_calc/<PDB>.json` - Only paired residues
- `data/json/ls_fitting/<PDB>.json` - Only paired residues

**Files that DON'T change**:
- `data/json/base_pair/<PDB>.json` - Same (all pairs found)
- `data/json/hbond_list/<PDB>.json` - Same
- `data/json/pdb_atoms/<PDB>.json` - Same

### ✅ AFFECTED: Validation

**ResidueTracker**:
- Only tracks paired residues
- Validates only those against legacy
- Results in exact count match with legacy

### ❌ NOT AFFECTED: In-Memory Calculations

**Important**: The Structure object still has:
- ✅ All nucleotides identified
- ✅ All frames calculated
- ✅ All residues accessible in memory
- ✅ All base pairs found

**Only the JSON recording is filtered!**

---

## Potential Side Effects

### 1. Downstream Tools Reading JSON ⚠️

**Impact**: Tools that read `base_frame_calc/*.json` will only see paired residues

**Example**:
```python
# Reading modern JSON with --only-paired
with open('data/json/base_frame_calc/7EH2.json') as f:
    frames = json.load(f)

print(len(frames))  # 48, not 88!
```

**Solution**: 
- Document that JSON is filtered
- Or run without --only-paired for comprehensive analysis
- Or read from the Structure object directly (if using C++ API)

### 2. Frame Comparisons 📊

**Impact**: Can only compare paired residues

**Not an issue because**:
- Legacy also only has frames for paired residues
- Comparison scripts already handle this
- Unpaired residues have no legacy reference anyway

### 3. Residue Index Lookups 🔍

**Potential issue**:
```python
# If you expect residue at modern index 0 to be G:2 (first residue)
# But in --only-paired mode, index 0 might be G:15 (first PAIRED residue)
```

**Solution**: Always use **legacy_residue_idx** for lookups, not array indices

---

## Does It Affect Other Calculations?

### ❌ NO - Base Pair Finding

**Why**: Base pairs are found BEFORE filtering
- All nucleotides have frames
- All pairs are detected
- JSON records all pairs

**Result**: Same base pairs found regardless of mode

### ❌ NO - Geometric Parameters

**Why**: Parameters are calculated from base pairs
- Uses frames stored in BasePair objects
- Not affected by JSON filtering

### ❌ NO - H-Bond Detection  

**Why**: Uses residue objects directly
- Not dependent on JSON output
- Works with in-memory data

### ✅ YES - JSON-Based Comparisons

**Impact**: Comparison scripts reading JSON will only see paired residues

**Affected**:
- `scripts/compare_frames.py` - Will only compare paired
- `scripts/compare_ref_frames.py` - Will only see paired
- Any tool reading `base_frame_calc/*.json`

**Not affected**:
- Tools reading `base_pair/*.json` - Same data
- Tools using C++ API directly - Has all data
- Calculations done in-memory - Use all residues

---

## When to Use Each Mode

### Use --only-paired Mode

**Purpose**: Match legacy behavior exactly

**Use when**:
- ✅ Validating against legacy
- ✅ Comparing calculations with legacy
- ✅ Debugging differences from legacy
- ✅ Need exact legacy compatibility

**Why**:
- Ensures apple-to-apple comparison
- Same residues processed
- No confusion about missing data

### Use Normal Mode (No Flag)

**Purpose**: Comprehensive analysis

**Use when**:
- ✅ Analyzing all nucleotides (including unpaired)
- ✅ Studying loop regions
- ✅ Examining terminal overhangs  
- ✅ No legacy comparison needed
- ✅ Production analysis of new structures

**Why**:
- More complete information
- Can analyze unpaired regions
- Better for structural analysis

---

## Concrete Example: 7EH2

### The Structure

```
Chain G: TDADAADTDGDGDGDADGDCDTDGDTDCDADCDCDGDCDGDDA
         ^^^^^^^^^^^^^^^^__________^^^^^^^^^^^^^^^^
         unpaired (15)   paired (10) unpaired (1)
         
Chain H: DCDT...DCDT...DCDC...DGDCDCDCDTDADG
               ____________
               paired (12)
```

### What Each Mode Records

**Normal Mode**:
```json
{
  "base_frame_calc": [
    {"legacy_residue_idx": 6979, "chain": "G", "seq": 1},   // Unpaired
    {"legacy_residue_idx": 6980, "chain": "G", "seq": 2},   // Unpaired
    ...
    {"legacy_residue_idx": 6993, "chain": "G", "seq": 15},  // PAIRED ✓
    ...
    {"legacy_residue_idx": 7002, "chain": "G", "seq": 24},  // PAIRED ✓
    {"legacy_residue_idx": 7003, "chain": "G", "seq": 25},  // Unpaired
    ...
  ]  // Total: 88 entries
}
```

**--only-paired Mode**:
```json
{
  "base_frame_calc": [
    // Skips G:1-14 (unpaired)
    {"legacy_residue_idx": 6993, "chain": "G", "seq": 15},  // PAIRED ✓
    ...
    {"legacy_residue_idx": 7002, "chain": "G", "seq": 24},  // PAIRED ✓
    // Skips G:25 (unpaired)
    ...
  ]  // Total: 48 entries (matches legacy!)
}
```

---

## Impact on Existing Code

### Comparison Scripts

**Already handle this correctly**:

```python
# scripts/compare_frames.py
def compare_frames(legacy, modern):
    # Matches by legacy_residue_idx, not array index
    for leg_rec in legacy:
        leg_idx = leg_rec['legacy_residue_idx']
        
        # Find in modern by legacy_residue_idx
        mod_rec = find_by_index(modern, leg_idx)
        
        if mod_rec:
            compare(leg_rec, mod_rec)  # Works fine!
```

**Why it works**:
- Uses `legacy_residue_idx` for matching
- Doesn't assume sequential array indices
- Handles missing entries gracefully

### Potential Issues

**If code does this** ❌:
```python
# WRONG: Assumes array indices match
for i in range(len(legacy)):
    compare(legacy[i], modern[i])  # Breaks if different counts!
```

**Should do this** ✅:
```python
# CORRECT: Match by legacy_residue_idx
legacy_by_idx = {r['legacy_residue_idx']: r for r in legacy}
modern_by_idx = {r['legacy_residue_idx']: r for r in modern}

for idx in legacy_by_idx:
    if idx in modern_by_idx:
        compare(legacy_by_idx[idx], modern_by_idx[idx])
```

---

## Recommendation

### For Validation and Comparison

**Use --only-paired**:
- Ensures exact legacy match
- Simplifies debugging
- Validates correctly
- Required for validation scripts

### For Production Analysis

**Decide based on needs**:
- Need legacy comparison? → Use --only-paired
- Analyzing new structures? → Use normal mode
- Studying unpaired regions? → Use normal mode
- Maximum compatibility? → Use --only-paired

### Configuration

Could add to config file:
```yaml
# config/comparison.yaml
mode: only_paired  # or 'comprehensive'
reason: "Matching legacy for validation"
```

---

## Summary

### What --only-paired Does

✅ Calculates frames for ALL nucleotides (no change)  
✅ Finds ALL base pairs (no change)  
⚠️ Records ONLY paired residues to JSON (filtered)  
⚠️ Validates ONLY paired residues (filtered)  

### What It Doesn't Affect

✅ In-memory Structure object (has all data)  
✅ Base pair finding (finds all pairs)  
✅ Geometric calculations (uses all data)  
✅ H-bond detection (uses all data)  

### Potential Issues

⚠️ Downstream tools reading JSON will only see paired residues  
⚠️ Array index assumptions will break  
✅ Tools using legacy_residue_idx for matching will work fine  

### Bottom Line

**--only-paired is SAFE for validation and comparison**, and won't affect calculations. It only filters what gets written to JSON output to match legacy exactly.

For production use after validation, you can choose either mode based on your needs!

