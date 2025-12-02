# Comparison Index Rules - CRITICAL

**Date**: 2025-01-XX  
**Status**: ⚠️ **MANDATORY** - All comparisons MUST follow these rules

---

## ⚠️ CRITICAL RULE: Always Use Legacy Indices

**In ALL comparisons between legacy and modern code, you MUST use `legacy_residue_idx` (1-based).**

**NEVER use:**
- `residue_idx` (may be 0-based or different numbering)
- Simple counters
- Array positions
- Any other index system

---

## Why This Matters

Legacy code assigns residue indices sequentially (1-based) as residues are encountered during PDB parsing. Modern code may assign different indices due to:
- Different parsing order
- Different residue grouping
- Different filtering logic

**Using the wrong index means comparing completely different residues!**

---

## Rules for All Comparison Code

### Rule 1: JSON Field Names

**ALWAYS prefer `legacy_residue_idx` over `residue_idx`:**

```python
# ❌ WRONG
idx = record.get('residue_idx')

# ✅ CORRECT
idx = record.get('legacy_residue_idx') or record.get('residue_idx')
# But prefer legacy_residue_idx if both exist
```

**Better yet, be explicit:**
```python
# ✅ BEST
if 'legacy_residue_idx' in record:
    idx = record['legacy_residue_idx']
elif 'residue_idx' in record:
    idx = record['residue_idx']
else:
    raise ValueError(f"Record missing residue index: {record}")
```

### Rule 2: Command Line Arguments

**All comparison tools that take residue indices as arguments MUST document they expect legacy indices:**

```python
# ✅ CORRECT
parser.add_argument('residue_indices', nargs='+', type=int,
                   help='Legacy residue indices (1-based) to compare')
```

### Rule 3: Pair Comparisons

**When comparing pairs, use legacy indices from both records:**

```python
# ✅ CORRECT
r1 = record.get('legacy_residue_idx') or record.get('base_i')
r2 = record.get('legacy_residue_idx') or record.get('base_j')
# Or better:
r1 = record.get('base_i')  # base_i/base_j are already legacy indices
r2 = record.get('base_j')
```

### Rule 4: Frame Comparisons

**Frame records MUST use `legacy_residue_idx`:**

```python
# ✅ CORRECT
def load_frame_data(json_file):
    frames = {}
    for record in data:
        # Prefer legacy_residue_idx
        idx = record.get('legacy_residue_idx')
        if idx is None:
            idx = record.get('residue_idx')
        if idx is not None:
            frames[idx] = record
    return frames
```

### Rule 5: Validation Before Comparison

**Always validate that indices match before comparing:**

```python
# ✅ CORRECT
def compare_records(legacy_record, modern_record):
    leg_idx = legacy_record.get('legacy_residue_idx')
    mod_idx = modern_record.get('legacy_residue_idx')
    
    if leg_idx != mod_idx:
        raise ValueError(f"Index mismatch: legacy={leg_idx}, modern={mod_idx}")
    
    # Now safe to compare
    ...
```

---

## Code Patterns to Use

### Pattern 1: Safe Index Extraction

```python
def get_legacy_residue_idx(record: dict) -> Optional[int]:
    """Extract legacy residue index from record, preferring legacy_residue_idx."""
    # Prefer legacy_residue_idx
    if 'legacy_residue_idx' in record:
        return record['legacy_residue_idx']
    # Fallback to residue_idx (but log warning)
    if 'residue_idx' in record:
        print(f"WARNING: Using residue_idx instead of legacy_residue_idx")
        return record['residue_idx']
    # Try base_i/base_j for pair records
    if 'base_i' in record:
        return record['base_i']
    return None
```

### Pattern 2: Index Validation

```python
def validate_indices_match(legacy_record: dict, modern_record: dict) -> bool:
    """Validate that legacy and modern records refer to same residue."""
    leg_idx = get_legacy_residue_idx(legacy_record)
    mod_idx = get_legacy_residue_idx(modern_record)
    
    if leg_idx is None or mod_idx is None:
        return False
    
    if leg_idx != mod_idx:
        print(f"ERROR: Index mismatch - legacy={leg_idx}, modern={mod_idx}")
        return False
    
    return True
```

### Pattern 3: Pair Index Extraction

```python
def get_pair_indices(record: dict) -> Optional[Tuple[int, int]]:
    """Extract pair indices, preferring legacy indices."""
    # Try base_i/base_j first (these are legacy indices)
    if 'base_i' in record and 'base_j' in record:
        return (record['base_i'], record['base_j'])
    # Try residue1_idx/residue2_idx
    if 'residue1_idx' in record and 'residue2_idx' in record:
        r1 = record['residue1_idx']
        r2 = record['residue2_idx']
        # Validate these are legacy indices (should be > 0, typically < 10000)
        if r1 > 0 and r2 > 0:
            return (r1, r2)
    return None
```

---

## Files That Need Updates

### Python Scripts
- [ ] `scripts/compare_frames.py` - Line 36: Prefer `legacy_residue_idx`
- [ ] `scripts/compare_json.py` - Check all index usage
- [ ] `scripts/compare_base_pairs.py` - Verify uses legacy indices
- [ ] `scripts/compare_validation_geometry.py` - Verify uses legacy indices
- [ ] `scripts/compare_best_partner.py` - Verify uses legacy indices
- [ ] `scripts/compare_mutual_best.py` - Verify uses legacy indices
- [ ] `x3dna_json_compare/*.py` - All comparison modules

### C++ Tools
- [ ] `tools/compare_quality_scores.cpp` - Verify uses legacy indices
- [ ] `tools/compare_hbond_detection.cpp` - Verify uses legacy indices
- [ ] `tools/compare_validation_discrepancy.cpp` - Verify uses legacy indices
- [ ] All other `tools/compare_*.cpp` files

---

## Testing Checklist

Before using any comparison tool:

1. [ ] Verify modern JSON was generated with `--fix-indices`
2. [ ] Check that JSON records contain `legacy_residue_idx` field
3. [ ] Verify tool documentation states it uses legacy indices
4. [ ] Test with known pair to verify correct residues are compared

---

## Common Mistakes

### ❌ Mistake 1: Using residue_idx without checking

```python
# WRONG
idx = record['residue_idx']  # May not be legacy index!
```

### ❌ Mistake 2: Assuming array position equals index

```python
# WRONG
for i, record in enumerate(records):
    idx = i + 1  # WRONG! Array position != legacy index
```

### ❌ Mistake 3: Not validating indices match

```python
# WRONG
legacy_data[idx] vs modern_data[idx]  # May be different residues!
```

### ✅ Correct Approach

```python
# CORRECT
leg_idx = get_legacy_residue_idx(legacy_record)
mod_idx = get_legacy_residue_idx(modern_record)
if leg_idx == mod_idx:
    # Now safe to compare
    compare(legacy_record, modern_record)
```

---

## Verification Commands

### Check if JSON has legacy indices:

```bash
# Check frame_calc records
jq '.[] | select(.legacy_residue_idx != null) | .legacy_residue_idx' data/json/frame_calc/<PDB_ID>.json | head -5

# Check pair records
jq '.[] | select(.base_i != null) | {base_i, base_j}' data/json/find_bestpair_selection/<PDB_ID>.json | head -5
```

### Verify indices match between legacy and modern:

```bash
# Compare residue indices
python3 scripts/compare_json.py compare <PDB_ID> --record-type residue_indices
```

---

## Summary

**Golden Rule**: In ALL comparisons, use `legacy_residue_idx` (1-based) from JSON records. Never assume indices match or use array positions.

**When in doubt**: Extract `legacy_residue_idx` explicitly, validate it exists, and verify it matches between legacy and modern records before comparing.

---

*This document must be followed for all comparison code to ensure accurate results.*

