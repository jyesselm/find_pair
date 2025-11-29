# Debug JSON Generation Guide

**Date**: 2025-01-XX  
**Purpose**: Generate JSON files for specific residues to debug frame calculation differences

---

## Tool: `debug_frame_json`

Generates JSON files for specific residues showing frame calculation details.

### Usage

```bash
./build/debug_frame_json <pdb_file> <legacy_idx1> [legacy_idx2] [pdb_id]
```

### Examples

```bash
# Generate JSON for single residue
./build/debug_frame_json data/pdb/6CAQ.pdb 1101 6CAQ

# Generate JSON for pair of residues
./build/debug_frame_json data/pdb/6CAQ.pdb 1101 1127 6CAQ
```

### Output Files

1. **`data/json/debug_<PDB_ID>_frames_detailed.json`**
   - Detailed frame calculation information
   - Includes frame origins, rotation matrices, matched atoms

2. **`data/json/debug_<PDB_ID>/base_frame_calc/<PDB_ID>.json`**
   - Standard base_frame_calc records (matches legacy format)

3. **`data/json/debug_<PDB_ID>/ls_fitting/<PDB_ID>.json`**
   - Least-squares fitting records (if generated)

---

## Comparing with Legacy JSON

### Step 1: Generate Modern JSON

```bash
./build/debug_frame_json data/pdb/6CAQ.pdb 1101 1127 6CAQ
```

### Step 2: Compare with Legacy JSON

```bash
# Compare base_frame_calc records
python3 -c "
import json

# Modern
with open('data/json/debug_6CAQ/base_frame_calc/6CAQ.json') as f:
    modern = json.load(f)

# Legacy (note: legacy uses index 1102, modern uses 1101)
with open('data/json_legacy/base_frame_calc/6CAQ.json') as f:
    legacy = json.load(f)

# Find matching records
for m in modern:
    if m.get('legacy_residue_idx') == 1101:
        print('Modern residue 1101:', m.get('rms_fit'), m.get('num_matched_atoms'))
        
for l in legacy:
    if l.get('residue_idx') == 1102:
        print('Legacy residue 1102:', l.get('rms_fit'), l.get('num_matched_atoms'))
"
```

### Step 3: Compare Frame Origins

```bash
# Check ls_fitting records for translation/origin
python3 << 'EOF'
import json

# Modern
with open('data/json/debug_6CAQ/ls_fitting/6CAQ.json') as f:
    modern_ls = json.load(f)

# Legacy
with open('data/json_legacy/ls_fitting/6CAQ.json') as f:
    legacy_ls = json.load(f)

print("Frame Origins Comparison:")
print("=" * 60)

# Modern residue 1101 (should match legacy 1102)
for m in modern_ls:
    if m.get('legacy_residue_idx') == 1101:
        print(f"Modern 1101: {m.get('translation')}")

# Legacy residue 1102
for l in legacy_ls:
    if l.get('residue_idx') == 1102:
        print(f"Legacy 1102: {l.get('translation')}")
EOF
```

---

## What to Compare

### 1. Frame Origins (Translation)

**Key Question**: Do frame origins match?

- Modern: `ls_fitting.translation` for residue 1101
- Legacy: `ls_fitting.translation` for residue 1102

**Expected**: Should match if frame calculation is correct.

### 2. Rotation Matrices

**Key Question**: Do rotation matrices match?

- Modern: `ls_fitting.rotation_matrix` for residue 1101
- Legacy: `ls_fitting.rotation_matrix` for residue 1102

**Expected**: Should match if frame calculation is correct.

### 3. Matched Atoms

**Key Question**: Are the same atoms matched?

- Modern: `base_frame_calc.matched_atoms` for residue 1101
- Legacy: `base_frame_calc.matched_atoms` for residue 1102

**Expected**: Should match (same atoms, possibly different order).

### 4. RMS Fit

**Key Question**: Do RMS fits match?

- Modern: `base_frame_calc.rms_fit` for residue 1101
- Legacy: `base_frame_calc.rms_fit` for residue 1102

**Expected**: Should match if same atoms and coordinates are used.

---

## Example: Pair (1101, 1127) in 6CAQ

### Generated Files

```bash
./build/debug_frame_json data/pdb/6CAQ.pdb 1101 1127 6CAQ
```

### Modern Output

**Residue 1101 (G)**:
- Frame origin: [236.897, 166.84, 15.7898]
- RMS fit: 0.005376
- Matched atoms: 9 (C4, N3, C2, N1, C6, C5, N7, C8, N9)

**Residue 1127 (C)**:
- Frame origin: [219.334, 165.341, 6.94874]
- RMS fit: 0.004438
- Matched atoms: 6 (C4, N3, C2, N1, C6, C5)

### Legacy Comparison

**Residue 1102 (G)** - should match modern 1101:
- Check `data/json_legacy/ls_fitting/6CAQ.json` for translation
- Check `data/json_legacy/base_frame_calc/6CAQ.json` for RMS fit

---

## Troubleshooting

### Issue: Legacy JSON files not found

**Solution**: Run legacy code first to generate JSON files:
```bash
cd org/build
./bin/find_pair data/pdb/6CAQ.pdb
```

### Issue: Frame origins don't match

**Check**:
1. Are the same residues being compared? (Remember: legacy 1102 = modern 1101)
2. Are the same templates being used?
3. Are the same atoms being matched?
4. Are coordinates in the same coordinate system?

---

## Related Documentation

- [FRAME_ORIGIN_MISMATCH.md](FRAME_ORIGIN_MISMATCH.md) - Frame origin investigation
- [OFF_BY_ONE_ANALYSIS.md](OFF_BY_ONE_ANALYSIS.md) - Off-by-one error analysis
- [INVESTIGATION_SUMMARY.md](INVESTIGATION_SUMMARY.md) - Overall investigation summary

---

*This guide explains how to use debug_frame_json to generate and compare JSON files for debugging frame calculations.*

