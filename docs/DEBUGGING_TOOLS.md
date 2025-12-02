# Debugging Tools Guide

**Created**: 2025-12-01  
**Purpose**: Guide to debugging tools for investigating legacy-modern code differences

---

## Available Tools

### 1. `compare_frames.py` - Frame Comparison Tool

**Purpose**: Compare reference frame origins and orientations between legacy and modern.

**Usage**:
```bash
# Compare specific residues
python3 scripts/compare_frames.py 9CF3 25 27

# Compare all residues
python3 scripts/compare_frames.py 9CF3 --all
```

**Output**: Shows frame origins, calculates distances between matching residues, highlights mismatches.

**Use When**:
- Investigating frame calculation differences
- Verifying frame origins match between legacy and modern
- Finding residues with incorrect frames

---

### 2. `analyze_validation_frame_bug.py` - Validation Frame Bug Detector

**Purpose**: Detect when validation uses wrong frames by comparing frame_calc frames vs validation dorg.

**Usage**:
```bash
python3 scripts/analyze_validation_frame_bug.py <PDB_ID> <res_i> <res_j>
python3 scripts/analyze_validation_frame_bug.py 9CF3 25 27
```

**Output**:
- Shows frame origins for both residues
- Compares validation dorg vs expected dorg from frame origins
- **Detects frame retrieval bugs** - when validation uses wrong frames

**Use When**:
- Validation geometry values seem wrong
- Validation rejects pairs that should be valid
- Investigating frame retrieval/setting bugs

**Example Output**:
```
🚨 FRAME RETRIEVAL BUG DETECTED!
   Validation uses dorg=18.490534
   But frame origins yield dorg=4.874563
   Difference: 13.615971 Å
```

---

### 3. `compare_validation_geometry.py` - Validation Geometry Comparison

**Purpose**: Side-by-side comparison of all validation geometry values and checks.

**Usage**:
```bash
python3 scripts/compare_validation_geometry.py <PDB_ID> <res_i> <res_j>
python3 scripts/compare_validation_geometry.py 9CF3 25 27
```

**Output**:
- All geometry values (dorg, d_v, plane_angle, dNN, quality_score)
- Direction vectors (dir_x, dir_y, dir_z)
- Validation status (is_valid, bp_type_id)
- Which validation checks pass/fail
- Differences highlighted

**Use When**:
- Comparing validation results for a specific pair
- Understanding why a pair is rejected
- Debugging validation threshold issues

---

### 4. `compare_best_partner.py` - Best Partner Candidate Comparison

**Purpose**: Compare best partner candidates for a specific residue.

**Usage**:
```bash
python3 scripts/compare_best_partner.py <PDB_ID> <res_i> [--verbose]
python3 scripts/compare_best_partner.py 9CF3 27 --verbose
```

**Output**:
- All candidate pairs
- Quality scores for each candidate
- Which candidate is selected as best
- Differences between legacy and modern

**Use When**:
- Investigating pair selection differences
- Understanding why a different partner is selected
- Comparing quality scores

---

### 5. `compare_mutual_best.py` - Mutual Best Decision Comparison

**Purpose**: Compare mutual best decisions.

**Usage**:
```bash
python3 scripts/compare_mutual_best.py <PDB_ID> [--verbose]
python3 scripts/compare_mutual_best.py 9CF3 --verbose
```

**Output**:
- Which pairs are mutual best
- Which pairs were selected
- Differences in mutual best decisions

---

### 6. `compare_iteration.py` - Iteration State Comparison

**Purpose**: Compare iteration states during pair selection.

**Usage**:
```bash
python3 scripts/compare_iteration.py <PDB_ID> [--verbose]
```

**Output**:
- Pairs found in each iteration
- Matched residues
- Iteration-by-iteration differences

---

## Debugging Workflows

### Workflow 1: Investigating Pair Selection Differences

```bash
# 1. Compare final pair selection
python3 scripts/compare_json.py compare <PDB_ID>

# 2. For a specific residue with different partner:
python3 scripts/compare_best_partner.py <PDB_ID> <res_i> --verbose

# 3. Check mutual best decisions
python3 scripts/compare_mutual_best.py <PDB_ID> --verbose
```

### Workflow 2: Investigating Validation Failures

```bash
# 1. Compare validation geometry
python3 scripts/compare_validation_geometry.py <PDB_ID> <res_i> <res_j>

# 2. Check for frame retrieval bug
python3 scripts/analyze_validation_frame_bug.py <PDB_ID> <res_i> <res_j>

# 3. Compare frames
python3 scripts/compare_frames.py <PDB_ID> <res_i> <res_j>
```

### Workflow 3: Investigating Frame Issues

```bash
# 1. Compare all frames
python3 scripts/compare_frames.py <PDB_ID> --all

# 2. For a specific pair with validation issues:
python3 scripts/analyze_validation_frame_bug.py <PDB_ID> <res_i> <res_j>

# 3. Compare validation geometry
python3 scripts/compare_validation_geometry.py <PDB_ID> <res_i> <res_j>
```

---

## Tool Status

| Tool | Status | Purpose |
|------|--------|---------|
| `compare_frames.py` | ✅ Created | Compare frame origins/orientations |
| `analyze_validation_frame_bug.py` | ✅ Created | Detect frame retrieval bugs |
| `compare_validation_geometry.py` | ✅ Created | Compare validation results |
| `compare_best_partner.py` | ✅ Exists | Compare best partner candidates |
| `compare_mutual_best.py` | ✅ Exists | Compare mutual best decisions |
| `compare_iteration.py` | ✅ Exists | Compare iteration states |
| `debug_ref_frames.py` | ✅ Exists | Debug ref_frames.dat output |
| `compare_json.py` | ✅ Exists | Comprehensive JSON comparison |

---

## Future Tools (Ideas)

1. **Frame Orientation Comparison**: Compare rotation matrices, not just origins
2. **Residue Mapping Inspector**: Verify legacy index to residue mapping
3. **H-Bond Visualizer**: Show which H-bonds are detected
4. **Quality Score Breakdown**: Show components of quality score calculation
5. **Automated Bug Reporter**: Auto-detect and report common issues

---

## Related Documentation

- `docs/DEBUGGING_WORKFLOW.md` - Step-by-step debugging workflow
- `docs/LEGACY_INDICES_GUIDE.md` - Understanding legacy indices
- `docs/TESTING_GUIDE.md` - Testing and comparison guide

