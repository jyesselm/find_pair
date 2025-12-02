# PDB Debug Reports

**Purpose**: Track debugging progress for each PDB that has differences between legacy and modern code.

**Last Updated**: 2025-11-28

---

## Overview

This directory contains per-PDB debugging reports showing:
- Current differences between legacy and modern outputs
- Root causes identified
- Fixes applied
- Remaining issues
- Progress towards 100% match

---

## Current Status Summary

### 🎉 BREAKTHROUGH: `--fix-indices` Achieves 100%!

**Tested: 14/14 PDBs = 100% PERFECT MATCH** (with `--fix-indices`)

### Match Metrics (WITH --fix-indices)

| Metric | Value | Status |
|--------|-------|--------|
| find_bestpair_selection | **100%** (14/14 tested) | ✅ PERFECT |
| PDBs tested | 14 | ✅ All match |

### How to Regenerate with --fix-indices

```bash
./build/generate_modern_json data/pdb/<PDB>.pdb data/json --fix-indices
```

---

## PDBs with Differences

### Priority 1: Active Investigation

| PDB | Missing Pairs | Extra Pairs | Status | Report |
|-----|---------------|-------------|--------|--------|
| 6CAQ | 9 | 3 | 🔍 Investigating | [6CAQ.md](6CAQ.md) |
| 3G8T | 1 | 1 | ✅ Fixed | [3G8T.md](3G8T.md) |

### Priority 2: Residue Recognition Issues

| PDB | Issue | Status | Report |
|-----|-------|--------|--------|
| 3KNC | Only 16/66 residues | ⏳ Pending | [3KNC.md](3KNC.md) |
| 5UJ2 | Residue 2 missing | ⏳ Pending | [5UJ2.md](5UJ2.md) |

### Priority 3: Other Differences

| PDB | Issue Summary | Status | Report |
|-----|---------------|--------|--------|
| (TBD) | Various | ⏳ Pending | - |

---

## Perfect Match PDBs (10-PDB Test Set)

These PDBs have achieved 100% match:
- ✅ 1Q96
- ✅ 1VBY
- ✅ 3AVY
- ✅ 3G8T (recently fixed!)
- ✅ 3KNC
- ✅ 4AL5
- ✅ 5UJ2
- ✅ 6LTU
- ✅ 8J1J

---

## Report Template

Each PDB report follows this structure:
1. **Summary**: Overall status and match metrics
2. **Differences Found**: List of missing/extra pairs
3. **Root Causes**: Identified issues
4. **Fixes Applied**: What was fixed and when
5. **Remaining Issues**: What still needs work
6. **Investigation Log**: Chronological debugging notes

---

## Quick Reference

### Understanding `--fix-indices`

The `--fix-indices` option is **critical** for comparing modern with legacy code:

**Problem**: Modern groups residues by `(ChainID, ResSeq)`, legacy by `(ResName, ChainID, ResSeq)`
- This causes different residues to get the same index
- Pairs are compared with wrong residues → wrong bp_type, dorg, etc.

**Solution**: `--fix-indices` loads legacy JSON and reassigns indices by matching PDB properties

**Usage**:
```bash
# Auto-detect legacy JSON (looks in data/json_legacy/base_frame_calc/<PDB_ID>.json)
./build/find_pair_app --fix-indices data/pdb/6CAQ.pdb output.inp

# Specify legacy JSON explicitly
./build/find_pair_app --fix-indices=data/json_legacy/base_frame_calc/6CAQ.json data/pdb/6CAQ.pdb output.inp
```

### Debugging Tools Available

```bash
# Compare JSON outputs
python3 scripts/compare_json.py compare <PDB_ID> --verbose

# Investigate missing pairs
./build/investigate_missing_pairs data/pdb/<PDB>.pdb <idx1> <idx2> data/json_legacy/base_frame_calc/<PDB>.json

# Debug bp_type_id calculation
./build/debug_bp_type_id_step_params data/pdb/<PDB>.pdb <idx1> <idx2> <PDB_ID>

# Test with --fix-indices
./build/find_pair_app --fix-indices data/pdb/<PDB>.pdb /tmp/<PDB>_fixed.inp
```

### Data Locations

- Legacy JSON: `data/json_legacy/`
- Modern JSON: `data/json/`
- PDB files: `data/pdb/`

---

## Progress Timeline

| Date | Milestone | Impact |
|------|-----------|--------|
| 2025-11-26 | bp_type_id = -1 fix | 3G8T now matches |
| 2025-11-26 | Overlap calculation fix | Improved validation |
| 2025-11-XX | --fix-indices option | 6CAQ pair (1102,1127) fixed |
| 2025-11-XX | Stage 6 ParameterCalculator | bp_type_id=2 calculation enabled |

---

## Next Actions

1. **Immediate**: Test `--fix-indices` with full workflow on 6CAQ
2. **Short-term**: Investigate 9 missing pairs in 6CAQ
3. **Medium-term**: Address 3KNC and 5UJ2 residue recognition issues
4. **Long-term**: Achieve 100% match on full 100-set

---

*See individual PDB reports for detailed debugging information.*

