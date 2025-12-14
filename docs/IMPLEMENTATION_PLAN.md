# Implementation Plan: Systematically Fix X3DNA Legacy/Modern Validation Differences

**Last Updated**: 2024-12-14
**Status**: Phase 1 Complete, Phase 2 In Progress

## Current Status Summary

### What Works (Completed)
- **Phase 1: Debug Infrastructure** - COMPLETE
  - `include/x3dna/debug/pair_validation_debugger.hpp` - Debug class for pair validation
  - `src/x3dna/debug/pair_validation_debugger.cpp` - Implementation
  - `tools/compare_single_pair.cpp` - Tool for comparing specific pairs

- **Phase 4: bp_type_id and WC Detection** - COMPLETE
  - Watson-Crick pair detection working correctly
  - -2.0 quality bonus applied when bp_type_id == 2

### Current Test Results (100-PDB Test Set)

| Stage | Pass Rate | Status |
|-------|-----------|--------|
| 1: Atoms | 94% (94/100) | 6 failures: 2CV2, 4RQF, 5VOE, 8UPT, 8Z1P, 8ZYC |
| 2-5: Frames | 100% | All pass |
| 6-7: Pairs | ~65% | 34 PDBs have extra pairs in modern |
| 8: H-bonds | 99% | 1 failure |
| 9: Base Pair | ~65% | Cascading from pair differences |
| 10: Selection | 99% | 1 failure |
| 11-12: Steps | ~5% | Cascading from pair differences |

### Key Finding
- Modern identifies MORE pairs than legacy (not fewer)
- 65/100 PDBs have matching pair_validation records
- 34/100 PDBs have extra pairs in modern that legacy doesn't have
- 0/100 PDBs are missing pairs that legacy has

## What Needs To Be Done

### Phase 2: Fix H-bond Counting (Current Focus)

**Problem**: Modern `count_simple()` is more lenient than legacy, causing extra pairs to be validated.

**Investigation Notes**:
- Example: 1ASY pair (29, 42) - C-C pair with O2-O2 H-bond at 3.76 Å
- Modern counts this as valid H-bond, legacy doesn't
- Need to trace exact difference in `is_baseatom()` + `good_hbatoms()` logic

**Files to investigate**:
- `src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_counter.cpp`
- `src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_utils.cpp`
- Legacy: `org/src/cmn_fncs.c` lines 4624-4636

**Debug command**:
```bash
./build/compare_single_pair data/pdb/1ASY.pdb 29 42 --json-dir data/json_legacy --verbose
```

### Phase 3: Fix Quality Score Calculation
- Verify MFACTOR (10000) scaling
- Match [2.5, 3.5] Å distance range for good H-bonds
- Verify -3.0 for 2+ good H-bonds, -1.0 for 1, 0.0 for none

### Phase 5: Fix is_nucleotide() Detection
- Compare ring atom matching
- Verify RMSD threshold (NT_CUTOFF = 0.2618)

### Phase 6: Fix Pair Validation Recording
- Ensure only valid pairs (bpid != 0) are recorded
- Ensure i < j ordering to avoid duplicates

### Phase 7: Fix Atom Indexing (6 failing PDBs)
- HETATM vs ATOM record handling
- alt_loc selection
- Insertion code handling

### Phase 8: Integration Testing
- Achieve 95%+ on all stages
- Document remaining edge cases

## Key Files

| File | Purpose |
|------|---------|
| `tools/compare_single_pair.cpp` | Compare specific pair between legacy/modern |
| `src/x3dna/algorithms/base_pair_validator.cpp` | H-bond counting, validation |
| `src/x3dna/algorithms/quality_score_calculator.cpp` | Quality score, bp_type_id |
| `src/x3dna/algorithms/hydrogen_bond/hydrogen_bond_counter.cpp` | Simple H-bond counting |

## Validation Commands

```bash
# Build
make release

# Test specific pair
./build/compare_single_pair data/pdb/1EHZ.pdb 1 72 --json-dir data/json_legacy --verbose

# Run stage validation
X3DNA=/Users/jyesselman2/local/installs/x3dna fp2-validate validate 6 --test-set 100

# Compare specific PDB
X3DNA=/Users/jyesselman2/local/installs/x3dna fp2-validate compare 1ASY --verbose
```

## User Preferences
- **Match Level**: 95%+ with documented exceptions
- **Debug Infrastructure**: Permanent debug layer via env var
- **Architecture**: Keep modern OOP design
