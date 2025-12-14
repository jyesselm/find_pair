# Implementation Plan: Systematically Fix X3DNA Legacy/Modern Validation Differences

**Last Updated**: 2025-12-14
**Status**: Phase 5 Complete - Stages 6-12 at 99-100%

## Current Status Summary

### What Works (Completed)
- **Phase 1: Debug Infrastructure** - COMPLETE
  - `include/x3dna/debug/pair_validation_debugger.hpp` - Debug class for pair validation
  - `src/x3dna/debug/pair_validation_debugger.cpp` - Implementation
  - `tools/compare_single_pair.cpp` - Tool for comparing specific pairs

- **Phase 2: Fix pair_validation Recording Timing** - COMPLETE
  - **Root Cause Found**: Modern was recording pair_validation during Phase 1 (all pairs validation),
    while legacy records during Phase 2 (greedy selection)
  - **Fix Applied**: Code already fixed to record during find_best_partner, matching legacy behavior
  - **Result**: Stages 6-7 now at 100%

- **Phase 3: Fix hbond_list and base_pair Recording** - COMPLETE
  - **Key Finding**: Legacy pair_validation and base_pair have IDENTICAL pairs
  - base_pair recording moved to greedy selection (same location as pair_validation)
  - hbond_list recording happens during greedy selection validation
  - **Result**: Stages 8-9 now at 100%

- **Phase 4: bp_type_id and WC Detection** - COMPLETE
  - Watson-Crick pair detection working correctly
  - -2.0 quality bonus applied when bp_type_id == 2

- **Phase 5: Fix Steps/Helical Comparison** - COMPLETE
  - **Key Finding**: Legacy and modern compute steps in different orders due to five2three algorithm
  - **Solution**: mst_org position matching in `x3dna_json_compare/step_comparison.py`
  - Steps matched by midstep origin position (0.01 Å threshold)
  - Sign-inverted parameter matches accepted (for steps computed A→B vs B→A)
  - **Result**: Stages 11-12 now at 100% with tests_python validation

### Current Test Results

#### 100-PDB Test Set
| Stage | Pass Rate | Status |
|-------|-----------|--------|
| 1: Atoms | 94% (94/100) | 6 failures: 2CV2, 4RQF, 5VOE, 8UPT, 8Z1P, 8ZYC |
| 2-5: Frames | 99% (99/100) | 1 failure: 6OZK (inosine recognition) |
| 6: pair_validation | **100%** | ✅ Fixed |
| 7: distance_checks | **100%** | ✅ Fixed |
| 8: hbond_list | **100%** | ✅ Fixed |
| 9: base_pair | **100%** | ✅ Fixed |
| 10: Selection | **100%** | ✅ Fixed |
| 11: Steps | **100%** | ✅ Fixed (via position matching) |
| 12: Helical | **100%** | ✅ Fixed (via position matching) |

#### Full fast_pdbs Dataset (3602 PDBs)
| Stage Group | Pass Rate | Status |
|-------------|-----------|--------|
| 6-9: Pairs (pair_validation, distance_checks, hbond_list, base_pair) | **99.9%** (3598/3602) | ✅ 4 edge cases |
| 11: Steps | **99.8%** (3594/3602) | ✅ 8 tolerance boundary cases |
| 12: Helical | **99.97%** (3601/3602) | ✅ 1 edge case (8H0I) |

Failed PDBs for pairs (4): 2XD0, 4E8R, 8ANE, 9D5J - minor numerical differences
Failed PDBs for steps (8): 1KOG, 4YVI, 7AFA, 7Z4I, 8D29, + 3 others - tolerance boundary (diff ~0.0001-0.0004)
Failed PDB for helical (1): 8H0I - rise=null due to extreme geometry (shift=-48Å, twist=180°)

### Key Findings from Phase 3 Investigation

1. **pair_validation Recording Issue (SOLVED)**
   - The issue was NOT H-bond counting logic
   - The issue was WHEN pair_validation was recorded
   - Legacy records during greedy selection (best_pair loop), modern was recording during Phase 1
   - After regenerating JSON with correct recording timing, stages 6-7 are 100%

2. **hbond_list and base_pair Recording Issue (SOLVED)**
   - Key insight: Legacy pair_validation and base_pair have IDENTICAL pairs
   - Legacy records both during greedy selection (find_best_partner)
   - hbond_list also recorded during validation phase via get_hbond_ij
   - Fix: Record base_pair at same location as pair_validation in modern code
   - Result: Stages 8-9 now at 100%

3. **Multiplet Pairs Clarification**
   - Some pairs appear in hbond_list but not in base_pair (e.g., 1QRU pair 37,38)
   - These come from "network mode" validation which records hbond_list but not base_pair
   - The comparison normalizes pair keys so (37,38) and (38,37) are considered identical

## What Needs To Be Done

### Phase 6: Fix is_nucleotide() Detection (Stage 2-5)
- Compare ring atom matching
- Verify RMSD threshold (NT_CUTOFF = 0.2618)
- Affected: 6OZK (inosine recognition issue)

### Phase 7: Fix Atom Indexing (Stage 1, 6 failing PDBs)
- HETATM vs ATOM record handling
- alt_loc selection
- Insertion code handling
- Affected PDBs: 2CV2, 4RQF, 5VOE, 8UPT, 8Z1P, 8ZYC

### Phase 8: Integration Testing
- Already achieved 99%+ on all stages!
- Remaining edge cases documented above

## Key Files

| File | Purpose |
|------|---------|
| `tools/compare_single_pair.cpp` | Compare specific pair between legacy/modern |
| `src/x3dna/algorithms/base_pair_finder.cpp` | Main pair finding logic, greedy selection |
| `src/x3dna/algorithms/base_pair_validator.cpp` | H-bond counting, validation |
| `src/x3dna/algorithms/quality_score_calculator.cpp` | Quality score, bp_type_id |

## Validation Commands

```bash
# Build
make release

# Regenerate JSON for test set
./build/generate_modern_json --pdb-list=data/test_set_100_list.txt --pdb-dir=data/pdb data/json --stage=validation

# Run stage validation
X3DNA=/Users/jyesselman2/local/installs/x3dna fp2-validate validate 6 --test-set 100
X3DNA=/Users/jyesselman2/local/installs/x3dna fp2-validate validate pairs --test-set 100

# Compare specific PDB
X3DNA=/Users/jyesselman2/local/installs/x3dna fp2-validate compare 1QRU --verbose
```

## User Preferences
- **Match Level**: 95%+ with documented exceptions
- **Debug Infrastructure**: Permanent debug layer via env var
- **Architecture**: Keep modern OOP design
