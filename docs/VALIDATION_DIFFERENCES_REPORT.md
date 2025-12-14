# Validation Differences Report

**Date**: 2025-12-14 (Updated)
**Test Set**: 100 PDBs from `data/valid_pdbs_fast.json`

## Executive Summary

| Stage | Name | Pass Rate | Notes |
|-------|------|-----------|-------|
| 1 | pdb_atoms | 94% (94/100) | 6 failures due to atom indexing differences |
| 2 | residue_indices | 100% (100/100) | All pass |
| 3 | base_frame_calc | 100% (100/100) | All pass |
| 4 | ls_fitting | 100% (100/100) | All pass |
| 5 | frame_calc | 100% (100/100) | All pass |
| 6 | pair_validation | 45% (45/100) | Different pair selection between legacy/modern |
| 7 | distance_checks | 45% (45/100) | Same pairs as stage 6 |
| 8 | hbond_list | 99% (99/100) | 1 failure |
| 9 | base_pair | 34% (34/100) | Different pairs pass through to base_pair output |
| 10 | find_bestpair_selection | 99% (99/100) | 1 failure |
| 11 | bpstep_params | 2% (2/100) | Cascading from pair selection differences |
| 12 | helical_params | 1% (1/100) | Cascading from pair selection differences |

## Critical Finding

**The main issue is different pair selection between legacy and modern code.** The greedy selection algorithm produces different results due to:
1. Different iteration order when identifying candidate pairs
2. Different tie-breaking behavior when multiple pairs have similar quality scores
3. Legacy code may identify additional pairs that modern code doesn't

---

## Stage 1: Atoms (94% pass)

### Failing PDBs
| PDB | Legacy Records | Modern Records | Issues |
|-----|----------------|----------------|--------|
| 2CV2 | 11,215 | 11,215 | 40 missing, 40 mismatched |
| 4RQF | 8,619 | 8,619 | 363 missing, 37 extra, 63 mismatched |
| 5VOE | 3,129 | 3,129 | 49 missing, 12 extra |
| 8UPT | - | - | Atom indexing differences |
| 8Z1P | - | - | Atom indexing differences |
| 8ZYC | - | - | Atom indexing differences |

### Root Cause
The atom failures appear to be related to atom indexing differences between legacy and modern code. The total counts match, but specific atom-to-index mappings differ in some cases.

---

## Stage 2: Residue Indices (100% pass)

All 100 PDBs pass residue index validation.

---

## Stages 3-5: Frames (100% pass)

All 100 PDBs pass frame calculation validation (base_frame_calc, ls_fitting, frame_calc).

---

## Stages 6-7: Pairs (45% pass)

### Failure Categories

#### 1. Different pair candidate identification
Modern and legacy code identify different candidate pairs due to:
- Different distance threshold implementations
- Different geometric check ordering

#### 2. Different final pair selection
The greedy selection algorithm is sensitive to:
- Iteration order through candidate pairs
- Quality score calculation differences
- Tie-breaking behavior

### Example: 1EHZ
- Legacy identifies more base pairs than modern
- Legacy base_pair: 35 records
- Modern base_pair: 30 records
- 5 pairs in legacy are not found in modern

---

## Stage 8: H-bond List (99% pass)

Only 1 failure in 100 PDBs. H-bond identification is relatively consistent between legacy and modern.

---

## Stage 9: Base Pair (34% pass)

### Root Cause
Stage 9 has the lowest pass rate (34%) because it checks the final base_pair records output. Differences here indicate:
1. Different pairs were selected in earlier stages
2. Different pairs pass the validation checks
3. Legacy outputs more base_pair records than modern for many structures

---

## Stage 10: Find Bestpair Selection (99% pass)

Only 1 failure. The final pair selection is mostly consistent once the candidate pairs are the same.

---

## Stages 11-12: Steps/Helical (1-2% pass)

### Root Cause
Step parameters are calculated between consecutive base pairs. If different pairs are selected:
1. Different steps are calculated
2. Parameters for those steps differ
3. Cascading failures from pair selection differences

### Regeneration Required
To properly compare step params:
1. First fix pair selection differences (stages 6-10)
2. Regenerate modern step params after pair selection is fixed

---

## File Counts Summary

| JSON Type | Modern | Legacy |
|-----------|--------|--------|
| base_pair | 3,707 | 4,128 |
| find_bestpair_selection | 4,126 | 3,921 |
| pair_validation | 4,124 | 3,606 |
| distance_checks | 4,126 | 3,607 |
| hbond_list | 4,126 | 3,607 |

Note: Legacy has MORE base_pair files (4,128 vs 3,707) but FEWER of the other file types. This suggests legacy code outputs more base_pair records per structure.

---

## Action Items

### High Priority

1. **Investigate pair selection algorithm differences** - The 45% pass rate on stages 6-7 indicates fundamental differences in how pairs are selected
2. **Fix base_pair output differences** - The 34% pass rate suggests the pair validation/selection flow differs between legacy and modern

### Medium Priority

3. **Investigate atom indexing differences** - 6 PDBs fail stage 1 due to atom index mismatches
4. **Regenerate step params** - After fixing pair selection, regenerate modern step params

### Low Priority

5. **Document known representation differences** - Some differences may be intentional design changes

---

## Validation Commands Reference

```bash
# Environment setup
export X3DNA=/Users/jyesselman2/local/installs/x3dna

# Validate individual stages
fp2-validate validate 1 --test-set 100        # Atoms
fp2-validate validate 2 --test-set 100        # Residue indices
fp2-validate validate frames --test-set 100   # Frames (3-5)
fp2-validate validate pairs --test-set 100    # Pairs (6-10)
fp2-validate validate steps --test-set 100    # Steps (11-12)

# Validate single PDB with verbose output
fp2-validate compare 1EHZ --verbose

# Validate all stages
fp2-validate validate all --test-set 100
```
