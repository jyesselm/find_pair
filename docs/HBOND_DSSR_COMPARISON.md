# H-Bond Detection: DSSR Comparison Analysis

## Overview

This document summarizes the investigation into matching our H-bond detection output with DSSR (Dissecting the Spatial Structure of RNA), the gold standard for RNA structural analysis.

## Final Results

**Match Rate: 98.1%** (78,989 / 80,515 DSSR H-bonds matched across 95 structures)

| Metric | Value |
|--------|-------|
| Total DSSR H-bonds | 80,515 |
| Matched | 78,989 (98.1%) |
| DSSR-only (unmatched) | 1,526 |
| Valid DSSR-only (non-questionable) | 464 |

## Key Findings

### 1. Distance Cutoff Differences

DSSR uses a 3.5Å cutoff for most H-bond types. Our initial implementation used 4.0Å for legacy compatibility with base-pair detection, which caused over-detection.

**Solution**: Created `dssr_like()` preset with context-specific 3.5Å thresholds.

### 2. Missing Sugar-Backbone Context Classification

**Bug Found**: The `determine_nucleotide_context()` function in `geometry.cpp` was missing a case for sugar-backbone interactions (e.g., O2'-OP1, O2'-OP2 bonds).

These bonds returned `UNKNOWN` context and were filtered out when using the `RNA_INTERNAL` interaction filter.

**Impact**: 96% of missing RNA-only H-bonds (3,296 → 140) were recovered after this fix.

**Fix Location**: `src/x3dna/algorithms/hydrogen_bond/geometry.cpp:166-168`
```cpp
// Sugar-backbone interactions (e.g., O2'-OP1, O2'-OP2)
if ((a1_sugar && a2_backbone) || (a1_backbone && a2_sugar)) {
    return HBondContext::BASE_BACKBONE;
}
```

### 3. Interaction Filter Not Applied

**Bug Found**: The `interaction_filter` parameter existed in `HBondDetectionParams` but wasn't being checked during bond detection.

**Fix Location**: `src/x3dna/algorithms/hydrogen_bond/detector.cpp`
- Added `context_to_interaction_type()` helper function
- Added `passes_interaction_filter()` helper function
- Added filter check in `find_candidate_bonds()`

### 4. Protein H-Bond Support

DSSR reports H-bonds across all molecule types: RNA internal, protein-protein, and RNA-protein interactions.

**Solution**: Updated `dssr_like()` preset to use `HBondInteractionType::ANY` and added protein-specific distance thresholds.

## Remaining Differences (464 bonds)

Analysis of the 464 valid DSSR-only bonds that we don't detect:

### Distance Distribution
| Distance Range | Count | Percentage |
|----------------|-------|------------|
| ≤3.0Å | 10 | 2.2% |
| 3.0-3.2Å | 36 | 7.8% |
| 3.2-3.4Å | 90 | 19.4% |
| 3.4-3.5Å | 73 | 15.7% |
| >3.5Å | 255 | **55.0%** |

**Key Insight**: 55% of remaining unmatched bonds are beyond our 3.5Å cutoff. DSSR appears to use slightly relaxed thresholds in certain cases.

### Top Missing Atom Pairs
| Atom Pair | Count |
|-----------|-------|
| N7-OP2 | 27 |
| N7-O4' | 17 |
| O2-O2' | 16 |
| N2-OP1 | 14 |
| N6-OP2 | 13 |

### DSSR donAcc_type Distribution
| Type | Count | Percentage |
|------|-------|------------|
| acceptable | 211 | 45.5% |
| standard | 177 | 38.1% |
| questionable | 76 | 16.4% |

## Implementation Details

### HBondDetectionParams::dssr_like() Preset

```cpp
HBondDetectionParams HBondDetectionParams::dssr_like() {
    HBondDetectionParams params;
    // All contexts use 3.5Å cutoff
    params.distances.base_base_max = 3.5;
    params.distances.base_backbone_max = 3.5;
    params.distances.backbone_backbone_max = 3.5;
    params.distances.base_sugar_max = 3.5;
    params.distances.sugar_sugar_max = 3.5;
    params.distances.protein_mainchain_max = 3.5;
    params.distances.protein_sidechain_max = 3.5;
    params.distances.base_protein_max = 3.5;
    params.distances.protein_ligand_max = 3.5;
    params.distances.min_distance = 2.0;
    params.distances.conflict_filter_distance = 4.5;
    params.allowed_elements = ".O.N.";
    // Include ALL interactions
    params.interaction_filter = core::HBondInteractionType::ANY;
    return params;
}
```

### HBondContext Types

| Context | Description |
|---------|-------------|
| BASE_BASE | Between nucleobase atoms |
| BASE_BACKBONE | Between base and phosphate backbone |
| BACKBONE_BACKBONE | Between backbone atoms |
| BASE_SUGAR | Between base and sugar ring |
| SUGAR_SUGAR | Between sugar atoms |
| PROTEIN_MAINCHAIN | Between protein backbone atoms |
| PROTEIN_SIDECHAIN | Involving protein side chains |
| BASE_PROTEIN | Between nucleobase and protein |
| SUGAR_PROTEIN | Between sugar and protein |
| BACKBONE_PROTEIN | Between NA backbone and protein |
| BASE_LIGAND | Between base and ligand |
| PROTEIN_LIGAND | Between protein and ligand |
| LIGAND_LIGAND | Between ligand atoms |

### Interaction Filters

| Filter | Includes |
|--------|----------|
| BASE_BASE | Only base-base H-bonds |
| RNA_INTERNAL | All RNA internal (base, backbone, sugar combinations) |
| ANY | All H-bond types including protein |

## Validation Tools

### DSSR Comparison Script

```bash
# Compare single structure
python x3dna_json_compare/dssr_comparison.py 1EHZ

# Systematic analysis across multiple structures
python x3dna_json_compare/dssr_comparison.py --analyze
```

### Generate All H-Bonds JSON

```bash
./build/generate_modern_json data/pdb/1EHZ.pdb data/json --stage=all_hbonds
```

## Conclusions

1. **98.1% match rate** is excellent agreement with DSSR
2. Most remaining differences are due to DSSR using slightly longer distance cutoffs (>3.5Å) in certain cases
3. Two bugs were fixed during this investigation:
   - Missing sugar-backbone context classification
   - Interaction filter not being applied
4. The `dssr_like()` preset provides comprehensive H-bond detection matching DSSR behavior
