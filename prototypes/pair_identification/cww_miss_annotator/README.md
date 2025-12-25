# cWW Miss Annotator

Comprehensive diagnostic system for analyzing Watson-Crick base pair classification failures by combining H-bond analysis, geometric analysis, and DSSR reference comparison.

## Overview

This module provides a complete pipeline to diagnose why Watson-Crick (cWW) base pairs are misclassified:

1. **Data Loading**: Load DSSR pairs, slot-based H-bonds, and PDB structures
2. **H-bond Analysis**: Compare detected H-bonds against expected canonical patterns
3. **Geometric Analysis**: Compute RMSD to templates and identify outliers
4. **Miss Categorization**: Assign reason codes explaining classification failures
5. **Reporting**: Generate detailed JSON reports with diagnostics

## Files

- `annotator.py` - Main orchestrator combining all analyses
- `loaders.py` - Data loading (DSSR, slot H-bonds)
- `hbond_analyzer.py` - H-bond pattern analysis
- `geometric_analyzer.py` - RMSD-based geometric analysis
- `diagnostics.py` - Dataclasses for diagnostic information
- `report.py` - Report generation and formatting
- `test_annotator.py` - Test script
- `example_usage.py` - Demonstration script
- `__init__.py` - Package initialization

## Quick Start

```python
from pathlib import Path
from cww_miss_annotator.annotator import MissAnnotator

# Initialize annotator
annotator = MissAnnotator(
    pdb_dir=Path("data/pdb"),
    hbond_dir=Path("data/json/slot_hbonds"),
    dssr_dir=Path("data/json_dssr"),
    idealized_dir=Path("basepair-idealized"),
    exemplar_dir=Path("basepair-exemplars"),
)

# Annotate a PDB
report = annotator.annotate_pdb("1EHZ")

# Print summary
print(f"Total canonical cWW: {report.total_canonical_cww}")
print(f"True positives: {report.true_positives}")
print(f"False negatives: {report.false_negatives}")

# Examine false negatives
for ann in report.fn_annotations:
    print(f"\n{ann.res_id1} - {ann.res_id2} ({ann.sequence})")
    print(f"  Reasons: {', '.join(ann.reasons)}")
    print(f"  Our prediction: {ann.our_prediction}")
    print(f"  RMSD gap: {ann.geometric_diagnostics.rmsd_gap:.3f}")
```

## Running the Test

```bash
python prototypes/pair_identification/cww_miss_annotator/test_annotator.py
```

Example output:
```
Testing annotator on 100D...

Results for 100D:
  Total canonical cWW: 10
  True positives: 6
  False negatives: 4
  False positives: 0

False negative examples:

  1. A-DG-3 - B-DC-18 (GC)
     Saenger: 19-XIX
     Our prediction: tHS
     Reasons: rmsd_prefers_other
     H-bonds found: 3
     H-bonds missing: 0
     RMSD cWW: 0.123
     RMSD best: 0.017
```

## Reason Codes

The annotator categorizes misclassifications with the following reason codes:

### H-bond Based Reasons
- **no_hbonds**: No hydrogen bonds detected between the pair
- **missing_hbonds**: Some expected canonical H-bonds are missing
- **wrong_atoms**: H-bonds detected but to incorrect acceptor atoms
- **extra_hbonds**: Additional H-bonds beyond the expected pattern
- **distance_issues**: H-bond distances outside normal range (2.5-3.5Å)
- **overloaded_acceptor**: Acceptor atom receiving more than 2 H-bonds

### Geometric Reasons
- **rmsd_prefers_other**: Better RMSD fit to a non-cWW template (>0.1Å gap)
- **geometric_outlier**: Interbase angle >15° or N1-N9 distance abnormal

### Classification Reasons
- **non_canonical**: DSSR Saenger classification is "--" (non-canonical)

## Architecture

### MissAnnotator
Main orchestrator that coordinates all analyses:
- Loads DSSR pairs, slot H-bonds, and PDB structures
- Runs H-bond and geometric analysis for each pair
- Categorizes misses and generates comprehensive reports

### HBondAnalyzer
Analyzes hydrogen bonding patterns:
- Compares detected H-bonds against expected canonical patterns
- Identifies missing, extra, and incorrect H-bonds
- Validates H-bond distances and acceptor saturation
- Parses DSSR H-bond descriptions for comparison

### GeometricAnalyzer
Performs RMSD-based geometric analysis:
- Aligns residue pairs to idealized templates for all LW classes
- Computes RMSD to cWW template and best-fitting template
- Calculates interbase angle and N1-N9 distance
- Identifies geometric outliers

## Data Formats

### DSSR Format
- Nucleotide IDs: `"A.G1"`, `"0.C530"`, `"D.DG10"`
- Pair data includes: LW class, Saenger type, H-bond descriptions, angles, distances

### Slot Format
- Residue IDs: `"A-G-1"`, `"0-C-530"`, `"D-DG-10"`
- H-bond data includes: donor/acceptor atoms, distance, context, slot indices, alignment

## Key Functions

### `dssr_to_slot_id(dssr_nt: str) -> Optional[str]`
Convert DSSR nucleotide format to slot format.

### `slot_to_dssr_id(slot_id: str) -> str`
Convert slot format to DSSR nucleotide format.

### `load_dssr_pairs(dssr_path: Path, lw_filter: str = "cWW") -> Dict[Tuple[str, str], DSSRPair]`
Load DSSR pairs filtered by Leontis-Westhof class. Returns dict keyed by (slot_id1, slot_id2).

### `load_slot_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[SlotHBond]]`
Load slot-based H-bond data. Returns dict keyed by (res_id_i, res_id_j) with bidirectional lookup.

## Data Classes

### `DSSRPair`
```python
@dataclass
class DSSRPair:
    nt1: str              # First nucleotide (slot format)
    nt2: str              # Second nucleotide (slot format)
    lw_class: str         # Leontis-Westhof classification (e.g., "cWW")
    saenger: str          # Saenger classification (e.g., "19-XIX")
    hbonds_desc: str      # H-bond description from DSSR
    interbase_angle: float
    n1n9_distance: float
    bp: str               # Base pair sequence (e.g., "G-C")
```

### `SlotHBond`
```python
@dataclass
class SlotHBond:
    donor_atom: str
    acceptor_atom: str
    distance: float
    context: str                  # "base_base", "base_sugar", "sugar_sugar"
    h_slot: Optional[int]         # Hydrogen slot index (0 or 1)
    lp_slot: Optional[int]        # Lone pair slot index (0 or 1)
    alignment: Optional[float]    # Geometric alignment score
```

## Example Output

```
DSSR Pair: Q-C-0 - Q-G-22
  Base pair: C-G
  Saenger: 19-XIX
  DSSR H-bonds: O2(carbonyl)-N2(amino)[2.46],N3-N1(imino)[2.75],N4(amino)-O6(carbonyl)[3.02]
  Interbase angle: 25.2°
  N1-N9 distance: 8.73Å
  Slot H-bonds (3 total):
    N2 -> O2: 2.46Å (h_slot=1, lp_slot=0)
    N1 -> N3: 2.75Å (h_slot=0, lp_slot=0)
    N4 -> O6: 3.02Å (h_slot=1, lp_slot=0)
```

## Running the Example

```bash
python3 prototypes/pair_identification/cww_miss_annotator/example_usage.py
```

This will:
1. Load DSSR pairs and slot H-bonds for 1A9N
2. Show H-bond evidence for each DSSR cWW pair
3. Identify residue pairs with H-bonds not classified as cWW by DSSR
