# cWW Classification Reason Codes

This document explains the reason codes used by the cWW miss annotator to diagnose why a Watson-Crick base pair might be misclassified.

## Overview

When analyzing a base pair that DSSR classifies as cWW (cis Watson-Crick/Watson-Crick), the annotator checks multiple criteria. If the pair fails any criteria AND has a BP score < 0.70, it is flagged with one or more reason codes explaining the issue.

**BP Score Threshold**: Pairs with BP score >= 0.70 (grade C or better) are considered good cWW pairs and are NOT flagged, regardless of minor issues.

---

## Reason Codes Summary

| Category | Reason | Threshold | Description |
|----------|--------|-----------|-------------|
| H-bond | `no_hbonds` | - | No H-bonds detected between bases |
| H-bond | `missing_hbonds` | - | Some expected H-bonds not found |
| H-bond | `wrong_hbonds` | - | H-bonds to incorrect acceptor atoms |
| H-bond | `extra_hbonds` | - | Unexpected H-bonds beyond canonical |
| H-bond | `long_hbonds` | > 4.0Å | Extremely stretched H-bonds |
| H-bond | `short_hbonds` | < 2.0Å | Clash-like distances |
| H-bond | `overloaded_acceptor` | > 2 | Too many H-bonds to one acceptor |
| Geometry | `poor_planarity` | angle ≥ 20° | Bases not coplanar |
| Geometry | `geometric_outlier` | varies | Parameters outside normal range |
| Geometry | `rmsd_prefers_other` | gap > 0.5Å | Different LW class fits better |
| Classification | `non_canonical` | - | DSSR Saenger is "--" |
| Classification | `dssr_questionable` | RMSD>1.0 + no H-bonds | Likely DSSR error |

---

## H-Bond Related Reasons

### `no_hbonds`
**Description**: No hydrogen bonds detected between the bases.

**Interpretation**: The H-bond detection algorithm found no valid H-bonds. This could mean:
- All expected H-bond distances are beyond the detection threshold (~3.5Å)
- The bases are not positioned to form H-bonds
- The pair geometry is incompatible with cWW

---

### `missing_hbonds`
**Description**: Some expected H-bonds are not found.

**Expected H-bonds for cWW**:
- **G-C / C-G**: 3 H-bonds (N1-N3, N2-O2, N4-O6)
- **A-U / U-A**: 2 H-bonds (N6-O4, N3-N1)

**Interpretation**: The pair has some but not all expected H-bonds. This could indicate a partial opening of the pair or alternative H-bonding pattern.

---

### `wrong_hbonds`
**Description**: H-bonds are formed to incorrect acceptor atoms.

**Example**: G.N2 donating to C.N3 instead of C.O2

**Interpretation**: The bases are forming H-bonds but not in the canonical Watson-Crick pattern. This could indicate:
- A different base pair type (wobble, Hoogsteen)
- Structural distortion causing non-canonical contacts
- Modified bases with altered H-bonding capability

---

### `extra_hbonds`
**Description**: Unexpected H-bonds detected beyond the canonical pattern.

**Interpretation**: Additional H-bonding interactions beyond Watson-Crick. Could indicate bifurcated H-bonds, water-mediated contacts, or non-canonical geometry.

---

### `long_hbonds`
**Description**: H-bonds detected but distances exceed 4.0Å (extreme stretching).

**Distance ranges**:
- Ideal: 2.7-3.2Å
- Acceptable: 2.0-3.5Å
- Stretched: 3.5-4.0Å (not flagged)
- **Extreme (flagged): > 4.0Å**

**Interpretation**: The H-bonds are extremely stretched - beyond what would typically be considered a real hydrogen bond. This indicates severe structural distortion.

---

### `short_hbonds`
**Description**: H-bonds detected but distances are less than 2.0Å.

**Interpretation**: Indicates a steric clash - atoms are too close together. Usually indicates poor structure refinement or incorrect atom positions.

---

### `overloaded_acceptor`
**Description**: An acceptor atom is receiving more than 2 H-bonds.

**Interpretation**: Chemically, sp2 oxygen/nitrogen can accept at most 2 H-bonds (one per lone pair). More than 2 suggests detection errors or unusual geometry.

---

## Geometry Related Reasons

### `poor_planarity`
**Description**: Bases are not coplanar.

**Threshold**: Interbase angle ≥ 20°

**Interpretation**: The two bases are significantly tilted relative to each other. Watson-Crick pairs should have nearly parallel bases (typically < 15°). Large angles indicate propeller twist, buckle deformation, or non-planar pair type.

---

### `geometric_outlier`
**Description**: Geometric parameters are outside normal ranges.

**Criteria** (from DSSR data):
- Interbase angle > 15°
- N1-N9 distance outside expected range for the sequence

**Interpretation**: The overall geometry deviates from ideal cWW. Bases may be tilted, shifted, or have unusual glycosidic bond positions.

---

### `rmsd_prefers_other`
**Description**: A non-Watson-Crick template fits better than cWW.

**Threshold**: RMSD gap > 0.5Å AND best template uses different edges (not WW)

**Interpretation**: The geometry matches a different Leontis-Westhof class better. For example, if cWS has RMSD 0.3Å and cWW has RMSD 0.9Å, the pair is likely misclassified.

**Note**: cWW and tWW both use Watson edges and are not distinguished by this check.

---

## Classification Related Reasons

### `non_canonical`
**Description**: DSSR reports Saenger class as "--" (non-canonical).

**Interpretation**: DSSR itself indicates this is not a standard Watson-Crick pair. It may be classified as cWW based on edge geometry but lacks the canonical H-bonding pattern (Saenger XIX for G-C, XX for A-U).

---

### `dssr_questionable`
**Description**: DSSR classification is likely incorrect.

**Criteria**:
- RMSD to cWW template > 1.0Å (very poor geometry)
- No H-bonds detected

**Interpretation**: This pair has neither good geometry nor H-bonds to support cWW classification. DSSR may have made an error. These cases should NOT count against our classification accuracy.

---

## BP Score

The BP (Base Pair) Score determines whether to flag a pair:

**Formula**: `Score = 0.30×RMSD + 0.40×Coverage + 0.30×Quality`

| Component | Description |
|-----------|-------------|
| RMSD Score | 1.0 if ≤0.3Å, 0.0 if ≥1.0Å, linear between |
| Coverage | (found H-bonds) / (expected H-bonds) |
| Quality | Average H-bond quality (distance + alignment) |

### Geometry-Adjusted H-Bond Leniency

If geometry is good (low RMSD), long H-bonds are penalized less. This handles "stretched but real" cWW pairs - good geometric fit but elongated H-bonds due to structural strain.

| RMSD | Geometry Leniency | Acceptable H-Bond Distance |
|------|------------------|---------------------------|
| ≤ 0.5Å | 100% (full) | Up to 4.2Å without penalty |
| 0.5-0.8Å | Partial (linear) | 3.2-4.2Å sliding scale |
| ≥ 0.8Å | 0% (none) | Penalized starting at 3.2Å |

**Rationale**: If a base pair fits the cWW template geometrically, it's likely a real Watson-Crick pair even if H-bonds are stretched. Poor geometry + long H-bonds = likely not a cWW pair.

### Extended H-Bond Search

For pairs with:
- Good geometry (RMSD < 1.0Å)
- Good planarity (interbase angle < 30°)
- Few/no H-bonds at normal threshold

The system recalculates H-bonds with:
- Extended distance threshold: 5.0Å (vs normal 3.5Å)
- Alignment-weighted scoring: as distance increases, alignment quality matters more

**Quality Scoring for Extended H-Bonds**:
| Distance | Distance Weight | Alignment Weight |
|----------|-----------------|------------------|
| < 3.2Å | 70% | 30% |
| 3.2-4.0Å | 50% | 50% |
| > 4.0Å | 30% | 70% |

**Rationale**: Poor refinement can stretch H-bonds beyond normal thresholds. If the geometry is correct, we look for H-bonds at longer distances and rely on alignment (angle between donor-H and acceptor lone pair) to validate them.

**Grades**:
| Grade | Score | Interpretation |
|-------|-------|----------------|
| A | ≥ 0.90 | Excellent cWW |
| B | ≥ 0.80 | Good cWW |
| C | ≥ 0.70 | Acceptable (not flagged) |
| D | ≥ 0.60 | Marginal |
| F | < 0.60 | Poor/questionable |

---

## Usage

```bash
# List all outliers for a PDB
python visualize_outlier.py --pdb 1FKA

# Filter by specific reason
python visualize_outlier.py --pdb 1FKA --reason no_hbonds

# Render a specific outlier
python visualize_outlier.py --pdb 1FKA --render 0
```
