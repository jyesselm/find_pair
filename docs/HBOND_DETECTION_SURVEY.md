# H-Bond Detection Survey: Barnaba, FR3D, HBPLUS, ClaRNA

A comprehensive survey of hydrogen bond detection implementations across four established RNA/protein structure analysis tools. This document summarizes key findings and identifies features worth adopting or avoiding.

## Executive Summary

| Tool | Language | Distance Cutoff | Angle Cutoff | Key Innovation |
|------|----------|-----------------|--------------|----------------|
| **Barnaba** | Python | 3.3 Å | 60° (plane) | Ellipsoidal distance metric, LCS geometry |
| **FR3D** | MATLAB | 4.0-4.5 Å | 110° D-H-A | Isostericity-based atom assignment, percentile scoring |
| **HBPLUS** | C | 3.9 Å D-A, 2.5 Å H-A | 90° DHA/DAAA/HAAA | Aromatic ring geometry, multi-angle validation |
| **ClaRNA** | Python | 4.0 Å (C/O), 3.5 Å (N/H) | 130° | ML-based classification, graph consistency |

---

## 1. Distance Thresholds Comparison

| Tool | Donor-Acceptor | H-Acceptor | Notes |
|------|----------------|------------|-------|
| Barnaba | 3.3 Å | N/A (heavy atoms only) | Squared distance check (0.1089 nm²) |
| FR3D | 4.5 Å (heavy), 4.0 Å (with H) | Implicit in angle | Different thresholds based on H availability |
| HBPLUS | 3.9 Å | 2.5 Å | Most conservative; requires both checks |
| ClaRNA | 4.0 Å (C/O), 3.5 Å (N/H) | N/A | Atom-type dependent thresholds |
| **Our current** | 3.5-4.0 Å | N/A | Context-dependent (dssr_like preset) |

### Recommendations
- **ADOPT**: FR3D's dual-threshold approach (looser for heavy-atom only, stricter with H)
- **ADOPT**: ClaRNA's atom-type dependent cutoffs (N vs O donors)
- **REJECT**: Single fixed cutoffs (too rigid for diverse geometries)

---

## 2. Angle Criteria Comparison

| Tool | Donor Angle (D-H-A) | Acceptor Angle | Special Handling |
|------|---------------------|----------------|------------------|
| Barnaba | None | None | Only checks plane angles (60°) |
| FR3D | ≥ 110° | Via coplanarity | Requires hydrogen to be present |
| HBPLUS | ≥ 90° | ≥ 90° (DAAA, HAAA) | Aromatic: ≤ 20° to ring axis |
| ClaRNA | ≥ 130° | N/A | Stricter than others |
| **Our current** | ≥ 90° (optional) | ≥ 70° (optional) | Disabled by default |

### Recommendations
- **ADOPT**: HBPLUS's multi-angle validation (DHA + DAAA + HAAA)
- **ADOPT**: HBPLUS's aromatic ring special handling (DAAX/HAAX ≤ 20°)
- **CONSIDER**: ClaRNA's stricter 130° threshold for high-confidence detection
- **REJECT**: Barnaba's no-angle approach (misses geometric validation entirely)

---

## 3. Donor/Acceptor Classification

### Approaches Compared

| Tool | Method | Flexibility | Modified Base Support |
|------|--------|-------------|----------------------|
| Barnaba | Static per-residue dict | Low | None (hardcoded lists) |
| FR3D | Isostericity table lookup | Medium | Parent base fallback |
| HBPLUS | Residue template + runtime | High | Runtime addition possible |
| ClaRNA | Per-base atom pairs | Low | None |
| **Our current** | Element-based fallback | Medium | Limited |

### Key Donor/Acceptor Atoms by Tool

```
              Barnaba          FR3D            HBPLUS          ClaRNA
Adenine
  Donors:     N6,C2,C8,O2'    (from table)    N6,C2,C8        C2,C8,N6
  Acceptors:  N1,N3,N7,O2'    (from table)    N1,N3,N7        (implicit)

Guanine
  Donors:     N1,N2,C8,O2'    (from table)    N1,N2,C8        N1,C8,N2
  Acceptors:  O6,N3,N7,O2'    (from table)    O6,N3,N7        (implicit)
```

### Recommendations
- **ADOPT**: FR3D's isostericity-based assignment (grounded in structural biology)
- **ADOPT**: HBPLUS's runtime residue addition for unknown modifications
- **REJECT**: Static hardcoded lists (cannot handle modified bases gracefully)

---

## 4. Quality Scoring Approaches

| Tool | Scoring Method | Output | Continuous? |
|------|---------------|--------|-------------|
| Barnaba | H-bond count only | Binary (2+ for AU, 3+ for CG) | No |
| FR3D | Percentile-based coplanarity | 0-1 score | Yes |
| HBPLUS | Pass/fail filters | Binary | No |
| ClaRNA | ML classifier confidence | Fuzzy (?) marking | Semi |
| **Our current** | Gaussian distance + linear angle | 0-100 score | Yes |

### FR3D Coplanarity Score (Worth Adopting)
```
Coplanar = min(Gap1Val, Gap2Val, dot1Val, dot2Val, dot3Val, MinDistVal)
- Each component: percentile from known RNA structures (70th/90th/97th)
- Final score: 0-1 (minimum of all components)
```

### Recommendations
- **ADOPT**: FR3D's percentile-based normalization (grounded in real structure statistics)
- **ADOPT**: Our Gaussian + linear scoring (modern, continuous)
- **CONSIDER**: ClaRNA's fuzzy classification concept for borderline cases
- **REJECT**: Binary-only classification (loses nuance)

---

## 5. Geometric Innovations Worth Adopting

### 5.1 Ellipsoidal Distance (Barnaba)
```python
# Scale factors for flattened base pair geometry
f_factors = [0.5, 0.5, 0.3]  # x, y, z weights
scaled_dist = sqrt((dx/0.5)² + (dy/0.5)² + (dz/0.3)²)
```
**Verdict**: CONSIDER - accounts for anisotropic base geometry, but adds complexity.

### 5.2 Aromatic Ring Geometry (HBPLUS)
```c
// For aromatic acceptors, check angle to ring plane axis
if (daax_ang < cosMAX_DAAX)  // 20° threshold
    reject();
```
**Verdict**: ADOPT - critical for proper pi-stacking and aromatic acceptor handling.

### 5.3 LCS-Based Edge Classification (Barnaba)
```python
# Local coordinate system defines Watson/Hoogsteen/Sugar edges
edge_angles = arctan2(vectors[:,1], vectors[:,0]) - theta1
bins = [0, 1.84, 3.84, 2*pi]  # Edge boundaries in radians
edge = digitize(edge_angles, bins) - 1
```
**Verdict**: ADOPT - enables Leontis-Westhof edge classification.

### 5.4 SVD Superimposition (ClaRNA)
```python
# Fit actual base to reference structure
R, t = svd_superimpose(actual_coords, reference_coords)
rms = compute_rms(actual_coords @ R + t, reference_coords)
```
**Verdict**: CONSIDER - useful for base normalization but we already have frame fitting.

---

## 6. H-Bond Requirements by Pair Type

### Barnaba (Simple Count)
| Pair Type | Minimum H-bonds |
|-----------|-----------------|
| AU | 2 |
| CG | 3 |
| GU | 2 |

### FR3D (Flexible)
| Expected H-bonds | Required Good | Notes |
|------------------|---------------|-------|
| 4 | 3 | Allow 1 missing |
| 3 | 2 | Allow 1 missing |
| 1-2 | All | Strict for small counts |

**Verdict**: ADOPT FR3D's flexible approach - better handles real structure variation.

---

## 7. Features to ADOPT

### High Priority (Clear Benefit)

1. **Multi-angle validation (HBPLUS)**
   - Add DAAA angle check (donor-acceptor-acceptor-antecedent)
   - Add HAAA angle check (H-acceptor-acceptor-antecedent)
   - Validates acceptor orbital geometry

2. **Aromatic ring handling (HBPLUS)**
   - Detect aromatic acceptors (Phe, Tyr, Trp, His)
   - Calculate perpendicular ring axis
   - Apply 20° maximum angle to axis

3. **Flexible H-bond requirements (FR3D)**
   - For 3+ expected H-bonds, allow 1 missing
   - For 1-2 expected, require all
   - Handles crystal structure imperfections

4. **Leontis-Westhof edge classification (Barnaba/FR3D)**
   - Classify edges: Watson, Hoogsteen, Sugar
   - Enable proper non-canonical pair annotation

### Medium Priority (Good to Have)

5. **Percentile-based scoring (FR3D)**
   - Normalize scores against known structure statistics
   - Report quality as percentile (e.g., "top 10%")

6. **Bidirectional H-bond checking (Barnaba)**
   - Check D1→A2 and D2→A1
   - Ensures symmetric detection

7. **Atom-type dependent thresholds (ClaRNA)**
   - Different cutoffs for N vs O donors
   - More chemically accurate

### Low Priority (Consider Later)

8. **Ellipsoidal distance metric (Barnaba)**
   - Account for flattened base geometry
   - May improve edge cases

9. **Graph consistency (ClaRNA)**
   - Post-process to enforce bidirectional agreement
   - Useful for complex structures

---

## 8. Features to REJECT

1. **No angle validation (Barnaba)**
   - Only checks plane angles, not H-bond angles
   - Misses many invalid H-bonds

2. **Fixed 90° threshold for all angles (HBPLUS)**
   - Too strict for modern analysis
   - Reject many valid H-bonds

3. **Black-box ML classification (ClaRNA)**
   - Not interpretable
   - Hard to debug/tune

4. **Binary-only classification (HBPLUS)**
   - Loses quality information
   - Can't distinguish strong vs weak bonds

5. **Python 2 code patterns (ClaRNA)**
   - `.has_key()`, `xrange`
   - Not maintainable

6. **MATLAB-only implementation (FR3D)**
   - Not portable
   - Expensive licensing

---

## 9. Implementation Priorities

### Phase 1: Core Improvements
1. Add DAAA angle calculation and validation
2. Add aromatic ring special handling
3. Implement flexible H-bond count requirements

### Phase 2: Edge Classification
4. Add LCS-based edge detection (W/H/S)
5. Enable Leontis-Westhof annotation output

### Phase 3: Scoring Refinements
6. Add percentile normalization to quality scores
7. Add bidirectional consistency check

### Phase 4: Advanced Features
8. Consider ellipsoidal distance option
9. Add graph consistency post-processing

---

## References

- **Barnaba**: Bottaro et al. (2019) "Barnaba: software for analysis of nucleic acid structures and trajectories"
- **FR3D**: Sarver et al. (2008) "FR3D: Finding local and composite recurrent structural motifs in RNA 3D structures"
- **HBPLUS**: McDonald & Thornton (1994) "Satisfying Hydrogen Bonding Potential in Proteins"
- **ClaRNA**: Waleń et al. (2014) "ClaRNA: a classifier of contacts in RNA 3D structures"
- **Leontis-Westhof**: Leontis & Westhof (2001) "Geometric nomenclature and classification of RNA base pairs"
