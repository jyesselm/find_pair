# H-Bond Advanced Features Implementation Plan

## Overview

Implement advanced H-bond detection features from industry survey (Barnaba, FR3D, HBPLUS, ClaRNA) into the existing H-bond detection system while maintaining legacy compatibility.

**Reference**: See `docs/HBOND_DETECTION_SURVEY.md` for full analysis of external tools.

---

## Phase 1: Core Improvements (High Priority)

### 1.1 DAAA Angle Support

Add Donor-Acceptor-Acceptor-Antecedent angle calculation (from HBPLUS).

**Rationale**: Validates acceptor orbital geometry. The antecedent is the atom bonded to the acceptor, and the angle D...A-AA should be reasonable for proper orbital alignment.

**Files to modify:**
- `include/x3dna/core/hbond.hpp` - Add `daaa_angle`, `antecedent_atom` fields
- `include/x3dna/core/hbond_quality.hpp` - Add `daaa_angle_score` to score struct
- `include/x3dna/algorithms/hydrogen_bond/geometry.hpp` - Add antecedent lookup
- `src/x3dna/algorithms/hydrogen_bond/geometry.cpp` - Implement antecedent lookup
- `src/x3dna/algorithms/hydrogen_bond/quality_scorer.cpp` - Add DAAA scoring

**Key data structure:**
```cpp
// Antecedent lookup: acceptor atom -> bonded neighbor for DAAA angle
const std::unordered_map<std::string, std::string> ANTECEDENT_LOOKUP = {
    // Adenine acceptors
    {"N1", "C6"}, {"N3", "C4"}, {"N7", "C8"},
    // Guanine acceptors
    {"O6", "C5"}, {"N3", "C4"}, {"N7", "C8"},
    // Cytosine acceptors
    {"O2", "N3"}, {"N3", "C4"},
    // Uracil/Thymine acceptors
    {"O2", "N3"}, {"O4", "C5"},
};
```

**Updated scoring weights:**
```cpp
struct HBondScoringParams {
    double weight_distance = 0.40;        // was 0.45
    double weight_donor_angle = 0.25;     // was 0.30
    double weight_acceptor_angle = 0.20;  // was 0.25
    double weight_daaa_angle = 0.15;      // NEW
};
```

---

### 1.2 Aromatic Ring Special Handling

Add special treatment for aromatic acceptors (from HBPLUS).

**Rationale**: Aromatic rings have different H-bonding geometry. The H-bond should be roughly perpendicular to the ring plane (within 20°).

**New files:**
- `include/x3dna/algorithms/hydrogen_bond/aromatic_classifier.hpp`
- `src/x3dna/algorithms/hydrogen_bond/aromatic_classifier.cpp`

**Interface:**
```cpp
class AromaticAcceptorClassifier {
public:
    // Returns true if atom is in/adjacent to aromatic ring
    [[nodiscard]] static bool is_aromatic_acceptor(
        const std::string& atom_name,
        const std::string& residue_name);

    // Get atoms defining the aromatic ring
    [[nodiscard]] static std::vector<std::string> get_ring_atoms(
        const std::string& residue_name);

    // Calculate perpendicular axis to ring plane
    [[nodiscard]] static std::optional<geometry::Vector3D> calculate_ring_perpendicular(
        const core::Residue& residue);

    // Calculate angle between D-A vector and ring perpendicular (DAAX)
    [[nodiscard]] static double calculate_daax_angle(
        const geometry::Vector3D& donor_pos,
        const geometry::Vector3D& acceptor_pos,
        const geometry::Vector3D& ring_perpendicular);
};
```

**Aromatic atoms:**
```cpp
// Nucleobase ring atoms (all aromatic)
static const std::unordered_set<std::string> NA_AROMATIC = {
    "N1", "N3", "N7", "O6", "O2", "O4", "N6", "N2", "N4"
};

// Protein aromatic sidechains
static const std::unordered_map<std::string, std::unordered_set<std::string>> PROTEIN_AROMATIC = {
    {"HIS", {"ND1", "NE2"}},
    {"TYR", {"OH"}},
    {"TRP", {"NE1"}},
};
```

**Scoring modification:**
- For aromatic acceptors, apply 20° maximum DAAX angle threshold
- Penalize (or reject) H-bonds with DAAX angle > 20°

---

### 1.3 Flexible H-Bond Count Requirements

Allow partial H-bond satisfaction for pair validation (from FR3D).

**Rationale**: Real crystal structures often have 1 missing H-bond due to disorder, resolution limits, or genuine structural variation. FR3D allows this flexibility.

**New struct in detection_params.hpp:**
```cpp
struct HBondCountRequirements {
    int min_for_expected_4 = 3;  // 4 expected, require 3
    int min_for_expected_3 = 2;  // 3 expected, require 2
    int min_for_expected_2 = 2;  // 1-2 expected, require all
    int min_for_expected_1 = 1;

    [[nodiscard]] int min_required(int expected) const;
};
```

**Expected H-bonds by pair type:**
```cpp
// Watson-Crick pairs
A-U: 2 expected
G-C: 3 expected
G-U: 2 expected (wobble)

// Non-canonical
A-A: 2 expected
G-G: 2-3 expected
// etc.
```

**Location**: `src/x3dna/algorithms/pair_identification/base_pair_validator.cpp`

---

## Phase 2: Edge Classification (Medium Priority)

### 2.1 Leontis-Westhof Edge Classification

Add Watson/Hoogsteen/Sugar edge classification (from Barnaba/FR3D).

**Rationale**: Standard nomenclature for RNA base pair annotation. Enables proper non-canonical pair classification.

**New enum in hbond_types.hpp:**
```cpp
enum class BaseEdge {
    WATSON,     // W edge (N1/N6 side of purine, N3/N4 side of pyrimidine)
    HOOGSTEEN,  // H edge (N7/N6 side of purine)
    SUGAR,      // S edge (N3/O2'/C2' side)
    UNKNOWN
};
```

**New files:**
- `include/x3dna/algorithms/hydrogen_bond/edge_classifier.hpp`
- `src/x3dna/algorithms/hydrogen_bond/edge_classifier.cpp`

**Interface:**
```cpp
class EdgeClassifier {
public:
    // Classify which edge of the base the atom is on
    [[nodiscard]] static core::BaseEdge classify_edge(
        const core::Residue& residue,
        const std::string& atom_name);

    // Get angle from base center to atom (for debugging)
    [[nodiscard]] static double get_edge_angle(
        const core::Residue& residue,
        const std::string& atom_name);

private:
    // Angle ranges for each edge (degrees from reference axis)
    static constexpr double WATSON_MIN = -30.0;
    static constexpr double WATSON_MAX = 60.0;
    static constexpr double HOOGSTEEN_MIN = 60.0;
    static constexpr double HOOGSTEEN_MAX = 150.0;
    // SUGAR is the remainder
};
```

**Add to HBond struct:**
```cpp
BaseEdge donor_edge = BaseEdge::UNKNOWN;
BaseEdge acceptor_edge = BaseEdge::UNKNOWN;
```

---

## Phase 3: Scoring Refinements (Medium Priority)

### 3.1 Atom-Type Dependent Thresholds

Use different distance cutoffs for N vs O donors (from ClaRNA).

**Rationale**: Nitrogen donors form tighter H-bonds than oxygen donors.

**Add to HBondDistanceThresholds:**
```cpp
double n_donor_max = 3.5;  // Nitrogen donor max distance (tighter)
double o_donor_max = 4.0;  // Oxygen donor max distance (standard)

[[nodiscard]] double max_for_donor_element(const std::string& element) const;
```

**Modify find_candidate_bonds()** to use element-specific thresholds when `enable_element_thresholds` is true.

---

## File Summary

| File | Action | What |
|------|--------|------|
| `include/x3dna/core/hbond.hpp` | modify | Add `daaa_angle`, `antecedent_atom`, `is_aromatic_acceptor`, `donor_edge`, `acceptor_edge` |
| `include/x3dna/core/hbond_types.hpp` | modify | Add `BaseEdge` enum |
| `include/x3dna/core/hbond_quality.hpp` | modify | Add `daaa_angle_score` to `HBondQualityScore` |
| `include/x3dna/algorithms/hydrogen_bond/geometry.hpp` | modify | Add antecedent lookup |
| `src/x3dna/algorithms/hydrogen_bond/geometry.cpp` | modify | Implement antecedent lookup and DAAA angle |
| `include/x3dna/algorithms/hydrogen_bond/aromatic_classifier.hpp` | **create** | AromaticAcceptorClassifier |
| `src/x3dna/algorithms/hydrogen_bond/aromatic_classifier.cpp` | **create** | Aromatic detection implementation |
| `include/x3dna/algorithms/hydrogen_bond/edge_classifier.hpp` | **create** | EdgeClassifier |
| `src/x3dna/algorithms/hydrogen_bond/edge_classifier.cpp` | **create** | Edge classification implementation |
| `include/x3dna/algorithms/hydrogen_bond/detection_params.hpp` | modify | Add `HBondCountRequirements`, element thresholds |
| `src/x3dna/algorithms/hydrogen_bond/detection_params.cpp` | modify | Implement new presets |
| `include/x3dna/algorithms/hydrogen_bond/quality_scorer.hpp` | modify | Add DAAA scoring, aromatic handling |
| `src/x3dna/algorithms/hydrogen_bond/quality_scorer.cpp` | modify | Implement DAAA scoring |
| `src/x3dna/algorithms/hydrogen_bond/detector.cpp` | modify | Call classifiers |
| `src/x3dna/algorithms/pair_identification/base_pair_validator.cpp` | modify | Flexible H-bond count |

---

## Implementation Order

### Step 1: DAAA Angle Infrastructure
1. Add fields to `HBond` struct
2. Add `daaa_angle_score` to quality score
3. Add `ANTECEDENT_LOOKUP` map to geometry.cpp
4. Add `get_antecedent_atom_name()` and `find_antecedent_position()` functions
5. Integrate into `calculate_angles()` in detector.cpp
6. Add `score_daaa_angle()` to quality scorer
7. Update scoring weights

### Step 2: Aromatic Acceptor Classification
1. Create `aromatic_classifier.hpp/cpp`
2. Implement aromatic atom identification
3. Implement ring perpendicular calculation
4. Implement DAAX angle calculation
5. Add `is_aromatic_acceptor` to HBond
6. Modify scoring to use 20° threshold for aromatic

### Step 3: Flexible H-Bond Count
1. Add `HBondCountRequirements` struct
2. Add expected H-bond count lookup by pair type
3. Modify `base_pair_validator.cpp` to use flexible requirements

### Step 4: Edge Classification
1. Add `BaseEdge` enum
2. Create `edge_classifier.hpp/cpp`
3. Implement angle-based edge classification
4. Add `donor_edge`, `acceptor_edge` to HBond
5. Integrate into detection pipeline

### Step 5: Element-Specific Thresholds
1. Add `n_donor_max`, `o_donor_max` to distance thresholds
2. Add `max_for_donor_element()` method
3. Add `enable_element_thresholds` flag
4. Modify `find_candidate_bonds()` to use element thresholds

### Step 6: New Presets
1. Create `HBondDetectionParams::enhanced()` preset with all new features
2. Keep `legacy_compatible()` unchanged
3. Update `modern()` to optionally enable features

---

## Testing Strategy

### Unit Tests
- `test_daaa_angle.cpp` - Verify DAAA calculation on known structures
- `test_aromatic_classifier.cpp` - Verify aromatic atom identification
- `test_edge_classifier.cpp` - Verify edge angle calculation and classification

### Integration Tests
- Run `fp2-validate validate core --test-set fast` after each step
- Verify `legacy_compatible()` preset maintains 99%+ pass rate
- Test `enhanced()` preset on DSSR comparison set

### Validation Criteria
- Legacy compatibility: No regression in existing tests
- DSSR match rate: Should improve or maintain 98%+ with enhanced preset
- Performance: New features add <5% overhead (they're optional)

---

## Risks and Mitigations

1. **Performance impact**: DAAA and aromatic classification add atom lookups
   - *Mitigation*: Disabled by default in legacy preset

2. **Legacy compatibility**: New HBond fields affect serialization
   - *Mitigation*: Default values, only populate when features enabled

3. **Edge angle boundaries**: LW classification varies in literature
   - *Mitigation*: Use published Leontis-Westhof 2001 ranges, make configurable

4. **Aromatic detection failures**: Modified bases may lack aromaticity
   - *Mitigation*: Limit to standard residues, fallback to non-aromatic

5. **Flexible count ambiguity**: Different pair types have different expectations
   - *Mitigation*: Build lookup from standard pair definitions

---

## References

- Leontis & Westhof (2001) "Geometric nomenclature and classification of RNA base pairs"
- McDonald & Thornton (1994) "Satisfying Hydrogen Bonding Potential in Proteins" (HBPLUS)
- Sarver et al. (2008) "FR3D: Finding local and composite recurrent structural motifs"
- Bottaro et al. (2019) "Barnaba: software for analysis of nucleic acid structures"
