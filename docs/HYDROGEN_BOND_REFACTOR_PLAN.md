# Hydrogen Bond Refactoring Plan

## Goal

Extend the hydrogen bond detection system to capture richer information about hydrogen bonds including molecular context (RNA-RNA, RNA-protein, etc.), structural context (base-base, base-sugar, etc.), and geometric angles (X-D-A, D-A-X, dihedral) while maintaining full backward compatibility with legacy JSON output.

---

## Current State Analysis

### Existing Data Structures

#### 1. `HydrogenBond` struct (`include/x3dna/core/hydrogen_bond.hpp`)
```cpp
struct HydrogenBond {
    std::string donor_atom;          // Donor atom name (e.g., "N6")
    std::string acceptor_atom;       // Acceptor atom name (e.g., "O4")
    double distance = 0.0;           // Bond distance in Angstroms
    char type = ' ';                 // '-' for standard, ' ' for non-standard
    std::optional<size_t> hbond_idx; // Optional index for tracking
};
```

#### 2. `HydrogenBondResult` struct (`include/x3dna/algorithms/hydrogen_bond_finder.hpp`)
```cpp
struct HydrogenBondResult {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;        // '-' for standard, '*' for non-standard, ' ' for invalid
    int linkage_type; // Conflict resolution marker (0, 1, 2, or 18)
};
```

#### 3. `DetailedHBondResult` struct (`include/x3dna/algorithms/hydrogen_bond_finder.hpp`)
```cpp
struct DetailedHBondResult {
    std::vector<HydrogenBondResult> initial_hbonds;
    std::vector<HydrogenBondResult> after_conflict_resolution;
    std::vector<HydrogenBondResult> after_validation;
    std::vector<HydrogenBondResult> final_hbonds;
    int num_good_hb;
};
```

### Current Detection Flow

The `HydrogenBondFinder::find_hydrogen_bonds_detailed()` function implements:

1. **Initial Detection**: Loop through all atom pairs, check `good_hb_atoms()` and distance within `[hb_lower, hb_dist1]`
2. **Conflict Resolution**: `resolve_conflicts()` handles atoms with multiple potential H-bonds
3. **Validation**: `validate_hbonds()` calls `donor_acceptor()` to determine if standard H-bond
4. **Filtering**: Remove invalid bonds and count good H-bonds

### Available Residue Type Information

From `typing::ResidueClassification`:
- `molecule_type`: NucleicAcid, Protein, Water, Ion, Ligand, Unknown
- `nucleic_acid_type`: RNA, DNA, Unknown
- `base_type`: Adenine, Guanine, Cytosine, Thymine, Uracil, etc.
- `base_category`: Purine, Pyrimidine, Unknown
- `is_modified_nucleotide`, `is_modified_amino_acid`

From `INucleotide` interface:
- `is_rna()`, `is_dna()`, `is_purine()`, `is_pyrimidine()`
- `one_letter_code()`, `base_type()`

### Legacy JSON Output Format

From `data/json_legacy/hbond_list/100D.json`:
```json
{
  "type": "hbond_list",
  "base_i": 1,
  "base_j": 20,
  "num_hbonds": 7,
  "hb_info_string": "[3]  O2 - N2  2.77  N3 - N1  2.93  N4 - O6  2.93",
  "hbonds": [
    {
      "hbond_idx": 1,
      "donor_atom": " O2 ",
      "acceptor_atom": " N1 ",
      "distance": 3.716618,
      "type": " "
    },
    ...
  ]
}
```

---

## Proposed Data Structures

### Phase 1: New Enums for H-Bond Classification

**File**: `include/x3dna/core/hydrogen_bond_types.hpp` (new)

```cpp
namespace x3dna {
namespace core {

/**
 * @enum HBondMoleculeType
 * @brief Classification of H-bond by molecule types involved
 */
enum class HBondMoleculeType {
    RNA_RNA,        // Both residues are RNA
    DNA_DNA,        // Both residues are DNA
    RNA_DNA,        // One RNA, one DNA (hybrid)
    RNA_PROTEIN,    // RNA to protein
    DNA_PROTEIN,    // DNA to protein
    PROTEIN_PROTEIN,// Both residues are protein
    NUCLEIC_WATER,  // Nucleic acid to water
    PROTEIN_WATER,  // Protein to water
    NUCLEIC_ION,    // Nucleic acid to ion
    NUCLEIC_LIGAND, // Nucleic acid to ligand
    OTHER,          // Other combinations
    UNKNOWN         // Unable to classify
};

/**
 * @enum HBondAtomContext
 * @brief Classification of atom location within residue
 */
enum class HBondAtomContext {
    BASE,           // Base atom (N1, N2, N3, N4, N6, N7, O2, O4, O6, etc.)
    SUGAR,          // Sugar atom (O2', O3', O4', O5', etc.)
    BACKBONE,       // Backbone phosphate (O1P, O2P, O3', O5')
    SIDECHAIN,      // Amino acid sidechain
    MAINCHAIN,      // Amino acid mainchain (N, O, etc.)
    WATER,          // Water oxygen
    UNKNOWN         // Unable to classify
};

/**
 * @enum HBondStructuralType
 * @brief Combined structural classification (for nucleic acids)
 */
enum class HBondStructuralType {
    BASE_BASE,       // Both atoms in base region
    BASE_SUGAR,      // Base to sugar (e.g., N3 to O2')
    BASE_BACKBONE,   // Base to phosphate backbone
    SUGAR_SUGAR,     // Sugar to sugar (e.g., O2' to O2')
    SUGAR_BACKBONE,  // Sugar to backbone
    BACKBONE_BACKBONE, // Backbone to backbone
    NUCLEIC_PROTEIN, // Any nucleic acid to protein
    OTHER,           // Other combinations
    UNKNOWN          // Unable to classify
};

/**
 * @struct HBondAngles
 * @brief Geometric angles characterizing the H-bond
 */
struct HBondAngles {
    double xda_angle = 0.0;    // X-D-A angle (degrees), where X is neighbor of D
    double dax_angle = 0.0;    // D-A-X angle (degrees), where X is neighbor of A
    double xdax_dihedral = 0.0; // X-D-A-X dihedral angle (degrees)
    bool has_donor_neighbor = false;    // Whether donor neighbor was found
    bool has_acceptor_neighbor = false; // Whether acceptor neighbor was found
};

} // namespace core
} // namespace x3dna
```

### Phase 1-4: Extended HydrogenBond Struct

**File**: `include/x3dna/core/hydrogen_bond.hpp` (modified)

```cpp
struct HydrogenBond {
    // === Existing fields (unchanged) ===
    std::string donor_atom;
    std::string acceptor_atom;
    double distance = 0.0;
    char type = ' ';
    std::optional<size_t> hbond_idx;
    
    // === Phase 1: Classification fields (added) ===
    HBondMoleculeType molecule_type = HBondMoleculeType::UNKNOWN;
    HBondAtomContext donor_context = HBondAtomContext::UNKNOWN;
    HBondAtomContext acceptor_context = HBondAtomContext::UNKNOWN;
    HBondStructuralType structural_type = HBondStructuralType::UNKNOWN;
    
    // === Phase 2-3: Geometric angle fields (added) ===
    std::optional<HBondAngles> angles;
    
    // === Existing method (unchanged) ===
    [[nodiscard]] bool is_standard() const { return type == '-'; }
    
    // === Phase 1: New classification methods ===
    [[nodiscard]] bool is_base_base() const {
        return structural_type == HBondStructuralType::BASE_BASE;
    }
    
    [[nodiscard]] bool is_rna_rna() const {
        return molecule_type == HBondMoleculeType::RNA_RNA;
    }
    
    [[nodiscard]] bool involves_sugar() const {
        return donor_context == HBondAtomContext::SUGAR ||
               acceptor_context == HBondAtomContext::SUGAR;
    }
    
    [[nodiscard]] bool involves_protein() const {
        return molecule_type == HBondMoleculeType::RNA_PROTEIN ||
               molecule_type == HBondMoleculeType::DNA_PROTEIN ||
               molecule_type == HBondMoleculeType::PROTEIN_PROTEIN;
    }
};
```

---

## Implementation Plan

### Phase 1: Add New Fields Without Changing Detection Logic

**Goal**: Add classification fields to `HydrogenBond` and populate them after detection.

#### 1.1 Create New Header File

**File**: `include/x3dna/core/hydrogen_bond_types.hpp`
- Define `HBondMoleculeType`, `HBondAtomContext`, `HBondStructuralType`, `HBondAngles`
- Include from `hydrogen_bond.hpp`

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/core/hydrogen_bond_types.hpp` | create | New enums and HBondAngles struct |
| `include/x3dna/core/hydrogen_bond.hpp` | modify | Add new optional fields, include types header |

#### 1.2 Create Atom Context Classifier

**File**: `include/x3dna/algorithms/hydrogen_bond/atom_context_classifier.hpp`
**File**: `src/x3dna/algorithms/hydrogen_bond/atom_context_classifier.cpp`

```cpp
class AtomContextClassifier {
public:
    /**
     * @brief Classify atom by its position in residue structure
     * @param atom_name Atom name (e.g., "N6", "O2'", "O1P")
     * @param residue Optional residue for additional context
     * @return HBondAtomContext classification
     */
    [[nodiscard]] static HBondAtomContext classify(
        const std::string& atom_name,
        const core::poly::IResidue* residue = nullptr);
    
    /**
     * @brief Check if atom is a base atom (for nucleotides)
     */
    [[nodiscard]] static bool is_base_atom(const std::string& atom_name);
    
    /**
     * @brief Check if atom is a sugar atom
     */
    [[nodiscard]] static bool is_sugar_atom(const std::string& atom_name);
    
    /**
     * @brief Check if atom is a backbone atom
     */
    [[nodiscard]] static bool is_backbone_atom(const std::string& atom_name);
};
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/atom_context_classifier.hpp` | create | Atom context classification interface |
| `src/x3dna/algorithms/hydrogen_bond/atom_context_classifier.cpp` | create | Implementation with atom name patterns |

#### 1.3 Create H-Bond Type Classifier

**File**: `include/x3dna/algorithms/hydrogen_bond/hbond_classifier.hpp`
**File**: `src/x3dna/algorithms/hydrogen_bond/hbond_classifier.cpp`

```cpp
class HBondClassifier {
public:
    /**
     * @brief Classify H-bond by molecule types of both residues
     * @param res1 First residue (donor side)
     * @param res2 Second residue (acceptor side)
     * @return HBondMoleculeType classification
     */
    [[nodiscard]] static HBondMoleculeType classify_molecule_type(
        const core::poly::IResidue& res1,
        const core::poly::IResidue& res2);
    
    /**
     * @brief Classify H-bond structural type from atom contexts
     * @param donor_context Donor atom context
     * @param acceptor_context Acceptor atom context
     * @param mol_type Molecule type for context
     * @return HBondStructuralType classification
     */
    [[nodiscard]] static HBondStructuralType classify_structural_type(
        HBondAtomContext donor_context,
        HBondAtomContext acceptor_context,
        HBondMoleculeType mol_type);
    
    /**
     * @brief Fully classify an H-bond given residue information
     */
    static void classify_hbond(
        core::HydrogenBond& hbond,
        const core::poly::IResidue& res1,
        const core::poly::IResidue& res2);
};
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/hbond_classifier.hpp` | create | H-bond classification interface |
| `src/x3dna/algorithms/hydrogen_bond/hbond_classifier.cpp` | create | Classification logic implementation |

#### 1.4 Integrate Classification into HydrogenBondFinder

Modify `HydrogenBondFinder::find_hydrogen_bonds_detailed()` to optionally populate classification fields after detection is complete.

```cpp
// In hydrogen_bond_finder.hpp, add new method:
[[nodiscard]] static DetailedHBondResult find_hydrogen_bonds_detailed_classified(
    const core::poly::IResidue& res1,
    const core::poly::IResidue& res2,
    double hb_lower, double hb_dist1,
    double hb_dist2 = 4.5,
    bool classify = true);
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond_finder.hpp` | modify | Add classified version of find method |
| `src/x3dna/algorithms/hydrogen_bond_finder.cpp` | modify | Call classifier after validation |

#### Phase 1 Testing

- **Compile test**: All existing code compiles without changes
- **Unit tests**: New tests for `AtomContextClassifier` and `HBondClassifier`
  - `tests/unit/algorithms/hydrogen_bond/test_atom_context_classifier.cpp`
  - `tests/unit/algorithms/hydrogen_bond/test_hbond_classifier.cpp`
- **Regression test**: `fp2-validate validate hbonds --test-set 100` still passes
- **New fields test**: Verify new fields are populated correctly for known structures

| Test File | Action | Description |
|-----------|--------|-------------|
| `tests/unit/algorithms/hydrogen_bond/test_atom_context_classifier.cpp` | create | Atom context classification tests |
| `tests/unit/algorithms/hydrogen_bond/test_hbond_classifier.cpp` | create | H-bond classification tests |

---

### Phase 2: Add X-D-A and D-A-X Angle Calculations

**Goal**: Calculate angles at donor and acceptor atoms using bonded neighbors.

#### 2.1 Create Neighbor Finder Utility

**File**: `include/x3dna/algorithms/hydrogen_bond/neighbor_finder.hpp`
**File**: `src/x3dna/algorithms/hydrogen_bond/neighbor_finder.cpp`

```cpp
class NeighborFinder {
public:
    /**
     * @brief Find the covalently bonded neighbor of an atom
     * @param atom The atom to find neighbor for
     * @param residue The residue containing the atom
     * @param exclude_atom Optional atom to exclude (e.g., the H-bond partner)
     * @return Optional neighbor atom position
     *
     * Uses known bonding patterns for nucleotides/proteins or distance-based heuristics.
     */
    [[nodiscard]] static std::optional<geometry::Vector3D> find_neighbor(
        const core::Atom& atom,
        const core::poly::IResidue& residue,
        const std::string& exclude_atom = "");
    
    /**
     * @brief Get known bonded neighbors for common H-bond atoms
     * @param atom_name Atom name (e.g., "N6", "O4")
     * @return Vector of atom names that are typically bonded to this atom
     */
    [[nodiscard]] static std::vector<std::string> get_bonded_atoms(
        const std::string& atom_name);
};
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/neighbor_finder.hpp` | create | Neighbor finding interface |
| `src/x3dna/algorithms/hydrogen_bond/neighbor_finder.cpp` | create | Bonding pattern tables and distance fallback |

#### 2.2 Create Angle Calculator

**File**: `include/x3dna/algorithms/hydrogen_bond/angle_calculator.hpp`
**File**: `src/x3dna/algorithms/hydrogen_bond/angle_calculator.cpp`

```cpp
class HBondAngleCalculator {
public:
    /**
     * @brief Calculate X-D-A angle
     * @param neighbor_pos Position of donor's neighbor (X)
     * @param donor_pos Position of donor (D)
     * @param acceptor_pos Position of acceptor (A)
     * @return Angle in degrees
     */
    [[nodiscard]] static double calculate_xda_angle(
        const geometry::Vector3D& neighbor_pos,
        const geometry::Vector3D& donor_pos,
        const geometry::Vector3D& acceptor_pos);
    
    /**
     * @brief Calculate D-A-X angle
     * @param donor_pos Position of donor (D)
     * @param acceptor_pos Position of acceptor (A)
     * @param neighbor_pos Position of acceptor's neighbor (X)
     * @return Angle in degrees
     */
    [[nodiscard]] static double calculate_dax_angle(
        const geometry::Vector3D& donor_pos,
        const geometry::Vector3D& acceptor_pos,
        const geometry::Vector3D& neighbor_pos);
    
    /**
     * @brief Calculate all angles for an H-bond
     * @param donor_atom Donor atom
     * @param acceptor_atom Acceptor atom
     * @param res1 Residue containing donor
     * @param res2 Residue containing acceptor
     * @return HBondAngles struct with available angles
     */
    [[nodiscard]] static HBondAngles calculate_angles(
        const core::Atom& donor_atom,
        const core::Atom& acceptor_atom,
        const core::poly::IResidue& res1,
        const core::poly::IResidue& res2);
};
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/angle_calculator.hpp` | create | Angle calculation interface |
| `src/x3dna/algorithms/hydrogen_bond/angle_calculator.cpp` | create | Vector angle calculations |

#### 2.3 Integrate Angle Calculation

Modify `HBondClassifier::classify_hbond()` to also calculate angles:

```cpp
static void classify_hbond(
    core::HydrogenBond& hbond,
    const core::Atom& donor_atom,
    const core::Atom& acceptor_atom,
    const core::poly::IResidue& res1,
    const core::poly::IResidue& res2,
    bool calculate_angles = true);
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/hbond_classifier.hpp` | modify | Add angle calculation parameter |
| `src/x3dna/algorithms/hydrogen_bond/hbond_classifier.cpp` | modify | Call angle calculator |

#### Phase 2 Testing

- **Unit tests**: Test angle calculations with known geometries
  - `tests/unit/algorithms/hydrogen_bond/test_angle_calculator.cpp`
  - `tests/unit/algorithms/hydrogen_bond/test_neighbor_finder.cpp`
- **Integration test**: Verify angles match expected values for Watson-Crick pairs
- **Regression test**: `fp2-validate validate core --test-set 100` still passes

| Test File | Action | Description |
|-----------|--------|-------------|
| `tests/unit/algorithms/hydrogen_bond/test_angle_calculator.cpp` | create | Angle calculation tests |
| `tests/unit/algorithms/hydrogen_bond/test_neighbor_finder.cpp` | create | Neighbor finding tests |

---

### Phase 3: Add Dihedral Angle Calculation

**Goal**: Calculate X-D-A-X dihedral angle.

#### 3.1 Extend Angle Calculator

Add dihedral calculation to `HBondAngleCalculator`:

```cpp
/**
 * @brief Calculate X-D-A-X dihedral angle
 * @param donor_neighbor_pos Position of donor's neighbor (X1)
 * @param donor_pos Position of donor (D)
 * @param acceptor_pos Position of acceptor (A)
 * @param acceptor_neighbor_pos Position of acceptor's neighbor (X2)
 * @return Dihedral angle in degrees [-180, 180]
 */
[[nodiscard]] static double calculate_xdax_dihedral(
    const geometry::Vector3D& donor_neighbor_pos,
    const geometry::Vector3D& donor_pos,
    const geometry::Vector3D& acceptor_pos,
    const geometry::Vector3D& acceptor_neighbor_pos);
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/angle_calculator.hpp` | modify | Add dihedral method |
| `src/x3dna/algorithms/hydrogen_bond/angle_calculator.cpp` | modify | Implement dihedral calculation |

#### Phase 3 Testing

- **Unit tests**: Test dihedral calculations with known geometries
- **Integration test**: Verify dihedral angles for standard base pairs

---

### Phase 4: Add Full H-Bond Type Classification

**Goal**: Populate all classification fields for all H-bonds detected.

#### 4.1 Create Classification Post-Processor

Create a utility to classify all H-bonds in a base pair after detection:

```cpp
class HBondPostProcessor {
public:
    /**
     * @brief Classify all H-bonds in a base pair
     * @param hbonds H-bonds to classify (modified in place)
     * @param res1 First residue
     * @param res2 Second residue
     * @param calculate_angles Whether to calculate geometric angles
     */
    static void classify_all(
        std::vector<core::HydrogenBond>& hbonds,
        const core::poly::IResidue& res1,
        const core::poly::IResidue& res2,
        bool calculate_angles = true);
};
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/algorithms/hydrogen_bond/hbond_post_processor.hpp` | create | Post-processing interface |
| `src/x3dna/algorithms/hydrogen_bond/hbond_post_processor.cpp` | create | Batch classification |

#### 4.2 Optional Extended JSON Output

Add optional extended JSON format that includes new fields (only when requested):

```cpp
// In HydrogenBond struct
[[nodiscard]] nlohmann::json to_extended_json() const {
    auto j = nlohmann::json{
        {"donor_atom", donor_atom},
        {"acceptor_atom", acceptor_atom},
        {"distance", distance},
        {"type", std::string(1, type)}
    };
    
    // Add classification fields
    j["molecule_type"] = to_string(molecule_type);
    j["donor_context"] = to_string(donor_context);
    j["acceptor_context"] = to_string(acceptor_context);
    j["structural_type"] = to_string(structural_type);
    
    // Add angles if available
    if (angles) {
        j["angles"] = {
            {"xda_angle", angles->xda_angle},
            {"dax_angle", angles->dax_angle},
            {"xdax_dihedral", angles->xdax_dihedral},
            {"has_donor_neighbor", angles->has_donor_neighbor},
            {"has_acceptor_neighbor", angles->has_acceptor_neighbor}
        };
    }
    
    return j;
}
```

| File | Action | Description |
|------|--------|-------------|
| `include/x3dna/core/hydrogen_bond.hpp` | modify | Add to_extended_json() method |
| `include/x3dna/core/hydrogen_bond_types.hpp` | modify | Add to_string() functions for enums |

#### Phase 4 Testing

- **Full integration test**: Verify all H-bonds are classified correctly
- **Regression test**: Legacy JSON output unchanged
- **New output test**: Extended JSON format correct

---

## Files Summary

| File | Phase | Action | Description |
|------|-------|--------|-------------|
| `include/x3dna/core/hydrogen_bond_types.hpp` | 1 | create | Enums and HBondAngles struct |
| `include/x3dna/core/hydrogen_bond.hpp` | 1,4 | modify | Add new fields and methods |
| `include/x3dna/algorithms/hydrogen_bond/atom_context_classifier.hpp` | 1 | create | Atom context classification |
| `src/x3dna/algorithms/hydrogen_bond/atom_context_classifier.cpp` | 1 | create | Implementation |
| `include/x3dna/algorithms/hydrogen_bond/hbond_classifier.hpp` | 1 | create | H-bond type classification |
| `src/x3dna/algorithms/hydrogen_bond/hbond_classifier.cpp` | 1 | create | Implementation |
| `include/x3dna/algorithms/hydrogen_bond/neighbor_finder.hpp` | 2 | create | Neighbor finding utility |
| `src/x3dna/algorithms/hydrogen_bond/neighbor_finder.cpp` | 2 | create | Implementation |
| `include/x3dna/algorithms/hydrogen_bond/angle_calculator.hpp` | 2,3 | create | Angle calculations |
| `src/x3dna/algorithms/hydrogen_bond/angle_calculator.cpp` | 2,3 | create | Implementation |
| `include/x3dna/algorithms/hydrogen_bond/hbond_post_processor.hpp` | 4 | create | Batch post-processing |
| `src/x3dna/algorithms/hydrogen_bond/hbond_post_processor.cpp` | 4 | create | Implementation |
| `include/x3dna/algorithms/hydrogen_bond_finder.hpp` | 1 | modify | Add classified find method |
| `src/x3dna/algorithms/hydrogen_bond_finder.cpp` | 1 | modify | Integration |
| `tests/unit/algorithms/hydrogen_bond/test_atom_context_classifier.cpp` | 1 | create | Unit tests |
| `tests/unit/algorithms/hydrogen_bond/test_hbond_classifier.cpp` | 1 | create | Unit tests |
| `tests/unit/algorithms/hydrogen_bond/test_neighbor_finder.cpp` | 2 | create | Unit tests |
| `tests/unit/algorithms/hydrogen_bond/test_angle_calculator.cpp` | 2,3 | create | Unit tests |

---

## Testing Strategy Per Phase

### Phase 1 Tests

1. **Compile Test**: `make release` succeeds
2. **Unit Test**: AtomContextClassifier correctly classifies:
   - Base atoms: N1, N2, N3, N4, N6, N7, O2, O4, O6 -> BASE
   - Sugar atoms: O2', O3', O4', C1', C2', etc. -> SUGAR
   - Backbone atoms: O1P, O2P, O5' -> BACKBONE
3. **Unit Test**: HBondClassifier correctly classifies:
   - RNA+RNA -> RNA_RNA
   - RNA+Protein -> RNA_PROTEIN
   - BASE+BASE -> BASE_BASE
   - BASE+SUGAR -> BASE_SUGAR
4. **Regression Test**: `fp2-validate validate hbonds --test-set 100`
5. **Integration Test**: New fields populated for 1EHZ (known tRNA)

### Phase 2 Tests

1. **Unit Test**: X-D-A angle calculation matches geometric expectation
   - Test with known triangle coordinates
2. **Unit Test**: D-A-X angle calculation correct
3. **Unit Test**: Neighbor finder returns correct bonded atoms for:
   - N6 (bonded to C6)
   - O4 (bonded to C4)
   - O2' (bonded to C2')
4. **Integration Test**: Angles for Watson-Crick G-C pair in range [150-180] degrees

### Phase 3 Tests

1. **Unit Test**: Dihedral angle calculation matches known geometry
2. **Integration Test**: Dihedral angles for standard base pairs

### Phase 4 Tests

1. **Full Regression**: `fp2-validate validate core --test-set fast`
2. **Extended JSON**: New format includes all fields
3. **Legacy JSON**: Format unchanged

---

## Risk Mitigation

### Risk 1: Breaking Legacy JSON Output

**Mitigation**: 
- Keep all new fields as `optional` or with default values
- Never modify existing `to_json_legacy()` method
- Add new `to_extended_json()` method for new format
- Run regression tests after every change

### Risk 2: Performance Impact

**Mitigation**:
- Classification is optional (default off for legacy compatibility)
- Angle calculation is optional (can be enabled separately)
- Use lazy evaluation where possible
- Profile before and after with large structures (1VQ5)

### Risk 3: Neighbor Finding Accuracy

**Mitigation**:
- Use known bonding patterns from nucleotide/amino acid chemistry
- Fall back to distance-based heuristics (< 2.0 Angstroms) when patterns unknown
- Mark angles as unavailable rather than guess incorrectly

### Risk 4: Atom Name Format Variations

**Mitigation**:
- Use existing `trim()` utility for atom name normalization
- Handle both PDB v2 and v3 naming conventions
- Test with diverse PDB files (modified nucleotides, unusual names)

---

## Implementation Steps (Detailed)

### Step 1: Create hydrogen_bond_types.hpp
- [ ] Define HBondMoleculeType enum
- [ ] Define HBondAtomContext enum
- [ ] Define HBondStructuralType enum
- [ ] Define HBondAngles struct
- **Test**: Compiles

### Step 2: Extend HydrogenBond struct
- [ ] Add #include for hydrogen_bond_types.hpp
- [ ] Add optional classification fields with defaults
- [ ] Add convenience query methods
- **Test**: Existing code compiles, regression passes

### Step 3: Create AtomContextClassifier
- [ ] Create header with static methods
- [ ] Implement is_base_atom() using existing hydrogen_bond_utils
- [ ] Implement is_sugar_atom() with pattern matching
- [ ] Implement is_backbone_atom() with pattern matching
- [ ] Implement classify() combining the above
- **Test**: Unit test passes

### Step 4: Create HBondClassifier
- [ ] Create header with static methods
- [ ] Implement classify_molecule_type() using IResidue interface
- [ ] Implement classify_structural_type() combining contexts
- [ ] Implement classify_hbond() as convenience wrapper
- **Test**: Unit test passes

### Step 5: Integrate classification into HydrogenBondFinder
- [ ] Add find_hydrogen_bonds_detailed_classified() method
- [ ] Call HBondClassifier after validation step
- [ ] Ensure legacy method unchanged
- **Test**: Regression passes, new method returns classified H-bonds

### Step 6: Create NeighborFinder
- [ ] Define bonding pattern tables for nucleotides
- [ ] Implement get_bonded_atoms() lookup
- [ ] Implement find_neighbor() with pattern + distance fallback
- **Test**: Unit test passes

### Step 7: Create HBondAngleCalculator (XDA, DAX)
- [ ] Implement calculate_xda_angle() using vector math
- [ ] Implement calculate_dax_angle() using vector math
- [ ] Implement calculate_angles() combining above with neighbor finding
- **Test**: Unit test passes

### Step 8: Add dihedral calculation
- [ ] Implement calculate_xdax_dihedral()
- [ ] Integrate into calculate_angles()
- **Test**: Unit test passes

### Step 9: Create HBondPostProcessor
- [ ] Implement classify_all() for batch processing
- [ ] Integrate into pipeline (optional step)
- **Test**: Integration test passes

### Step 10: Add extended JSON output
- [ ] Add to_string() for each enum
- [ ] Add to_extended_json() method
- [ ] Ensure to_json_legacy() unchanged
- **Test**: Legacy regression passes, extended format correct

---

## Success Criteria

1. All existing tests pass (`make test`, `fp2-validate validate core --test-set fast`)
2. New classification fields populated for classified H-bonds
3. Angles calculated and within expected ranges for Watson-Crick pairs
4. Legacy JSON output byte-for-byte identical
5. No significant performance regression (< 5% slowdown)
6. New unit tests achieve > 90% coverage of new code
