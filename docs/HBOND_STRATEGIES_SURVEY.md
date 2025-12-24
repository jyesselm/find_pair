# H-Bond Tracking Strategies: Survey of Existing Programs

This document summarizes how four different programs handle hydrogen bond detection and optimal bond selection for nucleic acids. The goal is to inform the design of a new H-bond tracking system that respects chemical capacity limits (e.g., amines can donate 2, sp2 oxygens can accept 2).

## Executive Summary

| Program | Saturation Tracking | Best Bond Selection | Approach |
|---------|---------------------|---------------------|----------|
| **HBPLUS** | Capacity tables defined, but NOT enforced during detection | None - reports ALL valid H-bonds | Chemistry-aware capacity tables |
| **FR3D** | Edge-based (W/H/S), not per-atom | Exemplar distance comparison | Pre-classification, edge exclusivity |
| **barnaba** | None | None - simple counting | Distance-only, no geometry |
| **ClaRNA** | None | Classifier score | Pattern matching, not explicit H-bonds |

**Key Finding**: None of these programs implement the optimal H-bond selection system you're looking for. HBPLUS comes closest with its capacity tables but doesn't enforce them during detection.

---

## 1. HBPLUS

**Location**: `/Users/jyesselman2/local/installs/hbplus/`

### Approach
HBPLUS is the most sophisticated for general H-bond detection. It uses explicit **capacity tables** that define how many H-bonds each atom type can form.

### Donor/Acceptor Identification

Defined in lookup tables by residue type (`hbp_inpdb.c` lines 336-547):

```c
/* Acceptor capacity per atom type */
short accepts[MAXNAA][TOTNATM] = {
  {0, 0, 0, 2, 0, 2, 0, ...},  /* ALA: backbone O=2, OXT=2 */
  {0, 0, 0, 2, 0, 1, 2, ...},  /* CYS: O=2, SG=1, OXT=2 */
  ...
};

/* Donor capacity per atom type */
short donors[MAXNAA][TOTNATM] = {
  {1, 0, 0, 0, 0, 0, ...},     /* ALA: backbone N=1 */
  {1, 0, 0, 0, 0, 0, 0, 0, 3}, /* LYS: N=1, NZ=3 (primary amine) */
  ...
};
```

**Chemical accuracy examples**:
- Carbonyl O (sp2): accepts = 2
- Primary amine NH3 (LYS NZ): donors = 3
- Amide NH2 (ASN ND2): donors = 2
- Hydroxyl (SER OG): donors = 1, accepts = 2

### Saturation Tracking

**Tracked but NOT enforced.** HBPLUS counts bonds per atom:

```c
struct pdbatm {
    int ndonhb, nacchb;  /* count of H-bonds as donor/acceptor */
};
```

But these counters are updated AFTER detection, not used to limit during detection. All valid H-bonds are reported regardless of capacity.

### Geometry Validation

Comprehensive angle/distance checks (`hbp_hhb.c` lines 271-298):

| Criterion | Threshold |
|-----------|-----------|
| H...A distance | ≤ 2.5 Å |
| D...A distance | ≤ 3.9 Å (initial filter) |
| D-H-A angle | ≥ 90° |
| D-A-AA angle | ≥ 90° |
| H-A-AA angle | ≥ 90° |
| Aromatic perpendicular | ≤ 20° from normal |

### Best Bond Selection

**None.** Reports all valid H-bonds. No mutual exclusivity or greedy selection.

### Relevance for Your Implementation

HBPLUS provides the **best template for capacity tables**. You could use similar per-residue/per-atom capacity definitions, but add enforcement during detection.

---

## 2. FR3D (Find RNA 3D)

**Location**: `/Users/jyesselman2/local/installs/FR3D/`

### Approach
FR3D uses a **classification-first** approach. It first classifies the base pair type based on geometry, then validates that expected H-bonds exist.

### Donor/Acceptor Identification

**Pre-defined by pair classification**, not by scanning atoms. The function `zCheckHydrogen.m` hard-codes which atoms form H-bonds for each pair type:

```matlab
function [Hydrogen] = zCheckHydrogen(NT1,NT2,Class)
Paircode = 4*(Code2-1) + Code1;  % AA=1, CA=2, etc.

switch Paircode
  case 1,  % A-A
    switch Class
      case 1,  % cWW
        Hydrogen(1).Angle = zAngle(N1.Fit(6,:),N1.Fit(14,:),N2.Fit(4,:));
        Hydrogen(1).Distance = zDistance(N1.Fit(14,:),N2.Fit(4,:));
        ...
```

### Saturation Tracking

**Edge-based, not per-atom.** Each nucleotide has 3 edges (Watson-Crick, Hoogsteen, Sugar), each can participate in at most ONE base pair.

From `zEdgeMakesMultiplePairs.m`:
```matlab
for b = 1:3,                     % which edge (W, H, S)
  i = find(sum(w'==b) > 1);      % nucleotides using an edge twice
  % ... remove pairs with larger distance to exemplar
```

### Best Bond Selection

When edge conflicts occur, FR3D keeps the pair with **smallest distance to exemplar structure**:

```matlab
d(c) = zDistanceToExemplar(Exemplar,NT1,NT2,Class);
% ... keep pair with min(d)
```

### Geometry Validation

| Condition | Distance | Angle |
|-----------|----------|-------|
| With explicit H | ≤ 4.0 Å | ≥ 110° |
| Heavy atoms only | ≤ 4.5 Å | N/A |
| Gap cutoff | ≤ 2.0 Å | - |

### Relaxed Requirements

FR3D allows missing H-bonds:
```matlab
if (length(Pair.Hydrogen) == 4),  % 3 out of 4 is good enough
  gh = 3;
elseif (length(Pair.Hydrogen) == 3),  % 2 out of 3 is good enough
  gh = 2;
```

### Relevance for Your Implementation

FR3D's **edge-based exclusivity** is interesting - it's a higher-level constraint than per-atom saturation. The **exemplar distance for conflict resolution** is elegant but requires reference structures.

---

## 3. barnaba

**Location**: `/Users/jyesselman2/local/installs/barnaba/`

### Approach
barnaba uses the **simplest approach** - distance-only filtering with static donor/acceptor definitions.

### Donor/Acceptor Identification

Static per-base definitions in `definitions.py`:

```python
donors = {"A": ["N6","C2","C8","O2'"],  # Note: includes C-H donors
          "G": ["N1","N2","C8","O2'"],
          ...}
acceptors = {"A": ["N1","N3","N7","O2'"],
             "G": ["O6","N3","N7","O2'"],
             ...}
```

**Note**: barnaba treats C-H groups as donors (weak interactions).

### Saturation Tracking

**None.** Each residue pair is evaluated independently. The same donor/acceptor can be counted multiple times.

### Best Bond Selection

**None.** Simply counts all donor-acceptor pairs within distance threshold:

```python
combo_list = list(itertools.product(r1_donor, r2_acceptor)) + \
             list(itertools.product(r1_acceptor, r2_donor))
dist_sq = np.sum(delta**2, axis=2)
n_hbonds = (dist_sq < 0.1089).sum()  # 3.3 Å squared
```

### Geometry Validation

**Distance only** - 3.3 Å cutoff, no angle validation for H-bonds.

Base pairing has plane angle check (60° coplanarity).

### Relevance for Your Implementation

barnaba shows what NOT to do if you want accurate H-bond tracking. It's optimized for speed over accuracy.

---

## 4. ClaRNA

**Location**: `/Users/jyesselman2/local/installs/ClaRNA_download/`

### Approach
ClaRNA is fundamentally different - it's a **structural pattern matcher** using machine learning classifiers, not an explicit H-bond detector.

### Donor/Acceptor Identification

Pre-defined interaction pairs in `distances.py`:

```python
PH_INTERACTIONS = {
    'A': (('C2','H2'), ('C8','H8'), ('N6','1H6'), ('N6','2H6')),
    'G': (('N1','H1'), ('C8','H8'), ('N2','1H2'), ('N2','2H2')),
    ...
}
PH_OXYGENS = ["OP1","OP2","O5'","NEXT:O3'"]  # Phosphate oxygens
BR_OXYGENS = ["O2'","O3'","O4'"]              # Ribose oxygens
```

### Saturation Tracking

**None.** No mechanism to track or limit bonds per atom.

### Best Bond Selection

**Classifier score** - pattern matching to known structures:

```python
def aggregate_results(res):
    for r in res:
        if not ares.has_key(d) or r.score > ares[d].score:
            ares[d] = r  # Keep highest scoring match
```

### Geometry Validation

| Donor Type | Distance | Angle |
|------------|----------|-------|
| C-H | ≤ 4.0 Å | > 130° |
| N-H | ≤ 3.5 Å | > 130° |

### Relevance for Your Implementation

ClaRNA's approach is incompatible with explicit H-bond tracking. However, the interaction definitions could be useful for defining which atom pairs to consider.

---

## Recommendations for Your Implementation

Based on this survey, here's a proposed approach:

### 1. Define Capacity Tables (like HBPLUS, but for nucleic acids)

```python
DONOR_CAPACITY = {
    # Primary amines (NH2) - can donate 2
    ('A', 'N6'): 2,   # Adenine amino
    ('C', 'N4'): 2,   # Cytosine amino
    ('G', 'N2'): 2,   # Guanine amino

    # Imino (NH) - can donate 1
    ('G', 'N1'): 1,   # Guanine N1-H
    ('U', 'N3'): 1,   # Uracil N3-H
    ('T', 'N3'): 1,   # Thymine N3-H

    # C-H donors (weak) - can donate 1
    ('A', 'C2'): 1, ('A', 'C8'): 1,
    ('G', 'C8'): 1,
    ('C', 'C5'): 1, ('C', 'C6'): 1,
    ('U', 'C5'): 1, ('U', 'C6'): 1,
}

ACCEPTOR_CAPACITY = {
    # sp2 oxygens - can accept 2
    ('G', 'O6'): 2,   # Guanine O6
    ('U', 'O2'): 2, ('U', 'O4'): 2,
    ('C', 'O2'): 2,

    # Ring nitrogens - can accept 1-2
    ('A', 'N1'): 2, ('A', 'N3'): 2, ('A', 'N7'): 2,
    ('G', 'N3'): 2, ('G', 'N7'): 2,
    ('C', 'N3'): 2,
}
```

### 2. Track Usage During Detection

```python
class HBondTracker:
    def __init__(self):
        self.donor_usage = {}    # (res_id, atom_name) -> count
        self.acceptor_usage = {} # (res_id, atom_name) -> count

    def can_donate(self, res_id: str, atom_name: str, base_type: str) -> bool:
        key = (res_id, atom_name)
        used = self.donor_usage.get(key, 0)
        capacity = DONOR_CAPACITY.get((base_type, atom_name), 0)
        return used < capacity

    def can_accept(self, res_id: str, atom_name: str, base_type: str) -> bool:
        key = (res_id, atom_name)
        used = self.acceptor_usage.get(key, 0)
        capacity = ACCEPTOR_CAPACITY.get((base_type, atom_name), 0)
        return used < capacity

    def record_hbond(self, donor_res: str, donor_atom: str,
                     acceptor_res: str, acceptor_atom: str):
        self.donor_usage[(donor_res, donor_atom)] = \
            self.donor_usage.get((donor_res, donor_atom), 0) + 1
        self.acceptor_usage[(acceptor_res, acceptor_atom)] = \
            self.acceptor_usage.get((acceptor_res, acceptor_atom), 0) + 1
```

### 3. Greedy Best-First Selection

Process candidate H-bonds in quality order:
1. Rank all potential H-bonds by quality score (distance, angle, energy)
2. Take best candidate
3. Check if donor AND acceptor have remaining capacity
4. If yes, record the bond and update usage counts
5. Repeat until no valid candidates remain

```python
def select_optimal_hbonds(candidates: List[HBondCandidate],
                          tracker: HBondTracker) -> List[HBond]:
    # Sort by quality (best first)
    candidates.sort(key=lambda c: c.quality_score, reverse=True)

    selected = []
    for cand in candidates:
        if (tracker.can_donate(cand.donor_res, cand.donor_atom, cand.donor_base) and
            tracker.can_accept(cand.acceptor_res, cand.acceptor_atom, cand.acceptor_base)):

            tracker.record_hbond(cand.donor_res, cand.donor_atom,
                                cand.acceptor_res, cand.acceptor_atom)
            selected.append(cand.to_hbond())

    return selected
```

### 4. Geometry Validation (combine best practices)

| Criterion | Threshold | Source |
|-----------|-----------|--------|
| H...A distance | ≤ 2.5 Å (with H) | HBPLUS |
| D...A distance | ≤ 3.5 Å (heavy only) | DSSR-like |
| D-H-A angle | ≥ 120° | FR3D-inspired |
| D-A-AA angle | ≥ 90° | HBPLUS |

### 5. Angular Separation for Multi-Bond Atoms

**Critical insight**: When an atom can donate/accept multiple H-bonds, the bonds must be angularly separated based on hybridization geometry. Two H-bonds to the same atom can't be "stacked" on each other.

#### For Donors (multiple H atoms)

| Hybridization | Example | Expected Separation |
|---------------|---------|---------------------|
| sp2 NH2 | Adenine N6, Cytosine N4, Guanine N2 | ~120° between H atoms |
| sp3 NH3 | Lysine NZ | ~109.5° tetrahedral |

When checking if a second donation from an NH2 is valid:
```python
def check_donor_separation(donor_atom, existing_acceptor, new_acceptor, hybridization):
    """Check that two H-bonds from same donor are properly separated."""
    # Vector from donor to existing acceptor
    v1 = existing_acceptor.pos - donor_atom.pos
    # Vector from donor to new acceptor
    v2 = new_acceptor.pos - donor_atom.pos

    angle = angle_between(v1, v2)

    if hybridization == 'sp2':
        # Two H's are ~120° apart, so acceptors should be roughly opposite
        # Allow some tolerance: 60° - 180° is reasonable
        return angle >= 60.0
    elif hybridization == 'sp3':
        # Tetrahedral: ~109.5° separation
        return angle >= 70.0

    return True
```

#### For Acceptors (lone pairs)

| Hybridization | Example | Lone Pair Geometry |
|---------------|---------|-------------------|
| sp2 carbonyl O | Guanine O6, Uracil O2/O4 | Two lone pairs in plane, ~120° apart from C=O |
| sp2 ring N | Adenine N1/N3/N7 | One lone pair in plane |
| sp3 hydroxyl O | Ribose O2' | Two lone pairs tetrahedral |

When checking if a second acceptance is valid:
```python
def check_acceptor_separation(acceptor_atom, antecedent_atom,
                               existing_donor, new_donor, hybridization):
    """Check that two H-bonds to same acceptor use different lone pairs."""
    # For sp2 oxygen, lone pairs are in the plane perpendicular to C=O
    # and roughly 120° apart from each other

    v_antecedent = antecedent_atom.pos - acceptor_atom.pos  # e.g., C=O direction
    v1 = existing_donor.pos - acceptor_atom.pos
    v2 = new_donor.pos - acceptor_atom.pos

    # Angle between the two donor directions
    donor_angle = angle_between(v1, v2)

    if hybridization == 'sp2':
        # Lone pairs are ~120° apart, donors should approach from different directions
        # Minimum separation of ~60° is reasonable
        return donor_angle >= 60.0
    elif hybridization == 'sp3':
        # Tetrahedral lone pairs
        return donor_angle >= 70.0

    return True
```

### 6. Predicting H and Lone Pair Positions from Geometry

**Better approach**: Instead of just checking angles, compute the actual positions of hydrogens and lone pairs from known geometry, then check alignment.

#### Hydrogen Position Prediction

For donors, we can compute where the H atoms should be:

```python
import numpy as np

def predict_hydrogen_positions(donor_atom, antecedent_atoms, hybridization, base_normal=None):
    """
    Predict hydrogen positions based on hybridization and attached atoms.

    Args:
        donor_atom: Position of the donor (e.g., N6)
        antecedent_atoms: List of atoms bonded to donor (e.g., [C6] for N6)
        hybridization: 'sp2' or 'sp3'
        base_normal: Normal vector to base plane (for sp2 atoms)

    Returns:
        List of predicted H position unit vectors (relative to donor)
    """
    H_BOND_LENGTH = 1.01  # N-H bond length in Angstroms

    if hybridization == 'sp2' and len(antecedent_atoms) == 1:
        # sp2 NH2: Two H atoms in plane, 120° from C-N bond
        # Example: Adenine N6 bonded to C6

        c_to_n = normalize(donor_atom - antecedent_atoms[0])  # C6 -> N6 direction

        # Rotate 120° and -120° around base normal to get H directions
        h1_dir = rotate_vector(c_to_n, base_normal, 120.0)
        h2_dir = rotate_vector(c_to_n, base_normal, -120.0)

        return [h1_dir * H_BOND_LENGTH, h2_dir * H_BOND_LENGTH]

    elif hybridization == 'sp2' and len(antecedent_atoms) == 2:
        # sp2 imino NH: Single H opposite to the two ring atoms
        # Example: Guanine N1 bonded to C2 and C6

        # Average of the two antecedent directions, then negate
        avg_dir = normalize(antecedent_atoms[0] + antecedent_atoms[1] - 2*donor_atom)
        h_dir = -avg_dir  # H points away from ring

        return [h_dir * H_BOND_LENGTH]

    elif hybridization == 'sp2' and len(antecedent_atoms) == 1:
        # C-H: Single H opposite to the attached atom, in plane
        c_to_h = -normalize(antecedent_atoms[0] - donor_atom)
        return [c_to_h * H_BOND_LENGTH]

    return []


def rotate_vector(v, axis, angle_deg):
    """Rotate vector v around axis by angle_deg degrees."""
    angle_rad = np.radians(angle_deg)
    axis = normalize(axis)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)

    # Rodrigues rotation formula
    return v * cos_a + np.cross(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)
```

#### Lone Pair Position Prediction

For acceptors, compute where the lone pairs point:

```python
def predict_lone_pair_directions(acceptor_atom, antecedent_atoms, hybridization, base_normal=None):
    """
    Predict lone pair directions based on hybridization.

    Args:
        acceptor_atom: Position of acceptor (e.g., O6)
        antecedent_atoms: Atoms bonded to acceptor (e.g., [C6] for O6)
        hybridization: 'sp2' or 'sp3'
        base_normal: Normal to base plane

    Returns:
        List of lone pair direction unit vectors
    """
    if hybridization == 'sp2' and len(antecedent_atoms) == 1:
        # sp2 carbonyl oxygen: Two lone pairs in plane, 120° from C=O
        # Example: Guanine O6 bonded to C6

        c_to_o = normalize(acceptor_atom - antecedent_atoms[0])  # C6 -> O6

        # Lone pairs are 120° from C=O bond, in the base plane
        lp1_dir = rotate_vector(c_to_o, base_normal, 120.0)
        lp2_dir = rotate_vector(c_to_o, base_normal, -120.0)

        return [lp1_dir, lp2_dir]

    elif hybridization == 'sp2' and len(antecedent_atoms) == 2:
        # sp2 ring nitrogen: One lone pair pointing out of ring
        # Example: Adenine N1 bonded to C2 and C6

        # Lone pair points away from the ring (opposite to average of bonds)
        avg_bond = normalize(antecedent_atoms[0] + antecedent_atoms[1] - 2*acceptor_atom)
        lp_dir = -avg_bond

        return [lp_dir]

    elif hybridization == 'sp3':
        # sp3 oxygen: Two lone pairs in tetrahedral arrangement
        # More complex - need full geometry
        # For ribose O2', the lone pairs point roughly perpendicular to C-O-H plane
        pass

    return []
```

#### Scoring H-Bond Alignment

Now we can score how well an H-bond aligns with predicted positions:

```python
def score_hbond_alignment(donor_atom, acceptor_atom,
                          h_positions, lp_directions):
    """
    Score how well an H-bond aligns with predicted H and lone pair positions.

    Returns:
        (best_h_idx, best_lp_idx, alignment_score)
        where alignment_score is 0.0 (perfect) to 1.0 (worst)
    """
    donor_to_acceptor = normalize(acceptor_atom - donor_atom)
    acceptor_to_donor = -donor_to_acceptor

    # Find best matching hydrogen
    best_h_score = -1.0
    best_h_idx = 0
    for i, h_pos in enumerate(h_positions):
        h_dir = normalize(h_pos)  # Direction from donor to H
        # Acceptor should be roughly along H direction
        alignment = np.dot(h_dir, donor_to_acceptor)
        if alignment > best_h_score:
            best_h_score = alignment
            best_h_idx = i

    # Find best matching lone pair
    best_lp_score = -1.0
    best_lp_idx = 0
    for i, lp_dir in enumerate(lp_directions):
        # Donor should approach along lone pair direction
        alignment = np.dot(lp_dir, acceptor_to_donor)
        if alignment > best_lp_score:
            best_lp_score = alignment
            best_lp_idx = i

    # Combined score (higher is better, max = 2.0 for perfect alignment)
    total_score = best_h_score + best_lp_score

    return best_h_idx, best_lp_idx, total_score


def is_valid_hbond(donor_atom, acceptor_atom, h_positions, lp_directions,
                   min_alignment=0.5):
    """
    Check if H-bond has valid geometry by aligning with predicted H and LP.

    Args:
        min_alignment: Minimum dot product (0.5 = 60° max deviation)
    """
    best_h_idx, best_lp_idx, score = score_hbond_alignment(
        donor_atom, acceptor_atom, h_positions, lp_directions
    )

    # Both alignments should be reasonable
    # score of 1.0 means average alignment of 0.5 (60°) which is acceptable
    return score >= 2 * min_alignment
```

#### Tracking Which H/LP Slots Are Used

```python
class HBondSlotTracker:
    """Track which hydrogen and lone pair slots have been used."""

    def __init__(self):
        # Key: (res_id, atom_name), Value: set of used slot indices
        self.used_h_slots = {}      # Which H atoms have donated
        self.used_lp_slots = {}     # Which lone pairs have accepted

    def get_available_h_slots(self, res_id, atom_name, total_slots):
        """Return indices of H slots not yet used."""
        key = (res_id, atom_name)
        used = self.used_h_slots.get(key, set())
        return [i for i in range(total_slots) if i not in used]

    def get_available_lp_slots(self, res_id, atom_name, total_slots):
        """Return indices of LP slots not yet used."""
        key = (res_id, atom_name)
        used = self.used_lp_slots.get(key, set())
        return [i for i in range(total_slots) if i not in used]

    def mark_h_used(self, res_id, atom_name, slot_idx):
        key = (res_id, atom_name)
        if key not in self.used_h_slots:
            self.used_h_slots[key] = set()
        self.used_h_slots[key].add(slot_idx)

    def mark_lp_used(self, res_id, atom_name, slot_idx):
        key = (res_id, atom_name)
        if key not in self.used_lp_slots:
            self.used_lp_slots[key] = set()
        self.used_lp_slots[key].add(slot_idx)


def select_optimal_hbonds_with_slots(candidates, tracker, residue_data):
    """
    Select optimal H-bonds considering slot availability.

    Each NH2 has 2 H slots, each sp2 O has 2 LP slots, etc.
    """
    # Sort by quality (best first)
    candidates.sort(key=lambda c: c.distance)  # Shorter = better

    selected = []
    slot_tracker = HBondSlotTracker()

    for cand in candidates:
        # Get predicted H positions for donor
        donor_info = residue_data[cand.donor_res]
        h_positions = predict_hydrogen_positions(
            cand.donor_pos,
            donor_info.antecedents[cand.donor_atom],
            donor_info.hybridization[cand.donor_atom],
            donor_info.base_normal
        )

        # Get predicted LP directions for acceptor
        acceptor_info = residue_data[cand.acceptor_res]
        lp_directions = predict_lone_pair_directions(
            cand.acceptor_pos,
            acceptor_info.antecedents[cand.acceptor_atom],
            acceptor_info.hybridization[cand.acceptor_atom],
            acceptor_info.base_normal
        )

        # Find best alignment
        best_h, best_lp, score = score_hbond_alignment(
            cand.donor_pos, cand.acceptor_pos, h_positions, lp_directions
        )

        # Check if those slots are available
        available_h = slot_tracker.get_available_h_slots(
            cand.donor_res, cand.donor_atom, len(h_positions)
        )
        available_lp = slot_tracker.get_available_lp_slots(
            cand.acceptor_res, cand.acceptor_atom, len(lp_directions)
        )

        if best_h in available_h and best_lp in available_lp and score >= 1.0:
            # Valid H-bond! Mark slots as used
            slot_tracker.mark_h_used(cand.donor_res, cand.donor_atom, best_h)
            slot_tracker.mark_lp_used(cand.acceptor_res, cand.acceptor_atom, best_lp)
            selected.append(cand)

    return selected
```

#### Visual Example

```
Adenine N6 (NH2) donating to two acceptors:

         A1 (acceptor 1)
          \
           \ H1 slot (120° from C6-N6)
            \
    C6 ---- N6
            /
           / H2 slot (120° from C6-N6, opposite side)
          /
         A2 (acceptor 2)

- First H-bond: A1 aligns best with H1 slot → mark H1 used
- Second H-bond: A2 must use H2 slot (H1 already used)
- Third candidate: REJECTED (no available H slots)
```

```
Guanine O6 accepting from two donors:

    D1 (donor 1)
      \
       \ LP1 slot (120° from C6=O6)
        \
         O6 ==== C6
        /
       / LP2 slot (120° from C6=O6, opposite side)
      /
    D2 (donor 2)

- First H-bond: D1 aligns best with LP1 → mark LP1 used
- Second H-bond: D2 must use LP2 slot
- Third candidate: REJECTED (no available LP slots)
```

---

## Summary of Key Insights

1. **HBPLUS has the right chemistry** with per-atom capacity tables, but doesn't enforce them during detection

2. **FR3D uses edge-based exclusivity** which is a higher-level constraint appropriate for base pairs but not individual H-bonds

3. **barnaba and ClaRNA are too simple** for accurate H-bond tracking

4. **No existing tool does exactly what you need** - optimal selection with capacity enforcement

5. **The greedy approach with capacity tracking** is the recommended design pattern
