# Bifurcated H-Bond Support Plan

## What is Bifurcation?

Bifurcated H-bonds occur when:
1. **Bifurcated donor**: One H atom donates to TWO acceptors simultaneously
2. **Bifurcated acceptor**: One lone pair accepts from TWO donors simultaneously

This is NOT oversaturation - it's a specific geometric arrangement where partners share.

## Example: G-U Wobble Pair

```
      G                    U
      |                    |
   N1-H -----> O2 <----- H-N2
      |         ^          |
   O6 <-------- H-N3
```

DSSR reports 3 H-bonds:
- N1-H → O2 (standard)
- N3-H → O6 (standard)
- N2-H → O2 (bifurcated with N1-H → O2, sharing O2's lone pair)

## Current Problem

Our slot system marks O2's LP0 as "used" after N1→O2, then rejects N2→O2 because:
- LP0 is used
- LP1 has poor alignment (0.23 < 0.3 threshold)

## Proposed Solution

### 1. Allow LP Slot Sharing with Angular Check

Instead of binary "used/free", track the direction of incoming H-bonds per slot:

```python
@dataclass
class LPSlot:
    direction: np.ndarray
    bonds: List[Tuple[str, str, np.ndarray]]  # List of (donor_res, donor_atom, H_direction)
    max_bonds: int = 2  # Bifurcation limit
```

### 2. Bifurcation Geometry Criteria

For an LP slot to accept a second H-bond:
- Angle between the two donor directions must be >= 60° (not stacking on same vector)
- Each individual alignment must still be reasonable (>= 0.0)
- Combined alignment score has lower threshold for bifurcated bonds

```python
def can_bifurcate(self, existing_h_dir: np.ndarray, new_h_dir: np.ndarray) -> bool:
    """Check if two H-bonds can share this LP slot."""
    angle = angle_between(existing_h_dir, new_h_dir)
    return 60.0 <= angle <= 180.0  # Must be angularly separated
```

### 3. Updated Selection Algorithm

```python
def select_optimal(self, candidates):
    # Sort by distance (shortest first)
    candidates.sort(key=lambda c: c.distance)

    selected = []
    for c in candidates:
        h_slots = get_h_slots(c.donor_atom)
        lp_slots = get_lp_slots(c.acceptor_atom)

        # Try to find available or bifurcatable slot pair
        for h_slot in h_slots:
            if h_slot.is_full():
                continue

            for lp_slot in lp_slots:
                # Check if slot is available OR can bifurcate
                if lp_slot.is_available():
                    if alignment_ok(h_slot, lp_slot, c):
                        accept_bond(c, h_slot, lp_slot)
                        break
                elif lp_slot.can_add_bifurcated(c.donor_direction):
                    if bifurcated_alignment_ok(h_slot, lp_slot, c):
                        accept_bifurcated_bond(c, h_slot, lp_slot)
                        break
```

### 4. Slot Capacity Rules

| Atom Type | H Slots | Max per H | LP Slots | Max per LP |
|-----------|---------|-----------|----------|------------|
| NH2 (N6, N4, N2) | 2 | 1 each | - | - |
| Imino NH (N1, N3) | 1 | 1 | - | - |
| Carbonyl O (O6, O2, O4) | - | - | 2 | 2 each (bifurcation) |
| Ring N acceptor | - | - | 1 | 2 (bifurcation) |
| Phosphate O | - | - | 3 | 2 each |

### 5. Implementation Steps

1. **Update LPSlot class**:
   - Add `bonds` list to track current H-bonds using this slot
   - Add `can_add_bond()` method with bifurcation check

2. **Update HSlot class**:
   - Similar changes for donor-side bifurcation (less common)

3. **Update score_hbond_alignment()**:
   - Add bifurcation-aware scoring
   - Lower threshold for second bond to same slot

4. **Update select_optimal()**:
   - Check bifurcation possibility when slot is "used"
   - Track angular separation of bonds

5. **Add configuration**:
   - `allow_bifurcation: bool = True`
   - `min_bifurcation_angle: float = 60.0`
   - `bifurcation_alignment_threshold: float = 0.0`

### 6. Testing

Test on G-U wobble pairs in 1GID:
- Should now find N2-O2 as bifurcated with N1-O2
- Verify angular separation is correct
- Check that truly oversaturated cases are still rejected

## Expected Outcome

- Recall should increase (catch bifurcated bonds DSSR finds)
- Precision should stay high (still reject true oversaturation)
- Better match with DSSR while maintaining chemical validity
