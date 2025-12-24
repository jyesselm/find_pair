"""
H-Bond Optimizer with slot tracking.

Implements greedy best-first selection of H-bonds while respecting
chemical capacity limits (H slots for donors, LP slots for acceptors).
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict

try:
    from .geometry import (
        normalize, angle_between, compute_base_normal,
        predict_h_slots, predict_lp_slots, score_hbond_alignment,
        HSlot, LPSlot, DONOR_CAPACITY, ACCEPTOR_CAPACITY
    )
except ImportError:
    from geometry import (
        normalize, angle_between, compute_base_normal,
        predict_h_slots, predict_lp_slots, score_hbond_alignment,
        HSlot, LPSlot, DONOR_CAPACITY, ACCEPTOR_CAPACITY
    )


@dataclass
class Residue:
    """A nucleotide residue with its atoms."""
    res_id: str           # e.g., "A.G1" (DSSR format)
    base_type: str        # A, G, C, U, T, P, I (parent base type)
    atoms: Dict[str, np.ndarray] = field(default_factory=dict)
    base_normal: Optional[np.ndarray] = None
    residue_code: str = ""  # 3-letter code, e.g., "5MC", "H2U", "PSU"

    # Cached slots (computed on first access)
    _h_slots: Dict[str, List[HSlot]] = field(default_factory=dict)
    _lp_slots: Dict[str, List[LPSlot]] = field(default_factory=dict)

    def get_h_slots(self, atom_name: str) -> List[HSlot]:
        """Get H slots for a donor atom, computing if needed."""
        atom_name = atom_name.strip()
        if atom_name not in self._h_slots:
            if self.base_normal is None:
                self.base_normal = compute_base_normal(self.atoms)
            self._h_slots[atom_name] = predict_h_slots(
                self.base_type, atom_name, self.atoms, self.base_normal
            )
        return self._h_slots[atom_name]

    def get_lp_slots(self, atom_name: str) -> List[LPSlot]:
        """Get LP slots for an acceptor atom, computing if needed."""
        atom_name = atom_name.strip()
        if atom_name not in self._lp_slots:
            if self.base_normal is None:
                self.base_normal = compute_base_normal(self.atoms)
            self._lp_slots[atom_name] = predict_lp_slots(
                self.base_type, atom_name, self.atoms, self.base_normal
            )
        return self._lp_slots[atom_name]


@dataclass
class HBondCandidate:
    """A potential hydrogen bond between two residues."""
    donor_res_id: str
    donor_atom: str
    acceptor_res_id: str
    acceptor_atom: str
    distance: float
    donor_pos: np.ndarray
    acceptor_pos: np.ndarray

    # Alignment info (filled in during selection)
    h_slot_idx: int = -1
    lp_slot_idx: int = -1
    alignment_score: float = 0.0

    def quality_score(self) -> float:
        """Combined quality score (lower distance + higher alignment = better)."""
        # Prioritize distance - shorter is better
        # Alignment is used as a filter, not for scoring
        return -self.distance


@dataclass
class HBond:
    """A selected hydrogen bond."""
    donor_res_id: str
    donor_atom: str
    acceptor_res_id: str
    acceptor_atom: str
    distance: float
    h_slot_idx: int
    lp_slot_idx: int
    alignment_score: float


class HBondOptimizer:
    """
    Greedy H-bond optimizer with slot tracking.

    Selects optimal H-bonds by:
    1. Computing all candidate H-bonds within distance threshold
    2. Predicting H and LP positions from geometry
    3. Scoring alignment of each candidate with available slots
    4. Greedily selecting best candidates while respecting capacity
    """

    def __init__(self,
                 max_distance: float = 4.0,
                 min_alignment: float = 0.3,
                 min_bifurcation_angle: float = 45.0,
                 min_bifurcation_alignment: float = 0.5):
        """
        Args:
            max_distance: Maximum donor-acceptor distance in Angstroms
            min_alignment: Minimum alignment score (0-2) to accept H-bond
            min_bifurcation_angle: Minimum angle between bonds sharing a slot (degrees)
            min_bifurcation_alignment: Stricter alignment for bifurcated bonds (0-2)
        """
        self.max_distance = max_distance
        self.min_alignment = min_alignment
        self.min_bifurcation_angle = min_bifurcation_angle
        self.min_bifurcation_alignment = min_bifurcation_alignment

        self.residues: Dict[str, Residue] = {}
        self.candidates: List[HBondCandidate] = []
        self.selected: List[HBond] = []

    def add_residue(self, residue: Residue):
        """Add a residue to the optimizer."""
        self.residues[residue.res_id] = residue

    # Base atoms (cannot form intra-residue H-bonds with each other)
    BASE_ATOMS = {'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6'}

    def _is_intra_base_pair(self, res1_id: str, res2_id: str,
                            atom1: str, atom2: str) -> bool:
        """Check if this is an intra-residue base-to-base contact (covalent, not H-bond)."""
        if res1_id != res2_id:
            return False
        # Only filter out base-to-base contacts within same residue
        return atom1.strip() in self.BASE_ATOMS and atom2.strip() in self.BASE_ATOMS

    def find_candidates(self, res1_id: str, res2_id: str) -> List[HBondCandidate]:
        """
        Find all potential H-bond candidates between two residues.

        Checks all donor-acceptor pairs within distance threshold.
        """
        res1 = self.residues.get(res1_id)
        res2 = self.residues.get(res2_id)

        if not res1 or not res2:
            return []

        candidates = []

        # Check res1 donors -> res2 acceptors
        for donor_atom, donor_pos in res1.atoms.items():
            donor_key = (res1.base_type.upper(), donor_atom.strip())
            if donor_key not in DONOR_CAPACITY:
                continue

            for acceptor_atom, acceptor_pos in res2.atoms.items():
                acceptor_key = (res2.base_type.upper(), acceptor_atom.strip())
                if acceptor_key not in ACCEPTOR_CAPACITY:
                    continue

                # Skip intra-residue base-to-base contacts (covalent bonds)
                if self._is_intra_base_pair(res1_id, res2_id, donor_atom, acceptor_atom):
                    continue

                dist = np.linalg.norm(acceptor_pos - donor_pos)
                if dist <= self.max_distance:
                    candidates.append(HBondCandidate(
                        donor_res_id=res1_id,
                        donor_atom=donor_atom.strip(),
                        acceptor_res_id=res2_id,
                        acceptor_atom=acceptor_atom.strip(),
                        distance=dist,
                        donor_pos=donor_pos,
                        acceptor_pos=acceptor_pos
                    ))

        # Check res2 donors -> res1 acceptors
        for donor_atom, donor_pos in res2.atoms.items():
            donor_key = (res2.base_type.upper(), donor_atom.strip())
            if donor_key not in DONOR_CAPACITY:
                continue

            for acceptor_atom, acceptor_pos in res1.atoms.items():
                acceptor_key = (res1.base_type.upper(), acceptor_atom.strip())
                if acceptor_key not in ACCEPTOR_CAPACITY:
                    continue

                # Skip intra-residue base-to-base contacts (covalent bonds)
                if self._is_intra_base_pair(res2_id, res1_id, donor_atom, acceptor_atom):
                    continue

                dist = np.linalg.norm(acceptor_pos - donor_pos)
                if dist <= self.max_distance:
                    candidates.append(HBondCandidate(
                        donor_res_id=res2_id,
                        donor_atom=donor_atom.strip(),
                        acceptor_res_id=res1_id,
                        acceptor_atom=acceptor_atom.strip(),
                        distance=dist,
                        donor_pos=donor_pos,
                        acceptor_pos=acceptor_pos
                    ))

        return candidates

    def select_optimal(self, candidates: List[HBondCandidate]) -> List[HBond]:
        """
        Select optimal H-bonds using greedy slot-aware algorithm with bifurcation support.

        1. Reset all slots to unused
        2. Compute alignment scores for each candidate
        3. Sort by quality (distance)
        4. Greedily select best candidates, allowing bifurcation when angularly separated
        """
        if not candidates:
            return []

        # Reset slots for involved residues
        involved_res_ids = set()
        for c in candidates:
            involved_res_ids.add(c.donor_res_id)
            involved_res_ids.add(c.acceptor_res_id)

        for res_id in involved_res_ids:
            res = self.residues.get(res_id)
            if res:
                res._h_slots.clear()
                res._lp_slots.clear()

        # Compute alignment scores
        for c in candidates:
            donor_res = self.residues.get(c.donor_res_id)
            acceptor_res = self.residues.get(c.acceptor_res_id)

            if not donor_res or not acceptor_res:
                continue

            h_slots = donor_res.get_h_slots(c.donor_atom)
            lp_slots = acceptor_res.get_lp_slots(c.acceptor_atom)

            if h_slots and lp_slots:
                h_idx, lp_idx, score = score_hbond_alignment(
                    c.donor_pos, c.acceptor_pos, h_slots, lp_slots
                )
                c.h_slot_idx = h_idx
                c.lp_slot_idx = lp_idx
                c.alignment_score = score

        # Sort by quality (higher is better - shorter distance)
        candidates.sort(key=lambda c: c.quality_score(), reverse=True)

        selected = []
        for c in candidates:
            donor_res = self.residues.get(c.donor_res_id)
            acceptor_res = self.residues.get(c.acceptor_res_id)

            if not donor_res or not acceptor_res:
                continue

            h_slots = donor_res.get_h_slots(c.donor_atom)
            lp_slots = acceptor_res.get_lp_slots(c.acceptor_atom)

            if not h_slots or not lp_slots:
                continue

            # Check if best slots are still available
            if c.h_slot_idx >= len(h_slots) or c.lp_slot_idx >= len(lp_slots):
                continue

            h_slot = h_slots[c.h_slot_idx]
            lp_slot = lp_slots[c.lp_slot_idx]

            # Compute directions for bifurcation check
            donor_to_acceptor = normalize(c.acceptor_pos - c.donor_pos)
            acceptor_to_donor = -donor_to_acceptor

            # Check if we can use these slots (available or bifurcatable)
            h_can_use = h_slot.can_add_bond(donor_to_acceptor, self.min_bifurcation_angle)
            lp_can_use = lp_slot.can_add_bond(acceptor_to_donor, self.min_bifurcation_angle)

            # Determine if this is a bifurcated bond (slot already has bonds)
            is_bifurcated = len(h_slot.bond_directions) > 0 or len(lp_slot.bond_directions) > 0

            if not (h_can_use and lp_can_use):
                # Try to find alternative slots
                found_alt = False
                for hi, hs in enumerate(h_slots):
                    if not hs.can_add_bond(donor_to_acceptor, self.min_bifurcation_angle):
                        continue
                    for li, ls in enumerate(lp_slots):
                        if not ls.can_add_bond(acceptor_to_donor, self.min_bifurcation_angle):
                            continue
                        # Recompute alignment for this slot pair
                        h_align = np.dot(hs.direction, donor_to_acceptor)
                        lp_align = np.dot(ls.direction, acceptor_to_donor)
                        alt_score = h_align + lp_align
                        # Check if this would be bifurcated
                        alt_is_bifurcated = len(hs.bond_directions) > 0 or len(ls.bond_directions) > 0
                        min_align_required = self.min_bifurcation_alignment if alt_is_bifurcated else self.min_alignment
                        if alt_score >= min_align_required:
                            c.h_slot_idx = hi
                            c.lp_slot_idx = li
                            c.alignment_score = alt_score
                            h_slot = hs
                            lp_slot = ls
                            is_bifurcated = alt_is_bifurcated
                            found_alt = True
                            break
                    if found_alt:
                        break

                if not found_alt:
                    continue

            # Check minimum alignment (stricter for bifurcated bonds)
            min_align_required = self.min_bifurcation_alignment if is_bifurcated else self.min_alignment
            if c.alignment_score < min_align_required:
                continue

            # Accept this H-bond - record the bond direction for bifurcation tracking
            h_slot.add_bond(donor_to_acceptor)
            lp_slot.add_bond(acceptor_to_donor)

            selected.append(HBond(
                donor_res_id=c.donor_res_id,
                donor_atom=c.donor_atom,
                acceptor_res_id=c.acceptor_res_id,
                acceptor_atom=c.acceptor_atom,
                distance=c.distance,
                h_slot_idx=c.h_slot_idx,
                lp_slot_idx=c.lp_slot_idx,
                alignment_score=c.alignment_score
            ))

        return selected

    def optimize_pair(self, res1_id: str, res2_id: str) -> List[HBond]:
        """Find and select optimal H-bonds between two residues."""
        candidates = self.find_candidates(res1_id, res2_id)
        return self.select_optimal(candidates)


def format_hbond(hb: HBond) -> str:
    """Format an H-bond for display."""
    return (f"{hb.donor_atom:>4s} -> {hb.acceptor_atom:<4s}  "
            f"dist={hb.distance:.2f}  align={hb.alignment_score:.2f}  "
            f"slots=H{hb.h_slot_idx}/LP{hb.lp_slot_idx}")
