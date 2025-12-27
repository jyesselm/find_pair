"""Slot-based hydrogen bond finder with saturation tracking.

Implements greedy best-first H-bond selection respecting chemical capacity.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

import numpy as np

from .geometry import (
    ACCEPTOR_CAPACITY,
    DONOR_CAPACITY,
    HSlot,
    LPSlot,
    compute_base_normal,
    normalize,
    predict_h_slots,
    predict_lp_slots,
)

if TYPE_CHECKING:
    from ..residue import Residue


@dataclass
class HBondCandidate:
    """Potential H-bond before selection."""

    donor_res_id: str
    acceptor_res_id: str
    donor_atom: str
    acceptor_atom: str
    distance: float
    donor_pos: np.ndarray
    acceptor_pos: np.ndarray
    h_slot_idx: int = -1
    lp_slot_idx: int = -1
    alignment_score: float = 0.0

    def quality_score(self) -> float:
        """Combined quality score (shorter distance = better).

        Returns:
            Negative distance (higher is better).
        """
        return -self.distance


@dataclass
class HBond:
    """Selected H-bond with slot information."""

    donor_res_id: str
    acceptor_res_id: str
    donor_atom: str
    acceptor_atom: str
    distance: float
    h_slot_idx: int
    lp_slot_idx: int
    alignment_score: float


@dataclass
class ResidueSlotCache:
    """Cache for H and LP slots per residue."""

    res_id: str
    base_type: str
    atoms: Dict[str, np.ndarray]
    base_normal: Optional[np.ndarray] = None
    h_slots: Dict[str, List[HSlot]] = field(default_factory=dict)
    lp_slots: Dict[str, List[LPSlot]] = field(default_factory=dict)

    def get_h_slots(self, atom_name: str) -> List[HSlot]:
        """Get H slots for a donor atom, computing if needed.

        Args:
            atom_name: Name of the donor atom.

        Returns:
            List of HSlot objects.
        """
        atom_name = atom_name.strip()
        if atom_name not in self.h_slots:
            if self.base_normal is None:
                self.base_normal = compute_base_normal(self.atoms)
            self.h_slots[atom_name] = predict_h_slots(
                self.base_type, atom_name, self.atoms, self.base_normal
            )
        return self.h_slots[atom_name]

    def get_lp_slots(self, atom_name: str) -> List[LPSlot]:
        """Get LP slots for an acceptor atom, computing if needed.

        Args:
            atom_name: Name of the acceptor atom.

        Returns:
            List of LPSlot objects.
        """
        atom_name = atom_name.strip()
        if atom_name not in self.lp_slots:
            if self.base_normal is None:
                self.base_normal = compute_base_normal(self.atoms)
            self.lp_slots[atom_name] = predict_lp_slots(
                self.base_type, atom_name, self.atoms, self.base_normal
            )
        return self.lp_slots[atom_name]

    def clear_slots(self):
        """Clear all cached slots."""
        self.h_slots.clear()
        self.lp_slots.clear()


BASE_ATOMS = {"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"}


class HBondFinder:
    """Find H-bonds between residues using slot-based selection."""

    def __init__(
        self,
        max_distance: float = 4.0,
        min_alignment: float = 0.3,
        min_bifurcation_angle: float = 45.0,
        min_bifurcation_alignment: float = 0.5,
        short_distance_threshold: float = 3.2,
    ):
        """Initialize the H-bond finder.

        Args:
            max_distance: Maximum donor-acceptor distance (Angstroms).
            min_alignment: Minimum alignment score (0-2) to accept H-bond.
            min_bifurcation_angle: Min angle between bonds sharing a slot (deg).
            min_bifurcation_alignment: Stricter alignment for bifurcated bonds.
            short_distance_threshold: Below this, skip alignment check.
        """
        self.max_distance = max_distance
        self.min_alignment = min_alignment
        self.min_bifurcation_angle = min_bifurcation_angle
        self.min_bifurcation_alignment = min_bifurcation_alignment
        self.short_distance_threshold = short_distance_threshold
        self.residue_cache: Dict[str, ResidueSlotCache] = {}

    def _get_cache(self, residue: "Residue") -> ResidueSlotCache:
        """Get or create slot cache for a residue.

        Args:
            residue: Residue object.

        Returns:
            ResidueSlotCache for this residue.
        """
        if residue.res_id not in self.residue_cache:
            self.residue_cache[residue.res_id] = ResidueSlotCache(
                res_id=residue.res_id,
                base_type=residue.base_type,
                atoms=residue.atoms,
            )
        return self.residue_cache[residue.res_id]

    def _is_intra_base_pair(
        self, res1_id: str, res2_id: str, atom1: str, atom2: str
    ) -> bool:
        """Check if this is an intra-residue base-to-base contact.

        Args:
            res1_id: First residue ID.
            res2_id: Second residue ID.
            atom1: First atom name.
            atom2: Second atom name.

        Returns:
            True if same residue and both atoms are base atoms.
        """
        if res1_id != res2_id:
            return False
        return atom1.strip() in BASE_ATOMS and atom2.strip() in BASE_ATOMS

    def _add_directional_candidates(
        self,
        donor_res: "Residue",
        acceptor_res: "Residue",
        candidates: List[HBondCandidate],
    ):
        """Add candidates in one direction (donor -> acceptor).

        Args:
            donor_res: Residue with donor atoms.
            acceptor_res: Residue with acceptor atoms.
            candidates: List to append candidates to.
        """
        for donor_atom, donor_pos in donor_res.atoms.items():
            donor_key = (donor_res.base_type.upper(), donor_atom.strip())
            if donor_key not in DONOR_CAPACITY:
                continue

            for acceptor_atom, acceptor_pos in acceptor_res.atoms.items():
                acceptor_key = (acceptor_res.base_type.upper(), acceptor_atom.strip())
                if acceptor_key not in ACCEPTOR_CAPACITY:
                    continue

                if self._is_intra_base_pair(
                    donor_res.res_id, acceptor_res.res_id, donor_atom, acceptor_atom
                ):
                    continue

                dist = np.linalg.norm(acceptor_pos - donor_pos)
                if dist <= self.max_distance:
                    candidates.append(
                        HBondCandidate(
                            donor_res_id=donor_res.res_id,
                            donor_atom=donor_atom.strip(),
                            acceptor_res_id=acceptor_res.res_id,
                            acceptor_atom=acceptor_atom.strip(),
                            distance=dist,
                            donor_pos=donor_pos,
                            acceptor_pos=acceptor_pos,
                        )
                    )

    def find_candidates(self, res1: "Residue", res2: "Residue") -> List[HBondCandidate]:
        """Find all potential H-bond candidates between two residues.

        Args:
            res1: First residue.
            res2: Second residue.

        Returns:
            List of HBondCandidate objects within distance threshold.
        """
        candidates = []
        self._add_directional_candidates(res1, res2, candidates)
        self._add_directional_candidates(res2, res1, candidates)
        return candidates

    def _score_alignment(
        self,
        donor_pos: np.ndarray,
        acceptor_pos: np.ndarray,
        h_slots: List[HSlot],
        lp_slots: List[LPSlot],
    ) -> Tuple[int, int, float]:
        """Score how well an H-bond aligns with H and LP slots.

        Args:
            donor_pos: Donor atom position.
            acceptor_pos: Acceptor atom position.
            h_slots: H slots for the donor atom.
            lp_slots: LP slots for the acceptor atom.

        Returns:
            (best_h_idx, best_lp_idx, alignment_score)
        """
        if not h_slots or not lp_slots:
            return 0, 0, 0.0

        donor_to_acceptor = normalize(acceptor_pos - donor_pos)
        acceptor_to_donor = -donor_to_acceptor

        # Find best matching H slot
        best_h_score = -2.0
        best_h_idx = 0
        for i, h_slot in enumerate(h_slots):
            alignment = np.dot(h_slot.direction, donor_to_acceptor)
            if alignment > best_h_score:
                best_h_score = alignment
                best_h_idx = i

        # Find best matching LP slot
        best_lp_score = -2.0
        best_lp_idx = 0
        for i, lp_slot in enumerate(lp_slots):
            alignment = np.dot(lp_slot.direction, acceptor_to_donor)
            if alignment > best_lp_score:
                best_lp_score = alignment
                best_lp_idx = i

        total_score = best_h_score + best_lp_score
        return best_h_idx, best_lp_idx, total_score

    def _reset_slots(self, candidates: List[HBondCandidate]):
        """Reset slots for all residues involved in candidates.

        Args:
            candidates: List of candidate H-bonds.
        """
        involved_res_ids = set()
        for c in candidates:
            involved_res_ids.add(c.donor_res_id)
            involved_res_ids.add(c.acceptor_res_id)

        for res_id in involved_res_ids:
            if res_id in self.residue_cache:
                self.residue_cache[res_id].clear_slots()

    def _compute_alignments(
        self, candidates: List[HBondCandidate], residues: Dict[str, "Residue"]
    ):
        """Compute alignment scores for all candidates.

        Args:
            candidates: List of candidate H-bonds.
            residues: Dict mapping res_id to Residue objects.
        """
        for c in candidates:
            donor_res = residues.get(c.donor_res_id)
            acceptor_res = residues.get(c.acceptor_res_id)

            if not donor_res or not acceptor_res:
                continue

            donor_cache = self._get_cache(donor_res)
            acceptor_cache = self._get_cache(acceptor_res)

            h_slots = donor_cache.get_h_slots(c.donor_atom)
            lp_slots = acceptor_cache.get_lp_slots(c.acceptor_atom)

            if h_slots and lp_slots:
                h_idx, lp_idx, score = self._score_alignment(
                    c.donor_pos, c.acceptor_pos, h_slots, lp_slots
                )
                c.h_slot_idx = h_idx
                c.lp_slot_idx = lp_idx
                c.alignment_score = score

    def select_optimal(
        self, candidates: List[HBondCandidate], residues: Dict[str, "Residue"]
    ) -> List[HBond]:
        """Select optimal H-bonds using greedy slot-aware algorithm.

        Args:
            candidates: List of candidate H-bonds.
            residues: Dict mapping res_id to Residue objects.

        Returns:
            List of selected HBond objects.
        """
        if not candidates:
            return []

        self._reset_slots(candidates)
        self._compute_alignments(candidates, residues)
        candidates.sort(key=lambda c: c.quality_score(), reverse=True)

        selected = []
        for c in candidates:
            donor_res = residues.get(c.donor_res_id)
            acceptor_res = residues.get(c.acceptor_res_id)

            if not donor_res or not acceptor_res:
                continue

            donor_cache = self._get_cache(donor_res)
            acceptor_cache = self._get_cache(acceptor_res)

            h_slots = donor_cache.get_h_slots(c.donor_atom)
            lp_slots = acceptor_cache.get_lp_slots(c.acceptor_atom)

            if not self._try_select_bond(c, h_slots, lp_slots):
                continue

            selected.append(
                HBond(
                    donor_res_id=c.donor_res_id,
                    acceptor_res_id=c.acceptor_res_id,
                    donor_atom=c.donor_atom,
                    acceptor_atom=c.acceptor_atom,
                    distance=c.distance,
                    h_slot_idx=c.h_slot_idx,
                    lp_slot_idx=c.lp_slot_idx,
                    alignment_score=c.alignment_score,
                )
            )

        return selected

    def _try_select_bond(
        self,
        candidate: HBondCandidate,
        h_slots: List[HSlot],
        lp_slots: List[LPSlot],
    ) -> bool:
        """Try to select a bond with available slots.

        Args:
            candidate: The candidate H-bond.
            h_slots: H slots for the donor atom.
            lp_slots: LP slots for the acceptor atom.

        Returns:
            True if bond was selected, False otherwise.
        """
        if not h_slots or not lp_slots:
            return False

        if candidate.h_slot_idx >= len(h_slots):
            return False
        if candidate.lp_slot_idx >= len(lp_slots):
            return False

        h_slot = h_slots[candidate.h_slot_idx]
        lp_slot = lp_slots[candidate.lp_slot_idx]

        donor_to_acceptor = normalize(candidate.acceptor_pos - candidate.donor_pos)
        acceptor_to_donor = -donor_to_acceptor

        # Check if we can use these slots
        h_can_use = h_slot.can_add_bond(donor_to_acceptor, self.min_bifurcation_angle)
        lp_can_use = lp_slot.can_add_bond(acceptor_to_donor, self.min_bifurcation_angle)

        is_bifurcated = (
            len(h_slot.bond_directions) > 0 or len(lp_slot.bond_directions) > 0
        )

        if not (h_can_use and lp_can_use):
            # Try alternative slots
            if not self._find_alternative_slots(
                candidate, h_slots, lp_slots, donor_to_acceptor, acceptor_to_donor
            ):
                return False
            # Update h_slot, lp_slot, and is_bifurcated
            h_slot = h_slots[candidate.h_slot_idx]
            lp_slot = lp_slots[candidate.lp_slot_idx]
            is_bifurcated = (
                len(h_slot.bond_directions) > 0 or len(lp_slot.bond_directions) > 0
            )

        # Check minimum alignment
        if candidate.distance >= self.short_distance_threshold:
            min_align = (
                self.min_bifurcation_alignment if is_bifurcated else self.min_alignment
            )
            if candidate.alignment_score < min_align:
                return False

        # Accept this H-bond
        h_slot.add_bond(donor_to_acceptor)
        lp_slot.add_bond(acceptor_to_donor)
        return True

    def _find_alternative_slots(
        self,
        candidate: HBondCandidate,
        h_slots: List[HSlot],
        lp_slots: List[LPSlot],
        donor_to_acceptor: np.ndarray,
        acceptor_to_donor: np.ndarray,
    ) -> bool:
        """Try to find alternative slots for a candidate.

        Args:
            candidate: The candidate H-bond.
            h_slots: H slots for the donor atom.
            lp_slots: LP slots for the acceptor atom.
            donor_to_acceptor: Direction vector from donor to acceptor.
            acceptor_to_donor: Direction vector from acceptor to donor.

        Returns:
            True if alternative slots found, False otherwise.
        """
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

                # Check if bifurcated
                alt_is_bifurcated = (
                    len(hs.bond_directions) > 0 or len(ls.bond_directions) > 0
                )
                min_align_required = (
                    self.min_bifurcation_alignment
                    if alt_is_bifurcated
                    else self.min_alignment
                )

                if alt_score >= min_align_required:
                    candidate.h_slot_idx = hi
                    candidate.lp_slot_idx = li
                    candidate.alignment_score = alt_score
                    return True

        return False

    def find_between(self, res1: "Residue", res2: "Residue") -> List[HBond]:
        """Find H-bonds between two residues.

        Args:
            res1: First residue.
            res2: Second residue.

        Returns:
            List of selected HBond objects.
        """
        candidates = self.find_candidates(res1, res2)
        residues = {res1.res_id: res1, res2.res_id: res2}
        return self.select_optimal(candidates, residues)
