"""Selection strategies for pair finding.

This module provides strategies for selecting pairs from a pool of candidates.
The main strategy is MutualBestStrategy which implements greedy selection with
mutual best criterion - a pair is selected only if each residue's best partner
is the other residue.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Protocol, Set

from prototypes.pair_identification.finder.cache import CandidateInfo


class SelectionStrategy(Protocol):
    """Protocol for pair selection strategies."""

    def select(self, candidates: List[CandidateInfo]) -> List[CandidateInfo]:
        """Select pairs from candidates.

        Args:
            candidates: List of candidate pairs.

        Returns:
            List of selected pairs.
        """
        ...


@dataclass
class SelectionResult:
    """Result of pair selection.

    Attributes:
        selected: List of selected pairs.
        rejected: List of rejected pairs with reasons.
        used_residues: Set of residue IDs used in selected pairs.
    """

    selected: List[CandidateInfo]
    rejected: List[tuple]  # (CandidateInfo, reason)
    used_residues: Set[str]


class MutualBestStrategy:
    """Greedy selection with mutual best criterion.

    A pair is selected only if:
    1. It has valid geometry (passes validation)
    2. Each residue's best-scoring partner is the other residue
    3. Neither residue is already used in a selected pair

    This matches the C++ BasePairFinder mutual selection logic.

    Attributes:
        min_score: Minimum quality score to consider (default: 0.0).
        require_mutual: If True, require mutual best (default: True).
    """

    def __init__(self, min_score: float = 0.0, require_mutual: bool = True):
        """Initialize strategy.

        Args:
            min_score: Minimum quality score to consider.
            require_mutual: If True, require mutual best criterion.
        """
        self.min_score = min_score
        self.require_mutual = require_mutual

    def select(self, candidates: List[CandidateInfo]) -> List[CandidateInfo]:
        """Select pairs using mutual best criterion.

        Args:
            candidates: List of candidate pairs with quality scores.

        Returns:
            List of selected pairs.
        """
        result = self.select_with_details(candidates)
        return result.selected

    def select_with_details(self, candidates: List[CandidateInfo]) -> SelectionResult:
        """Select pairs with detailed results.

        Args:
            candidates: List of candidate pairs.

        Returns:
            SelectionResult with selected, rejected, and used residues.
        """
        valid = self._filter_valid(candidates)
        sorted_candidates = sorted(valid, key=lambda c: -c.quality_score)

        selected: List[CandidateInfo] = []
        rejected: List[tuple] = []
        used: Set[str] = set()

        best_for_residue = self._build_best_map(sorted_candidates)

        for candidate in sorted_candidates:
            res1, res2 = candidate.res_id1, candidate.res_id2

            if res1 in used or res2 in used:
                rejected.append((candidate, "residue_already_used"))
                continue

            if self.require_mutual:
                if not self._is_mutual_best(candidate, best_for_residue):
                    rejected.append((candidate, "not_mutual_best"))
                    continue

            selected.append(candidate)
            used.add(res1)
            used.add(res2)

            self._update_best_map(best_for_residue, res1, res2)

        return SelectionResult(selected=selected, rejected=rejected, used_residues=used)

    def _filter_valid(self, candidates: List[CandidateInfo]) -> List[CandidateInfo]:
        """Filter to valid candidates above minimum score."""
        return [
            c
            for c in candidates
            if c.validation.is_valid and c.quality_score >= self.min_score
        ]

    def _build_best_map(
        self, candidates: List[CandidateInfo]
    ) -> Dict[str, CandidateInfo]:
        """Build map of residue -> best candidate.

        Args:
            candidates: Sorted candidates (best first).

        Returns:
            Dict mapping residue ID to best candidate containing it.
        """
        best: Dict[str, CandidateInfo] = {}

        for candidate in candidates:
            for res_id in [candidate.res_id1, candidate.res_id2]:
                if res_id not in best:
                    best[res_id] = candidate

        return best

    def _is_mutual_best(
        self, candidate: CandidateInfo, best_map: Dict[str, CandidateInfo]
    ) -> bool:
        """Check if candidate is mutual best for both residues.

        Args:
            candidate: Candidate to check.
            best_map: Map of residue ID to best candidate.

        Returns:
            True if both residues' best candidate is this pair.
        """
        best1 = best_map.get(candidate.res_id1)
        best2 = best_map.get(candidate.res_id2)

        if best1 is None or best2 is None:
            return False

        return self._same_pair(candidate, best1) and self._same_pair(candidate, best2)

    @staticmethod
    def _same_pair(c1: CandidateInfo, c2: CandidateInfo) -> bool:
        """Check if two candidates represent the same pair."""
        return (c1.res_id1 == c2.res_id1 and c1.res_id2 == c2.res_id2) or (
            c1.res_id1 == c2.res_id2 and c1.res_id2 == c2.res_id1
        )

    def _update_best_map(
        self, best_map: Dict[str, CandidateInfo], res1: str, res2: str
    ) -> None:
        """Update best map after selecting a pair.

        Removes entries for used residues and updates any candidates
        that had these residues as their best option.

        Args:
            best_map: Map to update in place.
            res1: First residue ID.
            res2: Second residue ID.
        """
        best_map.pop(res1, None)
        best_map.pop(res2, None)


class GreedyBestStrategy:
    """Simple greedy selection by score without mutual criterion.

    Selects pairs in descending score order, skipping any pair where
    either residue is already used. Simpler but may miss some pairs
    that mutual best would catch.
    """

    def __init__(self, min_score: float = 0.0):
        """Initialize strategy.

        Args:
            min_score: Minimum quality score to consider.
        """
        self.min_score = min_score

    def select(self, candidates: List[CandidateInfo]) -> List[CandidateInfo]:
        """Select pairs greedily by score.

        Args:
            candidates: List of candidate pairs.

        Returns:
            List of selected pairs.
        """
        valid = [
            c
            for c in candidates
            if c.validation.is_valid and c.quality_score >= self.min_score
        ]
        sorted_candidates = sorted(valid, key=lambda c: -c.quality_score)

        selected: List[CandidateInfo] = []
        used: Set[str] = set()

        for candidate in sorted_candidates:
            if candidate.res_id1 in used or candidate.res_id2 in used:
                continue

            selected.append(candidate)
            used.add(candidate.res_id1)
            used.add(candidate.res_id2)

        return selected
