"""H-bond analyzer for cWW miss annotation.

This module analyzes hydrogen bonding patterns between Watson-Crick base pairs,
comparing detected H-bonds against expected canonical patterns and DSSR reference.
"""

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

from .diagnostics import ExpectedHBond, HBondDiagnostics
from .loaders import SlotHBond


# Expected H-bond patterns for canonical Watson-Crick base pairs
CWW_EXPECTED_PATTERNS: Dict[str, List[Tuple[str, str]]] = {
    "GC": [("N1", "N3"), ("N2", "O2"), ("O6", "N4")],
    "CG": [("N4", "O6"), ("N3", "N1"), ("O2", "N2")],
    "AU": [("N6", "O4"), ("N1", "N3")],
    "UA": [("O4", "N6"), ("N3", "N1")],
}


def parse_dssr_hbonds(hbonds_desc: str) -> List[Tuple[str, str, float]]:
    """Parse DSSR hbonds_desc string.

    Args:
        hbonds_desc: DSSR format like "O6(carbonyl)-N4(amino)[2.83],N1(imino)-N3[2.88]"

    Returns:
        List of (donor_atom, acceptor_atom, distance) tuples.
        Note: DSSR format doesn't explicitly label donor vs acceptor, so we
        return both directions and let the analyzer handle matching.

    Examples:
        >>> parse_dssr_hbonds("O6(carbonyl)-N4(amino)[2.83],N1(imino)-N3[2.88]")
        [('O6', 'N4', 2.83), ('N1', 'N3', 2.88)]
        >>> parse_dssr_hbonds("N6(amino)-O4(carbonyl)[2.95]")
        [('N6', 'O4', 2.95)]
        >>> parse_dssr_hbonds("")
        []
    """
    if not hbonds_desc:
        return []

    hbonds = []

    # Pattern: atom1[(type)]-atom2[(type)][distance]
    # Example: "O6(carbonyl)-N4(amino)[2.83]"
    pattern = r'([A-Z]\d+(?:\*)?)\([^)]+\)-([A-Z]\d+(?:\*)?)\([^)]+\)\[(\d+\.\d+)\]'

    for match in re.finditer(pattern, hbonds_desc):
        atom1 = match.group(1).rstrip('*')  # Remove asterisk if present
        atom2 = match.group(2).rstrip('*')
        distance = float(match.group(3))
        hbonds.append((atom1, atom2, distance))

    # Also handle simpler format without parentheses
    # Example: "N1-N3[2.88]"
    simple_pattern = r'([A-Z]\d+(?:\*)?)-([A-Z]\d+(?:\*)?)\[(\d+\.\d+)\]'

    for match in re.finditer(simple_pattern, hbonds_desc):
        atom1 = match.group(1).rstrip('*')
        atom2 = match.group(2).rstrip('*')
        distance = float(match.group(3))
        # Only add if not already added by the detailed pattern
        if (atom1, atom2, distance) not in hbonds:
            hbonds.append((atom1, atom2, distance))

    return hbonds


class HBondAnalyzer:
    """Analyzes H-bond patterns for Watson-Crick base pairs.

    Compares detected H-bonds against expected canonical patterns,
    identifies missing/extra/wrong H-bonds, and provides diagnostics
    for classification failures.
    """

    def __init__(self):
        """Initialize the H-bond analyzer with expected patterns."""
        self.expected_patterns = CWW_EXPECTED_PATTERNS

    def analyze(
        self,
        sequence: str,
        found_hbonds: List[SlotHBond],
        dssr_hbonds_desc: str,
    ) -> HBondDiagnostics:
        """Analyze H-bond pattern for a base pair.

        Args:
            sequence: Two-letter sequence like "GC", "CG", "AU", "UA"
            found_hbonds: List of detected SlotHBond objects
            dssr_hbonds_desc: DSSR hbonds description string for reference

        Returns:
            HBondDiagnostics with detailed analysis of H-bond patterns
        """
        # Get expected pattern for this sequence
        expected = self.expected_patterns.get(sequence, [])
        expected_hbonds = [
            ExpectedHBond(donor_atom=d, acceptor_atom=a, source="canonical")
            for d, a in expected
        ]

        # Parse DSSR reference
        dssr_hbonds = parse_dssr_hbonds(dssr_hbonds_desc)

        # Convert found_hbonds to dict format for diagnostics
        found_hbond_dicts = [
            {
                "donor_atom": hb.donor_atom,
                "acceptor_atom": hb.acceptor_atom,
                "distance": hb.distance,
                "context": hb.context,
                "h_slot": hb.h_slot,
                "lp_slot": hb.lp_slot,
                "alignment": hb.alignment,
            }
            for hb in found_hbonds
        ]

        # Initialize diagnostics
        diagnostics = HBondDiagnostics(
            expected_hbonds=expected_hbonds,
            found_hbonds=found_hbond_dicts,
        )

        # Find matched, missing, and extra H-bonds
        matched = []
        missing = []
        extra = []

        for exp_hb in expected_hbonds:
            match = self._match_hbond(found_hbonds, exp_hb.donor_atom, exp_hb.acceptor_atom)
            if match:
                matched.append(match)
            else:
                missing.append(exp_hb)

        # Find extra H-bonds (not in expected)
        for hb in found_hbonds:
            if not self._is_expected(hb, expected):
                extra.append({
                    "donor_atom": hb.donor_atom,
                    "acceptor_atom": hb.acceptor_atom,
                    "distance": hb.distance,
                    "context": hb.context,
                })

        diagnostics.missing_hbonds = missing
        diagnostics.extra_hbonds = extra

        # Identify wrong atoms (donor matches but acceptor differs)
        wrong_atoms = {}
        for hb in found_hbonds:
            wrong_desc = self._find_wrong_atom(hb, expected)
            if wrong_desc:
                key = f"{hb.donor_atom}->{hb.acceptor_atom}"
                wrong_atoms[key] = wrong_desc

        diagnostics.wrong_atoms = wrong_atoms

        # Check distance issues (outside normal range 2.5-3.5A)
        distance_issues = []
        for hb in found_hbonds:
            if hb.distance < 2.5:
                distance_issues.append((
                    f"{hb.donor_atom}-{hb.acceptor_atom}",
                    hb.distance,
                    "too_short"
                ))
            elif hb.distance > 3.5:
                distance_issues.append((
                    f"{hb.donor_atom}-{hb.acceptor_atom}",
                    hb.distance,
                    "too_long"
                ))

        diagnostics.distance_issues = distance_issues

        # Check for overloaded acceptors (>2 H-bonds to same acceptor)
        acceptor_counts = Counter(hb.acceptor_atom for hb in found_hbonds)
        overloaded = [atom for atom, count in acceptor_counts.items() if count > 2]
        diagnostics.overloaded_acceptors = overloaded

        return diagnostics

    def _match_hbond(
        self,
        found_hbonds: List[SlotHBond],
        donor: str,
        acceptor: str,
    ) -> Optional[SlotHBond]:
        """Find a matching H-bond in the found list.

        Args:
            found_hbonds: List of detected H-bonds
            donor: Expected donor atom name
            acceptor: Expected acceptor atom name

        Returns:
            Matching SlotHBond if found, None otherwise.
            Checks both directions (donor->acceptor and acceptor->donor).
        """
        for hb in found_hbonds:
            # Direct match
            if hb.donor_atom == donor and hb.acceptor_atom == acceptor:
                return hb
            # Reverse match (some H-bonds may have swapped donor/acceptor labels)
            if hb.donor_atom == acceptor and hb.acceptor_atom == donor:
                return hb
        return None

    def _is_expected(
        self,
        hb: SlotHBond,
        expected: List[Tuple[str, str]],
    ) -> bool:
        """Check if an H-bond is in the expected pattern.

        Args:
            hb: SlotHBond to check
            expected: List of (donor, acceptor) tuples

        Returns:
            True if H-bond matches any expected pattern (either direction)
        """
        for donor, acceptor in expected:
            if hb.donor_atom == donor and hb.acceptor_atom == acceptor:
                return True
            if hb.donor_atom == acceptor and hb.acceptor_atom == donor:
                return True
        return False

    def _find_wrong_atom(
        self,
        hb: SlotHBond,
        expected: List[Tuple[str, str]],
    ) -> Optional[str]:
        """Check if donor matches but acceptor differs from expected.

        Args:
            hb: SlotHBond to check
            expected: List of (donor, acceptor) tuples

        Returns:
            Description string if wrong acceptor, None otherwise.
            Example: "Expected N1->N3, got N1->O2"
        """
        for donor, acceptor in expected:
            # Donor matches but acceptor differs
            if hb.donor_atom == donor and hb.acceptor_atom != acceptor:
                return f"Expected {donor}->{acceptor}, got {donor}->{hb.acceptor_atom}"
            # Reverse: acceptor matches (as donor) but donor differs
            if hb.donor_atom == acceptor and hb.acceptor_atom != donor:
                return f"Expected {acceptor}->{donor}, got {acceptor}->{hb.acceptor_atom}"
        return None
