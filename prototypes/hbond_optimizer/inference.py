"""
H-bond donor/acceptor inference module.

Infers H-bond properties from molecular geometry and validates against
reference data from DSSR/literature.
"""

import json
import os
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass

try:
    from .hybridization import (
        analyze_residue_hbond_properties,
        classify_nitrogen,
        classify_oxygen,
        detect_atom_hybridization,
        get_element,
        Hybridization
    )
except ImportError:
    from hybridization import (
        analyze_residue_hbond_properties,
        classify_nitrogen,
        classify_oxygen,
        detect_atom_hybridization,
        get_element,
        Hybridization
    )


# Path to reference data
DATA_DIR = Path(__file__).parent.parent.parent / "data"
DONORS_FILE = DATA_DIR / "hbond_donors.json"
ACCEPTORS_FILE = DATA_DIR / "hbond_acceptors.json"


def load_reference_data() -> Tuple[Dict, Dict]:
    """Load reference donor and acceptor data from JSON files."""
    donors = {}
    acceptors = {}

    if DONORS_FILE.exists():
        with open(DONORS_FILE) as f:
            donors = json.load(f)

    if ACCEPTORS_FILE.exists():
        with open(ACCEPTORS_FILE) as f:
            acceptors = json.load(f)

    return donors, acceptors


@dataclass
class ValidationResult:
    """Result of validating inferred vs reference H-bond properties."""
    residue_name: str
    inferred_donors: Dict[str, int]
    reference_donors: Dict[str, int]
    inferred_acceptors: Dict[str, int]
    reference_acceptors: Dict[str, int]
    donor_matches: int = 0
    donor_mismatches: int = 0
    donor_missing: int = 0  # In reference but not inferred
    donor_extra: int = 0    # Inferred but not in reference
    acceptor_matches: int = 0
    acceptor_mismatches: int = 0
    acceptor_missing: int = 0
    acceptor_extra: int = 0

    def calculate_stats(self):
        """Calculate match statistics."""
        # Donors
        all_donor_atoms = set(self.inferred_donors.keys()) | set(self.reference_donors.keys())
        for atom in all_donor_atoms:
            inferred = self.inferred_donors.get(atom, 0)
            reference = self.reference_donors.get(atom, 0)

            if reference is None:
                # Reference says null (unknown residue)
                continue

            if inferred == reference:
                if inferred > 0:
                    self.donor_matches += 1
            elif inferred == 0 and reference > 0:
                self.donor_missing += 1
            elif inferred > 0 and reference == 0:
                self.donor_extra += 1
            else:
                self.donor_mismatches += 1

        # Acceptors
        all_acceptor_atoms = set(self.inferred_acceptors.keys()) | set(self.reference_acceptors.keys())
        for atom in all_acceptor_atoms:
            inferred = self.inferred_acceptors.get(atom, 0)
            reference = self.reference_acceptors.get(atom, 0)

            if reference is None:
                continue

            if inferred == reference:
                if inferred > 0:
                    self.acceptor_matches += 1
            elif inferred == 0 and reference > 0:
                self.acceptor_missing += 1
            elif inferred > 0 and reference == 0:
                self.acceptor_extra += 1
            else:
                self.acceptor_mismatches += 1

    @property
    def donor_accuracy(self) -> float:
        """Accuracy for donor prediction."""
        total = self.donor_matches + self.donor_mismatches + self.donor_missing
        if total == 0:
            return 1.0
        return self.donor_matches / total

    @property
    def acceptor_accuracy(self) -> float:
        """Accuracy for acceptor prediction."""
        total = self.acceptor_matches + self.acceptor_mismatches + self.acceptor_missing
        if total == 0:
            return 1.0
        return self.acceptor_matches / total


class HBondInferenceEngine:
    """
    Engine for inferring H-bond properties from molecular geometry.

    Can use reference data as fallback for known residues.
    """

    def __init__(self, use_reference_fallback: bool = True):
        """
        Args:
            use_reference_fallback: If True, use reference data for known residues
        """
        self.use_reference_fallback = use_reference_fallback
        self.ref_donors, self.ref_acceptors = load_reference_data()

        # Atoms that should NOT be considered donors based on chemistry
        # (even if geometry suggests otherwise)
        self.excluded_donors = {
            # Phosphate oxygens - rarely donors
            'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P',
            # Carbonyl oxygens
            'O2', 'O4', 'O6',  # Base carbonyls
        }

        # Atoms that are known special cases
        self.known_donors = {
            # Amino groups
            'N6': 2,  # Adenine
            'N4': 2,  # Cytosine
            'N2': 2,  # Guanine
            # Imino
            'N1': 1,  # Guanine
            'N3': 1,  # Uracil/Thymine
            # Backbone
            'N': 2,   # Amino acid backbone (N-H + optionally N-H in some cases)
        }

        # Known acceptors with capacities
        self.known_acceptors = {
            # Carbonyl oxygens - 2 lone pairs each
            'O': 2,    # Backbone carbonyl
            'O2': 2,
            'O4': 2,
            'O6': 2,
            # Ring nitrogens - 1 lone pair
            'N1': 1,   # When not H-bonding (A, C)
            'N3': 1,
            'N7': 1,
            # Ribose oxygens
            "O2'": 2,
            "O3'": 2,
            "O4'": 1,  # Ring oxygen - less accessible
            "O5'": 2,
        }

    def infer_from_geometry(
        self,
        residue_name: str,
        atoms: Dict[str, np.ndarray]
    ) -> Tuple[Dict[str, int], Dict[str, int]]:
        """
        Infer donor and acceptor properties from atomic geometry.

        Args:
            residue_name: 3-letter residue code
            atoms: Dict mapping atom names to positions

        Returns:
            (donors, acceptors) - dicts mapping atom names to capacities
        """
        donors = {}
        acceptors = {}

        for atom_name, atom_pos in atoms.items():
            element = get_element(atom_name)

            # Skip non-heteroatoms
            if element not in ('N', 'O', 'S'):
                continue

            # Use specialized classifiers
            if element == 'N':
                d_cap, a_cap = classify_nitrogen(atom_name, atoms)
            elif element == 'O':
                d_cap, a_cap = classify_oxygen(atom_name, atoms)
            else:  # S
                # Simple heuristic for sulfur
                hybridization, num_bonds, _ = detect_atom_hybridization(atom_name, atoms)
                if hybridization == Hybridization.SP3 and num_bonds == 1:
                    d_cap, a_cap = 1, 2  # Thiol
                else:
                    d_cap, a_cap = 0, 2

            # Apply exclusions
            stripped_name = atom_name.strip()
            if stripped_name in self.excluded_donors:
                d_cap = 0

            if d_cap > 0:
                donors[stripped_name] = d_cap
            if a_cap > 0:
                acceptors[stripped_name] = a_cap

        return donors, acceptors

    def get_properties(
        self,
        residue_name: str,
        atoms: Optional[Dict[str, np.ndarray]] = None
    ) -> Tuple[Dict[str, int], Dict[str, int]]:
        """
        Get H-bond donor/acceptor properties for a residue.

        Uses reference data if available and fallback enabled,
        otherwise infers from geometry.

        Args:
            residue_name: 3-letter residue code
            atoms: Dict mapping atom names to positions (required for inference)

        Returns:
            (donors, acceptors) - dicts mapping atom names to capacities
        """
        residue_name = residue_name.upper().strip()

        # Try reference data first
        if self.use_reference_fallback:
            ref_donors = self.ref_donors.get(residue_name)
            ref_acceptors = self.ref_acceptors.get(residue_name)

            if ref_donors is not None and ref_acceptors is not None:
                return ref_donors, ref_acceptors

        # Fall back to geometry-based inference
        if atoms is not None:
            return self.infer_from_geometry(residue_name, atoms)

        return {}, {}

    def validate_inference(
        self,
        residue_name: str,
        atoms: Dict[str, np.ndarray]
    ) -> ValidationResult:
        """
        Validate inferred properties against reference data.

        Args:
            residue_name: 3-letter residue code
            atoms: Dict mapping atom names to positions

        Returns:
            ValidationResult with comparison details
        """
        residue_name = residue_name.upper().strip()

        # Get inferred values
        inferred_donors, inferred_acceptors = self.infer_from_geometry(residue_name, atoms)

        # Get reference values
        ref_donors = self.ref_donors.get(residue_name, {})
        ref_acceptors = self.ref_acceptors.get(residue_name, {})

        # Handle null (unknown residue)
        if ref_donors is None:
            ref_donors = {}
        if ref_acceptors is None:
            ref_acceptors = {}

        result = ValidationResult(
            residue_name=residue_name,
            inferred_donors=inferred_donors,
            reference_donors=ref_donors,
            inferred_acceptors=inferred_acceptors,
            reference_acceptors=ref_acceptors
        )
        result.calculate_stats()

        return result


def validate_amino_acids(engine: HBondInferenceEngine) -> Dict[str, ValidationResult]:
    """
    Validate inference for all 20 standard amino acids.

    Uses idealized coordinates if available, otherwise skips geometry check.
    """
    AMINO_ACIDS = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL'
    ]

    results = {}
    for aa in AMINO_ACIDS:
        # Check if in reference
        ref_donors = engine.ref_donors.get(aa)
        ref_acceptors = engine.ref_acceptors.get(aa)

        if ref_donors is None or ref_acceptors is None:
            print(f"Warning: {aa} not in reference data")
            continue

        # For now, just report reference data
        # (geometry validation requires actual coordinates)
        results[aa] = {
            'donors': ref_donors,
            'acceptors': ref_acceptors
        }

    return results


def print_amino_acid_summary():
    """Print a summary of amino acid H-bond properties from reference data."""
    engine = HBondInferenceEngine()

    print("Amino Acid H-Bond Properties (from reference data)")
    print("=" * 70)
    print(f"{'Residue':<8} {'Donors':<35} {'Acceptors':<35}")
    print("-" * 70)

    AMINO_ACIDS = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL'
    ]

    for aa in AMINO_ACIDS:
        donors = engine.ref_donors.get(aa, {})
        acceptors = engine.ref_acceptors.get(aa, {})

        if donors is None:
            donors = {}
        if acceptors is None:
            acceptors = {}

        # Format as compact strings
        donor_str = ", ".join(f"{k}:{v}" for k, v in sorted(donors.items()))
        acceptor_str = ", ".join(f"{k}:{v}" for k, v in sorted(acceptors.items()))

        # Truncate if too long
        if len(donor_str) > 33:
            donor_str = donor_str[:30] + "..."
        if len(acceptor_str) > 33:
            acceptor_str = acceptor_str[:30] + "..."

        print(f"{aa:<8} {donor_str:<35} {acceptor_str:<35}")


def analyze_questionable_donors(ref_donors: Dict) -> List[Tuple[str, str, int]]:
    """
    Find questionable donor entries in reference data.

    These are atoms that are rarely donors in practice:
    - Phosphate oxygens (O1P, O2P, OP1, OP2)
    - Carbonyl oxygens in typical contexts

    Returns:
        List of (residue, atom, capacity) tuples
    """
    questionable = []

    # Atoms that are rarely/never donors
    rarely_donors = {
        'O1P', 'O2P', 'OP1', 'OP2', 'O3P', 'OP3',  # Phosphate
        # Some oxygen atoms in specific contexts
    }

    for residue, atoms in ref_donors.items():
        if atoms is None:
            continue
        for atom, capacity in atoms.items():
            if atom in rarely_donors and capacity > 0:
                questionable.append((residue, atom, capacity))

    return questionable


if __name__ == "__main__":
    print("H-Bond Inference Module")
    print("=" * 50)

    # Print amino acid summary
    print("\n")
    print_amino_acid_summary()

    # Check for questionable donors in reference data
    print("\n\nQuestionable Donors in Reference Data")
    print("=" * 50)

    engine = HBondInferenceEngine()
    questionable = analyze_questionable_donors(engine.ref_donors)

    if questionable:
        print(f"Found {len(questionable)} questionable donor entries:")
        for res, atom, cap in questionable[:20]:  # Show first 20
            print(f"  {res}: {atom} = {cap}")
        if len(questionable) > 20:
            print(f"  ... and {len(questionable) - 20} more")
    else:
        print("No questionable entries found.")
