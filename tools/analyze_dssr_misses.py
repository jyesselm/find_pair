#!/usr/bin/env python3
"""
Analyze why we miss H-bonds that DSSR detects using slot-aware analysis.

For each DSSR H-bond we don't detect, explains:
- DONOR_OVERSATURATED: All H slots on donor atom already used
- ACCEPTOR_OVERSATURATED: All LP slots on acceptor atom already used
- POOR_ALIGNMENT: Alignment score below threshold
- DISTANCE_TOO_FAR: Distance > max_distance cutoff
- BIFURCATION_ANGLE_SMALL: Would bifurcate but angle too small
- NOT_VALID_DONOR_ACCEPTOR: Atoms not in donor/acceptor capacity tables
- NOT_BASE_PAIRED: Residues not base-paired (C++ only checks base pairs)

Usage:
    python tools/analyze_dssr_misses.py 1EHZ
    python tools/analyze_dssr_misses.py 1EHZ --verbose
    python tools/analyze_dssr_misses.py 1EHZ --use-cpp  # Use C++ output instead of prototype
"""

import argparse
import json
import sys
import numpy as np
from enum import Enum
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass, field

# Add project root and prototype to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class MissReason(Enum):
    """Reasons why an H-bond was missed."""
    DONOR_OVERSATURATED = "DONOR_OVERSATURATED"
    ACCEPTOR_OVERSATURATED = "ACCEPTOR_OVERSATURATED"
    POOR_ALIGNMENT = "POOR_ALIGNMENT"
    DISTANCE_TOO_FAR = "DISTANCE_TOO_FAR"
    BIFURCATION_ANGLE_SMALL = "BIFURCATION_ANGLE_SMALL"
    NOT_VALID_DONOR_ACCEPTOR = "NOT_VALID_DONOR_ACCEPTOR"
    NOT_BASE_PAIRED = "NOT_BASE_PAIRED"
    CONFLICT_LOST = "CONFLICT_LOST"
    UNKNOWN = "UNKNOWN"


@dataclass
class SlotUsageInfo:
    """Information about slot usage for a donor or acceptor."""
    atom_name: str
    res_id: str
    slot_type: str  # 'H' or 'LP'
    slots_total: int
    slots_used: int
    bond_directions: List[np.ndarray] = field(default_factory=list)
    competing_bonds: List[Dict] = field(default_factory=list)


@dataclass
class MissAnalysis:
    """Detailed analysis of why we missed a DSSR H-bond."""
    dssr_hb: Dict
    reason: MissReason
    details: str
    category: str  # base-base, base-sugar, etc.

    # Enhanced info from prototype analysis
    donor_slot_info: Optional[SlotUsageInfo] = None
    acceptor_slot_info: Optional[SlotUsageInfo] = None
    alignment_score: Optional[float] = None
    computed_distance: Optional[float] = None
    bifurcation_angle: Optional[float] = None
    competing_bond: Optional[Dict] = None


def get_atom_category(atom_name: str) -> str:
    """Categorize atom as base, sugar, or backbone."""
    name = atom_name.strip().upper()
    if "'" in name:
        return 'sugar'
    if name in ['P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P']:
        return 'backbone'
    return 'base'


def parse_dssr_atom_id(atom_id: str) -> Tuple[str, str, str]:
    """Parse DSSR atom ID into (atom_name, chain, resnum)."""
    if '@' not in atom_id:
        return atom_id, '', ''
    atom, dssr_res = atom_id.split('@', 1)
    if '.' not in dssr_res:
        return atom, '', dssr_res
    chain, rest = dssr_res.split('.', 1)
    resnum = ''
    for i in range(len(rest) - 1, -1, -1):
        if rest[i].isdigit():
            continue
        else:
            resnum = rest[i+1:]
            break
    if not resnum:
        for i, c in enumerate(rest):
            if c.isdigit():
                resnum = rest[i:]
                break
    return atom, chain, resnum


def dssr_res_id_to_internal(dssr_res: str) -> str:
    """Convert DSSR res ID to internal format.

    Note: The prototype uses DSSR format directly (A.G103), so we just
    normalize the format (handle insertion codes).
    """
    # Handle insertion code (^A format) -> convert to inline
    if '^' in dssr_res:
        dssr_res = dssr_res.replace('^', '')
    return dssr_res


def to_dssr_res_id(res_id: str) -> str:
    """Convert our res_id to DSSR format for matching."""
    parts = res_id.split('-')
    if len(parts) < 3:
        return res_id
    chain = parts[0]
    name = '-'.join(parts[1:-1])
    num_part = parts[-1]
    num = ''
    insertion = ''
    for i, c in enumerate(num_part):
        if c.isdigit() or c == '-':
            num += c
        else:
            insertion = num_part[i:]
            break
    if not num:
        num = num_part
    result = f"{chain}.{name}{num}"
    if insertion:
        result += f"^{insertion}"
    return result


def normalize_atom_name(name: str) -> str:
    """Normalize atom names for matching."""
    name = name.strip()
    atom_map = {'O1P': 'OP1', 'O2P': 'OP2', 'O3P': 'OP3'}
    return atom_map.get(name, name)


def make_key(atom1_id: str, atom2_id: str) -> Tuple[str, str]:
    """Create sorted key for matching."""
    return tuple(sorted([atom1_id, atom2_id]))


def load_dssr_hbonds(pdb_id: str, project_root: Path) -> List[Dict]:
    """Load DSSR H-bonds, excluding questionable."""
    dssr_file = project_root / "data" / "json_dssr" / f"{pdb_id}.json"
    if not dssr_file.exists():
        raise FileNotFoundError(f"DSSR file not found: {dssr_file}")

    with open(dssr_file) as f:
        data = json.load(f)

    hbonds = []
    for hb in data.get('hbonds', []):
        if hb.get('donAcc_type') != 'questionable':
            hbonds.append(hb)
    return hbonds


class PrototypeMissAnalyzer:
    """Analyze DSSR misses using the Python H-bond optimizer prototype."""

    def __init__(self, pdb_id: str, project_root: Path,
                 max_distance: float = 4.0,
                 min_alignment: float = 0.3,
                 min_bifurcation_angle: float = 45.0,
                 min_bifurcation_alignment: float = 0.5):
        self.pdb_id = pdb_id
        self.project_root = project_root
        self.max_distance = max_distance
        self.min_alignment = min_alignment
        self.min_bifurcation_angle = min_bifurcation_angle
        self.min_bifurcation_alignment = min_bifurcation_alignment

        # Import prototype modules
        from prototypes.hbond_optimizer import (
            HBondOptimizer, parse_pdb_residues,
            DONOR_CAPACITY, ACCEPTOR_CAPACITY,
            predict_h_slots, predict_lp_slots, score_hbond_alignment,
            normalize, angle_between
        )

        self.HBondOptimizer = HBondOptimizer
        self.DONOR_CAPACITY = DONOR_CAPACITY
        self.ACCEPTOR_CAPACITY = ACCEPTOR_CAPACITY
        self.predict_h_slots = predict_h_slots
        self.predict_lp_slots = predict_lp_slots
        self.score_hbond_alignment = score_hbond_alignment
        self.normalize = normalize
        self.angle_between = angle_between

        # Load PDB residues
        pdb_path = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        self.residues = parse_pdb_residues(pdb_path)

        # Load DSSR H-bonds
        self.dssr_hbonds = load_dssr_hbonds(pdb_id, project_root)

        # Run optimizer on all DSSR residue pairs to establish slot state
        self._run_optimizer()

    def _get_dssr_residue_pairs(self) -> Set[Tuple[str, str]]:
        """Get unique residue pairs from DSSR H-bonds."""
        pairs = set()
        for hb in self.dssr_hbonds:
            atom1_id = hb.get('atom1_id', '')
            atom2_id = hb.get('atom2_id', '')

            if '@' in atom1_id and '@' in atom2_id:
                res1 = dssr_res_id_to_internal(atom1_id.split('@')[1])
                res2 = dssr_res_id_to_internal(atom2_id.split('@')[1])
                pairs.add(tuple(sorted([res1, res2])))
        return pairs

    def _run_optimizer(self):
        """Run optimizer on all DSSR residue pairs."""
        self.optimizer = self.HBondOptimizer(
            max_distance=self.max_distance,
            min_alignment=self.min_alignment,
            min_bifurcation_angle=self.min_bifurcation_angle,
            min_bifurcation_alignment=self.min_bifurcation_alignment
        )

        # Add all residues
        for res_id, residue in self.residues.items():
            self.optimizer.add_residue(residue)

        # Run on all DSSR pairs
        pairs = self._get_dssr_residue_pairs()
        self.selected_hbonds = {}  # (res1, res2) -> [HBond, ...]

        for res1, res2 in pairs:
            if res1 in self.residues and res2 in self.residues:
                hbonds = self.optimizer.optimize_pair(res1, res2)
                self.selected_hbonds[(res1, res2)] = hbonds

    def _get_selected_key(self, donor_res: str, donor_atom: str,
                          acceptor_res: str, acceptor_atom: str) -> Optional[str]:
        """Check if this H-bond was selected by the optimizer."""
        pair_key = tuple(sorted([donor_res, acceptor_res]))
        if pair_key not in self.selected_hbonds:
            return None

        for hb in self.selected_hbonds[pair_key]:
            if ((hb.donor_res_id == donor_res and hb.donor_atom == donor_atom and
                 hb.acceptor_res_id == acceptor_res and hb.acceptor_atom == acceptor_atom) or
                (hb.donor_res_id == acceptor_res and hb.donor_atom == acceptor_atom and
                 hb.acceptor_res_id == donor_res and hb.acceptor_atom == donor_atom)):
                return "MATCHED"
        return None

    def analyze_miss(self, dssr_hb: Dict) -> MissAnalysis:
        """Analyze why a specific DSSR H-bond was missed."""
        atom1_id = dssr_hb.get('atom1_id', '')
        atom2_id = dssr_hb.get('atom2_id', '')
        distance = dssr_hb.get('distance', 0.0)

        # Parse atom IDs
        atom1, _, _ = parse_dssr_atom_id(atom1_id)
        atom2, _, _ = parse_dssr_atom_id(atom2_id)

        res1 = dssr_res_id_to_internal(atom1_id.split('@')[1]) if '@' in atom1_id else ''
        res2 = dssr_res_id_to_internal(atom2_id.split('@')[1]) if '@' in atom2_id else ''

        # Get category
        cat1 = get_atom_category(atom1)
        cat2 = get_atom_category(atom2)
        category = '-'.join(sorted([cat1, cat2]))

        # Check if we have the residues
        if res1 not in self.residues or res2 not in self.residues:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.NOT_BASE_PAIRED,
                details=f"Residue not found in PDB parsing",
                category=category
            )

        residue1 = self.residues[res1]
        residue2 = self.residues[res2]

        # Determine donor/acceptor (try both directions)
        analysis = self._analyze_donor_acceptor(
            dssr_hb, residue1, atom1, residue2, atom2, category
        )
        if analysis.reason != MissReason.NOT_VALID_DONOR_ACCEPTOR:
            return analysis

        # Try reverse direction
        return self._analyze_donor_acceptor(
            dssr_hb, residue2, atom2, residue1, atom1, category
        )

    def _analyze_donor_acceptor(self, dssr_hb: Dict,
                                 donor_res, donor_atom: str,
                                 acceptor_res, acceptor_atom: str,
                                 category: str) -> MissAnalysis:
        """Analyze a specific donor->acceptor direction."""
        distance = dssr_hb.get('distance', 0.0)

        # Check if valid donor/acceptor
        donor_key = (donor_res.base_type.upper(), donor_atom.strip())
        acceptor_key = (acceptor_res.base_type.upper(), acceptor_atom.strip())

        if donor_key not in self.DONOR_CAPACITY:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.NOT_VALID_DONOR_ACCEPTOR,
                details=f"{donor_atom} on {donor_res.base_type} is not a known donor",
                category=category
            )

        if acceptor_key not in self.ACCEPTOR_CAPACITY:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.NOT_VALID_DONOR_ACCEPTOR,
                details=f"{acceptor_atom} on {acceptor_res.base_type} is not a known acceptor",
                category=category
            )

        # Check distance
        if distance > self.max_distance:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.DISTANCE_TOO_FAR,
                details=f"Distance {distance:.2f}A > {self.max_distance}A cutoff",
                category=category,
                computed_distance=distance
            )

        # Get atom positions
        if donor_atom not in donor_res.atoms or acceptor_atom not in acceptor_res.atoms:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.UNKNOWN,
                details=f"Atom positions not found in residue",
                category=category
            )

        donor_pos = donor_res.atoms[donor_atom]
        acceptor_pos = acceptor_res.atoms[acceptor_atom]

        # Get slots
        h_slots = donor_res.get_h_slots(donor_atom)
        lp_slots = acceptor_res.get_lp_slots(acceptor_atom)

        if not h_slots:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.DONOR_OVERSATURATED,
                details=f"No H slots available for {donor_atom}",
                category=category,
                donor_slot_info=SlotUsageInfo(
                    atom_name=donor_atom,
                    res_id=donor_res.res_id,
                    slot_type='H',
                    slots_total=0,
                    slots_used=0
                )
            )

        if not lp_slots:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.ACCEPTOR_OVERSATURATED,
                details=f"No LP slots available for {acceptor_atom}",
                category=category,
                acceptor_slot_info=SlotUsageInfo(
                    atom_name=acceptor_atom,
                    res_id=acceptor_res.res_id,
                    slot_type='LP',
                    slots_total=0,
                    slots_used=0
                )
            )

        # Compute alignment
        h_idx, lp_idx, alignment_score = self.score_hbond_alignment(
            donor_pos, acceptor_pos, h_slots, lp_slots
        )

        # Check slot availability
        h_slot = h_slots[h_idx]
        lp_slot = lp_slots[lp_idx]

        donor_to_acceptor = self.normalize(acceptor_pos - donor_pos)
        acceptor_to_donor = -donor_to_acceptor

        # Check if slots are available or can bifurcate
        h_slots_used = len(h_slot.bond_directions)
        lp_slots_used = len(lp_slot.bond_directions)

        h_can_add = h_slot.can_add_bond(donor_to_acceptor, self.min_bifurcation_angle)
        lp_can_add = lp_slot.can_add_bond(acceptor_to_donor, self.min_bifurcation_angle)

        # Build slot info
        donor_slot_info = SlotUsageInfo(
            atom_name=donor_atom,
            res_id=donor_res.res_id,
            slot_type='H',
            slots_total=len(h_slots),
            slots_used=h_slots_used,
            bond_directions=[d.tolist() for d in h_slot.bond_directions]
        )

        acceptor_slot_info = SlotUsageInfo(
            atom_name=acceptor_atom,
            res_id=acceptor_res.res_id,
            slot_type='LP',
            slots_total=len(lp_slots),
            slots_used=lp_slots_used,
            bond_directions=[d.tolist() for d in lp_slot.bond_directions]
        )

        # Check if donor slot is full
        if not h_can_add:
            # Calculate angle to existing bonds
            if h_slot.bond_directions:
                angles = [self.angle_between(donor_to_acceptor, d)
                          for d in h_slot.bond_directions]
                min_angle = min(angles)
                return MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason=MissReason.BIFURCATION_ANGLE_SMALL,
                    details=f"H slot used, bifurcation angle {min_angle:.1f}° < {self.min_bifurcation_angle}°",
                    category=category,
                    donor_slot_info=donor_slot_info,
                    acceptor_slot_info=acceptor_slot_info,
                    bifurcation_angle=min_angle
                )
            else:
                return MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason=MissReason.DONOR_OVERSATURATED,
                    details=f"All {len(h_slots)} H slots on {donor_atom} fully used",
                    category=category,
                    donor_slot_info=donor_slot_info
                )

        # Check if acceptor slot is full
        if not lp_can_add:
            if lp_slot.bond_directions:
                angles = [self.angle_between(acceptor_to_donor, d)
                          for d in lp_slot.bond_directions]
                min_angle = min(angles)
                return MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason=MissReason.BIFURCATION_ANGLE_SMALL,
                    details=f"LP slot used, bifurcation angle {min_angle:.1f}° < {self.min_bifurcation_angle}°",
                    category=category,
                    donor_slot_info=donor_slot_info,
                    acceptor_slot_info=acceptor_slot_info,
                    bifurcation_angle=min_angle
                )
            else:
                return MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason=MissReason.ACCEPTOR_OVERSATURATED,
                    details=f"All {len(lp_slots)} LP slots on {acceptor_atom} fully used",
                    category=category,
                    acceptor_slot_info=acceptor_slot_info
                )

        # Check alignment threshold
        is_bifurcated = h_slots_used > 0 or lp_slots_used > 0
        min_align = self.min_bifurcation_alignment if is_bifurcated else self.min_alignment

        if alignment_score < min_align:
            return MissAnalysis(
                dssr_hb=dssr_hb,
                reason=MissReason.POOR_ALIGNMENT,
                details=f"Alignment score {alignment_score:.2f} < {min_align} threshold" +
                        (" (bifurcated)" if is_bifurcated else ""),
                category=category,
                donor_slot_info=donor_slot_info,
                acceptor_slot_info=acceptor_slot_info,
                alignment_score=alignment_score
            )

        # If we get here, it should have been selected - check competing bonds
        pair_key = tuple(sorted([donor_res.res_id, acceptor_res.res_id]))
        if pair_key in self.selected_hbonds:
            selected = self.selected_hbonds[pair_key]
            if selected:
                winner = selected[0]  # First (best) bond
                return MissAnalysis(
                    dssr_hb=dssr_hb,
                    reason=MissReason.CONFLICT_LOST,
                    details=f"Lost to {winner.donor_atom}->{winner.acceptor_atom} (d={winner.distance:.2f}A)",
                    category=category,
                    donor_slot_info=donor_slot_info,
                    acceptor_slot_info=acceptor_slot_info,
                    competing_bond={'donor': winner.donor_atom, 'acceptor': winner.acceptor_atom,
                                    'distance': winner.distance}
                )

        return MissAnalysis(
            dssr_hb=dssr_hb,
            reason=MissReason.UNKNOWN,
            details=f"Passed all checks but not selected - needs investigation",
            category=category,
            donor_slot_info=donor_slot_info,
            acceptor_slot_info=acceptor_slot_info,
            alignment_score=alignment_score
        )

    def analyze_all_misses(self) -> List[MissAnalysis]:
        """Analyze all DSSR H-bonds we miss."""
        analyses = []

        for dssr_hb in self.dssr_hbonds:
            atom1_id = dssr_hb.get('atom1_id', '')
            atom2_id = dssr_hb.get('atom2_id', '')

            if '@' not in atom1_id or '@' not in atom2_id:
                continue

            # Convert to internal format
            atom1, _, _ = parse_dssr_atom_id(atom1_id)
            atom2, _, _ = parse_dssr_atom_id(atom2_id)
            res1 = dssr_res_id_to_internal(atom1_id.split('@')[1])
            res2 = dssr_res_id_to_internal(atom2_id.split('@')[1])

            # Check if we selected this bond
            if self._get_selected_key(res1, atom1, res2, atom2) == "MATCHED":
                continue

            # Analyze why we missed it
            analysis = self.analyze_miss(dssr_hb)
            analyses.append(analysis)

        return analyses


def print_analysis(analyses: List[MissAnalysis], verbose: bool = False):
    """Print analysis summary."""

    # Group by reason
    by_reason = {}
    for a in analyses:
        reason = a.reason.value
        if reason not in by_reason:
            by_reason[reason] = []
        by_reason[reason].append(a)

    print(f"\n{'='*70}")
    print(f"DSSR MISS ANALYSIS (Slot-Aware)")
    print(f"{'='*70}")
    print(f"Total misses: {len(analyses)}")
    print()

    # Summary by reason
    print("By reason:")
    for reason, items in sorted(by_reason.items(), key=lambda x: -len(x[1])):
        print(f"  {reason}: {len(items)}")
    print()

    # Summary by category
    by_category = {}
    for a in analyses:
        if a.category not in by_category:
            by_category[a.category] = []
        by_category[a.category].append(a)

    print("By atom category:")
    for cat, items in sorted(by_category.items(), key=lambda x: -len(x[1])):
        print(f"  {cat}: {len(items)}")
    print()

    # Details
    if verbose:
        for reason, items in sorted(by_reason.items(), key=lambda x: -len(x[1])):
            print(f"\n--- {reason} ({len(items)}) ---")
            for a in items[:10]:
                hb = a.dssr_hb
                print(f"  {hb['atom1_id']} -- {hb['atom2_id']} (d={hb['distance']:.2f}A, {hb.get('donAcc_type', '')})")
                print(f"    -> {a.details}")

                if a.donor_slot_info:
                    si = a.donor_slot_info
                    print(f"       Donor: {si.atom_name} H-slots {si.slots_used}/{si.slots_total}")
                if a.acceptor_slot_info:
                    si = a.acceptor_slot_info
                    print(f"       Acceptor: {si.atom_name} LP-slots {si.slots_used}/{si.slots_total}")
                if a.alignment_score is not None:
                    print(f"       Alignment: {a.alignment_score:.2f}")
                if a.bifurcation_angle is not None:
                    print(f"       Bifurcation angle: {a.bifurcation_angle:.1f}°")

            if len(items) > 10:
                print(f"  ... and {len(items) - 10} more")


def main():
    parser = argparse.ArgumentParser(description="Analyze why we miss DSSR H-bonds (slot-aware)")
    parser.add_argument("pdb_id", help="PDB ID to analyze")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed output")
    parser.add_argument("--max-dist", type=float, default=4.0, help="Max distance (default: 4.0)")
    parser.add_argument("--min-align", type=float, default=0.3, help="Min alignment (default: 0.3)")
    parser.add_argument("--min-bif-angle", type=float, default=45.0,
                        help="Min bifurcation angle (default: 45.0)")

    args = parser.parse_args()
    project_root = Path(__file__).parent.parent
    pdb_id = args.pdb_id.upper()

    try:
        print(f"Analyzing {pdb_id} with prototype H-bond optimizer...")
        print(f"Parameters: max_dist={args.max_dist}A, min_align={args.min_align}, "
              f"min_bif_angle={args.min_bif_angle}°")

        analyzer = PrototypeMissAnalyzer(
            pdb_id, project_root,
            max_distance=args.max_dist,
            min_alignment=args.min_align,
            min_bifurcation_angle=args.min_bif_angle
        )

        # Count matches
        total_dssr = len(analyzer.dssr_hbonds)
        analyses = analyzer.analyze_all_misses()
        matches = total_dssr - len(analyses)

        print(f"\nDSSR H-bonds: {total_dssr}")
        print(f"Matched: {matches} ({matches/total_dssr*100:.1f}%)")
        print(f"Missed: {len(analyses)} ({len(analyses)/total_dssr*100:.1f}%)")

        print_analysis(analyses, args.verbose)

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
