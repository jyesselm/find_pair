"""Data loaders for cWW miss annotator tool.

Loads DSSR pairs and slot-based H-bond data for comparison and analysis.
"""

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Import from parent directory core modules
sys.path.insert(0, str(Path(__file__).parent.parent))
from core.identifiers import dssr_to_res_id, res_id_to_dssr


def dssr_to_slot_id(dssr_nt: str) -> Optional[str]:
    """Convert DSSR nucleotide format to slot format.

    Args:
        dssr_nt: DSSR format like "A.G1", "0.C530", "D.DG10"

    Returns:
        Slot format like "A-G-1", "0-C-530", "D-DG-10", or None if parsing fails

    Examples:
        >>> dssr_to_slot_id("A.G1")
        "A-G-1"
        >>> dssr_to_slot_id("0.C530")
        "0-C-530"
        >>> dssr_to_slot_id("D.DG10")
        "D-DG-10"
    """
    # Use core module implementation
    return dssr_to_res_id(dssr_nt)


def slot_to_dssr_id(slot_id: str) -> str:
    """Convert slot format to DSSR nucleotide format.

    Args:
        slot_id: Slot format like "A-G-1", "0-C-530", "D-DG-10"

    Returns:
        DSSR format like "A.G1", "0.C530", "D.DG10"

    Examples:
        >>> slot_to_dssr_id("A-G-1")
        "A.G1"
        >>> slot_to_dssr_id("0-C-530")
        "0.C530"
        >>> slot_to_dssr_id("D-DG-10")
        "D.DG10"
    """
    # Use core module implementation
    result = res_id_to_dssr(slot_id)
    if result is None:
        raise ValueError(f"Invalid slot_id format: {slot_id}")
    return result


@dataclass
class DSSRPair:
    """DSSR base pair record."""
    nt1: str
    nt2: str
    lw_class: str
    saenger: str
    hbonds_desc: str
    interbase_angle: float
    n1n9_distance: float
    bp: str


def load_dssr_pairs(
    dssr_path: Path,
    lw_filter: str = "cWW"
) -> Dict[Tuple[str, str], DSSRPair]:
    """Load DSSR pairs and filter by Leontis-Westhof class.

    Args:
        dssr_path: Path to DSSR JSON file
        lw_filter: LW class to filter (default "cWW")

    Returns:
        Dict keyed by (slot_id1, slot_id2) tuple with DSSRPair values.
        Only includes standard sequences: GC, CG, AU, UA

    Examples:
        >>> pairs = load_dssr_pairs(Path("data/json_dssr/1EHZ.json"))
        >>> len(pairs)
        42
    """
    with open(dssr_path) as f:
        data = json.load(f)

    pairs = {}
    standard_seqs = {"GC", "CG", "AU", "UA"}

    for pair in data.get("pairs", []):
        lw_class = pair.get("LW", "")
        if lw_class != lw_filter:
            continue

        bp = pair.get("bp", "")
        # Extract base letters only (handles "C-G", "G-C", "U-A", etc.)
        bp_clean = bp.replace("-", "").replace("+", "")
        if bp_clean not in standard_seqs:
            continue

        nt1_dssr = pair.get("nt1", "")
        nt2_dssr = pair.get("nt2", "")

        nt1_slot = dssr_to_slot_id(nt1_dssr)
        nt2_slot = dssr_to_slot_id(nt2_dssr)

        if not nt1_slot or not nt2_slot:
            continue

        dssr_pair = DSSRPair(
            nt1=nt1_slot,
            nt2=nt2_slot,
            lw_class=lw_class,
            saenger=pair.get("Saenger", ""),
            hbonds_desc=pair.get("hbonds_desc", ""),
            interbase_angle=pair.get("interBase_angle", 0.0),
            n1n9_distance=pair.get("N1N9_dist", 0.0),
            bp=bp
        )

        pairs[(nt1_slot, nt2_slot)] = dssr_pair

    return pairs


@dataclass
class SlotHBond:
    """Slot-based hydrogen bond record."""
    donor_res_id: str
    donor_atom: str
    acceptor_res_id: str
    acceptor_atom: str
    distance: float
    context: str
    h_slot: Optional[int]
    lp_slot: Optional[int]
    alignment: Optional[float]


def load_slot_hbonds(hbond_path: Path) -> Dict[Tuple[str, str], List[SlotHBond]]:
    """Load slot-based H-bond data.

    Args:
        hbond_path: Path to slot hbonds JSON file

    Returns:
        Dict keyed by (res_id_i, res_id_j) tuple with list of SlotHBond values.
        Also includes reverse lookup (res_id_j, res_id_i) for bidirectional access.

    Examples:
        >>> hbonds = load_slot_hbonds(Path("data/json/slot_hbonds/1EHZ.json"))
        >>> len(hbonds)
        168
    """
    with open(hbond_path) as f:
        data = json.load(f)

    hbonds = {}

    for record in data:
        res_id_i = record.get("res_id_i", "")
        res_id_j = record.get("res_id_j", "")

        slot_hbonds = []
        for hbond in record.get("hbonds", []):
            slot_hbond = SlotHBond(
                donor_res_id=hbond.get("donor_res_id", res_id_i),
                donor_atom=hbond.get("donor_atom", ""),
                acceptor_res_id=hbond.get("acceptor_res_id", res_id_j),
                acceptor_atom=hbond.get("acceptor_atom", ""),
                distance=hbond.get("distance", 0.0),
                context=hbond.get("context", ""),
                h_slot=hbond.get("h_slot"),
                lp_slot=hbond.get("lp_slot"),
                alignment=hbond.get("alignment")
            )
            slot_hbonds.append(slot_hbond)

        # Store both directions for bidirectional lookup
        hbonds[(res_id_i, res_id_j)] = slot_hbonds
        hbonds[(res_id_j, res_id_i)] = slot_hbonds

    return hbonds
