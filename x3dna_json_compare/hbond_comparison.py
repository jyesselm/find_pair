"""
Hydrogen bond list comparison utilities.

Compares hbond_list records between legacy and modern JSON outputs.
Each hbond_list record contains detailed hydrogen bond information for a base pair.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field


@dataclass
class HBondComparison:
    """Result of hbond_list comparison for a single base pair."""

    missing_in_modern: List[Dict] = field(default_factory=list)
    extra_in_modern: List[Dict] = field(default_factory=list)
    mismatched_hbonds: List[Dict] = field(default_factory=list)
    num_hbonds_match: bool = True
    legacy_num_hbonds: int = 0
    modern_num_hbonds: int = 0


@dataclass
class HBondListComparison:
    """Result of hbond_list comparison across all pairs."""

    missing_in_modern: List[Dict] = field(default_factory=list)
    extra_in_modern: List[Dict] = field(default_factory=list)
    mismatched_pairs: List[Dict] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0
    pair_comparisons: Dict[Tuple[int, int], HBondComparison] = field(
        default_factory=dict
    )


def normalize_pair(base_i: int, base_j: int) -> Tuple[int, int]:
    """Normalize pair order to (min, max) for consistent comparison."""
    return (min(base_i, base_j), max(base_i, base_j))


def normalize_hbond(hbond: Dict) -> Dict:
    """Normalize H-bond for comparison (handle atom name variations)."""
    donor = hbond.get("donor_atom", "").strip()
    acceptor = hbond.get("acceptor_atom", "").strip()
    distance = round(hbond.get("distance", 0.0), 6)  # More precision for distance
    hb_type = hbond.get("type", " ")

    return {"donor": donor, "acceptor": acceptor, "distance": distance, "type": hb_type}


def compare_pair_hbonds(
    legacy_hbonds: List[Dict], modern_hbonds: List[Dict], tolerance: float = 1e-6
) -> HBondComparison:
    """
    Compare H-bonds for a single base pair.

    Args:
        legacy_hbonds: List of legacy H-bond records
        modern_hbonds: List of modern H-bond records
        tolerance: Tolerance for distance comparisons

    Returns:
        HBondComparison result
    """
    result = HBondComparison()

    # Normalize H-bonds for comparison
    legacy_norm = [normalize_hbond(hb) for hb in legacy_hbonds]
    modern_norm = [normalize_hbond(hb) for hb in modern_hbonds]

    result.legacy_num_hbonds = len(legacy_hbonds)
    result.modern_num_hbonds = len(modern_hbonds)
    result.num_hbonds_match = result.legacy_num_hbonds == result.modern_num_hbonds

    # Create sets for comparison using (donor, acceptor, distance, type)
    # Normalize atom order for donor/acceptor (some H-bonds may be recorded in reverse)
    legacy_set = set()
    modern_set = set()

    for hb in legacy_norm:
        # Create normalized key: sorted atoms + distance + type
        atoms_sorted = tuple(sorted([hb["donor"], hb["acceptor"]]))
        key = atoms_sorted + (hb["distance"], hb["type"])
        legacy_set.add(key)

    for hb in modern_norm:
        atoms_sorted = tuple(sorted([hb["donor"], hb["acceptor"]]))
        key = atoms_sorted + (hb["distance"], hb["type"])
        modern_set.add(key)

    missing_keys = legacy_set - modern_set
    extra_keys = modern_set - legacy_set
    common_keys = legacy_set & modern_set

    # Find missing H-bonds (with original atom order preserved)
    for hb in legacy_norm:
        atoms_sorted = tuple(sorted([hb["donor"], hb["acceptor"]]))
        key = atoms_sorted + (hb["distance"], hb["type"])
        if key in missing_keys:
            result.missing_in_modern.append(
                {
                    "donor_atom": hb["donor"],
                    "acceptor_atom": hb["acceptor"],
                    "distance": hb["distance"],
                    "type": hb["type"],
                }
            )

    # Find extra H-bonds
    for hb in modern_norm:
        atoms_sorted = tuple(sorted([hb["donor"], hb["acceptor"]]))
        key = atoms_sorted + (hb["distance"], hb["type"])
        if key in extra_keys:
            result.extra_in_modern.append(
                {
                    "donor_atom": hb["donor"],
                    "acceptor_atom": hb["acceptor"],
                    "distance": hb["distance"],
                    "type": hb["type"],
                }
            )

    # Check for mismatches in common H-bonds (distance differences within tolerance)
    # This handles cases where H-bonds match but distances differ slightly
    legacy_by_key = {}
    modern_by_key = {}

    for hb in legacy_norm:
        atoms_sorted = tuple(sorted([hb["donor"], hb["acceptor"]]))
        key = atoms_sorted + (hb["distance"], hb["type"])
        legacy_by_key[key] = hb

    for hb in modern_norm:
        atoms_sorted = tuple(sorted([hb["donor"], hb["acceptor"]]))
        key = atoms_sorted + (hb["distance"], hb["type"])
        modern_by_key[key] = hb

    # Check if there are H-bonds with same atoms and type but different distances
    for leg_key in legacy_set:
        if leg_key in common_keys:
            continue  # Exact match

        # Check if there's a modern H-bond with same atoms and type but different distance
        atoms_sorted = leg_key[:2]
        hb_type = leg_key[3]
        leg_dist = leg_key[2]

        for mod_key in modern_set:
            if (
                mod_key[:2] == atoms_sorted
                and mod_key[3] == hb_type
                and abs(mod_key[2] - leg_dist) > tolerance
            ):
                # Found a mismatch - same atoms/type but different distance
                result.mismatched_hbonds.append(
                    {
                        "donor_atom": (
                            atoms_sorted[0]
                            if atoms_sorted[0] < atoms_sorted[1]
                            else atoms_sorted[1]
                        ),
                        "acceptor_atom": (
                            atoms_sorted[1]
                            if atoms_sorted[0] < atoms_sorted[1]
                            else atoms_sorted[0]
                        ),
                        "legacy_distance": leg_dist,
                        "modern_distance": mod_key[2],
                        "type": hb_type,
                        "distance_diff": abs(mod_key[2] - leg_dist),
                    }
                )
                break

    return result


def compare_hbond_lists(
    legacy_records: List[Dict], modern_records: List[Dict], tolerance: float = 1e-6
) -> HBondListComparison:
    """
    Compare hbond_list records between legacy and modern JSON.

    Args:
        legacy_records: List of legacy hbond_list records
        modern_records: List of modern hbond_list records
        tolerance: Tolerance for floating point comparisons

    Returns:
        HBondListComparison result
    """
    # Deduplicate legacy records (legacy has duplicate record bug - generates extra copies)
    seen_pairs = set()
    legacy_dedup = []
    for rec in legacy_records:
        norm_pair = (min(rec['base_i'], rec['base_j']), max(rec['base_i'], rec['base_j']))
        if norm_pair in seen_pairs:
            continue
        seen_pairs.add(norm_pair)
        legacy_dedup.append(rec)
    
    # Deduplicate modern records (modern generates both (i,j) and (j,i))
    seen_pairs_modern = set()
    modern_dedup = []
    for rec in modern_records:
        norm_pair = (min(rec['base_i'], rec['base_j']), max(rec['base_i'], rec['base_j']))
        if norm_pair in seen_pairs_modern:
            continue
        seen_pairs_modern.add(norm_pair)
        modern_dedup.append(rec)
    
    result = HBondListComparison()

    # Build maps by normalized (base_i, base_j)
    legacy_map = {}
    modern_map = {}

    for rec in legacy_dedup:
        if rec.get("type") != "hbond_list":
            continue
        base_i = rec.get("base_i")
        base_j = rec.get("base_j")
        if base_i is not None and base_j is not None:
            key = normalize_pair(base_i, base_j)
            legacy_map[key] = rec

    # Modern now uses 1-based indices (same as legacy)
    for rec in modern_records:
        if rec.get("type") != "hbond_list":
            continue
        base_i = rec.get("base_i")
        base_j = rec.get("base_j")
        if base_i is not None and base_j is not None:
            key = normalize_pair(base_i, base_j)
            modern_map[key] = rec

    result.total_legacy = len(legacy_map)
    result.total_modern = len(modern_map)

    legacy_keys = set(legacy_map.keys())
    modern_keys = set(modern_map.keys())
    common_keys = legacy_keys & modern_keys

    result.common_count = len(common_keys)

    # Find missing in modern
    for key in legacy_keys - modern_keys:
        leg_rec = legacy_map[key]
        result.missing_in_modern.append(
            {
                "base_i": key[0],
                "base_j": key[1],
                "legacy_record": leg_rec,
                "num_hbonds": leg_rec.get("num_hbonds", 0),
            }
        )

    # Find extra in modern
    for key in modern_keys - legacy_keys:
        mod_rec = modern_map[key]
        result.extra_in_modern.append(
            {
                "base_i": key[0],
                "base_j": key[1],
                "modern_record": mod_rec,
                "num_hbonds": mod_rec.get("num_hbonds", 0),
            }
        )

    # Compare common pairs
    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]

        leg_hbonds = leg_rec.get("hbonds", [])
        mod_hbonds = mod_rec.get("hbonds", [])

        pair_comp = compare_pair_hbonds(leg_hbonds, mod_hbonds, tolerance)
        result.pair_comparisons[key] = pair_comp

        # Check for mismatches
        if (
            pair_comp.missing_in_modern
            or pair_comp.extra_in_modern
            or pair_comp.mismatched_hbonds
            or not pair_comp.num_hbonds_match
        ):
            result.mismatched_pairs.append(
                {
                    "base_i": key[0],
                    "base_j": key[1],
                    "legacy_num_hbonds": leg_rec.get("num_hbonds", 0),
                    "modern_num_hbonds": mod_rec.get("num_hbonds", 0),
                    "comparison": pair_comp,
                    "legacy_record": leg_rec,
                    "modern_record": mod_rec,
                }
            )

    return result
