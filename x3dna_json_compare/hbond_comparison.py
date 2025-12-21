"""
Hydrogen bond list comparison utilities.

Compares hbond_list records between legacy and modern JSON outputs.
Each hbond_list record contains detailed hydrogen bond information for a base pair.

Uses res_id-based matching for stable comparison that doesn't depend on index ordering.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
from .res_id_utils import get_res_id_i, get_res_id_j, make_pair_key
from .step_comparison import build_residue_idx_to_res_id_map


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
    matched_count: int = 0
    pair_comparisons: Dict[Tuple, HBondComparison] = field(
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

        atoms_sorted = leg_key[:2]
        hb_type = leg_key[3]
        leg_dist = leg_key[2]

        for mod_key in modern_set:
            if (
                mod_key[:2] == atoms_sorted
                and mod_key[3] == hb_type
                and abs(mod_key[2] - leg_dist) > tolerance
            ):
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


def get_hbond_pair_key(rec: Dict, residue_map: Optional[Dict[int, str]] = None) -> Optional[Tuple[str, str]]:
    """
    Get normalized pair key for an hbond_list record.

    Uses res_id_i/res_id_j if available (modern JSON), otherwise falls back
    to constructing from base_i/base_j using the residue map.

    Args:
        rec: hbond_list record
        residue_map: Optional mapping from residue_idx -> res_id

    Returns:
        Normalized (res_id_1, res_id_2) tuple or None
    """
    # Try res_id fields first (modern JSON)
    res_id_i = get_res_id_i(rec)
    res_id_j = get_res_id_j(rec)

    if res_id_i and res_id_j:
        return make_pair_key(res_id_i, res_id_j)

    # Fall back to constructing from base_i/base_j using residue map
    if residue_map:
        base_i = rec.get('base_i')
        base_j = rec.get('base_j')
        if base_i is not None and base_j is not None:
            res_id_i = residue_map.get(base_i)
            res_id_j = residue_map.get(base_j)
            if res_id_i and res_id_j:
                return make_pair_key(res_id_i, res_id_j)

    return None


def compare_hbond_lists(
    legacy_records: List[Dict],
    modern_records: List[Dict],
    tolerance: float = 1e-6,
    ignore_count_mismatch: bool = True,
    legacy_atoms: Optional[List[Dict]] = None,
) -> HBondListComparison:
    """
    Compare hbond_list records between legacy and modern JSON.

    Uses res_id-based matching for stable comparison when available.
    Falls back to index-based matching when res_id is not available.

    Args:
        legacy_records: List of legacy hbond_list records
        modern_records: List of modern hbond_list records
        tolerance: Tolerance for floating point comparisons
        ignore_count_mismatch: If True, don't flag num_hbonds differences as mismatches
        legacy_atoms: Optional legacy pdb_atoms records for res_id resolution

    Returns:
        HBondListComparison result
    """
    result = HBondListComparison()

    # Build residue map from legacy atoms
    residue_map = {}
    if legacy_atoms:
        residue_map = build_residue_idx_to_res_id_map(legacy_atoms)

    # Deduplicate records
    def dedup_records(records: List[Dict], use_res_id: bool = False) -> List[Dict]:
        seen = set()
        deduped = []
        for rec in records:
            if rec.get("type") != "hbond_list":
                continue

            if use_res_id:
                key = get_hbond_pair_key(rec, residue_map if not use_res_id else None)
            else:
                base_i = rec.get('base_i')
                base_j = rec.get('base_j')
                if base_i is not None and base_j is not None:
                    key = normalize_pair(base_i, base_j)
                else:
                    key = None

            if key and key not in seen:
                seen.add(key)
                deduped.append(rec)

        return deduped

    legacy_dedup = dedup_records(legacy_records, use_res_id=False)
    modern_dedup = dedup_records(modern_records, use_res_id=True)

    # Build maps - try res_id-based matching first
    legacy_map = {}
    modern_map = {}

    # For legacy, use residue_map to construct res_id keys
    for rec in legacy_dedup:
        key = get_hbond_pair_key(rec, residue_map)
        if key:
            legacy_map[key] = rec
        else:
            # Fall back to index-based key
            base_i = rec.get('base_i')
            base_j = rec.get('base_j')
            if base_i is not None and base_j is not None:
                legacy_map[('idx', base_i, base_j)] = rec

    # For modern, use res_id directly
    for rec in modern_dedup:
        key = get_hbond_pair_key(rec, None)
        if key:
            modern_map[key] = rec
        else:
            # Fall back to index-based key
            base_i = rec.get('base_i')
            base_j = rec.get('base_j')
            if base_i is not None and base_j is not None:
                modern_map[('idx', base_i, base_j)] = rec

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
                "pair_key": key,
                "base_i": leg_rec.get('base_i'),
                "base_j": leg_rec.get('base_j'),
                "legacy_record": leg_rec,
                "num_hbonds": leg_rec.get("num_hbonds", 0),
            }
        )

    # Find extra in modern
    for key in modern_keys - legacy_keys:
        mod_rec = modern_map[key]
        result.extra_in_modern.append(
            {
                "pair_key": key,
                "base_i": mod_rec.get('base_i'),
                "base_j": mod_rec.get('base_j'),
                "modern_record": mod_rec,
                "num_hbonds": mod_rec.get("num_hbonds", 0),
            }
        )

    # Compare common pairs
    matched_count = 0
    for key in common_keys:
        leg_rec = legacy_map[key]
        mod_rec = modern_map[key]

        leg_hbonds = leg_rec.get("hbonds", [])
        mod_hbonds = mod_rec.get("hbonds", [])

        pair_comp = compare_pair_hbonds(leg_hbonds, mod_hbonds, tolerance)
        result.pair_comparisons[key] = pair_comp

        # Check for mismatches
        if ignore_count_mismatch:
            has_mismatch = False
            matched_count += 1
        else:
            has_mismatch = (
                pair_comp.missing_in_modern
                or pair_comp.extra_in_modern
                or pair_comp.mismatched_hbonds
                or not pair_comp.num_hbonds_match
            )
            if not has_mismatch:
                matched_count += 1

        if has_mismatch:
            result.mismatched_pairs.append(
                {
                    "pair_key": key,
                    "base_i": leg_rec.get('base_i'),
                    "base_j": leg_rec.get('base_j'),
                    "legacy_num_hbonds": leg_rec.get("num_hbonds", 0),
                    "modern_num_hbonds": mod_rec.get("num_hbonds", 0),
                    "comparison": pair_comp,
                    "legacy_record": leg_rec,
                    "modern_record": mod_rec,
                }
            )

    result.matched_count = matched_count

    return result
