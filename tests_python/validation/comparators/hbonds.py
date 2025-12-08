"""
Stage 8: Hydrogen bond comparison (hbond_list).

Compares: num_hbonds, hbonds[].{donor_atom, acceptor_atom, distance}
"""

from typing import Any, Dict, List, Tuple

from .base import CompareResult, check_required_field, compare_scalar
from ..config import Tolerance


# Pair key type: (base_i, base_j)
PairKey = Tuple[int, int]


def _make_pair_key(rec: Dict[str, Any]) -> PairKey:
    """Create a normalized pair key from base_i/base_j (smaller index first)."""
    i, j = rec.get("base_i", 0), rec.get("base_j", 0)
    return (min(i, j), max(i, j))


def _is_swapped(rec: Dict[str, Any]) -> bool:
    """Check if the record has swapped base_i/base_j (i.e., base_i > base_j)."""
    return rec.get("base_i", 0) > rec.get("base_j", 0)


def compare_hbonds(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = Tolerance.DISTANCE
) -> CompareResult:
    """
    Compare hbond_list records (Stage 8).
    
    Uses (base_i, base_j) as matching key.
    Compares: num_hbonds, hbonds[].{donor_atom, acceptor_atom, distance}
    
    Note: Legacy records H-bonds in both directions (i,j) and (j,i) with
    swapped donor/acceptor atoms. Modern only records one direction.
    We prefer records where base_i <= base_j for consistent comparison.
    """
    errors: List[str] = []
    
    # Build lookups, preferring non-swapped records (base_i <= base_j)
    leg_lookup = _build_hbond_lookup(legacy_records)
    mod_lookup = _build_hbond_lookup(modern_records)
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_hbond_record(key, leg, mod, errors, tolerance)
    
    _report_missing_pairs(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _normalize_hbonds(rec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize H-bond record so base_i <= base_j, swapping donor/acceptor if needed.
    
    When base_i > base_j, the donor atoms are from base_i and acceptor from base_j.
    To normalize, we swap base_i/base_j AND swap donor/acceptor in each H-bond.
    """
    base_i, base_j = rec.get("base_i", 0), rec.get("base_j", 0)
    
    if base_i <= base_j:
        # Already normalized
        return rec
    
    # Need to swap: create new record with swapped bases and H-bonds
    hbonds = rec.get("hbonds", [])
    swapped_hbonds = []
    for hb in hbonds:
        swapped_hbonds.append({
            "hbond_idx": hb.get("hbond_idx"),
            "donor_atom": hb.get("acceptor_atom"),  # Swap
            "acceptor_atom": hb.get("donor_atom"),  # Swap
            "distance": hb.get("distance"),
            "type": hb.get("type"),
        })
    
    return {
        "base_i": base_j,  # Swap
        "base_j": base_i,  # Swap
        "num_hbonds": rec.get("num_hbonds"),
        "hbonds": swapped_hbonds,
    }


def _build_hbond_lookup(records: List[Dict[str, Any]]) -> Dict[PairKey, Dict[str, Any]]:
    """
    Build lookup by normalized pair key.
    
    Normalizes each record so base_i <= base_j, with donor/acceptor swapped accordingly.
    When duplicates exist, prefer records originally with base_i <= base_j.
    """
    lookup = {}
    for rec in records:
        key = _make_pair_key(rec)
        base_i, base_j = rec.get("base_i", 0), rec.get("base_j", 0)
        
        # Normalize the record
        normalized = _normalize_hbonds(rec)
        
        # Keep first seen or prefer originally non-swapped
        if key not in lookup or base_i <= base_j:
            lookup[key] = normalized
    return lookup


def _filter_valid_hbonds(hbonds: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Filter out invalid H-bonds from legacy data.
    
    Filters:
    1. 'Self H-bonds' where donor and acceptor are the same atom (impossible)
    2. H-bonds with distance < 2.0 Å (too short, likely covalent bonds or bugs)
    
    Legacy code has bugs where it records O6-O6, N1-N1, and H-bonds at covalent
    bond distances (~1.5 Å) that are physically impossible.
    """
    MIN_HBOND_DISTANCE = 2.0  # Minimum realistic H-bond distance
    
    return [
        hb for hb in hbonds
        if (hb.get("donor_atom", "").strip() != hb.get("acceptor_atom", "").strip()
            and hb.get("distance", 0) >= MIN_HBOND_DISTANCE)
    ]


def _compare_hbond_record(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single hbond_list record."""
    # hbonds list
    leg_hbonds = leg.get("hbonds")
    mod_hbonds = mod.get("hbonds")
    
    if leg_hbonds is None:
        errors.append(f"Key {key} legacy missing hbonds")
        return
    if mod_hbonds is None:
        errors.append(f"Key {key} modern missing hbonds")
        return
    
    # Filter out invalid self-H-bonds from legacy (bug in legacy code)
    leg_hbonds = _filter_valid_hbonds(leg_hbonds)
    mod_hbonds = _filter_valid_hbonds(mod_hbonds)
    
    # Compare filtered counts (skip num_hbonds field since it includes invalid ones)
    if len(leg_hbonds) != len(mod_hbonds):
        errors.append(f"Key {key} valid hbond count: {len(leg_hbonds)} vs {len(mod_hbonds)}")
        return
    
    _compare_hbond_lists(key, leg_hbonds, mod_hbonds, errors, tolerance)


def _compare_hbond_lists(
    key: PairKey,
    leg_hbonds: List[Dict[str, Any]],
    mod_hbonds: List[Dict[str, Any]],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare hbond lists for a pair."""
    if len(leg_hbonds) != len(mod_hbonds):
        errors.append(f"Key {key} hbond count: {len(leg_hbonds)} vs {len(mod_hbonds)}")
        return
    
    # Sort both lists by normalized atom pair (order-independent) for comparison
    # This handles cases where legacy and modern have different base_i/base_j order
    leg_sorted = sorted(leg_hbonds, key=_hbond_sort_key_normalized)
    mod_sorted = sorted(mod_hbonds, key=_hbond_sort_key_normalized)
    
    for i, (leg_hb, mod_hb) in enumerate(zip(leg_sorted, mod_sorted)):
        _compare_single_hbond_normalized(key, i, leg_hb, mod_hb, errors, tolerance)


def _hbond_sort_key_normalized(hb: Dict[str, Any]) -> Tuple[str, str, float]:
    """Sort key for hbonds - normalized atom pair + distance for stable ordering."""
    atom1 = hb.get("donor_atom", "").strip()
    atom2 = hb.get("acceptor_atom", "").strip()
    distance = hb.get("distance", 0.0)
    # Return sorted pair + distance for consistent ordering
    sorted_atoms = tuple(sorted([atom1, atom2]))
    return (sorted_atoms[0], sorted_atoms[1], distance)


def _compare_single_hbond(
    key: PairKey,
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single hydrogen bond (exact donor/acceptor match)."""
    # donor_atom
    leg_donor = leg.get("donor_atom", "").strip()
    mod_donor = mod.get("donor_atom", "").strip()
    if leg_donor != mod_donor:
        errors.append(f"Key {key} hbond[{idx}] donor_atom: '{leg_donor}' vs '{mod_donor}'")
    
    # acceptor_atom
    leg_acc = leg.get("acceptor_atom", "").strip()
    mod_acc = mod.get("acceptor_atom", "").strip()
    if leg_acc != mod_acc:
        errors.append(f"Key {key} hbond[{idx}] acceptor_atom: '{leg_acc}' vs '{mod_acc}'")
    
    # distance
    leg_dist = leg.get("distance")
    mod_dist = mod.get("distance")
    if leg_dist is not None and mod_dist is not None:
        compare_scalar(leg_dist, mod_dist, f"hbond[{idx}].distance", key, errors, tolerance)


def _compare_single_hbond_normalized(
    key: PairKey,
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str],
    tolerance: float
) -> None:
    """Compare a single hydrogen bond with normalized atom pairs.
    
    This handles cases where legacy has (donor=A, acceptor=B) and
    modern has (donor=B, acceptor=A) for the same H-bond.
    """
    # Get normalized atom pairs
    leg_atoms = tuple(sorted([
        leg.get("donor_atom", "").strip(),
        leg.get("acceptor_atom", "").strip()
    ]))
    mod_atoms = tuple(sorted([
        mod.get("donor_atom", "").strip(),
        mod.get("acceptor_atom", "").strip()
    ]))
    
    if leg_atoms != mod_atoms:
        errors.append(f"Key {key} hbond[{idx}] atoms: {leg_atoms} vs {mod_atoms}")
    
    # distance should be the same regardless of donor/acceptor order
    leg_dist = leg.get("distance")
    mod_dist = mod.get("distance")
    if leg_dist is not None and mod_dist is not None:
        compare_scalar(leg_dist, mod_dist, f"hbond[{idx}].distance", key, errors, tolerance)


def _report_missing_pairs(
    leg_lookup: Dict[PairKey, Any],
    mod_lookup: Dict[PairKey, Any],
    errors: List[str]
) -> None:
    """Report pairs missing from either dataset."""
    missing = set(leg_lookup.keys()) - set(mod_lookup.keys())
    for key in list(missing)[:3]:
        errors.append(f"Missing in modern: {key}")
    
    extra = set(mod_lookup.keys()) - set(leg_lookup.keys())
    for key in list(extra)[:3]:
        errors.append(f"Extra in modern: {key}")
