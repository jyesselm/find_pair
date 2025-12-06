"""
Stage 8: Hydrogen bond comparison (hbond_list).

Compares: num_hbonds, hbonds (list of donor-acceptor pairs)
"""

from typing import Any, Dict, List, Tuple

from .base import CompareResult, check_required_field, compare_scalar


# Pair key type
PairKey = Tuple[Tuple[str, int, str], Tuple[str, int, str]]


def _make_pair_key(rec: Dict[str, Any]) -> PairKey:
    """Create a pair key from a record."""
    return (
        (rec.get("chain_id_i", ""), rec.get("residue_seq_i", 0), rec.get("insertion_i", " ")),
        (rec.get("chain_id_j", ""), rec.get("residue_seq_j", 0), rec.get("insertion_j", " "))
    )


def compare_hbonds(
    legacy_records: List[Dict[str, Any]],
    modern_records: List[Dict[str, Any]],
    tolerance: float = 0.0
) -> CompareResult:
    """
    Compare hbond_list records (Stage 8).
    
    Args:
        legacy_records: Legacy hbond_list records
        modern_records: Modern hbond_list records
        tolerance: Unused
        
    Returns:
        (passed, errors) tuple
    """
    errors: List[str] = []
    
    leg_lookup = {_make_pair_key(r): r for r in legacy_records}
    mod_lookup = {_make_pair_key(r): r for r in modern_records}
    
    for key in leg_lookup.keys() & mod_lookup.keys():
        leg, mod = leg_lookup[key], mod_lookup[key]
        _compare_hbond_record(key, leg, mod, errors)
    
    _report_missing_pairs(leg_lookup, mod_lookup, errors)
    return len(errors) == 0, errors


def _compare_hbond_record(
    key: PairKey,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str]
) -> None:
    """Compare a single hbond_list record."""
    # num_hbonds
    leg_num, mod_num = check_required_field(leg, mod, "num_hbonds", key, errors)
    compare_scalar(leg_num, mod_num, "num_hbonds", key, errors)
    
    # hbonds list
    leg_hbonds = leg.get("hbonds")
    mod_hbonds = mod.get("hbonds")
    
    if leg_hbonds is None:
        errors.append(f"Key {key} legacy missing hbonds")
        return
    if mod_hbonds is None:
        errors.append(f"Key {key} modern missing hbonds")
        return
    
    _compare_hbond_lists(key, leg_hbonds, mod_hbonds, errors)


def _compare_hbond_lists(
    key: PairKey,
    leg_hbonds: List[Dict[str, Any]],
    mod_hbonds: List[Dict[str, Any]],
    errors: List[str]
) -> None:
    """Compare hbond lists for a pair."""
    if len(leg_hbonds) != len(mod_hbonds):
        errors.append(f"Key {key} hbond count: {len(leg_hbonds)} vs {len(mod_hbonds)}")
        return
    
    # Sort both lists by (donor_atom, acceptor_atom) for consistent comparison
    leg_sorted = sorted(leg_hbonds, key=lambda h: (h.get("donor_atom", ""), h.get("acceptor_atom", "")))
    mod_sorted = sorted(mod_hbonds, key=lambda h: (h.get("donor_atom", ""), h.get("acceptor_atom", "")))
    
    for i, (leg_hb, mod_hb) in enumerate(zip(leg_sorted, mod_sorted)):
        _compare_single_hbond(key, i, leg_hb, mod_hb, errors)


def _compare_single_hbond(
    key: PairKey,
    idx: int,
    leg: Dict[str, Any],
    mod: Dict[str, Any],
    errors: List[str]
) -> None:
    """Compare a single hydrogen bond."""
    leg_donor = leg.get("donor_atom")
    mod_donor = mod.get("donor_atom")
    if leg_donor != mod_donor:
        errors.append(f"Key {key} hbond[{idx}] donor: '{leg_donor}' vs '{mod_donor}'")
    
    leg_acc = leg.get("acceptor_atom")
    mod_acc = mod.get("acceptor_atom")
    if leg_acc != mod_acc:
        errors.append(f"Key {key} hbond[{idx}] acceptor: '{leg_acc}' vs '{mod_acc}'")


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

