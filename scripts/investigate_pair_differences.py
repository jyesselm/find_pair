#!/usr/bin/env python3
"""
Investigate PDBs with pair differences between legacy and modern.
Saves detailed information about missing/extra pairs.
"""

import json
import os
import csv
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Set
from datetime import datetime

# PDBs with differences (from verification run)
PDBS_WITH_DIFFERENCES = [
    '3F2T', '5V0O', '7O5H', '7VNV', '8BYV', '8TKK', '8VPV', '9CF3'
]


@dataclass
class PairDifference:
    """Details about a pair difference."""
    pdb_id: str
    pair: Tuple[int, int]
    diff_type: str  # 'missing' or 'extra'
    legacy_quality: float = 0.0
    modern_quality: float = 0.0
    notes: str = ""


def load_pairs(filepath: str) -> Set[Tuple[int, int]]:
    """Load pairs from JSON file."""
    if not os.path.exists(filepath):
        return set()
    
    with open(filepath) as f:
        data = json.load(f)
    
    pairs = set()
    for rec in data:
        if 'pairs' in rec:
            for pair in rec['pairs']:
                if len(pair) >= 2:
                    pairs.add((min(pair[0], pair[1]), max(pair[0], pair[1])))
    return pairs


def load_base_pair_details(filepath: str) -> Dict[Tuple[int, int], dict]:
    """Load base pair details from JSON."""
    if not os.path.exists(filepath):
        return {}
    
    with open(filepath) as f:
        data = json.load(f)
    
    details = {}
    for rec in data:
        # Try different field names
        i = rec.get('base_i', rec.get('residue_idx1', rec.get('residue_i')))
        j = rec.get('base_j', rec.get('residue_idx2', rec.get('residue_j')))
        if i is not None and j is not None:
            pair = (min(i, j), max(i, j))
            details[pair] = rec
    return details


def load_validation_details(filepath: str) -> Dict[Tuple[int, int], dict]:
    """Load validation details (quality scores) from JSON."""
    if not os.path.exists(filepath):
        return {}
    
    with open(filepath) as f:
        data = json.load(f)
    
    details = {}
    for rec in data:
        i = rec.get('residue_idx1', rec.get('residue_i', rec.get('base_i')))
        j = rec.get('residue_idx2', rec.get('residue_j', rec.get('base_j')))
        if i is not None and j is not None:
            pair = (min(i, j), max(i, j))
            # Keep the one with quality_score if available
            if pair not in details or 'quality_score' in rec:
                details[pair] = rec
    return details


def investigate_pdb(pdb_id: str, base_dir: str) -> List[PairDifference]:
    """Investigate a single PDB with differences."""
    legacy_pairs_path = os.path.join(base_dir, f'data/json_legacy/find_bestpair_selection/{pdb_id}.json')
    modern_pairs_path = os.path.join(base_dir, f'data/json/find_bestpair_selection/{pdb_id}.json')
    legacy_bp_path = os.path.join(base_dir, f'data/json_legacy/base_pair/{pdb_id}.json')
    modern_bp_path = os.path.join(base_dir, f'data/json/base_pair/{pdb_id}.json')
    legacy_val_path = os.path.join(base_dir, f'data/json_legacy/pair_validation/{pdb_id}.json')
    modern_val_path = os.path.join(base_dir, f'data/json/pair_validation/{pdb_id}.json')
    
    legacy_pairs = load_pairs(legacy_pairs_path)
    modern_pairs = load_pairs(modern_pairs_path)
    
    legacy_bp = load_base_pair_details(legacy_bp_path)
    modern_bp = load_base_pair_details(modern_bp_path)
    legacy_val = load_validation_details(legacy_val_path)
    modern_val = load_validation_details(modern_val_path)
    
    differences = []
    
    # Missing pairs (in legacy, not in modern)
    missing = legacy_pairs - modern_pairs
    for pair in sorted(missing):
        bp_info = legacy_bp.get(pair, {})
        val_info = legacy_val.get(pair, modern_val.get(pair, {}))
        diff = PairDifference(
            pdb_id=pdb_id,
            pair=pair,
            diff_type='missing',
            legacy_quality=val_info.get('quality_score', 0),
            modern_quality=modern_val.get(pair, {}).get('quality_score', 0),
            notes=f"bp_type={bp_info.get('bp_type', '?')}"
        )
        differences.append(diff)
    
    # Extra pairs (in modern, not in legacy)
    extra = modern_pairs - legacy_pairs
    for pair in sorted(extra):
        bp_info = modern_bp.get(pair, {})
        val_info = modern_val.get(pair, legacy_val.get(pair, {}))
        diff = PairDifference(
            pdb_id=pdb_id,
            pair=pair,
            diff_type='extra',
            legacy_quality=legacy_val.get(pair, {}).get('quality_score', 0),
            modern_quality=val_info.get('quality_score', 0),
            notes=f"bp_type={bp_info.get('bp_type', '?')}"
        )
        differences.append(diff)
    
    return differences


def main():
    base_dir = '/Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2'
    
    print("=" * 80)
    print("INVESTIGATING PDBs WITH PAIR DIFFERENCES")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"PDBs to investigate: {len(PDBS_WITH_DIFFERENCES)}")
    print("-" * 80)
    
    all_differences = []
    
    for pdb_id in PDBS_WITH_DIFFERENCES:
        print(f"\nInvestigating {pdb_id}...")
        diffs = investigate_pdb(pdb_id, base_dir)
        all_differences.extend(diffs)
        
        missing = [d for d in diffs if d.diff_type == 'missing']
        extra = [d for d in diffs if d.diff_type == 'extra']
        
        print(f"  Missing: {len(missing)}, Extra: {len(extra)}")
        for d in diffs:
            print(f"    {d.diff_type}: {d.pair}")
    
    # Save to CSV
    csv_path = os.path.join(base_dir, 'data/pair_differences.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pdb_id', 'pair_i', 'pair_j', 'diff_type', 'legacy_quality', 'modern_quality', 'notes'])
        for d in all_differences:
            writer.writerow([d.pdb_id, d.pair[0], d.pair[1], d.diff_type, d.legacy_quality, d.modern_quality, d.notes])
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total differences: {len(all_differences)}")
    print(f"  Missing in modern: {len([d for d in all_differences if d.diff_type == 'missing'])}")
    print(f"  Extra in modern: {len([d for d in all_differences if d.diff_type == 'extra'])}")
    print(f"\nResults saved to: {csv_path}")
    
    # Print detailed table
    print("\n" + "-" * 80)
    print("DETAILED DIFFERENCES")
    print("-" * 80)
    print(f"{'PDB':<8} {'Pair':<15} {'Type':<10} {'Legacy Q':<12} {'Modern Q':<12} {'Notes':<15}")
    print("-" * 80)
    for d in all_differences:
        print(f"{d.pdb_id:<8} ({d.pair[0]}, {d.pair[1]})  {d.diff_type:<10} {d.legacy_quality:<12.4f} {d.modern_quality:<12.4f} {d.notes}")


if __name__ == '__main__':
    main()

