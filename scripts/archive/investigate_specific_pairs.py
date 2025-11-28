#!/usr/bin/env python3
"""
Investigate specific pairs that differ between legacy and modern.

Uses available JSON data (base_pair, hbond_list) to analyze why pairs differ.
Since pair_validation was deleted, we work with what we have.

Usage:
    python3 scripts/investigate_specific_pairs.py 3CF5 3239 3680
    python3 scripts/investigate_specific_pairs.py --pdb 3CF5 --missing
    python3 scripts/investigate_specific_pairs.py --pdb 3CF5 --all-differences
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from x3dna_json_compare.json_file_finder import find_json_file


def load_json_records(pdb_id: str, record_type: str, is_legacy: bool) -> List[Dict]:
    """Load JSON records."""
    base_dir = project_root / "data" / ("json_legacy" if is_legacy else "json")
    file_path = find_json_file(base_dir, pdb_id, record_type)
    
    if not file_path or not file_path.exists():
        return []
    
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            if isinstance(data, list):
                return data
            elif isinstance(data, dict):
                return [data]
            return []
    except Exception:
        return []


def find_pair_record(records: List[Dict], pair: Tuple[int, int]) -> Optional[Dict]:
    """Find a record for a specific pair."""
    i, j = pair
    for record in records:
        if isinstance(record, dict):
            bi = record.get('base_i')
            bj = record.get('base_j')
            # Handle both 0-based and 1-based indices
            if bi is not None and bj is not None:
                # Normalize to 1-based for comparison
                bi_1based = bi + 1 if isinstance(bi, int) and bi >= 0 else bi
                bj_1based = bj + 1 if isinstance(bj, int) and bj >= 0 else bj
                
                if ((bi_1based == i and bj_1based == j) or 
                    (bi_1based == j and bj_1based == i)):
                    return record
    return None


def analyze_pair(pdb_id: str, pair: Tuple[int, int]) -> Dict:
    """Analyze a specific pair."""
    i, j = pair
    analysis = {
        'pdb_id': pdb_id,
        'pair': pair,
        'legacy': {},
        'modern': {},
    }
    
    # Check base_pair records
    legacy_bp = load_json_records(pdb_id, 'base_pair', True)
    modern_bp = load_json_records(pdb_id, 'base_pair', False)
    
    legacy_bp_record = find_pair_record(legacy_bp, pair)
    modern_bp_record = find_pair_record(modern_bp, pair)
    
    analysis['legacy']['base_pair'] = legacy_bp_record is not None
    analysis['modern']['base_pair'] = modern_bp_record is not None
    
    if legacy_bp_record:
        analysis['legacy']['bp_type'] = legacy_bp_record.get('bp_type', 'N/A')
        analysis['legacy']['has_hbonds'] = len(legacy_bp_record.get('hbonds', [])) > 0
        analysis['legacy']['num_hbonds'] = len(legacy_bp_record.get('hbonds', []))
    
    if modern_bp_record:
        analysis['modern']['bp_type'] = modern_bp_record.get('bp_type', 'N/A')
        analysis['modern']['has_hbonds'] = len(modern_bp_record.get('hbonds', [])) > 0
        analysis['modern']['num_hbonds'] = len(modern_bp_record.get('hbonds', []))
    
    # Check hbond_list records
    legacy_hb = load_json_records(pdb_id, 'hbond_list', True)
    modern_hb = load_json_records(pdb_id, 'hbond_list', False)
    
    legacy_hb_record = find_pair_record(legacy_hb, pair)
    modern_hb_record = find_pair_record(modern_hb, pair)
    
    analysis['legacy']['hbond_list'] = legacy_hb_record is not None
    analysis['modern']['hbond_list'] = modern_hb_record is not None
    
    if legacy_hb_record:
        hbonds = legacy_hb_record.get('hbonds', [])
        analysis['legacy']['hbond_count'] = len(hbonds)
        # Count good H-bonds (distance in [2.5, 3.5])
        good_hb = [hb for hb in hbonds 
                  if isinstance(hb, dict) and 
                  'distance' in hb and 
                  2.5 <= hb['distance'] <= 3.5]
        analysis['legacy']['good_hbond_count'] = len(good_hb)
    
    if modern_hb_record:
        hbonds = modern_hb_record.get('hbonds', [])
        analysis['modern']['hbond_count'] = len(hbonds)
        # Count good H-bonds
        good_hb = [hb for hb in hbonds 
                  if isinstance(hb, dict) and 
                  'distance' in hb and 
                  2.5 <= hb['distance'] <= 3.5]
        analysis['modern']['good_hbond_count'] = len(good_hb)
    
    # Check find_bestpair_selection
    legacy_selection = load_json_records(pdb_id, 'find_bestpair_selection', True)
    modern_selection = load_json_records(pdb_id, 'find_bestpair_selection', False)
    
    legacy_selected = False
    modern_selected = False
    
    for rec in legacy_selection:
        if isinstance(rec, dict) and 'pairs' in rec:
            for p in rec.get('pairs', []):
                if isinstance(p, list) and len(p) >= 2:
                    if (min(p[0], p[1]), max(p[0], p[1])) == pair:
                        legacy_selected = True
                        break
    
    for rec in modern_selection:
        if isinstance(rec, dict) and 'pairs' in rec:
            for p in rec.get('pairs', []):
                if isinstance(p, list) and len(p) >= 2:
                    if (min(p[0], p[1]), max(p[0], p[1])) == pair:
                        modern_selected = True
                        break
    
    analysis['legacy']['selected'] = legacy_selected
    analysis['modern']['selected'] = modern_selected
    
    return analysis


def find_alternative_partners(pdb_id: str, residue: int, is_legacy: bool) -> List[Tuple[int, int]]:
    """Find all pairs involving a specific residue."""
    selection = load_json_records(pdb_id, 'find_bestpair_selection', is_legacy)
    pairs = []
    
    for rec in selection:
        if isinstance(rec, dict) and 'pairs' in rec:
            for p in rec.get('pairs', []):
                if isinstance(p, list) and len(p) >= 2:
                    pair = (min(p[0], p[1]), max(p[0], p[1]))
                    if pair[0] == residue or pair[1] == residue:
                        pairs.append(pair)
    
    return pairs


def investigate_same_residue_case(pdb_id: str, shared_residue: int, 
                                  missing_pair: Tuple[int, int],
                                  extra_pair: Tuple[int, int]) -> Dict:
    """Investigate why same residue selected different partner."""
    result = {
        'pdb_id': pdb_id,
        'shared_residue': shared_residue,
        'missing_pair': missing_pair,
        'extra_pair': extra_pair,
        'missing_analysis': analyze_pair(pdb_id, missing_pair),
        'extra_analysis': analyze_pair(pdb_id, extra_pair),
    }
    
    # Find all partners for this residue in both legacy and modern
    legacy_partners = find_alternative_partners(pdb_id, shared_residue, True)
    modern_partners = find_alternative_partners(pdb_id, shared_residue, False)
    
    result['legacy_partners'] = legacy_partners
    result['modern_partners'] = modern_partners
    
    return result


def print_analysis(analysis: Dict):
    """Print formatted analysis."""
    print("=" * 80)
    print(f"PAIR ANALYSIS: {analysis['pdb_id']} - {analysis['pair']}")
    print("=" * 80)
    print()
    
    print("LEGACY:")
    print(f"  Selected: {analysis['legacy']['selected']}")
    print(f"  Base pair record: {analysis['legacy']['base_pair']}")
    if analysis['legacy'].get('bp_type'):
        print(f"  BP type: {analysis['legacy']['bp_type']}")
    if 'num_hbonds' in analysis['legacy']:
        print(f"  H-bonds: {analysis['legacy']['num_hbonds']}")
        print(f"  Good H-bonds: {analysis['legacy'].get('good_hbond_count', 'N/A')}")
    print()
    
    print("MODERN:")
    print(f"  Selected: {analysis['modern']['selected']}")
    print(f"  Base pair record: {analysis['modern']['base_pair']}")
    if analysis['modern'].get('bp_type'):
        print(f"  BP type: {analysis['modern']['bp_type']}")
    if 'num_hbonds' in analysis['modern']:
        print(f"  H-bonds: {analysis['modern']['num_hbonds']}")
        print(f"  Good H-bonds: {analysis['modern'].get('good_hbond_count', 'N/A')}")
    print()
    
    # Comparison
    if analysis['legacy']['selected'] != analysis['modern']['selected']:
        print("⚠️  SELECTION DIFFERENCE:")
        if analysis['legacy']['selected']:
            print("  Legacy selected this pair, modern did not")
        else:
            print("  Modern selected this pair, legacy did not")
    else:
        print("✅ Both agree on selection status")
    print()


def print_same_residue_investigation(investigation: Dict):
    """Print investigation of same residue case."""
    print("=" * 80)
    print(f"SAME RESIDUE INVESTIGATION: {investigation['pdb_id']}")
    print(f"Residue {investigation['shared_residue']} selected different partners")
    print("=" * 80)
    print()
    
    print(f"Missing pair: {investigation['missing_pair']}")
    print(f"Extra pair: {investigation['extra_pair']}")
    print()
    
    print("LEGACY PARTNERS:")
    for pair in investigation['legacy_partners']:
        print(f"  {pair}")
    print()
    
    print("MODERN PARTNERS:")
    for pair in investigation['modern_partners']:
        print(f"  {pair}")
    print()
    
    print("MISSING PAIR ANALYSIS:")
    print_analysis(investigation['missing_analysis'])
    
    print("EXTRA PAIR ANALYSIS:")
    print_analysis(investigation['extra_analysis'])


def main():
    parser = argparse.ArgumentParser(
        description="Investigate specific pairs that differ"
    )
    parser.add_argument(
        'pdb_id',
        nargs='?',
        help='PDB ID to investigate'
    )
    parser.add_argument(
        'res1',
        type=int,
        nargs='?',
        help='First residue index'
    )
    parser.add_argument(
        'res2',
        type=int,
        nargs='?',
        help='Second residue index'
    )
    parser.add_argument(
        '--pdb',
        help='PDB ID (alternative to positional)'
    )
    parser.add_argument(
        '--missing',
        action='store_true',
        help='Investigate all missing pairs for PDB'
    )
    parser.add_argument(
        '--all-differences',
        action='store_true',
        help='Investigate all differences for PDB'
    )
    parser.add_argument(
        '--same-residue',
        action='store_true',
        help='Investigate same residue cases'
    )
    
    args = parser.parse_args()
    
    pdb_id = args.pdb_id or args.pdb
    if not pdb_id:
        parser.print_help()
        sys.exit(1)
    
    if args.res1 and args.res2:
        # Analyze specific pair
        pair = (min(args.res1, args.res2), max(args.res1, args.res2))
        analysis = analyze_pair(pdb_id, pair)
        print_analysis(analysis)
    
    elif args.same_residue or args.all_differences:
        # Load differences
        from scripts.analyze_find_bestpair_differences import analyze_single_pdb
        result = analyze_single_pdb(pdb_id)
        
        if result.get('status') != 'analyzed':
            print(f"Error: Could not analyze {pdb_id}")
            sys.exit(1)
        
        missing_pairs = result.get('missing_in_modern', [])
        extra_pairs = result.get('extra_in_modern', [])
        
        # Find same residue cases
        same_residue_cases = []
        for missing_pair in missing_pairs:
            i, j = missing_pair
            for extra_pair in extra_pairs:
                ei, ej = extra_pair
                if i == ei or i == ej or j == ei or j == ej:
                    shared = i if i == ei or i == ej else j
                    same_residue_cases.append({
                        'shared': shared,
                        'missing': missing_pair,
                        'extra': extra_pair
                    })
                    break
        
        if same_residue_cases:
            print(f"Found {len(same_residue_cases)} same residue cases for {pdb_id}\n")
            for case in same_residue_cases[:5]:  # Show first 5
                investigation = investigate_same_residue_case(
                    pdb_id, case['shared'], case['missing'], case['extra']
                )
                print_same_residue_investigation(investigation)
                print("\n" + "="*80 + "\n")
        else:
            print(f"No same residue cases found for {pdb_id}")
    
    elif args.missing:
        # Investigate all missing pairs
        from scripts.analyze_find_bestpair_differences import analyze_single_pdb
        result = analyze_single_pdb(pdb_id)
        
        if result.get('status') != 'analyzed':
            print(f"Error: Could not analyze {pdb_id}")
            sys.exit(1)
        
        missing_pairs = result.get('missing_in_modern', [])
        print(f"Investigating {len(missing_pairs)} missing pairs for {pdb_id}\n")
        
        for pair in missing_pairs[:10]:  # Show first 10
            analysis = analyze_pair(pdb_id, tuple(pair))
            print_analysis(analysis)
            print()
    
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

