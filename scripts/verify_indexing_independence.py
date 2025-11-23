#!/usr/bin/env python3
"""
Verify that differences between legacy and modern JSON are NOT due to indexing issues.

This script:
1. Matches residues by chain_id + residue_seq (not by residue_idx)
2. Compares atom sets for matched residues
3. Verifies that the same physical residue is being compared
4. Shows detailed comparison for a specific PDB

Usage:
    python3 scripts/verify_indexing_independence.py 1H4S
    python3 scripts/verify_indexing_independence.py 1H4S --verbose
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict


def load_json_file(json_path: Path) -> Dict:
    """Load JSON file."""
    with open(json_path) as f:
        return json.load(f)


def extract_residue_key(record: Dict) -> Tuple[str, int, str]:
    """
    Extract unique key for residue matching.
    Uses chain_id + residue_seq + residue_name (not residue_idx).
    
    Returns: (chain_id, residue_seq, residue_name)
    """
    chain_id = record.get('chain_id', '')
    residue_seq = record.get('residue_seq', 0)
    residue_name = record.get('residue_name', '').strip()
    return (chain_id, residue_seq, residue_name)


def get_base_frame_calc_records(data: Dict) -> List[Dict]:
    """Extract base_frame_calc records from JSON."""
    records = []
    for calc in data.get('calculations', []):
        if calc.get('type') == 'base_frame_calc':
            records.append(calc)
    return records


def normalize_atom_set(atoms: List[str]) -> Set[str]:
    """Normalize atom names for comparison."""
    return {atom.strip() for atom in atoms}


def compare_residues(
    legacy_record: Dict,
    modern_record: Dict,
    pdb_id: str
) -> Dict:
    """
    Compare two residue records and return detailed comparison.
    
    Returns dict with:
    - match: bool
    - legacy_atoms: set
    - modern_atoms: set
    - only_in_legacy: set
    - only_in_modern: set
    - rms_diff: float
    - details: dict with full record info
    """
    legacy_atoms = normalize_atom_set(legacy_record.get('matched_atoms', []))
    modern_atoms = normalize_atom_set(modern_record.get('matched_atoms', []))
    
    only_in_legacy = legacy_atoms - modern_atoms
    only_in_modern = modern_atoms - legacy_atoms
    
    match = (legacy_atoms == modern_atoms)
    
    legacy_rms = legacy_record.get('rms_fit', 0.0)
    modern_rms = modern_record.get('rms_fit', 0.0)
    rms_diff = abs(legacy_rms - modern_rms)
    
    return {
        'match': match,
        'legacy_atoms': legacy_atoms,
        'modern_atoms': modern_atoms,
        'only_in_legacy': only_in_legacy,
        'only_in_modern': only_in_modern,
        'rms_diff': rms_diff,
        'details': {
            'pdb_id': pdb_id,
            'chain_id': legacy_record.get('chain_id', ''),
            'residue_seq': legacy_record.get('residue_seq', 0),
            'residue_name': legacy_record.get('residue_name', '').strip(),
            'base_type': legacy_record.get('base_type', ''),
            'legacy_residue_idx': legacy_record.get('residue_idx', -1),
            'modern_residue_idx': modern_record.get('residue_idx', -1),
            'legacy_rms': legacy_rms,
            'modern_rms': modern_rms,
            'legacy_num_atoms': legacy_record.get('num_matched_atoms', 0),
            'modern_num_atoms': modern_record.get('num_matched_atoms', 0),
        }
    }


def verify_indexing_independence(
    pdb_id: str,
    project_root: Path,
    verbose: bool = False
) -> Dict:
    """
    Verify that differences are NOT due to indexing.
    
    Returns summary statistics and detailed comparisons.
    """
    legacy_json_path = project_root / 'data' / 'json_legacy' / f'{pdb_id}.json'
    modern_json_path = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    if not legacy_json_path.exists():
        raise FileNotFoundError(f"Legacy JSON not found: {legacy_json_path}")
    if not modern_json_path.exists():
        raise FileNotFoundError(f"Modern JSON not found: {modern_json_path}")
    
    # Load JSON files
    legacy_data = load_json_file(legacy_json_path)
    modern_data = load_json_file(modern_json_path)
    
    # Extract base_frame_calc records
    legacy_records = get_base_frame_calc_records(legacy_data)
    modern_records = get_base_frame_calc_records(modern_data)
    
    # Build index by residue key (chain_id + residue_seq + residue_name)
    legacy_by_key: Dict[Tuple[str, int, str], Dict] = {}
    for record in legacy_records:
        key = extract_residue_key(record)
        legacy_by_key[key] = record
    
    modern_by_key: Dict[Tuple[str, int, str], Dict] = {}
    for record in modern_records:
        key = extract_residue_key(record)
        modern_by_key[key] = record
    
    # Find common residues (matched by key, not by index)
    legacy_keys = set(legacy_by_key.keys())
    modern_keys = set(modern_by_key.keys())
    common_keys = legacy_keys & modern_keys
    
    # Compare common residues
    comparisons = []
    matches = 0
    differences = 0
    
    for key in sorted(common_keys):
        legacy_record = legacy_by_key[key]
        modern_record = modern_by_key[key]
        
        comparison = compare_residues(legacy_record, modern_record, pdb_id)
        comparisons.append(comparison)
        
        if comparison['match']:
            matches += 1
        else:
            differences += 1
    
    # Find residues only in legacy or modern
    only_in_legacy = legacy_keys - modern_keys
    only_in_modern = modern_keys - legacy_keys
    
    # Group differences by pattern
    difference_patterns = defaultdict(int)
    for comp in comparisons:
        if not comp['match']:
            # Create pattern key from atom differences
            legacy_only = sorted(comp['only_in_legacy'])
            modern_only = sorted(comp['only_in_modern'])
            pattern = f"Legacy+{','.join(legacy_only)} | Modern+{','.join(modern_only)}"
            difference_patterns[pattern] += 1
    
    return {
        'pdb_id': pdb_id,
        'total_legacy': len(legacy_records),
        'total_modern': len(modern_records),
        'common_residues': len(common_keys),
        'only_in_legacy': len(only_in_legacy),
        'only_in_modern': len(only_in_modern),
        'matches': matches,
        'differences': differences,
        'match_rate': matches / len(common_keys) if common_keys else 0.0,
        'comparisons': comparisons,
        'difference_patterns': dict(difference_patterns),
        'only_in_legacy_keys': sorted(only_in_legacy),
        'only_in_modern_keys': sorted(only_in_modern),
    }


def print_report(result: Dict, verbose: bool = False):
    """Print detailed comparison report."""
    print("=" * 80)
    print(f"INDEXING INDEPENDENCE VERIFICATION: {result['pdb_id']}")
    print("=" * 80)
    print()
    
    print("SUMMARY:")
    print(f"  Total Legacy Residues: {result['total_legacy']}")
    print(f"  Total Modern Residues: {result['total_modern']}")
    print(f"  Common Residues (matched by chain_id + residue_seq): {result['common_residues']}")
    print(f"  Only in Legacy: {result['only_in_legacy']}")
    print(f"  Only in Modern: {result['only_in_modern']}")
    print()
    
    print("COMPARISON RESULTS (matched by chain_id + residue_seq, NOT by residue_idx):")
    print(f"  Exact Matches: {result['matches']}")
    print(f"  Differences: {result['differences']}")
    print(f"  Match Rate: {result['match_rate']:.1%}")
    print()
    
    if result['differences'] > 0:
        print("DIFFERENCE PATTERNS:")
        sorted_patterns = sorted(
            result['difference_patterns'].items(),
            key=lambda x: x[1],
            reverse=True
        )
        for pattern, count in sorted_patterns[:10]:
            print(f"  {count:3d}x: {pattern}")
        print()
        
        if verbose:
            print("DETAILED DIFFERENCES:")
            print("-" * 80)
            
            for comp in result['comparisons']:
                if not comp['match']:
                    d = comp['details']
                    print(f"\nResidue: {d['residue_name']} {d['chain_id']}:{d['residue_seq']}")
                    print(f"  Base Type: {d['base_type']}")
                    print(f"  Legacy residue_idx: {d['legacy_residue_idx']}")
                    print(f"  Modern residue_idx: {d['modern_residue_idx']}")
                    print(f"  Index Difference: {abs(d['legacy_residue_idx'] - d['modern_residue_idx'])}")
                    print(f"  Legacy RMS: {d['legacy_rms']:.6f}")
                    print(f"  Modern RMS: {d['modern_rms']:.6f}")
                    print(f"  RMS Difference: {comp['rms_diff']:.6f}")
                    print(f"  Legacy Atoms ({d['legacy_num_atoms']}): {sorted(comp['legacy_atoms'])}")
                    print(f"  Modern Atoms ({d['modern_num_atoms']}): {sorted(comp['modern_atoms'])}")
                    if comp['only_in_legacy']:
                        print(f"  Only in Legacy: {sorted(comp['only_in_legacy'])}")
                    if comp['only_in_modern']:
                        print(f"  Only in Modern: {sorted(comp['only_in_modern'])}")
                    print()
    
    if result['only_in_legacy'] > 0:
        print(f"RESIDUES ONLY IN LEGACY ({result['only_in_legacy']}):")
        for key in result['only_in_legacy_keys'][:10]:
            print(f"  {key[2]} {key[0]}:{key[1]}")
        if len(result['only_in_legacy_keys']) > 10:
            print(f"  ... and {len(result['only_in_legacy_keys']) - 10} more")
        print()
    
    if result['only_in_modern'] > 0:
        print(f"RESIDUES ONLY IN MODERN ({result['only_in_modern']}):")
        for key in result['only_in_modern_keys'][:10]:
            print(f"  {key[2]} {key[0]}:{key[1]}")
        if len(result['only_in_modern_keys']) > 10:
            print(f"  ... and {len(result['only_in_modern_keys']) - 10} more")
        print()
    
    print("=" * 80)
    print("VERIFICATION:")
    print("=" * 80)
    print()
    print("✓ Residues are matched by chain_id + residue_seq (NOT by residue_idx)")
    print("✓ This ensures we compare the same physical residue")
    print("✓ Any differences are REAL differences, not indexing artifacts")
    print()
    
    if result['differences'] > 0:
        print("⚠ DIFFERENCES FOUND:")
        print("  These differences are NOT due to indexing because:")
        print("  1. Residues are matched by chain_id + residue_seq")
        print("  2. The same physical residue is being compared")
        print("  3. Differences are in atom sets, not residue identification")
        print()
        print("  Root causes likely:")
        print("  - Different atom matching logic (C4 inclusion, N7 for pyrimidines, etc.)")
        print("  - Different template matching")
        print("  - Different filtering criteria")
    else:
        print("✓ NO DIFFERENCES FOUND - Perfect match!")


def main():
    parser = argparse.ArgumentParser(
        description="Verify that differences are NOT due to indexing issues"
    )
    parser.add_argument(
        'pdb_id',
        help='PDB ID to analyze (e.g., 1H4S)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show detailed differences'
    )
    parser.add_argument(
        '--project-root',
        type=Path,
        default=Path(__file__).parent.parent,
        help='Project root directory'
    )
    
    args = parser.parse_args()
    
    try:
        result = verify_indexing_independence(
            args.pdb_id,
            args.project_root,
            args.verbose
        )
        print_report(result, args.verbose)
        
        # Exit with error code if differences found
        if result['differences'] > 0:
            sys.exit(1)
        else:
            sys.exit(0)
            
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

