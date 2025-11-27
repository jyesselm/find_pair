#!/usr/bin/env python3
"""
Analyze differences in find_bestpair_selection between legacy and modern code.

This script identifies specific pairs that differ and investigates why they differ
by analyzing quality scores, validation results, and selection logic.

Usage:
    # Analyze all differences in test set
    python3 scripts/analyze_find_bestpair_differences.py --test-set 100

    # Analyze specific PDB
    python3 scripts/analyze_find_bestpair_differences.py 3G8T

    # Save detailed report
    python3 scripts/analyze_find_bestpair_differences.py --test-set 100 --output differences_report.md
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from x3dna_json_compare.json_file_finder import find_json_file
from x3dna_json_compare.find_bestpair_comparison import compare_find_bestpair_selections


def load_find_bestpair_selection(pdb_id: str, is_legacy: bool) -> List[Dict]:
    """Load find_bestpair_selection records for a PDB."""
    base_dir = project_root / "data" / ("json_legacy" if is_legacy else "json")
    file_path = find_json_file(base_dir, pdb_id, "find_bestpair_selection")
    
    if not file_path or not file_path.exists():
        return []
    
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            if isinstance(data, list):
                return data
            elif isinstance(data, dict) and 'pairs' in data:
                return data['pairs'] if isinstance(data['pairs'], list) else [data['pairs']]
            return []
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return []


def get_pairs_set(records: List[Dict]) -> Set[Tuple[int, int]]:
    """Extract pairs as set of (min, max) tuples."""
    pairs = set()
    for record in records:
        if isinstance(record, dict):
            if 'pairs' in record:
                # Format: {"type": "find_bestpair_selection", "pairs": [[i, j], ...]}
                for pair in record.get('pairs', []):
                    if isinstance(pair, list) and len(pair) >= 2:
                        i, j = pair[0], pair[1]
                        pairs.add((min(i, j), max(i, j)))
            elif 'base_i' in record and 'base_j' in record:
                # Format: {"type": "find_bestpair_selection", "base_i": i, "base_j": j}
                i, j = record['base_i'], record['base_j']
                pairs.add((min(i, j), max(i, j)))
    return pairs


def load_quality_scores(pdb_id: str, is_legacy: bool) -> Dict[Tuple[int, int], Dict]:
    """Load quality scores for pairs from pair_validation or base_pair records."""
    base_dir = project_root / "data" / ("json_legacy" if is_legacy else "json")
    quality_scores = {}
    
    # Try pair_validation first (has quality_score)
    validation_file = find_json_file(base_dir, pdb_id, "pair_validation")
    if validation_file and validation_file.exists():
        try:
            with open(validation_file, 'r') as f:
                records = json.load(f)
                if isinstance(records, list):
                    for record in records:
                        if isinstance(record, dict) and 'base_i' in record and 'base_j' in record:
                            i, j = record['base_i'], record['base_j']
                            pair_key = (min(i, j), max(i, j))
                            if 'calculated_values' in record:
                                quality_scores[pair_key] = record['calculated_values']
                            elif 'quality_score' in record:
                                quality_scores[pair_key] = {'quality_score': record['quality_score']}
        except Exception:
            pass
    
    # Try base_pair as fallback
    if not quality_scores:
        base_pair_file = find_json_file(base_dir, pdb_id, "base_pair")
        if base_pair_file and base_pair_file.exists():
            try:
                with open(base_pair_file, 'r') as f:
                    records = json.load(f)
                    if isinstance(records, list):
                        for record in records:
                            if isinstance(record, dict) and 'base_i' in record and 'base_j' in record:
                                i, j = record['base_i'], record['base_j']
                                pair_key = (min(i, j), max(i, j))
                                # Extract any quality-related fields
                                quality_scores[pair_key] = {k: v for k, v in record.items() 
                                                           if 'quality' in k.lower() or 'score' in k.lower()}
            except Exception:
                pass
    
    return quality_scores


def analyze_single_pdb(pdb_id: str) -> Dict:
    """Analyze differences for a single PDB."""
    legacy_records = load_find_bestpair_selection(pdb_id, True)
    modern_records = load_find_bestpair_selection(pdb_id, False)
    
    if not legacy_records and not modern_records:
        return {
            'pdb_id': pdb_id,
            'status': 'no_data',
            'error': 'No find_bestpair_selection records found'
        }
    
    legacy_pairs = get_pairs_set(legacy_records)
    modern_pairs = get_pairs_set(modern_records)
    
    missing_in_modern = legacy_pairs - modern_pairs
    extra_in_modern = modern_pairs - legacy_pairs
    common_pairs = legacy_pairs & modern_pairs
    
    result = {
        'pdb_id': pdb_id,
        'status': 'analyzed',
        'legacy_count': len(legacy_pairs),
        'modern_count': len(modern_pairs),
        'common_count': len(common_pairs),
        'missing_in_modern': sorted(missing_in_modern),
        'extra_in_modern': sorted(extra_in_modern),
        'missing_count': len(missing_in_modern),
        'extra_count': len(extra_in_modern),
    }
    
    # Try to load quality scores for analysis
    if missing_in_modern or extra_in_modern:
        legacy_quality = load_quality_scores(pdb_id, True)
        modern_quality = load_quality_scores(pdb_id, False)
        
        result['quality_analysis'] = {
            'missing_pairs_with_scores': [],
            'extra_pairs_with_scores': [],
        }
        
        for pair in missing_in_modern:
            pair_info = {
                'pair': pair,
                'legacy_quality': legacy_quality.get(pair, {}),
                'modern_quality': modern_quality.get(pair, {}),
            }
            result['quality_analysis']['missing_pairs_with_scores'].append(pair_info)
        
        for pair in extra_in_modern:
            pair_info = {
                'pair': pair,
                'legacy_quality': legacy_quality.get(pair, {}),
                'modern_quality': modern_quality.get(pair, {}),
            }
            result['quality_analysis']['extra_pairs_with_scores'].append(pair_info)
    
    return result


def analyze_test_set(test_set_size: int) -> Dict:
    """Analyze all PDBs in test set."""
    test_set_file = project_root / "data" / "test_sets" / f"test_set_{test_set_size}.json"
    
    if not test_set_file.exists():
        return {'error': f'Test set file not found: {test_set_file}'}
    
    with open(test_set_file, 'r') as f:
        test_data = json.load(f)
        pdb_ids = test_data.get('pdb_ids', [])
    
    results = {}
    all_missing = []
    all_extra = []
    pdbs_with_differences = []
    
    print(f"Analyzing {len(pdb_ids)} PDBs...")
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        if i % 10 == 0:
            print(f"  Progress: {i}/{len(pdb_ids)}")
        
        result = analyze_single_pdb(pdb_id)
        results[pdb_id] = result
        
        if result.get('status') == 'analyzed':
            if result.get('missing_count', 0) > 0 or result.get('extra_count', 0) > 0:
                pdbs_with_differences.append(pdb_id)
                all_missing.extend([(pdb_id, pair) for pair in result.get('missing_in_modern', [])])
                all_extra.extend([(pdb_id, pair) for pair in result.get('extra_in_modern', [])])
    
    return {
        'total_pdbs': len(pdb_ids),
        'pdbs_analyzed': len([r for r in results.values() if r.get('status') == 'analyzed']),
        'pdbs_with_differences': len(pdbs_with_differences),
        'total_missing': len(all_missing),
        'total_extra': len(all_extra),
        'all_missing_pairs': all_missing,
        'all_extra_pairs': all_extra,
        'results': results,
    }


def generate_report(analysis: Dict, output_file: Optional[Path] = None) -> str:
    """Generate detailed report."""
    lines = []
    
    lines.append("=" * 80)
    lines.append("FIND_BESTPAIR SELECTION DIFFERENCES ANALYSIS")
    lines.append("=" * 80)
    lines.append("")
    
    if 'error' in analysis:
        lines.append(f"Error: {analysis['error']}")
        return "\n".join(lines)
    
    lines.append("SUMMARY")
    lines.append("-" * 80)
    lines.append(f"Total PDBs: {analysis['total_pdbs']}")
    lines.append(f"PDBs Analyzed: {analysis['pdbs_analyzed']}")
    lines.append(f"PDBs with Differences: {analysis['pdbs_with_differences']}")
    lines.append(f"Total Missing Pairs: {analysis['total_missing']}")
    lines.append(f"Total Extra Pairs: {analysis['total_extra']}")
    lines.append("")
    
    # Group by PDB
    pdbs_with_diffs = analysis.get('pdbs_with_differences', [])
    if isinstance(pdbs_with_diffs, int):
        # If it's just a count, get the actual list from results
        pdbs_with_diffs = [pdb_id for pdb_id, result in analysis.get('results', {}).items()
                          if result.get('status') == 'analyzed' and 
                          (result.get('missing_count', 0) > 0 or result.get('extra_count', 0) > 0)]
    
    lines.append("=" * 80)
    lines.append("DIFFERENCES BY PDB")
    lines.append("=" * 80)
    lines.append("")
    
    for pdb_id in sorted(pdbs_with_diffs):
        result = analysis['results'][pdb_id]
        if result.get('status') != 'analyzed':
            continue
        
        lines.append(f"{pdb_id}:")
        lines.append(f"  Legacy pairs: {result['legacy_count']}")
        lines.append(f"  Modern pairs: {result['modern_count']}")
        lines.append(f"  Common pairs: {result['common_count']}")
        
        if result.get('missing_count', 0) > 0:
            lines.append(f"  Missing in modern ({result['missing_count']}): {result['missing_in_modern'][:10]}")
            if len(result['missing_in_modern']) > 10:
                lines.append(f"    ... and {len(result['missing_in_modern']) - 10} more")
        
        if result.get('extra_count', 0) > 0:
            lines.append(f"  Extra in modern ({result['extra_count']}): {result['extra_in_modern'][:10]}")
            if len(result['extra_in_modern']) > 10:
                lines.append(f"    ... and {len(result['extra_in_modern']) - 10} more")
        
        lines.append("")
    
    # List all missing pairs
    if analysis.get('all_missing_pairs'):
        lines.append("=" * 80)
        lines.append(f"ALL MISSING PAIRS ({len(analysis['all_missing_pairs'])})")
        lines.append("=" * 80)
        lines.append("")
        
        for pdb_id, pair in analysis['all_missing_pairs'][:50]:  # Show first 50
            lines.append(f"  {pdb_id}: {pair}")
        
        if len(analysis['all_missing_pairs']) > 50:
            lines.append(f"  ... and {len(analysis['all_missing_pairs']) - 50} more")
        lines.append("")
    
    # List all extra pairs
    if analysis.get('all_extra_pairs'):
        lines.append("=" * 80)
        lines.append(f"ALL EXTRA PAIRS ({len(analysis['all_extra_pairs'])})")
        lines.append("=" * 80)
        lines.append("")
        
        for pdb_id, pair in analysis['all_extra_pairs'][:50]:  # Show first 50
            lines.append(f"  {pdb_id}: {pair}")
        
        if len(analysis['all_extra_pairs']) > 50:
            lines.append(f"  ... and {len(analysis['all_extra_pairs']) - 50} more")
        lines.append("")
    
    report_text = "\n".join(lines)
    
    if output_file:
        output_file.write_text(report_text)
        print(f"Report saved to: {output_file}")
    
    return report_text


def main():
    parser = argparse.ArgumentParser(
        description="Analyze differences in find_bestpair_selection"
    )
    parser.add_argument(
        'pdb_ids',
        nargs='*',
        help='Specific PDB IDs to analyze (or use --test-set)'
    )
    parser.add_argument(
        '--test-set',
        type=int,
        help='Analyze all PDBs in test set of specified size'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        help='Save report to file'
    )
    
    args = parser.parse_args()
    
    if args.test_set:
        analysis = analyze_test_set(args.test_set)
        report = generate_report(analysis, args.output)
        if not args.output:
            print(report)
    elif args.pdb_ids:
        results = {}
        for pdb_id in args.pdb_ids:
            results[pdb_id] = analyze_single_pdb(pdb_id)
        
        # Create analysis dict
        analysis = {
            'total_pdbs': len(args.pdb_ids),
            'pdbs_analyzed': len([r for r in results.values() if r.get('status') == 'analyzed']),
            'pdbs_with_differences': len([r for r in results.values() 
                                         if r.get('missing_count', 0) > 0 or r.get('extra_count', 0) > 0]),
            'total_missing': sum(r.get('missing_count', 0) for r in results.values()),
            'total_extra': sum(r.get('extra_count', 0) for r in results.values()),
            'all_missing_pairs': [(pdb_id, pair) for pdb_id, result in results.items()
                                 for pair in result.get('missing_in_modern', [])],
            'all_extra_pairs': [(pdb_id, pair) for pdb_id, result in results.items()
                               for pair in result.get('extra_in_modern', [])],
            'results': results,
        }
        
        report = generate_report(analysis, args.output)
        if not args.output:
            print(report)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

