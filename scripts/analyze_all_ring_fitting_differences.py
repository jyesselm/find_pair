#!/usr/bin/env python3
"""
Analyze ring fitting differences for all PDBs with differences.

This script processes all PDBs and generates comprehensive reports.
"""

import json
import sys
from pathlib import Path
from analyze_ring_fitting_differences import analyze_ring_fitting_differences, format_report


def get_all_pdb_ids(project_root: Path) -> List[str]:
    """Get all PDB IDs from legacy JSON directory."""
    json_dir = project_root / 'data' / 'json_legacy'
    if not json_dir.exists():
        return []
    
    pdb_ids = []
    for json_file in json_dir.glob('*.json'):
        pdb_id = json_file.stem
        pdb_ids.append(pdb_id)
    
    return sorted(pdb_ids)


def main():
    project_root = Path(__file__).parent.parent
    output_dir = project_root / 'docs' / 'ring_fitting_analyses'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get all PDB IDs
    all_pdb_ids = get_all_pdb_ids(project_root)
    
    print(f"Found {len(all_pdb_ids)} PDB files to analyze")
    print(f"Output directory: {output_dir}")
    print()
    
    results_with_differences = []
    
    for i, pdb_id in enumerate(all_pdb_ids, 1):
        print(f"[{i}/{len(all_pdb_ids)}] Processing {pdb_id}...", end=' ', flush=True)
        
        try:
            analysis = analyze_ring_fitting_differences(pdb_id, project_root)
            
            if analysis['errors']:
                print(f"ERROR: {analysis['errors'][0]}")
                continue
            
            if analysis['residue_analyses']:
                # Has differences - write detailed report
                output_file = output_dir / f'{pdb_id}_ring_fitting.txt'
                report = format_report(analysis)
                with open(output_file, 'w') as f:
                    f.write(report)
                
                results_with_differences.append({
                    'pdb_id': pdb_id,
                    'num_differences': len(analysis['residue_analyses']),
                    'residues': [r['residue_key'] for r in analysis['residue_analyses']]
                })
                print(f"✓ {len(analysis['residue_analyses'])} residues with differences")
            else:
                print("✓ No differences")
                
        except Exception as e:
            print(f"ERROR: {e}")
            continue
    
    # Generate summary report
    summary_file = output_dir / 'SUMMARY.txt'
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("RING FITTING DIFFERENCES SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total PDBs analyzed: {len(all_pdb_ids)}\n")
        f.write(f"PDBs with differences: {len(results_with_differences)}\n")
        f.write(f"PDBs with no differences: {len(all_pdb_ids) - len(results_with_differences)}\n\n")
        
        if results_with_differences:
            f.write("PDBs WITH DIFFERENCES:\n")
            f.write("-" * 80 + "\n")
            for result in results_with_differences:
                f.write(f"  {result['pdb_id']}: {result['num_differences']} residues\n")
                for residue in result['residues']:
                    chain = residue['chain_id'] if residue['chain_id'] != ' ' else '_'
                    seq = residue['residue_seq']
                    ins = residue['insertion'] if residue['insertion'] != ' ' else ''
                    f.write(f"    - {chain}:{seq}{ins}\n")
            f.write("\n")
            
            f.write("DETAILED REPORTS:\n")
            f.write("-" * 80 + "\n")
            for result in results_with_differences:
                f.write(f"  {result['pdb_id']}: {output_dir / result['pdb_id']}_ring_fitting.txt\n")
    
    print()
    print(f"Summary written to: {summary_file}")
    print(f"Detailed reports in: {output_dir}")


if __name__ == '__main__':
    main()

