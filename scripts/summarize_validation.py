#!/usr/bin/env python3
"""
Summarize index validation results for all PDBs.

This script:
1. Reads the index_validation_status.csv file
2. Generates summary statistics
3. Identifies any actual failures
4. Creates a detailed report

Usage:
    python3 scripts/summarize_validation.py
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict


def main():
    project_root = Path(__file__).parent.parent
    status_csv = project_root / "data" / "index_validation_status.csv"
    
    if not status_csv.exists():
        print(f"Error: {status_csv} not found")
        return 1
    
    # Read all statuses
    statuses = defaultdict(list)
    total = 0
    
    with open(status_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            total += 1
            status = row['match_status']
            statuses[status].append(row)
    
    # Print summary
    print("=" * 60)
    print("INDEX VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Total PDBs: {total}")
    print()
    
    for status in sorted(statuses.keys()):
        count = len(statuses[status])
        pct = 100.0 * count / total if total > 0 else 0
        
        emoji = "✅" if status == "PASS" else "❌" if status == "FAIL" else "⏭️"
        print(f"{emoji} {status:10s}: {count:5d} ({pct:5.1f}%)")
    
    # Show failures in detail
    if 'FAIL' in statuses and statuses['FAIL']:
        print()
        print("=" * 60)
        print("FAILURES (Index Mismatches)")
        print("=" * 60)
        
        for row in statuses['FAIL']:
            pdb_id = row['pdb_id']
            num_modern = row.get('num_modern', 0)
            num_legacy = row.get('num_legacy', 0)
            notes = row.get('notes', '')
            
            print(f"\n{pdb_id}:")
            print(f"  Modern: {num_modern}, Legacy: {num_legacy}")
            print(f"  Notes: {notes}")
    
    # Show skip reasons
    if 'SKIP' in statuses and statuses['SKIP']:
        print()
        print("=" * 60)
        print("SKIPPED PDBs")
        print("=" * 60)
        
        skip_reasons = defaultdict(list)
        for row in statuses['SKIP']:
            reason = row.get('notes', 'Unknown')
            skip_reasons[reason].append(row['pdb_id'])
        
        for reason, pdb_list in sorted(skip_reasons.items()):
            print(f"\n{reason}: {len(pdb_list)} PDBs")
            if len(pdb_list) <= 10:
                print(f"  {', '.join(pdb_list)}")
            else:
                print(f"  {', '.join(pdb_list[:10])}, ... and {len(pdb_list)-10} more")
    
    # Show statistics for passed PDBs
    if 'PASS' in statuses and statuses['PASS']:
        print()
        print("=" * 60)
        print("PASSED PDBs - Statistics")
        print("=" * 60)
        
        nucleotide_counts = []
        for row in statuses['PASS']:
            try:
                num_modern = int(row.get('num_modern', 0))
                if num_modern > 0:
                    nucleotide_counts.append(num_modern)
            except:
                pass
        
        if nucleotide_counts:
            nucleotide_counts.sort()
            print(f"Total with perfect index match: {len(nucleotide_counts)}")
            print(f"Smallest structure: {min(nucleotide_counts)} nucleotides")
            print(f"Largest structure: {max(nucleotide_counts)} nucleotides")
            print(f"Average: {sum(nucleotide_counts)/len(nucleotide_counts):.1f} nucleotides")
            print(f"Median: {nucleotide_counts[len(nucleotide_counts)//2]} nucleotides")
            
            # Show distribution
            bins = [0, 20, 50, 100, 200, 500, 1000, 10000]
            bin_counts = defaultdict(int)
            for count in nucleotide_counts:
                for i in range(len(bins)-1):
                    if bins[i] <= count < bins[i+1]:
                        bin_counts[f"{bins[i]}-{bins[i+1]}"] += 1
                        break
            
            print("\nSize distribution:")
            for bin_name in sorted(bin_counts.keys(), key=lambda x: int(x.split('-')[0])):
                count = bin_counts[bin_name]
                print(f"  {bin_name:12s}: {count:4d}")
    
    print()
    print("=" * 60)
    
    # Return status
    num_failures = len(statuses.get('FAIL', []))
    if num_failures > 0:
        print(f"\n⚠️  {num_failures} PDB(s) have index mismatches!")
        return 1
    else:
        print("\n✅ All tested PDBs have matching indices!")
        return 0


if __name__ == '__main__':
    sys.exit(main())

