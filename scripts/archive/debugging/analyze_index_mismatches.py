#!/usr/bin/env python3
"""
Analyze index mismatches between legacy and modern code.

This script:
1. Analyzes existing index mapping files for failed PDBs
2. Identifies patterns in index mismatches
3. Checks all PDBs for potential issues
4. Generates a detailed report

Usage:
    python3 scripts/analyze_index_mismatches.py
    python3 scripts/analyze_index_mismatches.py --recheck-all
"""

import sys
import json
import csv
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
from dataclasses import dataclass
import argparse


@dataclass
class MismatchAnalysis:
    """Analysis of a single PDB's index mismatch."""
    pdb_id: str
    num_modern: int
    num_legacy: int
    num_matched: int
    num_unmatched: int
    num_filtered: int
    
    # Detailed breakdowns
    modern_only: List[Dict]  # Residues in modern but not legacy
    legacy_only: List[Dict]  # Residues in legacy but not modern
    matched: List[Dict]      # Residues that match
    
    def summary(self) -> str:
        """Generate a summary string."""
        lines = [
            f"\n{'='*60}",
            f"PDB: {self.pdb_id}",
            f"{'='*60}",
            f"Modern residues:    {self.num_modern}",
            f"Legacy residues:    {self.num_legacy}",
            f"Matched:            {self.num_matched}",
            f"Modern only:        {self.num_unmatched}",
            f"Filtered:           {self.num_filtered}",
        ]
        
        if self.modern_only:
            lines.append(f"\nResidues in MODERN but not LEGACY ({len(self.modern_only)}):")
            for res in self.modern_only[:10]:  # Show first 10
                lines.append(
                    f"  {res['chain_id']}{res['residue_seq']}{res['insertion']:1s} "
                    f"{res['residue_name']:3s} (modern_idx={res['modern_index']}, "
                    f"legacy_idx={res['legacy_index']})"
                )
            if len(self.modern_only) > 10:
                lines.append(f"  ... and {len(self.modern_only) - 10} more")
        
        if self.legacy_only:
            lines.append(f"\nResidues in LEGACY but not MODERN ({len(self.legacy_only)}):")
            for res in self.legacy_only[:10]:
                lines.append(
                    f"  {res['chain_id']}{res['residue_seq']}{res['insertion']:1s} "
                    f"{res['residue_name']:3s} (legacy_idx={res['legacy_index']}, "
                    f"modern_idx={res['modern_index']})"
                )
            if len(self.legacy_only) > 10:
                lines.append(f"  ... and {len(self.legacy_only) - 10} more")
        
        if self.matched:
            lines.append(f"\nMatched residues ({len(self.matched)}):")
            for res in self.matched[:5]:
                lines.append(
                    f"  {res['chain_id']}{res['residue_seq']}{res['insertion']:1s} "
                    f"{res['residue_name']:3s} (legacy_idx={res['legacy_index']} -> "
                    f"modern_idx={res['modern_index']})"
                )
            if len(self.matched) > 5:
                lines.append(f"  ... and {len(self.matched) - 5} more")
        
        return "\n".join(lines)


def analyze_mapping_file(mapping_file: Path) -> Optional[MismatchAnalysis]:
    """Analyze a single index mapping file."""
    if not mapping_file.exists():
        return None
    
    try:
        with open(mapping_file) as f:
            residues = json.load(f)
        
        pdb_id = mapping_file.stem
        
        modern_only = []
        legacy_only = []
        matched = []
        
        for res in residues:
            legacy_idx = res.get('legacy_index', -1)
            modern_idx = res.get('modern_index', -1)
            filtered = res.get('filtered', False)
            
            if modern_idx >= 0 and legacy_idx >= 0:
                # Both have indices - matched
                matched.append(res)
            elif modern_idx >= 0 and legacy_idx < 0:
                # Modern has it, legacy doesn't
                modern_only.append(res)
            elif modern_idx < 0 and legacy_idx >= 0:
                # Legacy has it, modern doesn't
                legacy_only.append(res)
            # If both are -1, it's filtered out completely
        
        num_filtered = sum(1 for r in residues if r.get('filtered', False))
        
        return MismatchAnalysis(
            pdb_id=pdb_id,
            num_modern=sum(1 for r in residues if r.get('modern_index', -1) >= 0),
            num_legacy=sum(1 for r in residues if r.get('legacy_index', -1) >= 0),
            num_matched=len(matched),
            num_unmatched=len(modern_only),
            num_filtered=num_filtered,
            modern_only=modern_only,
            legacy_only=legacy_only,
            matched=matched
        )
    except Exception as e:
        print(f"Error analyzing {mapping_file}: {e}")
        return None


def check_pdb_for_mismatch(pdb_id: str, project_root: Path) -> Optional[MismatchAnalysis]:
    """Force regenerate mapping file for a PDB and analyze it."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    json_dir = project_root / "data" / "json"
    mapping_dir = project_root / "data" / "index_mapping"
    mapping_file = mapping_dir / f"{pdb_id}.json"
    
    if not pdb_file.exists():
        return None
    
    try:
        # Run generate_modern_json to create/update mapping
        cmd = [
            str(project_root / "build" / "generate_modern_json"),
            str(pdb_file),
            str(json_dir)
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # Check if mapping file was created (happens on failure)
        if mapping_file.exists():
            return analyze_mapping_file(mapping_file)
        else:
            # No mapping file means validation passed
            return None
            
    except Exception as e:
        print(f"Error checking {pdb_id}: {e}")
        return None


def analyze_residue_patterns(analyses: List[MismatchAnalysis]) -> Dict:
    """Analyze patterns across all mismatches."""
    patterns = {
        'total_pdbs_with_mismatches': len(analyses),
        'total_modern_only': 0,
        'total_legacy_only': 0,
        'residue_types_modern_only': defaultdict(int),
        'residue_types_legacy_only': defaultdict(int),
        'chains_modern_only': defaultdict(int),
        'chains_legacy_only': defaultdict(int),
    }
    
    for analysis in analyses:
        patterns['total_modern_only'] += len(analysis.modern_only)
        patterns['total_legacy_only'] += len(analysis.legacy_only)
        
        for res in analysis.modern_only:
            res_name = res['residue_name'].strip()
            patterns['residue_types_modern_only'][res_name] += 1
            patterns['chains_modern_only'][res['chain_id']] += 1
        
        for res in analysis.legacy_only:
            res_name = res['residue_name'].strip()
            patterns['residue_types_legacy_only'][res_name] += 1
            patterns['chains_legacy_only'][res['chain_id']] += 1
    
    return patterns


def main():
    parser = argparse.ArgumentParser(
        description='Analyze index mismatches between legacy and modern code'
    )
    parser.add_argument(
        '--recheck-all',
        action='store_true',
        help='Recheck all PDBs that had FAIL or SKIP status'
    )
    parser.add_argument(
        '--pdb',
        type=str,
        help='Check a specific PDB ID'
    )
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    mapping_dir = project_root / "data" / "index_mapping"
    status_csv = project_root / "data" / "index_validation_status.csv"
    
    analyses = []
    
    if args.pdb:
        # Check specific PDB
        print(f"Checking {args.pdb}...")
        analysis = check_pdb_for_mismatch(args.pdb, project_root)
        if analysis:
            analyses.append(analysis)
            print(analysis.summary())
        else:
            print(f"No mismatches found for {args.pdb}")
        return 0
    
    # Analyze existing mapping files
    print("Analyzing existing index mapping files...")
    if mapping_dir.exists():
        for mapping_file in sorted(mapping_dir.glob("*.json")):
            analysis = analyze_mapping_file(mapping_file)
            if analysis:
                analyses.append(analysis)
                print(f"  Found mismatch: {analysis.pdb_id}")
    
    # Optionally recheck all failed/skipped PDBs
    if args.recheck_all and status_csv.exists():
        print("\nRechecking all FAIL/SKIP PDBs...")
        with open(status_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                status = row['match_status']
                if status in ['FAIL', 'SKIP']:
                    pdb_id = row['pdb_id']
                    print(f"  Rechecking {pdb_id}...")
                    analysis = check_pdb_for_mismatch(pdb_id, project_root)
                    if analysis and analysis.pdb_id not in [a.pdb_id for a in analyses]:
                        analyses.append(analysis)
    
    # Generate report
    print(f"\n{'='*60}")
    print("INDEX MISMATCH ANALYSIS REPORT")
    print(f"{'='*60}")
    
    if not analyses:
        print("\n✅ No index mismatches found!")
        return 0
    
    # Print detailed analysis for each PDB
    for analysis in analyses:
        print(analysis.summary())
    
    # Print patterns
    if len(analyses) > 1:
        patterns = analyze_residue_patterns(analyses)
        print(f"\n{'='*60}")
        print("PATTERNS ACROSS ALL MISMATCHES")
        print(f"{'='*60}")
        print(f"Total PDBs with mismatches: {patterns['total_pdbs_with_mismatches']}")
        print(f"Total residues in modern only: {patterns['total_modern_only']}")
        print(f"Total residues in legacy only: {patterns['total_legacy_only']}")
        
        if patterns['residue_types_modern_only']:
            print("\nResidue types in MODERN only:")
            for res_type, count in sorted(
                patterns['residue_types_modern_only'].items(),
                key=lambda x: x[1],
                reverse=True
            ):
                print(f"  {res_type}: {count}")
        
        if patterns['residue_types_legacy_only']:
            print("\nResidue types in LEGACY only:")
            for res_type, count in sorted(
                patterns['residue_types_legacy_only'].items(),
                key=lambda x: x[1],
                reverse=True
            ):
                print(f"  {res_type}: {count}")
    
    # Write detailed CSV report
    report_file = project_root / "data" / "index_mismatch_report.csv"
    with open(report_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'pdb_id', 'chain_id', 'residue_seq', 'insertion', 'residue_name',
            'modern_index', 'legacy_index', 'status'
        ])
        
        for analysis in analyses:
            # Write modern only
            for res in analysis.modern_only:
                writer.writerow([
                    analysis.pdb_id,
                    res['chain_id'],
                    res['residue_seq'],
                    res['insertion'],
                    res['residue_name'].strip(),
                    res['modern_index'],
                    res['legacy_index'],
                    'MODERN_ONLY'
                ])
            
            # Write legacy only
            for res in analysis.legacy_only:
                writer.writerow([
                    analysis.pdb_id,
                    res['chain_id'],
                    res['residue_seq'],
                    res['insertion'],
                    res['residue_name'].strip(),
                    res['modern_index'],
                    res['legacy_index'],
                    'LEGACY_ONLY'
                ])
    
    print(f"\n📝 Detailed report written to: {report_file}")
    
    return 1 if analyses else 0


if __name__ == '__main__':
    sys.exit(main())

