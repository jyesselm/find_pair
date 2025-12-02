#!/usr/bin/env python3
"""
Verify that residue indices correctly identify residues across JSON files.

This script checks that for each residue:
1. The legacy_residue_idx in base_frame_calc matches the residue identity (chain/seq)
2. The same residue appears with the same legacy_residue_idx in base_pair JSON
3. Frame data is consistent when the same residue appears in multiple places

Usage:
    python3 scripts/verify_frame_residue_consistency.py [PDB_ID]
    python3 scripts/verify_frame_residue_consistency.py --all
"""

import sys
import json
import csv
from pathlib import Path
from typing import Dict, Set, List
import argparse


def check_pdb_consistency(pdb_id: str, project_root: Path) -> dict:
    """Check frame-residue consistency for one PDB."""
    result = {
        'pdb_id': pdb_id,
        'status': 'UNKNOWN',
        'issues': [],
        'residues_checked': 0
    }
    
    # Load base_frame_calc (has frame calc info per residue)
    base_frame_file = project_root / "data" / "json" / "base_frame_calc" / f"{pdb_id}.json"
    if not base_frame_file.exists():
        result['status'] = 'SKIP'
        result['issues'].append("No modern JSON found")
        return result
    
    try:
        with open(base_frame_file) as f:
            base_frames = json.load(f)
    except Exception as e:
        result['status'] = 'ERROR'
        result['issues'].append(f"Error loading base_frame_calc: {e}")
        return result
    
    # Build map: legacy_residue_idx -> residue identity
    residue_identity = {}  # {legacy_idx: (chain, seq, name)}
    
    for rec in base_frames:
        legacy_idx = rec.get('legacy_residue_idx', rec.get('residue_idx'))
        chain = rec.get('chain_id', '?')
        seq = rec.get('residue_seq', 0)
        name = rec.get('residue_name', '???').strip()
        
        if legacy_idx:
            if legacy_idx in residue_identity:
                # Duplicate index!
                old = residue_identity[legacy_idx]
                new = (chain, seq, name)
                if old != new:
                    result['issues'].append(
                        f"Duplicate legacy_idx {legacy_idx}: {old} vs {new}"
                    )
            else:
                residue_identity[legacy_idx] = (chain, seq, name)
    
    result['residues_checked'] = len(residue_identity)
    
    # Load base_pair JSON
    base_pair_file = project_root / "data" / "json" / "base_pair" / f"{pdb_id}.json"
    if base_pair_file.exists():
        try:
            with open(base_pair_file) as f:
                base_pairs = json.load(f)
            
            # Check that residues in base_pair match their identity
            for pair in base_pairs:
                # Check base_i
                idx_i = pair.get('base_i') or pair.get('legacy_residue_idx_i')
                chain_i = pair.get('chain_id_i', '?')
                seq_i = pair.get('residue_seq_i', 0)
                
                if idx_i and idx_i in residue_identity:
                    expected = residue_identity[idx_i]
                    if (chain_i, seq_i) != (expected[0], expected[1]):
                        result['issues'].append(
                            f"Residue {idx_i} mismatch in base_pair: "
                            f"expected {expected[0]}:{expected[1]}, "
                            f"got {chain_i}:{seq_i}"
                        )
                
                # Check base_j
                idx_j = pair.get('base_j') or pair.get('legacy_residue_idx_j')
                chain_j = pair.get('chain_id_j', '?')
                seq_j = pair.get('residue_seq_j', 0)
                
                if idx_j and idx_j in residue_identity:
                    expected = residue_identity[idx_j]
                    if (chain_j, seq_j) != (expected[0], expected[1]):
                        result['issues'].append(
                            f"Residue {idx_j} mismatch in base_pair: "
                            f"expected {expected[0]}:{expected[1]}, "
                            f"got {chain_j}:{seq_j}"
                        )
                        
        except Exception as e:
            result['issues'].append(f"Error checking base_pair: {e}")
    
    # Determine status
    if result['issues']:
        result['status'] = 'FAIL'
    else:
        result['status'] = 'PASS'
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Verify frame-residue consistency'
    )
    parser.add_argument('pdb_id', nargs='?', help='PDB ID to check')
    parser.add_argument('--all', action='store_true', 
                       help='Check all validated PDBs')
    
    args = parser.parse_args()
    
    project_root = Path(__file__).parent.parent
    
    if args.pdb_id:
        # Check single PDB
        result = check_pdb_consistency(args.pdb_id, project_root)
        
        print(f"\n{'='*60}")
        print(f"PDB: {result['pdb_id']}")
        print(f"Status: {result['status']}")
        print(f"Residues checked: {result['residues_checked']}")
        
        if result['issues']:
            print(f"\nIssues found ({len(result['issues'])}):")
            for issue in result['issues']:
                print(f"  ❌ {issue}")
        else:
            print(f"\n✅ All residue indices consistent!")
        
        return 0 if result['status'] == 'PASS' else 1
    
    elif args.all:
        # Check all validated PDBs
        status_csv = project_root / "data" / "index_validation_status.csv"
        
        if not status_csv.exists():
            print("Error: index_validation_status.csv not found")
            return 1
        
        # Get all PASS PDBs
        pass_pdbs = []
        with open(status_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['match_status'] == 'PASS':
                    pass_pdbs.append(row['pdb_id'])
        
        print(f"Checking {len(pass_pdbs)} validated PDBs...")
        
        total_checked = 0
        total_pass = 0
        total_fail = 0
        total_skip = 0
        
        for i, pdb_id in enumerate(pass_pdbs):
            if (i + 1) % 100 == 0:
                print(f"  Progress: {i+1}/{len(pass_pdbs)}...")
            
            result = check_pdb_consistency(pdb_id, project_root)
            total_checked += 1
            
            if result['status'] == 'PASS':
                total_pass += 1
            elif result['status'] == 'FAIL':
                total_fail += 1
                print(f"  ❌ {pdb_id}: {len(result['issues'])} issues")
                for issue in result['issues'][:3]:  # Show first 3
                    print(f"     - {issue}")
            elif result['status'] == 'SKIP':
                total_skip += 1
        
        print(f"\n{'='*60}")
        print(f"CONSISTENCY CHECK COMPLETE")
        print(f"{'='*60}")
        print(f"Total checked: {total_checked}")
        print(f"✅ PASS: {total_pass}")
        print(f"❌ FAIL: {total_fail}")
        print(f"⏭️ SKIP: {total_skip}")
        
        if total_fail == 0:
            print(f"\n✅ All residue indices are consistent!")
            return 0
        else:
            print(f"\n❌ Found {total_fail} PDBs with consistency issues")
            return 1
    
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())

