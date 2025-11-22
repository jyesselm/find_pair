#!/usr/bin/env python3
"""
Analyze frame calculation failures from the detailed JSON report.
Focuses on actual numerical differences, not just atom list ordering.
"""

import json
import sys
from collections import defaultdict

def main():
    json_file = 'docs/frame_calculation_failures.json'
    
    try:
        with open(json_file) as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {json_file} not found. Run debug_all_failures first.")
        sys.exit(1)
    
    summary = data.get('summary', {})
    failures = data.get('failures', [])
    
    print("=== Frame Calculation Failure Analysis ===\n")
    print(f"Total failures recorded: {summary.get('total_failures', 0)}")
    print(f"Failures in JSON (first 1000): {len(failures)}\n")
    
    # Categorize failures
    atom_list_only = []  # Only atom list differs, calculations match
    real_failures = []   # Actual numerical differences
    missing_residues = []  # Residue not found or no frame
    
    for failure in failures:
        reason = failure.get('failure_reason', '')
        differences = failure.get('differences', {})
        
        # Check if this is just atom list ordering
        rot_diff = differences.get('max_rot_diff', 0)
        trans_diff = differences.get('max_trans_diff', 0)
        rms_diff = differences.get('rms_diff', 0)
        atoms_differ = differences.get('atoms_differ', False)
        
        # If calculations match but atoms differ, it's just ordering
        if (reason == 'ATOMS' and 
            rot_diff < 1e-6 and trans_diff < 1e-6 and rms_diff < 1e-6 and
            failure['our']['num_matched'] == failure['legacy']['num_matched']):
            atom_list_only.append(failure)
        elif reason in ['RESIDUE_NOT_FOUND', 'NO_FRAME_INVALID', 'NO_FRAME_NOT_STORED']:
            missing_residues.append(failure)
        else:
            real_failures.append(failure)
    
    print("=== Categorization ===")
    print(f"Atom list ordering only (calculations match): {len(atom_list_only)}")
    print(f"Missing residues (not found or no frame): {len(missing_residues)}")
    print(f"Real failures (numerical differences): {len(real_failures)}\n")
    
    # Analyze real failures
    if real_failures:
        print("=== Real Failures Analysis ===\n")
        
        # Group by failure reason
        by_reason = defaultdict(list)
        for f in real_failures:
            by_reason[f['failure_reason']].append(f)
        
        print("Failures by reason:")
        for reason, fails in sorted(by_reason.items(), key=lambda x: -len(x[1])):
            print(f"  {reason}: {len(fails)}")
        
        # Analyze num_matched differences
        num_matched_diff = defaultdict(int)
        for f in real_failures:
            our = f['our']['num_matched']
            leg = f['legacy']['num_matched']
            if our != leg:
                key = f"{our} vs {leg}"
                num_matched_diff[key] += 1
        
        if num_matched_diff:
            print("\nNum matched differences:")
            for key, count in sorted(num_matched_diff.items(), key=lambda x: -x[1])[:10]:
                print(f"  {key}: {count}")
        
        # Analyze numerical differences
        rot_diffs = [f['differences']['max_rot_diff'] for f in real_failures 
                     if 'max_rot_diff' in f['differences']]
        trans_diffs = [f['differences']['max_trans_diff'] for f in real_failures 
                       if 'max_trans_diff' in f['differences']]
        rms_diffs = [f['differences']['rms_diff'] for f in real_failures 
                     if 'rms_diff' in f['differences']]
        
        if rot_diffs:
            print(f"\nRotation differences: max={max(rot_diffs):.6f}, "
                  f"avg={sum(rot_diffs)/len(rot_diffs):.6f}, "
                  f"median={sorted(rot_diffs)[len(rot_diffs)//2]:.6f}")
        if trans_diffs:
            print(f"Translation differences: max={max(trans_diffs):.6f}, "
                  f"avg={sum(trans_diffs)/len(trans_diffs):.6f}, "
                  f"median={sorted(trans_diffs)[len(trans_diffs)//2]:.6f}")
        if rms_diffs:
            print(f"RMS differences: max={max(rms_diffs):.6f}, "
                  f"avg={sum(rms_diffs)/len(rms_diffs):.6f}, "
                  f"median={sorted(rms_diffs)[len(rms_diffs)//2]:.6f}")
        
        # Show examples of different failure types
        print("\n=== Example Failures ===")
        
        # Find examples with different atom counts
        atom_count_examples = [f for f in real_failures 
                               if f['our']['num_matched'] != f['legacy']['num_matched']][:3]
        if atom_count_examples:
            print("\n1. Different atom counts:")
            for i, f in enumerate(atom_count_examples[:3], 1):
                print(f"   {i}. {f['pdb_name']} {f['chain_id']}:{f['seq_num']} {f['residue_name']}")
                print(f"      Our: {f['our']['num_matched']} atoms, Legacy: {f['legacy']['num_matched']} atoms")
                print(f"      Our atoms: {f['our']['matched_atoms'][:5]}...")
                print(f"      Legacy atoms: {f['legacy']['matched_atoms'][:5]}...")
                print(f"      Our RMS: {f['our']['rms']:.6f}, Legacy RMS: {f['legacy']['rms']:.6f}")
                if 'max_rot_diff' in f['differences']:
                    print(f"      Rot diff: {f['differences']['max_rot_diff']:.6f}")
                if 'max_trans_diff' in f['differences']:
                    print(f"      Trans diff: {f['differences']['max_trans_diff']:.6f}")
        
        # Find examples with large numerical differences but same atom count
        large_diff_examples = [f for f in real_failures 
                               if f['our']['num_matched'] == f['legacy']['num_matched'] and
                               f['differences'].get('max_rot_diff', 0) > 0.1][:3]
        if large_diff_examples:
            print("\n2. Large numerical differences (same atom count):")
            for i, f in enumerate(large_diff_examples[:3], 1):
                print(f"   {i}. {f['pdb_name']} {f['chain_id']}:{f['seq_num']} {f['residue_name']}")
                print(f"      Both: {f['our']['num_matched']} atoms")
                print(f"      Our RMS: {f['our']['rms']:.6f}, Legacy RMS: {f['legacy']['rms']:.6f}")
                if 'max_rot_diff' in f['differences']:
                    print(f"      Rot diff: {f['differences']['max_rot_diff']:.6f}")
                if 'max_trans_diff' in f['differences']:
                    print(f"      Trans diff: {f['differences']['max_trans_diff']:.6f}")
                print(f"      Our atoms: {f['our']['matched_atoms']}")
                print(f"      Legacy atoms: {f['legacy']['matched_atoms']}")
    
    # Analyze missing residues
    if missing_residues:
        print("\n=== Missing Residues Analysis ===")
        by_reason = defaultdict(list)
        for f in missing_residues:
            by_reason[f['failure_reason']].append(f)
        
        for reason, fails in sorted(by_reason.items(), key=lambda x: -len(x[1])):
            print(f"\n{reason}: {len(fails)} residues")
            # Show unique residue types
            residue_types = set(f['residue_name'].strip() for f in fails)
            print(f"  Unique residue types: {sorted(residue_types)[:20]}")
            
            # Show examples
            print("  Examples:")
            for f in fails[:5]:
                print(f"    {f['pdb_name']} {f['chain_id']}:{f['seq_num']} {f['residue_name']}")
    
    print("\n=== Summary ===")
    print(f"Out of {len(failures)} failures analyzed:")
    print(f"  - {len(atom_list_only)} are just atom list ordering (calculations match)")
    print(f"  - {len(missing_residues)} are missing residues (expected)")
    print(f"  - {len(real_failures)} are real failures needing investigation")

if __name__ == '__main__':
    main()

