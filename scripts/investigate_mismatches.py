#!/usr/bin/env python3
"""
Investigate the 69 mismatched residues and test with legacy mode.

This script:
1. Loads the mismatched residues from JSON
2. Tests each one with legacy mode to see if it fixes the issue
3. Reports which residues are fixed by legacy mode
"""

import json
import sys
from pathlib import Path
from x3dna_json_compare import JsonComparator, ComparisonResult

def investigate_mismatches(project_root: Path):
    """Investigate mismatches and test with legacy mode."""
    mismatches_file = project_root / "docs" / "MISMATCHED_RESIDUES.json"
    
    if not mismatches_file.exists():
        print(f"Error: {mismatches_file} not found. Run compare_all_json_executed.py first.")
        return
    
    with open(mismatches_file) as f:
        data = json.load(f)
    
    mismatches = data['mismatches']
    print("=" * 80)
    print("INVESTIGATING MISMATCHED RESIDUES")
    print("=" * 80)
    print(f"\nTotal mismatches: {len(mismatches)}\n")
    
    # Categorize
    missing = [m for m in mismatches if m['type'] == 'missing']
    atom_mismatches = [m for m in mismatches if m['type'] == 'mismatch']
    
    print(f"Missing residues: {len(missing)}")
    print(f"Atom mismatches: {len(atom_mismatches)}\n")
    
    # Test atom mismatches with legacy mode
    print("=" * 80)
    print("TESTING WITH LEGACY MODE")
    print("=" * 80)
    print("\nChecking if legacy mode fixes atom mismatches...\n")
    
    comparator = JsonComparator(enable_cache=False)
    legacy_mode_results = []
    
    # Group by PDB ID for efficiency
    pdb_groups = {}
    for mismatch in atom_mismatches:
        pdb_id = mismatch['pdb_id']
        if pdb_id not in pdb_groups:
            pdb_groups[pdb_id] = []
        pdb_groups[pdb_id].append(mismatch)
    
    fixed_by_legacy = []
    still_mismatched = []
    
    for pdb_id, group_mismatches in pdb_groups.items():
        # Check if legacy mode JSON exists
        legacy_file = project_root / f"data/json_legacy/{pdb_id}.json"
        modern_legacy_file = project_root / f"data/json/{pdb_id}_legacy.json"
        pdb_file = project_root / f"data/pdb/{pdb_id}.pdb"
        
        if not all([legacy_file.exists(), pdb_file.exists()]):
            continue
        
        if not modern_legacy_file.exists():
            print(f"⚠️  {pdb_id}: No legacy mode JSON found. Generate with --legacy flag.")
            continue
        
        # Compare with legacy mode
        try:
            result = comparator.compare_files(legacy_file, modern_legacy_file, pdb_file, pdb_id)
            
            if result.frame_comparison:
                fc = result.frame_comparison
                
                # Check each residue in this PDB
                for mismatch in group_mismatches:
                    chain_id, residue_seq_ins = mismatch['residue'].split(':')
                    residue_seq = int(residue_seq_ins[:-1] if residue_seq_ins[-1].isspace() else residue_seq_ins)
                    insertion = residue_seq_ins[-1] if residue_seq_ins[-1].isspace() else ' '
                    
                    key = (chain_id, residue_seq, insertion)
                    
                    # Check if this residue is still mismatched in legacy mode
                    still_mismatched_here = False
                    for mm in fc.mismatched_calculations:
                        if mm.residue_key == key:
                            leg_atoms = sorted(mm.legacy_matched_atoms)
                            mod_atoms = sorted(mm.modern_matched_atoms)
                            if leg_atoms != mod_atoms:
                                still_mismatched_here = True
                                still_mismatched.append({
                                    **mismatch,
                                    'legacy_mode_legacy_atoms': leg_atoms,
                                    'legacy_mode_modern_atoms': mod_atoms,
                                })
                            else:
                                fixed_by_legacy.append(mismatch)
                    
                    # If not in mismatches, it's fixed
                    if not still_mismatched_here:
                        fixed_by_legacy.append(mismatch)
        
        except Exception as e:
            print(f"Error processing {pdb_id}: {e}")
    
    print(f"\n✅ Fixed by legacy mode: {len(fixed_by_legacy)}")
    print(f"❌ Still mismatched: {len(still_mismatched)}\n")
    
    # Analyze patterns
    print("=" * 80)
    print("MISMATCH PATTERNS")
    print("=" * 80)
    
    # Group by atom differences
    patterns = {}
    for mismatch in atom_mismatches:
        only_legacy = tuple(mismatch.get('only_legacy', []))
        only_modern = tuple(mismatch.get('only_modern', []))
        pattern = (only_legacy, only_modern)
        
        if pattern not in patterns:
            patterns[pattern] = []
        patterns[pattern].append(mismatch)
    
    print(f"\nUnique patterns: {len(patterns)}\n")
    
    # Show top patterns
    sorted_patterns = sorted(patterns.items(), key=lambda x: len(x[1]), reverse=True)
    for i, (pattern, mismatches) in enumerate(sorted_patterns[:5], 1):
        only_legacy, only_modern = pattern
        print(f"Pattern {i}: {len(mismatches)} occurrences")
        if only_legacy:
            print(f"  Only in legacy: {list(only_legacy)}")
        if only_modern:
            print(f"  Only in modern: {list(only_modern)}")
        print()
    
    # Save report
    report_file = project_root / "docs" / "MISMATCH_INVESTIGATION.md"
    with open(report_file, 'w') as f:
        f.write("# Mismatch Investigation Report\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"- Total mismatches: {len(mismatches)}\n")
        f.write(f"- Missing residues: {len(missing)}\n")
        f.write(f"- Atom mismatches: {len(atom_mismatches)}\n")
        f.write(f"- Fixed by legacy mode: {len(fixed_by_legacy)}\n")
        f.write(f"- Still mismatched: {len(still_mismatched)}\n\n")
        f.write(f"## Pattern Analysis\n\n")
        for i, (pattern, mismatches) in enumerate(sorted_patterns[:10], 1):
            only_legacy, only_modern = pattern
            f.write(f"### Pattern {i}: {len(mismatches)} occurrences\n\n")
            if only_legacy:
                f.write(f"- Only in legacy: {list(only_legacy)}\n")
            if only_modern:
                f.write(f"- Only in modern: {list(only_modern)}\n")
            f.write(f"- Examples: {', '.join([m['pdb_id'] + ' ' + m['residue'] for m in mismatches[:5]])}\n\n")
    
    print(f"📄 Report saved to: {report_file}")
    
    # Show sample mismatches
    if still_mismatched:
        print("\n" + "=" * 80)
        print("SAMPLE STILL-MISMATCHED RESIDUES (First 10)")
        print("=" * 80)
        for mismatch in still_mismatched[:10]:
            print(f"\n{mismatch['pdb_id']} {mismatch['residue']}: {mismatch['residue_name']} ({mismatch['base_type']})")
            if mismatch.get('only_legacy'):
                print(f"  Only in legacy: {mismatch['only_legacy']}")
            if mismatch.get('only_modern'):
                print(f"  Only in modern: {mismatch['only_modern']}")
            if mismatch.get('rms_diff', 0) > 0.001:
                print(f"  RMS diff: {mismatch['rms_diff']:.6f}")

if __name__ == '__main__':
    project_root = Path(__file__).parent.parent
    investigate_mismatches(project_root)

