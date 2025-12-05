#!/usr/bin/env python3
"""
Quick scan to find all modified nucleotides with template mismatches or large RMS differences.
"""
import json
from pathlib import Path
from collections import defaultdict

def main():
    data_dir = Path("data")
    
    # Load fast PDB list
    with open(data_dir / "valid_pdbs_fast.json") as f:
        fast_pdbs = json.load(f)["valid_pdbs_with_atoms_and_frames"]
    
    print(f"Scanning {len(fast_pdbs)} PDBs for modified nucleotides...")
    print("="*70)
    
    issues = defaultdict(list)
    
    for pdb_id in fast_pdbs:
        legacy_file = data_dir / "json_legacy" / "base_frame_calc" / f"{pdb_id}.json"
        modern_file = data_dir / "json" / "base_frame_calc" / f"{pdb_id}.json"
        
        if not legacy_file.exists() or not modern_file.exists():
            continue
        
        try:
            with open(legacy_file) as f:
                legacy = json.load(f)
            with open(modern_file) as f:
                modern = json.load(f)
        except json.JSONDecodeError:
            # Skip malformed JSONs
            continue
        
        # Build lookup
        modern_by_key = {}
        for rec in modern:
            key = (rec.get('chain_id'), rec.get('residue_seq'))
            modern_by_key[key] = rec
        
        # Check each residue
        seen = set()
        for leg in legacy:
            idx = leg.get('residue_idx')
            if idx and idx in seen:
                continue
            if idx:
                seen.add(idx)
            
            res_name = leg.get('residue_name', '').strip()
            # Skip standard nucleotides
            if res_name in ['A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT']:
                continue
            
            key = (leg.get('chain_id'), leg.get('residue_seq'))
            if key not in modern_by_key:
                continue
            
            mod = modern_by_key[key]
            
            # Check for template mismatch
            leg_template = leg.get('standard_template', '').split('/')[-1]
            mod_template = mod.get('standard_template', '').split('/')[-1]
            
            leg_rms = leg.get('rms_fit', 0)
            mod_rms = mod.get('rms_fit', 0)
            rms_diff = abs(leg_rms - mod_rms)
            
            # Track if template differs OR large RMS difference
            if leg_template != mod_template or rms_diff > 0.01:
                issues[res_name].append({
                    'pdb': pdb_id,
                    'residue': f"{key[0]}{key[1]}",
                    'leg_template': leg_template,
                    'mod_template': mod_template,
                    'leg_rms': leg_rms,
                    'mod_rms': mod_rms,
                    'diff': rms_diff,
                    'template_mismatch': leg_template != mod_template
                })
    
    # Report
    print(f"\nModified nucleotides with issues:")
    print("="*70)
    
    for res_name in sorted(issues.keys()):
        instances = issues[res_name]
        template_mismatches = sum(1 for i in instances if i['template_mismatch'])
        
        print(f"\n{res_name}: {len(instances)} instances")
        if template_mismatches > 0:
            print(f"  ⚠️  {template_mismatches} template mismatches")
        
        # Show first example
        ex = instances[0]
        print(f"  Example: {ex['pdb']} {ex['residue']}")
        print(f"    Legacy: {ex['leg_template']} → RMS {ex['leg_rms']:.6f}")
        print(f"    Modern: {ex['mod_template']} → RMS {ex['mod_rms']:.6f}")
        print(f"    Diff: {ex['diff']:.2e}")
    
    print(f"\n{'='*70}")
    print(f"Total modified nucleotide types with issues: {len(issues)}")

if __name__ == "__main__":
    main()

