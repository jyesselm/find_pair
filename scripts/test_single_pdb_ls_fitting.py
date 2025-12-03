#!/usr/bin/env python3
"""
Test a single PDB for ls_fitting by regenerating modern JSON and comparing.
"""

import json
import subprocess
import sys
from pathlib import Path


def test_pdb(pdb_id):
    """Test a single PDB."""
    
    pdb_id = pdb_id.upper()
    
    # Paths
    pdb_file = Path(f"data/pdb/{pdb_id}.pdb")
    legacy_json = Path(f"data/json_legacy/ls_fitting/{pdb_id}.json")
    output_dir = Path("tmp/test_single")
    modern_json = output_dir / "ls_fitting" / f"{pdb_id}.json"
    
    if not pdb_file.exists():
        print(f"ERROR: PDB file not found: {pdb_file}")
        return False
    
    if not legacy_json.exists():
        print(f"ERROR: Legacy JSON not found: {legacy_json}")
        return False
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate modern JSON
    print(f"Generating modern JSON for {pdb_id}...")
    cmd = [
        "./build/generate_modern_json",
        str(pdb_file),
        str(output_dir),
        "--stage=ls_fitting"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR generating modern JSON:")
        print(result.stderr)
        return False
    
    if not modern_json.exists():
        print(f"ERROR: Modern JSON was not created: {modern_json}")
        return False
    
    # Load both JSONs
    with open(legacy_json) as f:
        legacy = json.load(f)
    with open(modern_json) as f:
        modern = json.load(f)
    
    # Deduplicate legacy (same logic as validation script)
    seen = set()
    unique_legacy = []
    for rec in legacy:
        key = (
            rec.get('chain_id', ''),
            rec.get('residue_seq', 0),
            rec.get('insertion', ' '),
            rec.get('residue_name', '').strip()
        )
        if key not in seen:
            seen.add(key)
            unique_legacy.append(rec)
    
    # Create residue sets
    legacy_residues = set()
    modern_residues = set()
    
    for rec in unique_legacy:
        res_key = (
            rec.get('chain_id', ''),
            rec.get('residue_seq', 0),
            rec.get('insertion', ' '),
            rec.get('residue_name', '').strip()
        )
        legacy_residues.add(res_key)
    
    for rec in modern:
        res_key = (
            rec.get('chain_id', ''),
            rec.get('residue_seq', 0),
            rec.get('insertion', ' '),
            rec.get('residue_name', '').strip()
        )
        modern_residues.add(res_key)
    
    # Compare
    print(f"\nResults for {pdb_id}:")
    print(f"  Legacy count: {len(unique_legacy)}")
    print(f"  Modern count: {len(modern)}")
    
    missing_in_modern = legacy_residues - modern_residues
    extra_in_modern = modern_residues - legacy_residues
    
    if missing_in_modern:
        print(f"  Missing in modern ({len(missing_in_modern)}):")
        for chain, seq, ins, name in sorted(missing_in_modern):
            ins_str = ins if ins != ' ' else ''
            print(f"    {chain} {seq}{ins_str} {name}")
    
    if extra_in_modern:
        print(f"  Extra in modern ({len(extra_in_modern)}):")
        for chain, seq, ins, name in sorted(extra_in_modern):
            ins_str = ins if ins != ' ' else ''
            print(f"    {chain} {seq}{ins_str} {name}")
    
    if not missing_in_modern and not extra_in_modern:
        print(f"  ✅ PERFECT MATCH!")
        return True
    else:
        print(f"  ⚠️  Count mismatch")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python test_single_pdb_ls_fitting.py <PDB_ID>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    success = test_pdb(pdb_id)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

