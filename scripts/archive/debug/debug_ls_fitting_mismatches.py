#!/usr/bin/env python3
"""
Debug ls_fitting count mismatches by examining each PDB one by one.
Identifies exactly what residues are missing in modern vs legacy.
"""

import json
import sys
from pathlib import Path
from collections import defaultdict


def deduplicate_legacy(records):
    """Deduplicate legacy records using chain, seq, insertion, name."""
    seen = set()
    unique = []
    for rec in records:
        key = (
            rec.get('chain_id', ''),
            rec.get('residue_seq', 0),
            rec.get('insertion', ' '),
            rec.get('residue_name', '').strip()
        )
        if key not in seen:
            seen.add(key)
            unique.append(rec)
    return unique


def compare_residue_lists(pdb_id, legacy_file, modern_file):
    """Compare residue lists between legacy and modern for a single PDB."""
    
    if not legacy_file.exists():
        return None, f"Legacy file not found: {legacy_file}"
    if not modern_file.exists():
        return None, f"Modern file not found: {modern_file}"
    
    try:
        with open(legacy_file) as f:
            legacy = json.load(f)
        with open(modern_file) as f:
            modern = json.load(f)
    except Exception as e:
        return None, f"Error loading JSON: {e}"
    
    # Deduplicate legacy
    legacy = deduplicate_legacy(legacy)
    
    # Create residue sets
    legacy_residues = set()
    modern_residues = set()
    
    for rec in legacy:
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
    
    # Find differences
    missing_in_modern = legacy_residues - modern_residues
    extra_in_modern = modern_residues - legacy_residues
    
    result = {
        'pdb_id': pdb_id,
        'legacy_count': len(legacy),
        'modern_count': len(modern),
        'missing_in_modern': sorted(list(missing_in_modern)),
        'extra_in_modern': sorted(list(extra_in_modern)),
    }
    
    return result, None


def analyze_all_mismatches(detailed_json_file, data_dir, tmp_dir):
    """Analyze all count mismatches."""
    
    with open(detailed_json_file) as f:
        results = json.load(f)
    
    count_mismatch_pdbs = results.get('count_mismatch_pdbs', [])
    
    print(f"Analyzing {len(count_mismatch_pdbs)} PDBs with count mismatches...")
    print("=" * 80)
    
    missing_bases_count = defaultdict(int)
    missing_bases_pdbs = defaultdict(list)
    
    for i, pdb_info in enumerate(count_mismatch_pdbs, 1):
        pdb_id = pdb_info['pdb_id']
        
        legacy_file = Path(data_dir) / f"json_legacy/ls_fitting/{pdb_id}.json"
        modern_file = Path(tmp_dir) / f"validation/ls_fitting/{pdb_id}.json"
        
        result, error = compare_residue_lists(pdb_id, legacy_file, modern_file)
        
        if error:
            print(f"[{i}/{len(count_mismatch_pdbs)}] {pdb_id}: ERROR - {error}")
            continue
        
        print(f"\n[{i}/{len(count_mismatch_pdbs)}] {pdb_id}")
        print(f"  Legacy: {result['legacy_count']}, Modern: {result['modern_count']}")
        
        if result['missing_in_modern']:
            print(f"  Missing in modern ({len(result['missing_in_modern'])}):")
            for chain, seq, ins, name in result['missing_in_modern']:
                ins_str = ins if ins != ' ' else ''
                print(f"    {chain} {seq}{ins_str} {name}")
                missing_bases_count[name] += 1
                missing_bases_pdbs[name].append(pdb_id)
        
        if result['extra_in_modern']:
            print(f"  Extra in modern ({len(result['extra_in_modern'])}):")
            for chain, seq, ins, name in result['extra_in_modern']:
                ins_str = ins if ins != ' ' else ''
                print(f"    {chain} {seq}{ins_str} {name}")
    
    print("\n" + "=" * 80)
    print("SUMMARY: Missing bases in modern (sorted by frequency)")
    print("=" * 80)
    
    for base, count in sorted(missing_bases_count.items(), key=lambda x: -x[1]):
        pdbs = missing_bases_pdbs[base]
        print(f"  {base:6s}: {count:3d} occurrences in {len(pdbs):2d} PDBs")
        if len(pdbs) <= 5:
            print(f"           PDBs: {', '.join(pdbs)}")
    
    return missing_bases_count, missing_bases_pdbs


def analyze_single_pdb(pdb_id, data_dir, tmp_dir):
    """Analyze a single PDB in detail."""
    
    legacy_file = Path(data_dir) / f"json_legacy/ls_fitting/{pdb_id}.json"
    modern_file = Path(tmp_dir) / f"validation/ls_fitting/{pdb_id}.json"
    
    result, error = compare_residue_lists(pdb_id, legacy_file, modern_file)
    
    if error:
        print(f"ERROR: {error}")
        return
    
    print(f"PDB: {pdb_id}")
    print(f"Legacy count: {result['legacy_count']}")
    print(f"Modern count: {result['modern_count']}")
    print(f"Difference: {result['legacy_count'] - result['modern_count']}")
    print()
    
    if result['missing_in_modern']:
        print(f"Missing in modern ({len(result['missing_in_modern'])}):")
        for chain, seq, ins, name in result['missing_in_modern']:
            ins_str = ins if ins != ' ' else ''
            print(f"  {chain} {seq}{ins_str} {name}")
        print()
    
    if result['extra_in_modern']:
        print(f"Extra in modern ({len(result['extra_in_modern'])}):")
        for chain, seq, ins, name in result['extra_in_modern']:
            ins_str = ins if ins != ' ' else ''
            print(f"  {chain} {seq}{ins_str} {name}")


def main():
    data_dir = Path("data")
    tmp_dir = Path("tmp")
    detailed_json = data_dir / "validation_results/ls_fitting_validation_detailed.json"
    
    if len(sys.argv) > 1:
        # Analyze single PDB
        pdb_id = sys.argv[1].upper()
        analyze_single_pdb(pdb_id, data_dir, tmp_dir)
    else:
        # Analyze all mismatches
        analyze_all_mismatches(detailed_json, data_dir, tmp_dir)


if __name__ == "__main__":
    main()

