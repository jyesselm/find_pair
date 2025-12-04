#!/usr/bin/env python3
"""
Test ls_fitting JSON with de-duplication of legacy records.

Legacy code had a bug where ls_fitting was called twice (app_fncs.c and ana_fncs.c),
resulting in duplicate records. Modern code correctly calls it once.

This test de-duplicates legacy records before comparison.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

def load_json(filepath: Path) -> List[Dict]:
    """Load JSON file."""
    with open(filepath) as f:
        return json.load(f)

def deduplicate_legacy(records: List[Dict]) -> List[Dict]:
    """
    De-duplicate legacy ls_fitting records by residue_idx.
    
    Legacy called record_ls_fitting twice, so we keep only the first occurrence.
    """
    seen = set()
    deduped = []
    
    for record in records:
        idx = record.get('residue_idx')
        if idx not in seen:
            seen.add(idx)
            deduped.append(record)
    
    return deduped

def compare_ls_fitting(pdb_id: str, modern_dir: Path, legacy_dir: Path) -> Tuple[bool, str]:
    """
    Compare modern and legacy ls_fitting JSON with de-duplication.
    
    Returns (matches, message)
    """
    modern_file = modern_dir / "ls_fitting" / f"{pdb_id}.json"
    legacy_file = legacy_dir / "ls_fitting" / f"{pdb_id}.json"
    
    if not modern_file.exists():
        return False, f"Modern JSON missing: {modern_file}"
    
    if not legacy_file.exists():
        return False, f"Legacy JSON missing: {legacy_file}"
    
    try:
        modern = load_json(modern_file)
        legacy_raw = load_json(legacy_file)
        legacy = deduplicate_legacy(legacy_raw)
        
        # Check counts
        if len(modern) != len(legacy):
            return False, f"Count mismatch after dedup: modern={len(modern)}, legacy={len(legacy)} (raw legacy={len(legacy_raw)})"
        
        # Compare each record
        mismatches = []
        for i, (m, l) in enumerate(zip(modern, legacy)):
            # Check residue_idx
            if m.get('residue_idx') != l.get('residue_idx'):
                mismatches.append(f"Record {i}: residue_idx {m.get('residue_idx')} vs {l.get('residue_idx')}")
                continue
            
            # Check num_points
            if m.get('num_points') != l.get('num_points'):
                mismatches.append(f"Record {i}: num_points {m.get('num_points')} vs {l.get('num_points')}")
            
            # Check rms_fit (with tolerance)
            m_rms = m.get('rms_fit', 0.0)
            l_rms = l.get('rms_fit', 0.0)
            if abs(m_rms - l_rms) > 1e-5:
                mismatches.append(f"Record {i}: rms_fit {m_rms} vs {l_rms} (diff={abs(m_rms - l_rms)})")
            
            # Check chain_id, residue_seq, residue_name
            for key in ['chain_id', 'residue_seq', 'residue_name']:
                if m.get(key) != l.get(key):
                    mismatches.append(f"Record {i}: {key} {m.get(key)} vs {l.get(key)}")
        
        if mismatches:
            return False, f"{len(mismatches)} field mismatches: " + "; ".join(mismatches[:3])
        
        return True, f"{len(modern)} records match (legacy had {len(legacy_raw)} with duplicates)"
        
    except Exception as e:
        return False, f"Error: {e}"

def main():
    """Run ls_fitting validation with de-duplication."""
    if len(sys.argv) < 2:
        print("Usage: test_ls_fitting_deduplicated.py PDB1 [PDB2 ...]")
        print("   or: test_ls_fitting_deduplicated.py --all")
        sys.exit(1)
    
    root_dir = Path(__file__).parent.parent
    modern_dir = root_dir / "data" / "json"
    legacy_dir = root_dir / "data" / "json_legacy"
    
    # Get list of PDBs
    if sys.argv[1] == "--all":
        # Find all PDBs with modern ls_fitting JSON
        pdbs = sorted([f.stem for f in (modern_dir / "ls_fitting").glob("*.json")])
    else:
        pdbs = sys.argv[1:]
    
    if not pdbs:
        print("No PDBs found to test")
        sys.exit(1)
    
    print(f"Testing {len(pdbs)} PDBs with de-duplicated legacy comparison\n")
    
    passed = 0
    failed = 0
    
    for pdb_id in pdbs:
        matches, message = compare_ls_fitting(pdb_id, modern_dir, legacy_dir)
        
        if matches:
            print(f"✅ {pdb_id}: PASS - {message}")
            passed += 1
        else:
            print(f"❌ {pdb_id}: FAIL - {message}")
            failed += 1
    
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Total: {len(pdbs)}")
    print(f"✅ Passed: {passed} ({100*passed/len(pdbs):.1f}%)")
    print(f"❌ Failed: {failed} ({100*failed/len(pdbs):.1f}%)")
    print(f"{'='*60}")
    
    sys.exit(0 if failed == 0 else 1)

if __name__ == "__main__":
    main()

