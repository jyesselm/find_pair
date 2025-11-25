#!/usr/bin/env python3
"""
Fix hbond_idx numbering in modern JSON files to match legacy (1-based per-pair).

Legacy uses 1-based per-pair indexing where each pair's H-bonds are numbered
starting from 1 (1, 2, 3, ...). Modern currently uses global indexing.

This script fixes existing JSON files. The proper fix should be in C++ code
where hbond_list records are written (see docs/HBOND_IDX_FIX_PLAN.md).
"""

import json
import sys
from pathlib import Path


def fix_hbond_idx_in_file(json_file: Path):
    """Fix hbond_idx numbering to use 1-based per-pair indexing."""
    if not json_file.exists():
        print(f"File not found: {json_file}")
        return False
    
    with open(json_file) as f:
        data = json.load(f)
    
    if not isinstance(data, list):
        print(f"Expected list format, got {type(data)}")
        return False
    
    fixed_count = 0
    for rec in data:
        if rec.get('type') != 'hbond_list':
            continue
        
        hbonds = rec.get('hbonds', [])
        if not hbonds:
            continue
        
        # Reset hbond_idx to 1-based per-pair (legacy style)
        for i, hb in enumerate(hbonds):
            old_idx = hb.get('hbond_idx', 'missing')
            new_idx = i + 1  # 1-based per-pair
            hb['hbond_idx'] = new_idx
            if old_idx != new_idx:
                fixed_count += 1
    
    if fixed_count > 0:
        # Write back
        with open(json_file, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Fixed {fixed_count} hbond_idx values in {json_file.name}")
        return True
    else:
        print(f"No changes needed in {json_file.name}")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: fix_hbond_idx_numbering.py <pdb_id> [pdb_id2 ...]")
        print("       fix_hbond_idx_numbering.py --all")
        print("\nFixes hbond_idx numbering in modern JSON files to match legacy")
        print("(1-based per-pair indexing instead of global indexing)")
        sys.exit(1)
    
    project_root = Path(__file__).parent.parent
    json_dir = project_root / "data" / "json"
    
    if sys.argv[1] == '--all':
        # Fix all PDB JSON files
        pdb_ids = []
        for json_file in json_dir.glob("*_hbond_list.json"):
            pdb_id = json_file.stem.replace("_hbond_list", "")
            pdb_ids.append(pdb_id)
    else:
        pdb_ids = sys.argv[1:]
    
    if not pdb_ids:
        print("No files to process")
        sys.exit(0)
    
    total_fixed = 0
    for pdb_id in pdb_ids:
        json_file = json_dir / f"{pdb_id}_hbond_list.json"
        if fix_hbond_idx_in_file(json_file):
            total_fixed += 1
    
    print(f"\nSummary: Fixed {total_fixed} file(s)")
    if total_fixed > 0:
        print("\nNote: This fixes existing JSON files. For permanent fix,")
        print("      update C++ code that writes hbond_list records.")
        print("      See docs/HBOND_IDX_FIX_PLAN.md for details.")


if __name__ == "__main__":
    main()

