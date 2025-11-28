#!/usr/bin/env python3
"""
Convert modern hbond_list JSON files to legacy format.

The legacy format includes an `hb_info_string` field that contains:
- Number of "good" H-bonds (excluding type=' ' invalid ones)
- For each good H-bond: " atom1typeatom2 dist" (e.g., " N1-N3 2.95")

Format: "[2] N1-N3 2.95 N6-O4 2.87"
"""

import json
import sys
from pathlib import Path
from typing import List, Dict, Any


def create_hb_info_string(hbonds: List[Dict]) -> str:
    """Create hb_info_string from hbonds array (legacy format)."""
    # Count "good" H-bonds (type != ' ' and type is not empty)
    good_hbonds = []
    for hb in hbonds:
        hb_type = hb.get("type", "").strip()
        # Legacy excludes type=' ' (invalid/space)
        if hb_type and hb_type != ' ':
            good_hbonds.append(hb)
    
    # Format: "[num_good] atom1typeatom2 dist atom1typeatom2 dist ..."
    num_good = len(good_hbonds)
    
    if num_good == 0:
        return f"[{num_good}]"
    
    parts = [f"[{num_good}]"]
    
    for hb in good_hbonds:
        donor = hb.get("donor_atom", "").strip()
        acceptor = hb.get("acceptor_atom", "").strip()
        hb_type_char = hb.get("type", "-").strip()
        distance = hb.get("distance", 0.0)
        
        # Convert atom names to PDB v3 format (4 chars, no extra spaces)
        # Remove extra spaces and pad to 4 chars if needed
        donor_clean = donor[:4] if len(donor) >= 4 else donor.ljust(4)
        acceptor_clean = acceptor[:4] if len(acceptor) >= 4 else acceptor.ljust(4)
        
        # Format: " atom1typeatom2 dist" (e.g., " N1-N3 2.95")
        parts.append(f" {donor_clean}{hb_type_char}{acceptor_clean} {distance:.2f}")
    
    return "".join(parts)


def convert_to_legacy_format(modern_record: Dict) -> Dict:
    """Convert a modern hbond_list record to legacy format."""
    legacy_record = {
        "type": "hbond_list",
        "base_i": modern_record.get("base_i"),
        "base_j": modern_record.get("base_j"),
        "num_hbonds": modern_record.get("num_hbonds"),
        "hb_info_string": create_hb_info_string(modern_record.get("hbonds", [])),
        "hbonds": []
    }
    
    # Convert hbonds - ensure hbond_idx is 1-based sequential
    hbonds = modern_record.get("hbonds", [])
    for idx, hb in enumerate(hbonds, start=1):
        legacy_hb = {
            "hbond_idx": idx,  # 1-based per-pair indexing
            "donor_atom": hb.get("donor_atom", ""),
            "acceptor_atom": hb.get("acceptor_atom", ""),
            "distance": hb.get("distance", 0.0),
            "type": hb.get("type", " ")
        }
        legacy_record["hbonds"].append(legacy_hb)
    
    return legacy_record


def process_pdb(pdb_id: str, project_root: Path) -> bool:
    """Convert modern hbond_list to legacy format."""
    modern_file = project_root / "data" / "json" / "hbond_list" / f"{pdb_id}.json"
    legacy_dir = project_root / "data" / "json_legacy" / "hbond_list"
    
    if not modern_file.exists():
        print(f"Error: Modern hbond_list file not found: {modern_file}")
        return False
    
    try:
        # Load modern hbond_list records
        with open(modern_file, 'r') as f:
            modern_records = json.load(f)
        
        if not isinstance(modern_records, list):
            print(f"Error: Modern file is not an array: {modern_file}")
            return False
        
        # Convert each record to legacy format
        legacy_records = [convert_to_legacy_format(record) for record in modern_records]
        
        # Create output directory
        legacy_dir.mkdir(parents=True, exist_ok=True)
        
        # Write legacy format file
        legacy_file = legacy_dir / f"{pdb_id}.json"
        with open(legacy_file, 'w') as f:
            json.dump(legacy_records, f, indent=2)
        
        print(f"✓ Created {legacy_file} with {len(legacy_records)} hbond_list records")
        return True
        
    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/create_legacy_hbond_list.py <pdb_id1> [pdb_id2] ...")
        print("Example: python3 scripts/create_legacy_hbond_list.py 3G8T 6CAQ")
        sys.exit(1)
    
    project_root = Path(__file__).parent.parent
    pdb_ids = sys.argv[1:]
    
    success_count = 0
    for pdb_id in pdb_ids:
        if process_pdb(pdb_id, project_root):
            success_count += 1
    
    print(f"\nCompleted: {success_count}/{len(pdb_ids)} PDBs processed")
    sys.exit(0 if success_count == len(pdb_ids) else 1)


if __name__ == "__main__":
    main()

