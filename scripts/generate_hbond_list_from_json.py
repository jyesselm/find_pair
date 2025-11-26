#!/usr/bin/env python3
"""
Extract hbond_list records from legacy JSON files and create separate hbond_list JSON files.

This script extracts hbond_list records from the main legacy JSON files (which may contain
all records) and creates separate hbond_list/<PDB_ID>.json files in the segmented directory structure.
"""

import json
import sys
from pathlib import Path
from typing import List, Dict, Any


def extract_hbond_list_records(json_data: Dict[str, Any]) -> List[Dict]:
    """Extract hbond_list records from JSON data."""
    records = []
    
    # Check if calculations is a list
    if isinstance(json_data, list):
        # Direct array format
        for record in json_data:
            if isinstance(record, dict) and record.get('type') == 'hbond_list':
                records.append(record)
    elif isinstance(json_data, dict):
        # Check if there's a calculations field
        calculations = json_data.get('calculations', {})
        
        if isinstance(calculations, list):
            # Array format: calculations is a list
            for record in calculations:
                if isinstance(record, dict) and record.get('type') == 'hbond_list':
                    records.append(record)
        elif isinstance(calculations, dict):
            # Grouped format: calculations is a dict with type keys
            hbond_group = calculations.get('hbond_list', [])
            if isinstance(hbond_group, list):
                records.extend(hbond_group)
    
    return records


def process_pdb(pdb_id: str, project_root: Path) -> bool:
    """Extract hbond_list records for a PDB and write to segmented file."""
    legacy_dir = project_root / "data" / "json_legacy"
    
    # Try to load from main JSON file
    main_json_file = legacy_dir / f"{pdb_id}.json"
    
    if not main_json_file.exists():
        print(f"Warning: Main JSON file not found: {main_json_file}")
        return False
    
    try:
        # Load JSON (might be incomplete due to crash, so parse carefully)
        with open(main_json_file, 'r') as f:
            try:
                json_data = json.load(f)
            except json.JSONDecodeError:
                # Try to extract hbond_list records manually from file
                print(f"Warning: Main JSON file is not valid JSON (possibly incomplete)")
                return False
        
        # Extract hbond_list records
        hbond_records = extract_hbond_list_records(json_data)
        
        if not hbond_records:
            print(f"No hbond_list records found in {main_json_file}")
            # Try checking split files instead
            split_files = [
                legacy_dir / f"{pdb_id}_hbond_list.json",
                legacy_dir / "hbond_list" / f"{pdb_id}.json"
            ]
            
            for split_file in split_files:
                if split_file.exists():
                    try:
                        with open(split_file, 'r') as sf:
                            hbond_records = json.load(sf)
                            if isinstance(hbond_records, list):
                                print(f"Found {len(hbond_records)} hbond_list records in {split_file}")
                                break
                    except Exception as e:
                        continue
            
            if not hbond_records:
                return False
        
        # Create output directory
        output_dir = legacy_dir / "hbond_list"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Write to segmented file
        output_file = output_dir / f"{pdb_id}.json"
        with open(output_file, 'w') as f:
            json.dump(hbond_records, f, indent=2)
        
        print(f"✓ Created {output_file} with {len(hbond_records)} hbond_list records")
        return True
        
    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/generate_hbond_list_from_json.py <pdb_id1> [pdb_id2] ...")
        print("Example: python3 scripts/generate_hbond_list_from_json.py 3G8T 6CAQ")
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

