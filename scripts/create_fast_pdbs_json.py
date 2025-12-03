#!/usr/bin/env python3
"""
Create a new JSON file with valid PDBs excluding slow ones.

Reads slow_pdbs.json and valid_pdbs.json, then creates a new file with
only the fast PDBs for efficient batch testing.
"""

import json
from pathlib import Path
from typing import Set, List

def main():
    project_root = Path(".")
    
    # Input files
    valid_pdbs_file = project_root / "data" / "valid_pdbs.json"
    slow_pdbs_file = project_root / "data" / "slow_pdbs.json"
    
    # Output file
    fast_pdbs_file = project_root / "data" / "valid_pdbs_fast.json"
    
    # Load valid PDBs
    if not valid_pdbs_file.exists():
        print(f"❌ Valid PDBs file not found: {valid_pdbs_file}")
        return 1
    
    with open(valid_pdbs_file) as f:
        valid_data = json.load(f)
    
    # Handle different possible structures
    if isinstance(valid_data, list):
        valid_pdb_ids = set(valid_data)
        output_structure = {"valid_pdbs": list(valid_data)}
    elif isinstance(valid_data, dict):
        # Use the same structure as input
        valid_pdb_ids = set(valid_data.get("valid_pdbs_with_atoms_and_frames", []))
        output_structure = valid_data.copy()
    else:
        print("❌ Invalid format in valid_pdbs.json")
        return 1
    
    print(f"Loaded {len(valid_pdb_ids)} valid PDBs")
    
    # Load slow PDBs (if file exists)
    slow_pdb_ids: Set[str] = set()
    if slow_pdbs_file.exists():
        with open(slow_pdbs_file) as f:
            slow_data = json.load(f)
        
        slow_pdbs_list = slow_data.get("slow_pdbs", [])
        slow_pdb_ids = {pdb["pdb_id"] for pdb in slow_pdbs_list}
        print(f"Found {len(slow_pdb_ids)} slow PDBs to exclude")
    else:
        print(f"⚠️  Slow PDBs file not found: {slow_pdbs_file}")
        print("   Will create file with all valid PDBs (none excluded)")
    
    # Calculate fast PDBs
    fast_pdb_ids = valid_pdb_ids - slow_pdb_ids
    
    print(f"Fast PDBs (for testing): {len(fast_pdb_ids)}")
    print(f"Excluded (slow): {len(slow_pdb_ids)}")
    print()
    
    # Create output structure
    if isinstance(output_structure, dict):
        # Update the main list
        if "valid_pdbs_with_atoms_and_frames" in output_structure:
            output_structure["valid_pdbs_with_atoms_and_frames"] = sorted(list(fast_pdb_ids))
        elif "valid_pdbs" in output_structure:
            output_structure["valid_pdbs"] = sorted(list(fast_pdb_ids))
        
        # Add metadata
        output_structure["metadata"] = {
            "source": "valid_pdbs.json",
            "excluded_slow_pdbs": len(slow_pdb_ids),
            "total_fast_pdbs": len(fast_pdb_ids),
            "total_original": len(valid_pdb_ids)
        }
        
        # Add list of excluded PDBs for reference
        if slow_pdb_ids:
            output_structure["excluded_pdbs"] = sorted(list(slow_pdb_ids))
    else:
        output_structure = {
            "valid_pdbs": sorted(list(fast_pdb_ids)),
            "metadata": {
                "source": "valid_pdbs.json",
                "excluded_slow_pdbs": len(slow_pdb_ids),
                "total_fast_pdbs": len(fast_pdb_ids),
                "total_original": len(valid_pdb_ids)
            }
        }
        if slow_pdb_ids:
            output_structure["excluded_pdbs"] = sorted(list(slow_pdb_ids))
    
    # Save output file
    fast_pdbs_file.parent.mkdir(parents=True, exist_ok=True)
    with open(fast_pdbs_file, 'w') as f:
        json.dump(output_structure, f, indent=2)
    
    print(f"✅ Created: {fast_pdbs_file}")
    print(f"   Contains {len(fast_pdb_ids)} fast PDBs for testing")
    
    if slow_pdb_ids:
        print(f"\nExcluded slow PDBs (first 20):")
        for pdb_id in sorted(list(slow_pdb_ids))[:20]:
            print(f"  - {pdb_id}")
        if len(slow_pdb_ids) > 20:
            print(f"  ... and {len(slow_pdb_ids) - 20} more")
    
    return 0

if __name__ == "__main__":
    exit(main())

