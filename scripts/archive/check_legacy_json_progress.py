#!/usr/bin/env python3
"""
Check progress of legacy JSON generation

Usage:
    python3 scripts/check_legacy_json_progress.py
"""

from pathlib import Path

def main():
    project_root = Path(__file__).parent.parent
    
    # Count existing legacy JSON files
    legacy_dir = project_root / "data" / "json_legacy" / "find_bestpair_selection"
    if legacy_dir.exists():
        existing = len(list(legacy_dir.glob("*.json")))
    else:
        existing = 0
    
    # Count total PDB files
    pdb_dir = project_root / "data" / "pdb"
    total_pdbs = len(list(pdb_dir.glob("*.pdb")))
    
    # Count PDBs <= 10 MB
    pdb_sizes = []
    for pdb_file in pdb_dir.glob("*.pdb"):
        size_mb = pdb_file.stat().st_size / (1024 * 1024)
        if size_mb <= 10.0:
            pdb_sizes.append((pdb_file.stem, size_mb))
    
    target_count = len(pdb_sizes)
    
    print("=" * 70)
    print("Legacy JSON Generation Progress")
    print("=" * 70)
    print(f"Total PDB files: {total_pdbs}")
    print(f"PDBs <= 10 MB: {target_count}")
    print(f"PDBs with legacy JSON: {existing}")
    print(f"Remaining: {target_count - existing}")
    
    if target_count > 0:
        progress = (existing / target_count) * 100
        print(f"Progress: {progress:.1f}%")
        print()
        
        # Size distribution
        if existing > 0:
            print("Size distribution of PDBs with legacy JSON:")
            existing_sizes = [size for pdb_id, size in pdb_sizes 
                            if (legacy_dir / f"{pdb_id}.json").exists()]
            if existing_sizes:
                print(f"  < 1 MB: {sum(1 for s in existing_sizes if s < 1)}")
                print(f"  1-5 MB: {sum(1 for s in existing_sizes if 1 <= s < 5)}")
                print(f"  5-10 MB: {sum(1 for s in existing_sizes if 5 <= s < 10)}")
    
    print()
    print("To generate more legacy JSON:")
    print("  python3 scripts/generate_legacy_json_batch.py --max-size 10 --workers 8")

if __name__ == "__main__":
    main()

