#!/usr/bin/env python3
"""
Regenerate modern JSON files only for PDBs with known differences.
Uses the diff lists from docs/ directory.
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

def load_pdb_list_from_file(filepath):
    """Load PDB IDs from a file, skipping comments."""
    pdbs = []
    if not os.path.exists(filepath):
        return pdbs
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Extract PDB ID (first word)
            pdb_id = line.split()[0] if line.split() else None
            if pdb_id:
                pdbs.append(pdb_id)
    return pdbs

def get_problematic_pdbs():
    """Load all problematic PDBs from doc files."""
    docs_dir = Path('docs')
    
    all_pdbs = set()
    pdb_sources = defaultdict(list)
    
    # Load from all diff files
    diff_files = {
        'small': docs_dir / 'pdbs_small_diff.txt',
        'medium': docs_dir / 'pdbs_medium_diff.txt',
        'large': docs_dir / 'pdbs_large_diff.txt',
        'very_large': docs_dir / 'pdbs_very_large_diff.txt',
        'representative': docs_dir / 'representative_pdbs.txt'
    }
    
    for category, filepath in diff_files.items():
        pdbs = load_pdb_list_from_file(filepath)
        for pdb_id in pdbs:
            all_pdbs.add(pdb_id)
            pdb_sources[pdb_id].append(category)
    
    return sorted(all_pdbs), pdb_sources

def generate_json_for_pdb(pdb_id, output_dir='data/json'):
    """Generate JSON for a single PDB using the test framework."""
    pdb_file = Path('data/pdb') / f'{pdb_id}.pdb'
    output_file = Path(output_dir) / f'{pdb_id}.json'
    
    if not pdb_file.exists():
        return {'pdb_id': pdb_id, 'status': 'error', 'message': f'PDB file not found: {pdb_file}'}
    
    try:
        # Use Python to call the C++ parser via a simple script
        # For now, we'll use subprocess to call test_json_generation filtered
        # But that's not easy. Let's create a C++ wrapper or use Python directly
        
        # Actually, let's import and use the modern parser if available
        # For now, let's write a temporary C++ program or use the existing test
        
        # Check if we can import the parser directly
        # If not, we'll need to use subprocess to run a test
        
        return {
            'pdb_id': pdb_id,
            'status': 'skipped',
            'message': 'Direct generation not yet implemented - use test_json_generation'
        }
    except Exception as e:
        return {'pdb_id': pdb_id, 'status': 'error', 'message': str(e)}

def save_problematic_list(pdbs, pdb_sources, output_file='docs/problematic_pdbs.txt'):
    """Save the list of problematic PDBs with their categories."""
    with open(output_file, 'w') as f:
        f.write("# Problematic PDBs - Automatically generated list\n")
        f.write("# PDBs with differences between legacy and modern JSON\n\n")
        
        # Group by category
        categories = defaultdict(list)
        for pdb_id, sources in pdb_sources.items():
            main_category = sources[0] if sources else 'unknown'
            categories[main_category].append(pdb_id)
        
        for category, pdb_list in sorted(categories.items()):
            f.write(f"# {category.upper()} differences\n")
            for pdb_id in sorted(pdb_list):
                sources_str = ', '.join(pdb_sources[pdb_id])
                f.write(f"{pdb_id}  # from: {sources_str}\n")
            f.write("\n")
        
        f.write(f"\n# Total: {len(pdbs)} problematic PDBs\n")

def create_single_pdb_test_cpp(pdb_id):
    """Create a simple C++ test program for one PDB."""
    # This is a helper to generate a test that processes just one PDB
    # We'll use the existing test framework instead
    pass

def main():
    print("Loading problematic PDB lists from docs/...")
    problematic_pdbs, pdb_sources = get_problematic_pdbs()
    
    print(f"Found {len(problematic_pdbs)} problematic PDBs")
    print(f"Categories: {set(c for sources in pdb_sources.values() for c in sources)}")
    
    # Save the combined list
    save_problematic_list(problematic_pdbs, pdb_sources)
    print(f"\nSaved combined list to: docs/problematic_pdbs.txt")
    
    # Show breakdown
    print("\nBreakdown by category:")
    category_counts = defaultdict(int)
    for sources in pdb_sources.values():
        for cat in sources:
            category_counts[cat] += 1
    
    for cat, count in sorted(category_counts.items()):
        print(f"  {cat}: {count} PDBs")
    
    print(f"\nTo regenerate JSON for these {len(problematic_pdbs)} PDBs:")
    print("  1. Use the test_json_generation with filtering")
    print("  2. Or use: python3 scripts/regenerate_problematic_pdbs.py --run")
    
    if '--run' in sys.argv:
        print("\nRegenerating JSON files...")
        print(f"This will process {len(problematic_pdbs)} PDBs")
        # TODO: Implement actual regeneration
        print("(Regeneration not yet fully implemented - use test_json_generation)")
    else:
        print(f"\nRun with --run flag to regenerate JSON for these PDBs")

if __name__ == '__main__':
    main()

