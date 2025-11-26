#!/usr/bin/env python3
"""
Batch generate residue ordering JSON files for all PDBs in a test set,
then compare modern vs legacy.

Usage: python3 scripts/generate_and_compare_residue_ordering_batch.py [test_set_size]
Default: 1000
"""

import json
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Dict, Tuple

def load_test_set(test_set_file: Path) -> List[str]:
    """Load PDB IDs from test set JSON file."""
    with open(test_set_file, 'r') as f:
        data = json.load(f)
    return data.get('pdb_ids', [])

def generate_modern_json(pdb_id: str, pdb_dir: Path, output_dir: Path, project_root: Path) -> bool:
    """Generate modern residue ordering JSON."""
    pdb_file = pdb_dir / f"{pdb_id}.pdb"
    output_json = output_dir / f"{pdb_id}_modern.json"
    
    if not pdb_file.exists():
        print(f"  ⚠️  PDB file not found: {pdb_file}")
        return False
    
    try:
        result = subprocess.run(
            [str(project_root / 'build' / 'generate_residue_ordering_json'), str(pdb_file), str(output_json)],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode != 0:
            print(f"  ❌ Modern JSON generation failed: {result.stderr}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"  ⏱️  Modern JSON generation timed out")
        return False
    except Exception as e:
        print(f"  ❌ Modern JSON generation error: {e}")
        return False

def generate_legacy_json(pdb_id: str, pdb_dir: Path, output_dir: Path, project_root: Path) -> bool:
    """Generate legacy residue ordering JSON."""
    pdb_file = pdb_dir / f"{pdb_id}.pdb"
    output_json = output_dir / f"{pdb_id}_legacy.json"
    
    if not pdb_file.exists():
        print(f"  ⚠️  PDB file not found: {pdb_file}")
        return False
    
    try:
        result = subprocess.run(
            [str(project_root / 'org' / 'build' / 'bin' / 'generate_residue_ordering_json'), str(pdb_file), str(output_json)],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode != 0:
            print(f"  ❌ Legacy JSON generation failed: {result.stderr}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"  ⏱️  Legacy JSON generation timed out")
        return False
    except Exception as e:
        print(f"  ❌ Legacy JSON generation error: {e}")
        return False

def compare_json(pdb_id: str, output_dir: Path, project_root: Path) -> Tuple[bool, str]:
    """Compare modern and legacy JSON files."""
    modern_json = output_dir / f"{pdb_id}_modern.json"
    legacy_json = output_dir / f"{pdb_id}_legacy.json"
    
    if not modern_json.exists() or not legacy_json.exists():
        return False, "Missing JSON files"
    
    try:
        result = subprocess.run(
            [str(project_root / 'build' / 'compare_residue_ordering'), str(modern_json), str(legacy_json)],
            capture_output=True,
            text=True,
            timeout=10
        )
        output = result.stdout + result.stderr
        return result.returncode == 0, output
    except Exception as e:
        return False, f"Comparison error: {e}"

def main():
    # Parse arguments
    test_set_size = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    
    # Setup paths
    project_root = Path(__file__).parent.parent
    test_set_file = project_root / "data" / "test_sets" / f"test_set_{test_set_size}.json"
    pdb_dir = project_root / "data" / "pdb"
    output_dir = project_root / "data" / "residue_ordering"
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load test set
    print(f"📋 Loading test set: {test_set_file}")
    pdb_ids = load_test_set(test_set_file)
    print(f"   Found {len(pdb_ids)} PDB IDs\n")
    
    # Statistics
    stats = {
        'total': len(pdb_ids),
        'modern_success': 0,
        'legacy_success': 0,
        'both_success': 0,
        'comparisons': 0,
        'matches': 0,
        'mismatches': 0,
        'errors': 0
    }
    
    # Process each PDB
    print(f"🔄 Processing {len(pdb_ids)} PDBs...\n")
    for i, pdb_id in enumerate(pdb_ids, 1):
        print(f"[{i}/{len(pdb_ids)}] {pdb_id}")
        
        # Generate modern JSON
        modern_ok = generate_modern_json(pdb_id, pdb_dir, output_dir, project_root)
        if modern_ok:
            stats['modern_success'] += 1
        
        # Generate legacy JSON
        legacy_ok = generate_legacy_json(pdb_id, pdb_dir, output_dir, project_root)
        if legacy_ok:
            stats['legacy_success'] += 1
        
        # Compare if both succeeded
        if modern_ok and legacy_ok:
            stats['both_success'] += 1
            stats['comparisons'] += 1
            match, output = compare_json(pdb_id, output_dir, project_root)
            if match:
                stats['matches'] += 1
                print(f"  ✅ Match")
            else:
                stats['mismatches'] += 1
                print(f"  ❌ Mismatch")
                # Print first few lines of comparison output
                lines = output.split('\n')[:5]
                for line in lines:
                    if line.strip():
                        print(f"     {line}")
        else:
            stats['errors'] += 1
        
        print()
    
    # Print summary
    print("=" * 60)
    print("📊 Summary")
    print("=" * 60)
    print(f"Total PDBs:           {stats['total']}")
    print(f"Modern JSON success:  {stats['modern_success']} ({100*stats['modern_success']/stats['total']:.1f}%)")
    print(f"Legacy JSON success:  {stats['legacy_success']} ({100*stats['legacy_success']/stats['total']:.1f}%)")
    print(f"Both succeeded:       {stats['both_success']} ({100*stats['both_success']/stats['total']:.1f}%)")
    print(f"Comparisons:          {stats['comparisons']}")
    print(f"Matches:              {stats['matches']} ({100*stats['matches']/stats['comparisons']:.1f}%)" if stats['comparisons'] > 0 else "Matches:              0")
    print(f"Mismatches:           {stats['mismatches']} ({100*stats['mismatches']/stats['comparisons']:.1f}%)" if stats['comparisons'] > 0 else "Mismatches:           0")
    print(f"Errors:               {stats['errors']}")
    print("=" * 60)
    
    # Save summary to JSON
    summary_file = output_dir / f"summary_{test_set_size}.json"
    with open(summary_file, 'w') as f:
        json.dump({
            'test_set_size': test_set_size,
            'statistics': stats,
            'pdb_ids': pdb_ids
        }, f, indent=2)
    print(f"\n💾 Summary saved to: {summary_file}")

if __name__ == '__main__':
    main()

