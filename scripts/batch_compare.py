#!/usr/bin/env python3
"""Batch comparison script that summarizes results across multiple PDBs."""

import sys
import subprocess
from pathlib import Path
import json

def extract_base_pairs(inp_file, is_legacy=False):
    """Extract base pairs from .inp file."""
    if not Path(inp_file).exists():
        return []
    
    pairs = []
    with open(inp_file) as f:
        lines = f.readlines()
        # Skip header: line 1 (PDB), line 2 (outfile), line 3 (duplex), line 4 (num bp), line 5 (flags)
        for line in lines[5:]:
            if line.strip().startswith('#') or line.strip().startswith('#####'):
                continue
            parts = line.split()
            if is_legacy:
                if len(parts) >= 2:
                    try:
                        res1, res2 = int(parts[0]), int(parts[1])
                        if res1 > 10 and res2 > 10:  # Skip header artifacts
                            pairs.append((res1, res2))
                    except ValueError:
                        continue
            else:
                if len(parts) >= 3 and parts[0].isdigit():
                    try:
                        res1, res2 = int(parts[1]), int(parts[2])
                        pairs.append((res1, res2))
                    except ValueError:
                        continue
    return pairs

def normalize_pairs(pairs):
    """Normalize pairs (sort and handle reverse)."""
    normalized = []
    for p1, p2 in pairs:
        if p1 < p2:
            normalized.append((p1, p2))
        else:
            normalized.append((p2, p1))
    return sorted(set(normalized))

def compare_pdb(pdb_name, project_root):
    """Compare a single PDB and return results."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_name}.pdb"
    legacy_inp = project_root / f"{pdb_name}_legacy.inp"
    modern_inp = project_root / f"{pdb_name}_modern.inp"
    
    if not pdb_file.exists():
        return {"pdb": pdb_name, "error": "PDB file not found"}
    
    # Run find_pair if needed
    if not legacy_inp.exists():
        subprocess.run(
            f"org/build/bin/find_pair_original data/pdb/{pdb_name}.pdb {pdb_name}_legacy.inp",
            shell=True,
            cwd=project_root,
            capture_output=True
        )
    
    if not modern_inp.exists():
        subprocess.run(
            f"./build/find_pair_app data/pdb/{pdb_name}.pdb {pdb_name}_modern.inp",
            shell=True,
            cwd=project_root,
            capture_output=True
        )
    
    # Extract base pairs
    legacy_pairs = extract_base_pairs(legacy_inp, is_legacy=True) if legacy_inp.exists() else []
    modern_pairs = extract_base_pairs(modern_inp, is_legacy=False) if modern_inp.exists() else []
    
    legacy_norm = normalize_pairs(legacy_pairs)
    modern_norm = normalize_pairs(modern_pairs)
    
    # Compare
    match = set(legacy_norm) == set(modern_norm)
    only_legacy = set(legacy_norm) - set(modern_norm)
    only_modern = set(modern_norm) - set(legacy_norm)
    
    return {
        "pdb": pdb_name,
        "legacy_pairs": len(legacy_pairs),
        "modern_pairs": len(modern_pairs),
        "legacy_unique": len(legacy_norm),
        "modern_unique": len(modern_norm),
        "match": match,
        "only_legacy_count": len(only_legacy),
        "only_modern_count": len(only_modern),
        "has_modern": len(modern_pairs) > 0,
        "has_legacy": len(legacy_pairs) > 0
    }

def main():
    if len(sys.argv) < 2:
        # Get list of PDBs from data/pdb
        project_root = Path(__file__).parent.parent
        pdbs = sorted([p.stem for p in (project_root / "data" / "pdb").glob("*.pdb")])
        print(f"Found {len(pdbs)} PDB files. Testing first 20...")
        pdbs = pdbs[:20]
    else:
        pdbs = sys.argv[1:]
    
    project_root = Path(__file__).parent.parent
    
    results = []
    for pdb in pdbs:
        print(f"Processing {pdb}...", end=" ", flush=True)
        result = compare_pdb(pdb, project_root)
        results.append(result)
        status = "✅" if result.get("match") else "⚠️"
        print(f"{status} Legacy: {result.get('legacy_pairs', 0)}, Modern: {result.get('modern_pairs', 0)}")
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}\n")
    
    total = len(results)
    matches = sum(1 for r in results if r.get("match"))
    modern_only = sum(1 for r in results if r.get("only_modern_count", 0) > 0 and r.get("legacy_pairs", 0) == 0)
    legacy_only = sum(1 for r in results if r.get("only_legacy_count", 0) > 0 and r.get("modern_pairs", 0) == 0)
    both_found = sum(1 for r in results if r.get("legacy_pairs", 0) > 0 and r.get("modern_pairs", 0) > 0)
    
    print(f"Total PDBs tested: {total}")
    print(f"Perfect matches: {matches} ({matches*100/total:.1f}%)")
    print(f"Modern found pairs, legacy didn't: {modern_only}")
    print(f"Legacy found pairs, modern didn't: {legacy_only}")
    print(f"Both found pairs: {both_found}")
    
    # Show mismatches
    mismatches = [r for r in results if not r.get("match") and r.get("legacy_pairs", 0) > 0 and r.get("modern_pairs", 0) > 0]
    if mismatches:
        print(f"\nMismatches (both found pairs but differ):")
        for r in mismatches[:10]:
            print(f"  {r['pdb']}: Legacy={r['legacy_pairs']}, Modern={r['modern_pairs']}, "
                  f"diff_legacy={r['only_legacy_count']}, diff_modern={r['only_modern_count']}")
    
    # Show modern-only finds
    modern_only_list = [r for r in results if r.get("only_modern_count", 0) > 0 and r.get("legacy_pairs", 0) == 0]
    if modern_only_list:
        print(f"\nModern found pairs, legacy didn't:")
        for r in modern_only_list[:10]:
            print(f"  {r['pdb']}: Modern={r['modern_pairs']} pairs")

if __name__ == "__main__":
    main()

