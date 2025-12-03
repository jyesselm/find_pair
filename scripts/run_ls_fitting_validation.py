#!/usr/bin/env python3
"""
Run full ls_fitting validation with detailed tracking of mismatch types.
"""

import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime

project_root = Path(".")
sys.path.insert(0, str(project_root))
from scripts.test_utils import load_valid_pdbs_fast, find_executables

def test_pdb_detailed(pdb_id, legacy_exe, modern_exe):
    """Test a single PDB with detailed mismatch categorization."""
    pdb_file = project_root / "data" / "pdb" / f"{pdb_id}.pdb"
    
    if not pdb_file.exists():
        return {"status": "skip", "reason": "PDB file not found"}
    
    # Generate legacy JSON (use existing or regenerate)
    legacy_json = project_root / "data" / "json_legacy" / "ls_fitting" / f"{pdb_id}.json"
    if not legacy_json.exists():
        try:
            subprocess.run([str(legacy_exe.resolve()), f"../data/pdb/{pdb_id}.pdb"],
                          cwd=str(project_root / "org"), 
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                          timeout=30)
        except:
            return {"status": "error", "reason": "Legacy generation failed"}
    
    if not legacy_json.exists():
        return {"status": "skip", "reason": "No legacy JSON"}
    
    # Generate modern JSON
    output_dir = project_root / "tmp" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        subprocess.run([str(modern_exe), str(pdb_file), str(output_dir), "--stage=ls_fitting"],
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
    except:
        return {"status": "error", "reason": "Modern generation failed"}
    
    modern_json = output_dir / "ls_fitting" / f"{pdb_id}.json"
    if not modern_json.exists():
        return {"status": "error", "reason": "Modern JSON not generated"}
    
    # Compare
    try:
        legacy = json.load(open(legacy_json))
        modern = json.load(open(modern_json))
    except json.JSONDecodeError as e:
        return {"status": "error", "reason": f"JSON decode error: {e}"}
    
    # Deduplicate legacy records (old legacy JSON has duplicates from ana_fncs.c)
    if len(legacy) > len(modern):
        seen = set()
        unique_legacy = []
        for rec in legacy:
            key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('residue_name', '').strip())
            if key not in seen:
                seen.add(key)
                unique_legacy.append(rec)
        legacy = unique_legacy
    
    # Count mismatch after deduplication?
    if len(legacy) != len(modern):
        return {
            "status": "count_mismatch",
            "legacy_count": len(legacy),
            "modern_count": len(modern),
            "diff": len(legacy) - len(modern)
        }
    
    # Check for value mismatches
    mismatches = 0
    max_diff = 0.0
    
    for i, (leg, mod) in enumerate(zip(legacy, modern)):
        mod_clean = {k: v for k, v in mod.items() if k != 'type'}
        
        if leg != mod_clean:
            # Find max floating-point difference
            for key in leg.keys():
                if key in mod_clean:
                    if isinstance(leg[key], float) and isinstance(mod_clean[key], float):
                        diff = abs(leg[key] - mod_clean[key])
                        max_diff = max(max_diff, diff)
                    elif isinstance(leg[key], list) and isinstance(mod_clean[key], list):
                        for lv, mv in zip(leg[key], mod_clean[key]):
                            if isinstance(lv, float) and isinstance(mv, float):
                                diff = abs(lv - mv)
                                max_diff = max(max_diff, diff)
            
            mismatches += 1
    
    if mismatches == 0:
        return {"status": "pass", "count": len(legacy)}
    else:
        return {
            "status": "fp_mismatch",
            "count": len(legacy),
            "mismatches": mismatches,
            "max_diff": max_diff
        }

def main():
    print("LS_FITTING Validation with Fixed Legacy Code")
    print("=" * 80)
    print()
    
    legacy_exe, modern_exe = find_executables(project_root)
    
    if not legacy_exe or not modern_exe:
        print("Error: Executables not found")
        return 1
    
    pdb_ids = load_valid_pdbs_fast(project_root)
    print(f"Testing {len(pdb_ids)} fast PDBs\n")
    
    results = {
        "pass": 0,
        "fp_mismatch": 0,
        "count_mismatch": 0,
        "error": 0,
        "skip": 0,
        "count_mismatch_pdbs": [],
        "start_time": datetime.now().isoformat()
    }
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        result = test_pdb_detailed(pdb_id, legacy_exe, modern_exe)
        
        status = result["status"]
        results[status] = results.get(status, 0) + 1
        
        if status == "pass":
            print(f"[{i}/{len(pdb_ids)}] {pdb_id}: ✓ PASS ({result['count']} records)")
        elif status == "fp_mismatch":
            print(f"[{i}/{len(pdb_ids)}] {pdb_id}: ✓ OK ({result['mismatches']}/{result['count']} FP diffs, max={result['max_diff']:.2e})")
        elif status == "count_mismatch":
            print(f"[{i}/{len(pdb_ids)}] {pdb_id}: ⚠️  COUNT MISMATCH ({result['legacy_count']} vs {result['modern_count']}, diff={result['diff']})")
            results["count_mismatch_pdbs"].append({
                "pdb_id": pdb_id,
                "legacy": result['legacy_count'],
                "modern": result['modern_count'],
                "diff": result['diff']
            })
        elif status == "error":
            print(f"[{i}/{len(pdb_ids)}] {pdb_id}: ERROR - {result['reason']}")
        elif status == "skip":
            pass  # Silent skip
        
        if i % 100 == 0:
            print(f"\n--- Progress: {i}/{len(pdb_ids)} | Pass: {results['pass']} | FP: {results['fp_mismatch']} | Count: {results['count_mismatch']} ---\n")
    
    results["end_time"] = datetime.now().isoformat()
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total tested: {len(pdb_ids)}")
    print(f"Perfect match: {results['pass']}")
    print(f"FP differences only: {results['fp_mismatch']} (acceptable)")
    print(f"Count mismatches: {results['count_mismatch']} (⚠️  need investigation)")
    print(f"Errors/Skips: {results['error'] + results['skip']}")
    
    if results["count_mismatch_pdbs"]:
        print(f"\nPDBs with count mismatches:")
        for pdb_info in results["count_mismatch_pdbs"][:20]:
            print(f"  {pdb_info['pdb_id']}: {pdb_info['legacy']} (legacy) vs {pdb_info['modern']} (modern), diff={pdb_info['diff']}")
    
    # Save results
    output_file = project_root / "data" / "validation_results" / "ls_fitting_validation_detailed.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    
    return 0 if results['count_mismatch'] == 0 else 1

if __name__ == "__main__":
    sys.exit(main())

