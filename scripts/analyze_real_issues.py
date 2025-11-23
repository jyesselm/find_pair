#!/usr/bin/env python3
"""
Analyze the 69 real atom set differences and the 1 significant RMS difference.

Focus on understanding what's wrong and if legacy mode fixes it.
"""

import json
from pathlib import Path
from x3dna_json_compare import JsonComparator

def analyze_real_issues(project_root: Path):
    """Analyze the 69 real differences."""
    differences_file = project_root / "docs" / "REAL_DIFFERENCES.json"
    
    if not differences_file.exists():
        print("Run find_real_differences.py first!")
        return
    
    with open(differences_file) as f:
        data = json.load(f)
    
    atom_set_diffs = data['atom_set_differences']
    rms_diffs = data['rms_differences']
    
    print("=" * 80)
    print("ANALYZING 69 REAL ATOM SET DIFFERENCES")
    print("=" * 80)
    
    # Analyze patterns
    print(f"\nTotal atom set differences: {len(atom_set_diffs)}\n")
    
    # Pattern 1: Legacy has N7, modern has C1' and C4 (pyrimidines getting purine atoms)
    pattern1 = [d for d in atom_set_diffs if ' N7 ' in d['only_legacy'] and " C1'" in d['only_modern']]
    print(f"Pattern 1 - Legacy N7 vs Modern C1'+C4 (pyrimidine issue): {len(pattern1)}")
    if pattern1:
        print(f"  Example: {pattern1[0]['pdb_id']} {pattern1[0]['residue']} ({pattern1[0]['base_type']})")
        print(f"    Legacy: {pattern1[0]['legacy_atoms']}")
        print(f"    Modern: {pattern1[0]['modern_atoms']}")
        print()
    
    # Pattern 2: Modern has C4, legacy doesn't
    pattern2 = [d for d in atom_set_diffs if ' C4 ' in d['only_modern'] and ' C4 ' not in d['only_legacy'] and len(d['only_legacy']) == 0]
    print(f"Pattern 2 - Modern has C4, legacy doesn't: {len(pattern2)}")
    if pattern2:
        print(f"  Example: {pattern2[0]['pdb_id']} {pattern2[0]['residue']} ({pattern2[0]['base_type']})")
        print(f"    Legacy: {pattern2[0]['legacy_atoms']}")
        print(f"    Modern: {pattern2[0]['modern_atoms']}")
        print()
    
    # Pattern 3: Legacy has H, modern has C1' and C4
    pattern3 = [d for d in atom_set_diffs if ' H' in d['only_legacy']]
    print(f"Pattern 3 - Legacy has H atom: {len(pattern3)}")
    if pattern3:
        print(f"  Example: {pattern3[0]['pdb_id']} {pattern3[0]['residue']} ({pattern3[0]['base_type']})")
        print(f"    Legacy: {pattern3[0]['legacy_atoms']}")
        print(f"    Modern: {pattern3[0]['modern_atoms']}")
        print()
    
    # Check if legacy mode fixes these
    print("=" * 80)
    print("TESTING WITH LEGACY MODE")
    print("=" * 80)
    
    # Get unique PDBs with issues
    unique_pdbs = set(d['pdb_id'] for d in atom_set_diffs)
    print(f"\nPDBs with real differences: {len(unique_pdbs)}")
    print(f"PDB IDs: {sorted(list(unique_pdbs))[:20]}...\n")
    
    # Test 1H4S since we know it has differences
    test_pdb = "1H4S"
    legacy_file = project_root / f"data/json_legacy/{test_pdb}.json"
    modern_file = project_root / f"data/json/{test_pdb}.json"
    modern_legacy_file = project_root / f"data/json/{test_pdb}_legacy.json"
    
    if modern_legacy_file.exists():
        print(f"Testing {test_pdb} with legacy mode...")
        comparator = JsonComparator(enable_cache=False)
        result = comparator.compare_files(legacy_file, modern_legacy_file, None, test_pdb)
        
        if result.frame_comparison:
            fc = result.frame_comparison
            legacy_mode_diffs = []
            for mismatch in fc.mismatched_calculations:
                leg_atoms = set(mismatch.legacy_matched_atoms)
                mod_atoms = set(mismatch.modern_matched_atoms)
                only_legacy = leg_atoms - mod_atoms
                only_modern = mod_atoms - leg_atoms
                
                if only_legacy or only_modern:
                    chain_id, residue_seq, insertion = mismatch.residue_key
                    legacy_mode_diffs.append({
                        'residue': f"{chain_id}:{residue_seq}{insertion}",
                        'only_legacy': sorted(only_legacy),
                        'only_modern': sorted(only_modern),
                    })
            
            print(f"  Legacy mode differences: {len(legacy_mode_diffs)}")
            if legacy_mode_diffs:
                print("  Still has differences:")
                for diff in legacy_mode_diffs[:5]:
                    print(f"    {diff['residue']}: Legacy only: {diff['only_legacy']}, Modern only: {diff['only_modern']}")
            else:
                print("  ✅ Legacy mode fixes all differences!")
    else:
        print(f"⚠️  {test_pdb}_legacy.json not found. Generate with --legacy flag.")
    
    # Analyze the significant RMS difference
    print("\n" + "=" * 80)
    print("SIGNIFICANT RMS DIFFERENCE")
    print("=" * 80)
    
    if rms_diffs:
        for diff in rms_diffs:
            print(f"\n{diff['pdb_id']} {diff['residue']}: {diff['residue_name']} ({diff['base_type']})")
            print(f"  RMS Legacy: {diff['rms_legacy']:.6f}")
            print(f"  RMS Modern: {diff['rms_modern']:.6f}")
            print(f"  RMS Diff: {diff['rms_diff']:.6f}")
            print(f"  Legacy atoms: {diff['legacy_atoms']}")
            print(f"  Modern atoms: {diff['modern_atoms']}")
            print(f"  Only legacy: {diff['only_legacy']}")
            print(f"  Only modern: {diff['only_modern']}")
    
    # Summary and recommendations
    print("\n" + "=" * 80)
    print("ISSUES FOUND")
    print("=" * 80)
    
    print("\n1. **Pyrimidine/Purine Mix-up**: Legacy code includes N7 for pyrimidines")
    print("   - Pattern: Legacy has N7, modern has C1' + C4")
    print("   - Likely cause: Legacy code bug or different base type detection")
    
    print("\n2. **C4 Atom Differences**: Modern includes C4, legacy doesn't")
    print("   - Pattern: Modern has C4, legacy excludes it")
    print("   - Likely cause: Legacy mode flag or C4 exclusion bug")
    
    print("\n3. **H Atom in Legacy**: Legacy includes H atom")
    print("   - Pattern: Legacy has H, modern doesn't")
    print("   - Likely cause: Legacy code bug or different atom matching")
    
    print("\n4. **Large RMS Difference**: 1H4S T:55 (PSU) has 1.5Å RMS difference")
    print("   - Likely cause: Different atom sets or matching issue")
    
    print("\n**Recommendations:**")
    print("1. Check if legacy mode fixes Pattern 2 (C4 differences)")
    print("2. Investigate why legacy includes N7 for pyrimidines (Pattern 1)")
    print("3. Investigate why legacy includes H atom (Pattern 3)")
    print("4. Deep dive into 1H4S T:55 PSU residue")

if __name__ == '__main__':
    project_root = Path(__file__).parent.parent
    analyze_real_issues(project_root)

