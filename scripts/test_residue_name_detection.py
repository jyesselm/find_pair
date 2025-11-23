#!/usr/bin/env python3
"""
Test residue name detection - compare legacy vs modern logic.
"""

import json
from pathlib import Path


def test_residue_name_parsing():
    """Test how residue names are parsed and recognized."""
    
    # Test cases from actual PDB
    test_names = [
        "A",      # Single letter (from PDB)
        "  A",     # With leading spaces (legacy format)
        " C",      # With one leading space
        "  C",     # With two leading spaces
        "G",       # Single letter
        "  G",     # With spaces
        "T",       # Single letter
        "  T",     # With spaces
        "U",       # Single letter
        "  U",     # With spaces
        "ADE",     # Three letter
        "CYT",     # Three letter
        "GUA",     # Three letter
        "THY",     # Three letter
        "URA",     # Three letter
        "DA",      # DNA format
        "DC",      # DNA format
        "DG",      # DNA format
        "DT",      # DNA format
    ]
    
    print("Testing residue name recognition:")
    print("=" * 80)
    
    # Simulate modern code logic
    def modern_one_letter_code(name):
        """Simulate modern one_letter_code() logic."""
        # trim_name() removes leading/trailing spaces
        trimmed = name.strip()
        
        # Check matches
        if trimmed == "A" or trimmed == "ADE" or trimmed == "DA":
            return 'A'
        if trimmed == "C" or trimmed == "CYT" or trimmed == "DC":
            return 'C'
        if trimmed == "G" or trimmed == "GUA" or trimmed == "DG":
            return 'G'
        if trimmed == "T" or trimmed == "THY" or trimmed == "DT":
            return 'T'
        if trimmed == "U" or trimmed == "URA" or trimmed == "DU":
            return 'U'
        return '?'
    
    # Simulate legacy code logic (NT_LIST matching)
    legacy_nt_list = [
        "  A", "  C", "  G", "  I", "  T", "  U",
        "ADE", "CYT", "GUA", "INO", "THY", "URA",
        " +A", " +C", " +G", " +I", " +T", " +U"
    ]
    
    def legacy_matches(name):
        """Check if name matches NT_LIST."""
        # Legacy uses num_strmatch which does exact string matching
        return name in legacy_nt_list
    
    print("\nResidue Name | Modern Code | Legacy Match | Status")
    print("-" * 80)
    
    for name in test_names:
        modern_code = modern_one_letter_code(name)
        legacy_match = legacy_matches(name)
        
        status = "✓" if (modern_code != '?' or legacy_match) else "✗"
        if modern_code != '?' and not legacy_match:
            status = "⚠ MODERN ONLY"
        elif legacy_match and modern_code == '?':
            status = "⚠ LEGACY ONLY"
        
        print(f"{repr(name):15} | {modern_code:11} | {str(legacy_match):11} | {status}")
    
    print("\n" + "=" * 80)
    print("\nKey Findings:")
    print("1. Modern code uses trim_name() which removes spaces")
    print("2. Legacy code uses exact string matching with NT_LIST")
    print("3. Single-letter names like 'A' should work in modern (trimmed to 'A')")
    print("4. But legacy expects '  A' (with spaces) in NT_LIST")
    print("\nIssue: Modern code should recognize 'A' but may not match legacy behavior")
    print("        if residue names are stored as 'A' instead of '  A'")


def check_actual_residue_names(pdb_id: str, project_root: Path):
    """Check actual residue names in parsed structure."""
    modern_json_file = project_root / 'data' / 'json' / f'{pdb_id}.json'
    
    if not modern_json_file.exists():
        print(f"Modern JSON not found: {modern_json_file}")
        return
    
    with open(modern_json_file) as f:
        data = json.load(f)
    
    # Get pdb_atoms record
    atoms_rec = None
    for rec in data.get('calculations', []):
        if rec.get('type') == 'pdb_atoms':
            atoms_rec = rec
            break
    
    if not atoms_rec:
        print("No pdb_atoms record found")
        return
    
    atoms = atoms_rec.get('atoms', [])
    
    # Group by residue
    residues = {}
    for atom in atoms:
        rname = atom.get('residue_name', '')
        chain = atom.get('chain_id', '')
        seq = atom.get('residue_seq', 0)
        key = (rname, chain, seq)
        if key not in residues:
            residues[key] = {
                'name': rname,
                'name_repr': repr(rname),
                'name_len': len(rname),
                'atoms': []
            }
        residues[key]['atoms'].append(atom.get('atom_name', ''))
    
    # Find nucleotide-like residues
    print(f"\nChecking residue names in {pdb_id}:")
    print("=" * 80)
    
    nuc_like = []
    for key, data in residues.items():
        rname = data['name']
        trimmed = rname.strip()
        if trimmed in ['A', 'C', 'G', 'T', 'U', 'ADE', 'CYT', 'GUA', 'THY', 'URA', 'DA', 'DC', 'DG', 'DT']:
            nuc_like.append((key, data))
    
    print(f"\nFound {len(nuc_like)} nucleotide-like residues:")
    print("\nResidue Name (repr) | Length | Trimmed | Modern Code | Chain:Seq")
    print("-" * 80)
    
    for (rname, chain, seq), data in nuc_like[:30]:
        trimmed = rname.strip()
        modern_code = '?'
        if trimmed == "A" or trimmed == "ADE" or trimmed == "DA":
            modern_code = 'A'
        elif trimmed == "C" or trimmed == "CYT" or trimmed == "DC":
            modern_code = 'C'
        elif trimmed == "G" or trimmed == "GUA" or trimmed == "DG":
            modern_code = 'G'
        elif trimmed == "T" or trimmed == "THY" or trimmed == "DT":
            modern_code = 'T'
        elif trimmed == "U" or trimmed == "URA" or trimmed == "DU":
            modern_code = 'U'
        
        print(f"{data['name_repr']:20} | {data['name_len']:6} | {trimmed:7} | {modern_code:11} | {chain}:{seq}")
    
    # Check if any have frame calculations
    frame_calcs = [r for r in data.get('calculations', [])
                   if r.get('type') in ['base_frame_calc', 'frame_calc', 'ls_fitting']]
    
    print(f"\nFrame calculations found: {len(frame_calcs)}")
    
    if frame_calcs:
        print("\nFirst 5 frame calculations:")
        for rec in frame_calcs[:5]:
            rname = rec.get('residue_name', '')
            chain = rec.get('chain_id', '')
            seq = rec.get('residue_seq', 0)
            base = rec.get('base_type', '?')
            print(f"  {rname} {base} {chain}:{seq}")


if __name__ == '__main__':
    import sys
    
    test_residue_name_parsing()
    
    if len(sys.argv) > 1:
        pdb_id = sys.argv[1]
        project_root = Path(__file__).parent.parent
        check_actual_residue_names(pdb_id, project_root)
    else:
        print("\nUsage: python3 test_residue_name_detection.py [PDB_ID]")
        print("Example: python3 test_residue_name_detection.py 1H4S")

