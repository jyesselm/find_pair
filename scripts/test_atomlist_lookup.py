#!/usr/bin/env python3
"""
Test atomlist.dat lookup logic to verify pattern matching.
"""

import re
from pathlib import Path

def load_atomlist(atomlist_file: Path):
    """Load atomlist.dat into a dictionary."""
    atomlist = {}
    with open(atomlist_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                aname = parts[0].upper()
                asym = parts[1].upper()
                if len(asym) == 1:
                    asym = ' ' + asym
                atomlist[aname] = asym
    return atomlist

def create_pattern(atom_name: str) -> str:
    """Create pattern by replacing non-alphabetic with '.' (matches legacy)."""
    pattern = list(atom_name)
    for i in range(min(4, len(pattern))):
        if not pattern[i].isalpha():
            pattern[i] = '.'
    return ''.join(pattern[:4])

def str_pmatch(pattern: str, aname: str) -> bool:
    """Match legacy str_pmatch: checks if pattern starts with aname."""
    return pattern.startswith(aname)

def aname2asym_sim(atom_name: str, atomlist: dict) -> str:
    """Simulate legacy aname2asym logic."""
    # Step 1: Create pattern
    aname_pattern = create_pattern(atom_name)
    
    # Step 2: Try to match against atomlist
    my_asym = None
    found_match = False
    for pattern, sym in atomlist.items():
        if str_pmatch(pattern, aname_pattern):
            my_asym = sym
            found_match = True
            break
    
    # Step 3: Fallback logic
    if not found_match:
        unknown = (aname_pattern == ".UNK")
        
        # Fallback case 1
        if (len(aname_pattern) >= 4 and aname_pattern[0] != '.' and 
            aname_pattern[1] != '.' and aname_pattern[2] == '.' and aname_pattern[3] == '.'):
            my_asym = ' ' + aname_pattern[1]
        # Fallback case 2
        elif (len(aname_pattern) >= 2 and aname_pattern[0] == '.' and 
              aname_pattern[1] != '.' and not unknown):
            my_asym = ' ' + aname_pattern[1]
        # Fallback case 3
        elif len(aname_pattern) >= 1 and aname_pattern[0] == 'H':
            my_asym = ' H'
        else:
            my_asym = 'XX'  # UNKATM
    
    return my_asym

def asym_idx(sym: str) -> int:
    """Map atomic symbol to idx (matches legacy atoms_list)."""
    atoms_list = {' C': 1, ' O': 2, ' H': 3, ' N': 4, ' S': 5, ' P': 6}
    return atoms_list.get(sym, 0)

def main():
    atomlist_file = Path('/Users/jyesselman2/installs/x3dna/config/atomlist.dat')
    atomlist = load_atomlist(atomlist_file)
    
    print(f"Loaded {len(atomlist)} entries from atomlist.dat\n")
    
    # Test cases from mismatched pairs
    test_atoms = [
        (' O4 ', '1VBY pair (20,21)'),
        (' N3 ', '1VBY pair (20,21)'),
        (' N3 ', '3AVY pair (1204,1223)'),
        (' N1 ', '3AVY pair (1204,1223)'),
    ]
    
    print("=== Testing atom name to idx conversion ===\n")
    for atom_name, description in test_atoms:
        pattern = create_pattern(atom_name)
        asym = aname2asym_sim(atom_name, atomlist)
        idx = asym_idx(asym)
        
        print(f"{description}:")
        print(f"  Atom name: '{atom_name}'")
        print(f"  Pattern: '{pattern}'")
        print(f"  Atomic symbol: '{asym}'")
        print(f"  idx: {idx}")
        
        # Check if pattern matches any atomlist entry
        matches = [p for p in atomlist.keys() if str_pmatch(p, pattern)]
        if matches:
            print(f"  Atomlist matches: {matches[:3]}")  # Show first 3
        else:
            print(f"  Atomlist matches: None (using fallback)")
        print()

if __name__ == '__main__':
    main()

