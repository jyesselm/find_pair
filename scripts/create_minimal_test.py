#!/usr/bin/env python3
"""
Create minimal test cases for ref_frames debugging.

Extracts a small fragment from a PDB file to create a minimal test case
with just 1-2 base pairs for easier debugging.
"""

import sys
from pathlib import Path
from typing import List, Tuple

def extract_pdb_fragment(pdb_file: Path, start_res: int, end_res: int, 
                         output_file: Path, chain_id: str = 'A'):
    """Extract a fragment from a PDB file."""
    atoms = []
    in_fragment = False
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                res_seq = int(line[22:26].strip())
                chain = line[21]
                
                if chain == chain_id and start_res <= res_seq <= end_res:
                    atoms.append(line)
                    in_fragment = True
                elif in_fragment:
                    # We've passed the fragment
                    break
    
    if not atoms:
        print(f"ERROR: No atoms found for chain {chain_id}, residues {start_res}-{end_res}")
        return False
    
    with open(output_file, 'w') as f:
        f.write("REMARK   Minimal test case extracted from " + str(pdb_file.name) + "\n")
        f.write("REMARK   Chain: " + chain_id + ", Residues: " + 
                f"{start_res}-{end_res}\n")
        for atom in atoms:
            f.write(atom)
        f.write("END\n")
    
    print(f"Created minimal PDB: {output_file}")
    print(f"  Chain: {chain_id}, Residues: {start_res}-{end_res}")
    print(f"  Atoms: {len(atoms)}")
    return True

def find_base_pairs_in_pdb(pdb_file: Path) -> List[Tuple[int, int, str, str]]:
    """Find potential base pairs in a PDB file (simple heuristic)."""
    # This is a placeholder - in reality, you'd need to run find_pair
    # to find actual base pairs
    print("Note: Run find_pair to find actual base pairs")
    return []

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Create minimal test cases for ref_frames debugging')
    parser.add_argument('pdb_file', type=Path, help='Input PDB file')
    parser.add_argument('output_file', type=Path, help='Output minimal PDB file')
    parser.add_argument('--start-res', type=int, required=True,
                       help='Start residue number')
    parser.add_argument('--end-res', type=int, required=True,
                       help='End residue number')
    parser.add_argument('--chain', type=str, default='A',
                       help='Chain ID (default: A)')
    
    args = parser.parse_args()
    
    if not args.pdb_file.exists():
        print(f"ERROR: PDB file not found: {args.pdb_file}")
        sys.exit(1)
    
    success = extract_pdb_fragment(
        args.pdb_file, args.start_res, args.end_res, 
        args.output_file, args.chain)
    
    if success:
        print(f"\nNext steps:")
        print(f"1. Run legacy find_pair on: {args.output_file}")
        print(f"2. Run modern find_pair on: {args.output_file}")
        print(f"3. Compare ref_frames using: python3 scripts/debug_ref_frames.py ...")
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()

