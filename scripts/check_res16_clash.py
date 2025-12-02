#!/usr/bin/env python3
"""Check if residue 16 and 59 are clashing and if it affects identification."""

import math

def parse_pdb_line(line):
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom_name = line[12:16].strip()
        return (x, y, z, atom_name)
    except:
        return None

# Get C1' atoms for residues 16 and 59 (chain D)
res16_c1 = None
res59_c1 = None

with open('data/pdb/1TTT.pdb') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21] if len(line) > 21 else ' '
            res_seq_str = line[22:26].strip()
            try:
                res_seq = int(res_seq_str)
                if chain == 'D':
                    if res_seq == 16:
                        atom_info = parse_pdb_line(line)
                        if atom_info and atom_info[3] == "C1'":
                            res16_c1 = atom_info[:3]
                    elif res_seq == 59:
                        atom_info = parse_pdb_line(line)
                        if atom_info and atom_info[3] == "C1'":
                            res59_c1 = atom_info[:3]
            except:
                pass

if res16_c1 and res59_c1:
    dist = math.sqrt(sum((a - b)**2 for a, b in zip(res16_c1, res59_c1)))
    print(f'Distance between Chain D ResSeq 16 C1\' and ResSeq 59 C1\': {dist:.2f} Å')
    
    # Check if this is a normal base pair distance (typically 10-20 Å)
    if dist < 5:
        print('⚠️  VERY CLOSE - Possible clash or overlap!')
    elif dist < 10:
        print('⚠️  Unusually close for base pairing')
    elif dist < 15:
        print('✅ Normal base pair distance range')
    else:
        print('✅ Typical base pair distance')
else:
    print('Could not find C1\' atoms')

# Get all atoms for residues 16 and 59 (chain D)
res16_atoms = []
res59_atoms = []

with open('data/pdb/1TTT.pdb') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21] if len(line) > 21 else ' '
            res_seq_str = line[22:26].strip()
            try:
                res_seq = int(res_seq_str)
                if chain == 'D':
                    atom_info = parse_pdb_line(line)
                    if atom_info:
                        if res_seq == 16:
                            res16_atoms.append(atom_info)
                        elif res_seq == 59:
                            res59_atoms.append(atom_info)
            except:
                pass

# Calculate minimum distance between any atoms
min_dist = float('inf')
clashing_atoms = []
for a16 in res16_atoms:
    for a59 in res59_atoms:
        dist = math.sqrt(sum((a16[i] - a59[i])**2 for i in range(3)))
        if dist < min_dist:
            min_dist = dist
            clashing_atoms = [(a16[3], a59[3], dist)]

print(f'\n=== Clash Analysis: Residue 16 vs Residue 59 (Chain D) ===\n')
print(f'Residue 16 atoms: {len(res16_atoms)}')
print(f'Residue 59 atoms: {len(res59_atoms)}')
print(f'\nMinimum distance between any atoms: {min_dist:.2f} Å')

if min_dist < 2.0:
    print('⚠️  SEVERE CLASH - atoms are overlapping!')
elif min_dist < 3.0:
    print('⚠️  MODERATE CLASH - very close contact')
elif min_dist < 4.0:
    print('⚠️  CLOSE CONTACT - may affect geometry')
else:
    print('✅ No clash detected')

if clashing_atoms:
    print(f'\nClosest atoms:')
    for a16_name, a59_name, dist in clashing_atoms:
        print(f'  Res16 {a16_name} ↔ Res59 {a59_name}: {dist:.2f} Å')

# Check NT_CUTOFF value
print(f'\n=== RMSD Check Context ===')
print('Legacy uses RMSD <= NT_CUTOFF to recognize nucleotides')
print('If residues 16 and 59 are clashing:')
print('1. This could distort the geometry of residue 16')
print('2. The RMSD calculation might fail (geometry doesn\'t match standard)')
print('3. This would explain why legacy skips residue 16')

