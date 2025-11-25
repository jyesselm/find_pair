#!/usr/bin/env python3
"""Calculate overall match rate across all PDBs"""

import json
from pathlib import Path

pdbs = ['1Q96', '1VBY', '3AVY', '4AL5', '3KNC', '5UJ2', '6LTU', '8J1J', '3G8T', '6CAQ']

total_legacy = 0
total_modern = 0
total_common = 0

for pdb in pdbs:
    try:
        with open(f'data/json/{pdb}_find_bestpair_selection.json') as f1:
            d1 = json.load(f1)
        with open(f'data/json_legacy/{pdb}_find_bestpair_selection.json') as f2:
            d2 = json.load(f2)
        
        p1 = set(tuple(sorted(p)) for p in d1[0]['pairs'])
        p2 = set(tuple(sorted(p)) for p in d2[0]['pairs'])
        
        total_legacy += len(p2)
        total_modern += len(p1)
        total_common += len(p1 & p2)
        
        match_rate = 100 * len(p1 & p2) / len(p2) if len(p2) > 0 else 0
        print(f'{pdb}: {len(p1 & p2)}/{len(p2)} = {match_rate:.1f}%')
    except Exception as e:
        print(f'{pdb}: Error - {e}')

print(f'\nOverall: {total_common}/{total_legacy} = {100*total_common/total_legacy:.2f}% match')
print(f'Total legacy pairs: {total_legacy}')
print(f'Total modern pairs: {total_modern}')
print(f'Total common pairs: {total_common}')
print(f'Total missing: {total_legacy - total_common}')
print(f'Total extra: {total_modern - total_common}')

