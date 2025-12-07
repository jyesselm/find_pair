#!/usr/bin/env python3
"""
Document modified residues in each PDB.

This script scans the base_frame_calc JSON files and identifies residues
with modified nucleotides (lowercase base_type indicates modified bases).

Also queries the RCSB PDB API to get official names for each modified residue type.

Output: data/modified_residues.json
"""

import json
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import urllib.request
import urllib.error
import argparse
import time


# GraphQL query for PDB chemical component info
PDB_GRAPHQL_QUERY = """
query molecule ($id: String!) {
    chem_comp(comp_id:$id){
        chem_comp {
            id
            name
            formula
            type
        }
    }
}
"""


def query_pdb_api(comp_id: str) -> dict | None:
    """Query PDB API for chemical component info."""
    url = "https://data.rcsb.org/graphql"
    
    variables = {"id": comp_id}
    payload = {
        "query": PDB_GRAPHQL_QUERY,
        "variables": variables
    }
    
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(
        url,
        data=data,
        headers={'Content-Type': 'application/json'}
    )
    
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            result = json.loads(response.read().decode('utf-8'))
            if result and 'data' in result and result['data']['chem_comp']:
                comp_data = result['data']['chem_comp']['chem_comp']
                return {
                    'name': comp_data['name'],
                    'formula': comp_data['formula'],
                    'type': comp_data['type']
                }
    except (urllib.error.URLError, TimeoutError) as e:
        pass  # Silently fail, we'll use fallback
    
    return None


def fetch_pdb_names(residue_codes: list, max_workers: int = 20) -> dict:
    """Fetch official names for all residue codes from PDB API in parallel."""
    print(f"\nFetching official names from RCSB PDB API for {len(residue_codes)} residue types...")
    
    pdb_info = {}
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(query_pdb_api, code): code for code in residue_codes}
        
        completed = 0
        for future in as_completed(futures):
            completed += 1
            if completed % 50 == 0:
                print(f"  Queried {completed}/{len(residue_codes)} from PDB API...")
            
            code = futures[future]
            try:
                result = future.result()
                if result:
                    pdb_info[code] = result
            except Exception:
                pass  # Skip failures
    
    print(f"  Retrieved info for {len(pdb_info)}/{len(residue_codes)} residue types")
    return pdb_info


def is_modified(base_type: str) -> bool:
    """Check if base_type indicates a modified nucleotide."""
    if not base_type:
        return False
    # Lowercase letters indicate modified bases (a, g, c, t, u, p, i)
    # Uppercase are standard (A, G, C, T, U, P, I)
    return base_type.islower()


def scan_pdb_for_modified_residues(json_path: Path) -> dict:
    """
    Scan a base_frame_calc JSON file for modified residues.
    
    Returns dict with:
        - modified_residues: list of {residue_name, base_type, chain_id, residue_seq}
        - unique_modified_types: set of unique 3-letter codes
    """
    with open(json_path) as f:
        frames = json.load(f)
    
    modified = []
    unique_types = set()
    
    for r in frames:
        base_type = r.get('base_type', '')
        if is_modified(base_type):
            res_name = r.get('residue_name', '').strip()
            modified.append({
                'residue_name': res_name,
                'base_type': base_type,
                'chain_id': r.get('chain_id', ''),
                'residue_seq': r.get('residue_seq'),
                'legacy_residue_idx': r.get('legacy_residue_idx') or r.get('residue_idx')
            })
            unique_types.add(res_name)
    
    return {
        'modified_residues': modified,
        'unique_modified_types': list(sorted(unique_types))
    }


def process_single_pdb(json_path: Path) -> tuple:
    """Process a single PDB file. Returns (pdb_id, data, error)."""
    pdb_id = json_path.stem
    try:
        data = scan_pdb_for_modified_residues(json_path)
        return (pdb_id, data, None)
    except Exception as e:
        return (pdb_id, None, str(e))


def main():
    parser = argparse.ArgumentParser(description='Document modified residues in PDB files')
    parser.add_argument('--data-dir', type=Path, 
                        default=Path('data/json/base_frame_calc'),
                        help='Directory containing base_frame_calc JSON files')
    parser.add_argument('--output', type=Path,
                        default=Path('data/modified_residues.json'),
                        help='Output JSON file path')
    parser.add_argument('--summary', action='store_true',
                        help='Print summary statistics')
    parser.add_argument('--workers', type=int, default=20,
                        help='Number of parallel workers (default: 20)')
    args = parser.parse_args()
    
    # Scan all PDB files
    results = {}
    all_modified_types = defaultdict(int)
    pdbs_with_modified = 0
    total_modified_residues = 0
    
    json_files = sorted(args.data_dir.glob('*.json'))
    print(f"Scanning {len(json_files)} PDB files with {args.workers} workers...")
    
    # Process in parallel
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_single_pdb, jp): jp for jp in json_files}
        
        completed = 0
        for future in as_completed(futures):
            completed += 1
            if completed % 500 == 0:
                print(f"  Processed {completed}/{len(json_files)} files...")
            
            pdb_id, data, error = future.result()
            
            if error:
                print(f"Error processing {pdb_id}: {error}")
                continue
            
            if data and data['modified_residues']:
                pdbs_with_modified += 1
                total_modified_residues += len(data['modified_residues'])
                
                results[pdb_id] = {
                    'count': len(data['modified_residues']),
                    'unique_types': data['unique_modified_types'],
                    'residues': data['modified_residues']
                }
                
                for res_name in data['unique_modified_types']:
                    all_modified_types[res_name] += 1
    
    # Fetch official names from PDB API
    pdb_info = fetch_pdb_names(list(all_modified_types.keys()), max_workers=args.workers)
    
    # Build modified types info with PDB data
    modified_types_info = {}
    for res_name, count in sorted(all_modified_types.items(), key=lambda x: -x[1]):
        info = {
            'pdb_count': count
        }
        if res_name in pdb_info:
            info['official_name'] = pdb_info[res_name]['name']
            info['formula'] = pdb_info[res_name]['formula']
            info['component_type'] = pdb_info[res_name]['type']
        modified_types_info[res_name] = info
    
    # Create output structure
    output = {
        'summary': {
            'total_pdbs_scanned': len(json_files),
            'pdbs_with_modified_residues': pdbs_with_modified,
            'total_modified_residues': total_modified_residues,
            'unique_modified_types_count': len(all_modified_types),
        },
        'modified_residue_types': modified_types_info,
        'pdbs': results
    }
    
    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nOutput written to: {args.output}")
    
    # Print summary
    if args.summary or True:  # Always print summary
        print(f"\n{'='*60}")
        print("SUMMARY")
        print(f"{'='*60}")
        print(f"Total PDBs scanned:          {len(json_files)}")
        print(f"PDBs with modified residues: {pdbs_with_modified}")
        print(f"Total modified residues:     {total_modified_residues}")
        print(f"Unique modified types:       {len(all_modified_types)}")
        
        print(f"\n{'='*80}")
        print("TOP 30 MODIFIED RESIDUE TYPES (by frequency)")
        print(f"{'='*80}")
        print(f"{'Code':<6} {'# PDBs':<8} {'Official Name'}")
        print(f"{'-'*80}")
        
        for res_name, count in list(sorted(all_modified_types.items(), 
                                           key=lambda x: -x[1]))[:30]:
            official_name = pdb_info.get(res_name, {}).get('name', '')
            print(f"{res_name:<6} {count:<8} {official_name}")


if __name__ == '__main__':
    main()

