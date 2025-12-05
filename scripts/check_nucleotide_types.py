#!/usr/bin/env python3
"""
Use RCSB PDB GraphQL API to verify modified nucleotide base types
"""
import json
import urllib.request
import urllib.error

# GraphQL query
QUERY = """
query molecule ($id: String!) {
    chem_comp(comp_id:$id){
        chem_comp {
            id
            name
            formula
            type
        }
        rcsb_chem_comp_info {
            atom_count
        }
        rcsb_chem_comp_synonyms {
          name
          type
        }
    }
}
"""

def query_pdb_api(comp_id):
    """Query PDB API for chemical component info"""
    url = "https://data.rcsb.org/graphql"
    
    variables = {"id": comp_id}
    payload = {
        "query": QUERY,
        "variables": variables
    }
    
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(
        url,
        data=data,
        headers={'Content-Type': 'application/json'}
    )
    
    try:
        with urllib.request.urlopen(req) as response:
            result = json.loads(response.read().decode('utf-8'))
            return result
    except urllib.error.URLError as e:
        print(f"Error querying {comp_id}: {e}")
        return None

def determine_base_type(name, synonyms):
    """Determine base type from chemical name"""
    # Combine name and synonyms for searching
    all_names = [name.upper()] + [s['name'].upper() for s in synonyms if s.get('name')]
    
    for n in all_names:
        if 'URIDINE' in n or 'URACIL' in n:
            return 'URACIL', 'u'
        elif 'CYTIDINE' in n or 'CYTOSINE' in n:
            return 'CYTOSINE', 'c'
        elif 'GUANOSINE' in n or 'GUANINE' in n:
            return 'GUANINE', 'g'
        elif 'THYMIDINE' in n or 'THYMINE' in n:
            return 'THYMINE', 't'
        elif 'ADENOSINE' in n or 'ADENINE' in n:
            return 'ADENINE', 'a'
        elif 'INOSINE' in n or 'HYPOXANTHINE' in n:
            return 'INOSINE', 'i'
    
    return 'UNKNOWN', '?'

def main():
    # Nucleotides to check
    nucleotides = ['CM0', 'JSP', 'NCA', '9DG']
    
    print("=" * 80)
    print("Verifying Modified Nucleotide Base Types Using RCSB PDB API")
    print("=" * 80)
    print()
    
    results = {}
    
    for nuc_id in nucleotides:
        print(f"Querying {nuc_id}...")
        result = query_pdb_api(nuc_id)
        
        if result and 'data' in result and result['data']['chem_comp']:
            comp_data = result['data']['chem_comp']
            chem_comp = comp_data['chem_comp']
            synonyms = comp_data.get('rcsb_chem_comp_synonyms', []) or []
            
            name = chem_comp['name']
            formula = chem_comp['formula']
            comp_type = chem_comp['type']
            
            base_type, one_letter = determine_base_type(name, synonyms)
            
            results[nuc_id] = {
                'name': name,
                'formula': formula,
                'type': comp_type,
                'base_type': base_type,
                'one_letter': one_letter
            }
            
            print(f"  ✓ {nuc_id}: {name}")
            print(f"    Base Type: {base_type} (code: {one_letter})")
            print(f"    Formula: {formula}")
            print()
        else:
            print(f"  ✗ Failed to get data for {nuc_id}")
            print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY - Required Registry Updates:")
    print("=" * 80)
    print()
    
    # Current registry values (what we have now)
    current = {
        'CM0': {'code': 't', 'type': 'THYMINE'},
        'JSP': {'code': 't', 'type': 'THYMINE'},
        'NCA': {'code': 't', 'type': 'THYMINE'},
        '9DG': {'code': 'g', 'type': 'GUANINE'}
    }
    
    for nuc_id, data in results.items():
        curr = current.get(nuc_id, {})
        print(f"{nuc_id}: {data['name']}")
        print(f"  Current registry: code='{curr.get('code', '?')}', type={curr.get('type', '?')}")
        print(f"  Correct values:   code='{data['one_letter']}', type={data['base_type']}")
        
        if curr.get('code') == data['one_letter'] and curr.get('type') == data['base_type']:
            print(f"  Status: ✓ CORRECT - No change needed")
        else:
            print(f"  Status: ✗ WRONG - Needs fixing!")
        print()
    
    # Output JSON format for easy updating
    print("=" * 80)
    print("JSON Format for Registry Update:")
    print("=" * 80)
    print()
    
    for nuc_id, data in results.items():
        is_purine = data['base_type'] in ['ADENINE', 'GUANINE', 'INOSINE']
        print(f'  "{nuc_id}": {{')
        print(f'    "code": "{data["one_letter"]}",')
        print(f'    "type": "{data["base_type"]}",')
        print(f'    "is_purine": {str(is_purine).lower()},')
        print(f'    "description": "{data["name"]}"')
        print(f'  }},')
    
    print()

if __name__ == '__main__':
    main()

