#!/usr/bin/env python3
"""
Compare helical parameters between legacy and modern code.
"""

import json
import sys
from pathlib import Path

def load_helical_params_from_json(json_file):
    """Load helical parameters from legacy JSON format."""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Handle both array format and object with calculations array
        records = []
        if isinstance(data, list):
            records = data
        elif isinstance(data, dict) and 'calculations' in data:
            records = data['calculations']
        elif isinstance(data, dict):
            # Single record
            records = [data]
        
        helical_params = []
        for record in records:
            if record.get('type') == 'helical_params':
                if 'params' in record:
                    params = record['params']
                    if isinstance(params, list) and len(params) >= 6:
                        helical_params.append({
                            'bp_idx1': record.get('bp_idx1', 0),
                            'bp_idx2': record.get('bp_idx2', 0),
                            'x_disp': params[0],
                            'y_disp': params[1],
                            'h_rise': params[2],
                            'inclination': params[3],
                            'tip': params[4],
                            'h_twist': params[5],
                        })
        return helical_params
    except Exception as e:
        print(f"Error loading {json_file}: {e}")
        return []

def compare_helical_params(legacy_file, modern_values):
    """Compare legacy helical params with modern calculated values."""
    legacy_params = load_helical_params_from_json(legacy_file)
    
    if not legacy_params:
        print(f"No helical parameters found in {legacy_file}")
        return False
    
    print(f"Found {len(legacy_params)} legacy helical parameters")
    print(f"Found {len(modern_values)} modern helical parameters")
    print()
    
    if len(legacy_params) != len(modern_values):
        print(f"⚠️  Count mismatch: Legacy={len(legacy_params)}, Modern={len(modern_values)}")
        print()
    
    matches = 0
    differences = []
    tolerance = 0.01
    
    for i, (leg, mod) in enumerate(zip(legacy_params, modern_values)):
        step_num = i + 1
        
        # Compare each parameter
        params_match = True
        step_diffs = []
        
        for param_name in ['x_disp', 'y_disp', 'h_rise', 'inclination', 'tip', 'h_twist']:
            leg_val = leg.get(param_name, 0.0)
            mod_val = mod.get(param_name, 0.0)
            diff = abs(leg_val - mod_val)
            
            if diff > tolerance:
                params_match = False
                step_diffs.append(f"  {param_name:12s}: Legacy={leg_val:10.4f}, Modern={mod_val:10.4f}, Diff={diff:10.4f}")
        
        if params_match:
            matches += 1
        else:
            differences.append({
                'step': step_num,
                'bp_idx1': leg.get('bp_idx1', 0),
                'bp_idx2': leg.get('bp_idx2', 0),
                'diffs': step_diffs
            })
    
    # Print results
    print(f"✅ Matches: {matches}/{len(legacy_params)}")
    if differences:
        print(f"❌ Differences: {len(differences)}/{len(legacy_params)}")
        print()
        print("Differences found:")
        for diff_info in differences:
            print(f"\nStep {diff_info['step']} (BP {diff_info['bp_idx1']}-{diff_info['bp_idx2']}):")
            for line in diff_info['diffs']:
                print(line)
    else:
        print("✅ All parameters match within tolerance!")
    
    return len(differences) == 0

def main():
    if len(sys.argv) < 2:
        print("Usage: compare_helical_params.py <pdb_id>")
        print("Example: compare_helical_params.py 6V9Q")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    
    # Look for legacy JSON files
    legacy_paths = [
        f"data/json_legacy/helical_params/{pdb_id}.json",
        f"data/json_legacy/bp_helical/{pdb_id}.json",
        f"data/json_legacy/{pdb_id}_globals.json",
    ]
    
    legacy_file = None
    for path in legacy_paths:
        if Path(path).exists():
            legacy_file = path
            break
    
    if not legacy_file:
        print(f"⚠️  No legacy helical parameter JSON found for {pdb_id}")
        print(f"   Checked: {', '.join(legacy_paths)}")
        print()
        print("To generate legacy JSON, run legacy analyze on the PDB file.")
        sys.exit(1)
    
    print(f"Comparing helical parameters for {pdb_id}")
    print(f"Legacy file: {legacy_file}")
    print()
    
    # For now, we need to run modern code and extract values
    # This is a placeholder - in practice, you'd extract from modern JSON or output
    print("⚠️  Modern values extraction not yet implemented")
    print("   Need to either:")
    print("   1. Run analyze_app and parse stdout")
    print("   2. Generate modern JSON with helical parameters")
    print("   3. Use comparison framework")
    
    # Load legacy for inspection
    legacy_params = load_helical_params_from_json(legacy_file)
    if legacy_params:
        print()
        print(f"Legacy helical parameters found: {len(legacy_params)} steps")
        print("\nFirst few legacy parameters:")
        for i, params in enumerate(legacy_params[:3]):
            print(f"Step {i+1} (BP {params['bp_idx1']}-{params['bp_idx2']}):")
            print(f"  X-disp={params['x_disp']:8.2f}, Y-disp={params['y_disp']:8.2f}, "
                  f"h-Rise={params['h_rise']:7.2f}")
            print(f"  Incl.={params['inclination']:8.2f}, Tip={params['tip']:8.2f}, "
                  f"h-Twist={params['h_twist']:8.2f}")

if __name__ == "__main__":
    main()

