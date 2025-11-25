#!/usr/bin/env python3
"""
Batch test donor_acceptor function on multiple test cases
Compares legacy vs modern donor_acceptor behavior
"""

import json
import glob
import os
import subprocess
import sys

def extract_test_cases_from_common_patterns():
    """Create test cases from common H-bond patterns"""
    common_hbonds = [
        ("A", "T", " N6 ", " O4 "),  # Standard A-T
        ("A", "T", " N1 ", " N3 "),  # Standard A-T
        ("G", "C", " N2 ", " O2 "),  # Standard G-C
        ("G", "C", " N1 ", " N3 "),  # Standard G-C
        ("G", "C", " O6 ", " N4 "),  # Standard G-C
        ("C", "G", " N3 ", " N2 "),  # Standard C-G (our problem case)
        ("C", "G", " N1 ", " O2'"),  # Backbone
        ("T", "A", " O4 ", " N6 "),  # Standard T-A
        ("T", "A", " N3 ", " N1 "),  # Standard T-A
        ("U", "A", " O4 ", " N6 "),  # Standard U-A
        ("A", "U", " N6 ", " O4 "),  # Standard A-U
        ("G", "U", " O6 ", " N3 "),  # G-U wobble
        ("U", "G", " N3 ", " O6 "),  # U-G wobble
    ]
    
    test_cases = []
    for base1, base2, atom1, atom2 in common_hbonds:
        test_cases.append({
            'pdb': 'common',
            'res1': '?',
            'res2': '?',
            'base1': base1,
            'base2': base2,
            'atom1': atom1,
            'atom2': atom2,
            'expected_type': '?'  # Will be determined from legacy
        })
    
    return test_cases

def extract_test_cases_from_json(json_dir="data/json", max_pdbs=10):
    """Extract H-bond test cases from modern JSON files or use common patterns"""
    # Try to load from file first
    if os.path.exists('donor_acceptor_test_cases.json'):
        try:
            with open('donor_acceptor_test_cases.json') as f:
                return json.load(f)
        except:
            pass
    
    # Otherwise use common patterns
    return extract_test_cases_from_common_patterns()

def test_legacy_donor_acceptor(base1, base2, atom1, atom2):
    """Test legacy donor_acceptor function"""
    try:
        result = subprocess.run(
            ["org/build/bin/test_donor_acceptor", base1, base2, atom1, atom2],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode != 0:
            return None, result.stderr
        
        # Parse JSON from output
        lines = result.stdout.split('\n')
        json_start = False
        json_lines = []
        for line in lines:
            if line.strip().startswith('{'):
                json_start = True
            if json_start:
                json_lines.append(line)
            if line.strip().startswith('}'):
                break
        
        if json_lines:
            data = json.loads('\n'.join(json_lines))
            return data.get('hbond_type', '?'), None
        return None, "No JSON output"
    except subprocess.TimeoutExpired:
        return None, "Timeout"
    except Exception as e:
        return None, str(e)

def test_modern_donor_acceptor(base1, base2, atom1, atom2):
    """Test modern donor_acceptor function"""
    # Try multiple possible paths
    possible_paths = [
        "build/bin/debug_donor_acceptor",
        "build/debug_donor_acceptor",
        "tools/debug_donor_acceptor",
        "./debug_donor_acceptor"
    ]
    
    tool_path = None
    for path in possible_paths:
        if os.path.exists(path):
            tool_path = path
            break
    
    if not tool_path:
        return None, "Tool not found (tried: " + ", ".join(possible_paths) + ")"
    
    try:
        result = subprocess.run(
            [tool_path, base1, base2, atom1, atom2],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode != 0:
            return None, result.stderr
        
        # Parse output - check for "Type: X" or "H-bond type: X"
        for line in result.stdout.split('\n'):
            if 'Type:' in line or 'H-bond type:' in line or 'type:' in line.lower():
                # Extract the character after the colon
                parts = line.split(':')
                if len(parts) >= 2:
                    type_char = parts[1].strip().split()[0]
                    # Should be a single character
                    if len(type_char) == 1:
                        return type_char, None
        # If not found, check stderr
        if result.stderr:
            return None, result.stderr
        return None, f"Could not parse output: {result.stdout[:200]}"
    except subprocess.TimeoutExpired:
        return None, "Timeout"
    except Exception as e:
        return None, str(e)

def main():
    print("Extracting test cases from JSON files...")
    test_cases = extract_test_cases_from_json(max_pdbs=10)
    
    if not test_cases:
        print("No test cases found!")
        return 1
    
    print(f"Found {len(test_cases)} test cases\n")
    print("Testing donor_acceptor function (legacy vs modern)...\n")
    
    matches = 0
    mismatches = 0
    errors = 0
    
    results = []
    
    for i, tc in enumerate(test_cases, 1):
        base1 = tc['base1']
        base2 = tc['base2']
        atom1 = tc['atom1']
        atom2 = tc['atom2']
        
        # Test legacy
        legacy_type, legacy_error = test_legacy_donor_acceptor(base1, base2, atom1, atom2)
        
        # Test modern
        modern_type, modern_error = test_modern_donor_acceptor(base1, base2, atom1, atom2)
        
        # Compare
        match = (legacy_type == modern_type) if (legacy_type and modern_type) else False
        
        result = {
            'case': i,
            'pdb': tc['pdb'],
            'base1': base1,
            'base2': base2,
            'atom1': atom1,
            'atom2': atom2,
            'legacy_type': legacy_type,
            'modern_type': modern_type,
            'expected_type': tc['expected_type'],
            'match': match,
            'legacy_error': legacy_error,
            'modern_error': modern_error
        }
        results.append(result)
        
        if legacy_error or modern_error:
            errors += 1
            status = "ERROR"
            print(f"{i:3d}. {tc['pdb']}: {base1}-{base2} {atom1} -> {atom2}")
            print(f"     Legacy: {legacy_type if legacy_type else legacy_error}")
            print(f"     Modern: {modern_type if modern_type else modern_error}")
        elif match:
            matches += 1
            status = "MATCH"
        else:
            mismatches += 1
            status = "MISMATCH"
            print(f"{i:3d}. {tc['pdb']}: {base1}-{base2} {atom1} -> {atom2}")
            print(f"     Legacy: {legacy_type}, Modern: {modern_type}, Expected: {tc['expected_type']}")
            print(f"     Status: {status}")
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total test cases: {len(test_cases)}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Errors: {errors}")
    print(f"Match rate: {matches/(len(test_cases)-errors)*100:.1f}%" if (len(test_cases)-errors) > 0 else "N/A")
    
    # Save detailed results
    output_file = "donor_acceptor_test_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nDetailed results saved to: {output_file}")
    
    if mismatches > 0:
        print("\n⚠️  MISMATCHES FOUND - Need to investigate!")
        return 1
    elif errors > 0:
        print("\n⚠️  ERRORS FOUND - Need to fix test tools!")
        return 1
    else:
        print("\n✅ All tests passed!")
        return 0

if __name__ == "__main__":
    sys.exit(main())

