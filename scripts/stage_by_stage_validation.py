#!/usr/bin/env python3
"""
Stage-by-stage validation of modern vs legacy JSON.

For each PDB in valid_pdbs.json, validate each stage and record:
- Which PDBs have full matches
- Which don't and what issues they have

Creates one JSON file per stage with validation results.
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional
import argparse


# Stage definitions
STAGES = {
    "stage1_atoms": {
        "name": "Atoms",
        "directories": ["pdb_atoms"],
        "description": "PDB atom parsing"
    },
    "stage2_frames": {
        "name": "Reference Frames",
        "directories": ["base_frame_calc", "frame_calc"],
        "description": "Frame calculation and reference frames"
    },
    "stage3_distances": {
        "name": "Distance Checks",
        "directories": ["distance_checks"],
        "description": "Geometric distance measurements"
    },
    "stage4_hbonds": {
        "name": "H-Bonds",
        "directories": ["hbond_list"],
        "description": "Hydrogen bond detection"
    },
    "stage5_validation": {
        "name": "Pair Validation",
        "directories": ["pair_validation"],
        "description": "Validation results for tested pairs"
    },
    "stage6_selection": {
        "name": "Pair Selection",
        "directories": ["find_bestpair_selection", "base_pair"],
        "description": "Final selected base pairs (PRIMARY OUTPUT)"
    },
    "stage7_steps": {
        "name": "Step Parameters",
        "directories": ["bpstep_params"],
        "description": "Step parameters (Shift, Slide, Rise, Tilt, Roll, Twist)"
    },
    "stage8_helical": {
        "name": "Helical Parameters",
        "directories": ["helical_params"],
        "description": "Helical parameters"
    }
}


def load_valid_pdbs(valid_pdbs_file: Path) -> List[str]:
    """Load list of valid PDBs."""
    with open(valid_pdbs_file) as f:
        data = json.load(f)
    return data.get("valid_pdbs_with_atoms_and_frames", [])


def compare_record_counts(legacy_file: Path, modern_file: Path) -> Dict:
    """Compare record counts between legacy and modern JSON files."""
    result = {
        "legacy_exists": legacy_file.exists(),
        "modern_exists": modern_file.exists(),
        "match": False,
        "legacy_count": None,
        "modern_count": None,
        "issue": None
    }
    
    if not result["legacy_exists"]:
        result["issue"] = "No legacy JSON"
        return result
    
    if not result["modern_exists"]:
        result["issue"] = "No modern JSON"
        return result
    
    try:
        with open(legacy_file) as f:
            legacy_data = json.load(f)
        legacy_count = len(legacy_data) if isinstance(legacy_data, list) else 1
        result["legacy_count"] = legacy_count
    except Exception as e:
        result["issue"] = f"Legacy JSON error: {str(e)}"
        return result
    
    try:
        with open(modern_file) as f:
            modern_data = json.load(f)
        modern_count = len(modern_data) if isinstance(modern_data, list) else 1
        result["modern_count"] = modern_count
    except Exception as e:
        result["issue"] = f"Modern JSON error: {str(e)}"
        return result
    
    result["match"] = (legacy_count == modern_count)
    if not result["match"]:
        result["issue"] = f"Count mismatch: legacy={legacy_count}, modern={modern_count}"
    
    return result


def validate_pdb_for_stage(pdb_id: str, stage_info: Dict, legacy_dir: Path, modern_dir: Path) -> Dict:
    """Validate a single PDB for a specific stage."""
    result = {
        "pdb_id": pdb_id,
        "stage_match": True,
        "directories": {},
        "overall_issue": None
    }
    
    # Check each directory for this stage
    for dir_name in stage_info["directories"]:
        legacy_file = legacy_dir / dir_name / f"{pdb_id}.json"
        modern_file = modern_dir / dir_name / f"{pdb_id}.json"
        
        comparison = compare_record_counts(legacy_file, modern_file)
        result["directories"][dir_name] = comparison
        
        if not comparison["match"]:
            result["stage_match"] = False
            if result["overall_issue"] is None:
                result["overall_issue"] = comparison["issue"]
    
    return result


def validate_stage(stage_id: str, stage_info: Dict, pdbs: List[str], 
                   legacy_dir: Path, modern_dir: Path, output_dir: Path,
                   limit: Optional[int] = None) -> Dict:
    """Validate all PDBs for a specific stage."""
    print(f"\n{'='*60}")
    print(f"Stage: {stage_info['name']}")
    print(f"Description: {stage_info['description']}")
    print(f"Directories: {', '.join(stage_info['directories'])}")
    print(f"{'='*60}")
    
    results = {
        "stage_id": stage_id,
        "stage_name": stage_info["name"],
        "description": stage_info["description"],
        "directories": stage_info["directories"],
        "total_pdbs_tested": 0,
        "passed": [],
        "failed": [],
        "summary": {
            "passed_count": 0,
            "failed_count": 0,
            "pass_rate": 0.0
        }
    }
    
    # Test PDBs
    test_pdbs = pdbs[:limit] if limit else pdbs
    results["total_pdbs_tested"] = len(test_pdbs)
    
    for i, pdb_id in enumerate(test_pdbs, 1):
        if i % 100 == 0:
            print(f"  Tested {i}/{len(test_pdbs)} PDBs...")
        
        validation = validate_pdb_for_stage(pdb_id, stage_info, legacy_dir, modern_dir)
        
        if validation["stage_match"]:
            results["passed"].append(pdb_id)
        else:
            results["failed"].append({
                "pdb_id": pdb_id,
                "issue": validation["overall_issue"],
                "details": validation["directories"]
            })
    
    # Calculate summary
    results["summary"]["passed_count"] = len(results["passed"])
    results["summary"]["failed_count"] = len(results["failed"])
    if results["total_pdbs_tested"] > 0:
        results["summary"]["pass_rate"] = (results["summary"]["passed_count"] / 
                                           results["total_pdbs_tested"]) * 100
    
    # Save results
    output_file = output_dir / f"{stage_id}_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Print summary
    print(f"\n  Results:")
    print(f"    Passed: {results['summary']['passed_count']}")
    print(f"    Failed: {results['summary']['failed_count']}")
    print(f"    Pass Rate: {results['summary']['pass_rate']:.2f}%")
    print(f"  Saved to: {output_file}")
    
    # Show first few failures
    if results["failed"]:
        print(f"\n  First 5 failures:")
        for fail in results["failed"][:5]:
            print(f"    - {fail['pdb_id']}: {fail['issue']}")
    
    return results


def main():
    parser = argparse.ArgumentParser(description="Stage-by-stage JSON validation")
    parser.add_argument("--stage", type=str, help="Specific stage to test (stage1_atoms, stage2_frames, etc.)")
    parser.add_argument("--limit", type=int, help="Limit number of PDBs to test (for quick testing)")
    parser.add_argument("--legacy-dir", type=str, default="data/json_legacy", 
                       help="Legacy JSON directory")
    parser.add_argument("--modern-dir", type=str, default="data/json", 
                       help="Modern JSON directory")
    parser.add_argument("--output-dir", type=str, default="data/validation_results",
                       help="Output directory for result JSON files")
    parser.add_argument("--valid-pdbs", type=str, default="data/valid_pdbs.json",
                       help="Path to valid_pdbs.json file")
    
    args = parser.parse_args()
    
    # Setup paths
    legacy_dir = Path(args.legacy_dir)
    modern_dir = Path(args.modern_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Load valid PDBs
    valid_pdbs_file = Path(args.valid_pdbs)
    pdbs = load_valid_pdbs(valid_pdbs_file)
    print(f"Loaded {len(pdbs)} valid PDBs from {valid_pdbs_file}")
    
    if args.limit:
        print(f"Limiting to first {args.limit} PDBs")
    
    # Run validation
    if args.stage:
        # Validate specific stage
        if args.stage not in STAGES:
            print(f"Error: Unknown stage '{args.stage}'")
            print(f"Available stages: {', '.join(STAGES.keys())}")
            return
        
        stage_info = STAGES[args.stage]
        validate_stage(args.stage, stage_info, pdbs, legacy_dir, modern_dir, output_dir, args.limit)
    else:
        # Validate all stages
        print(f"\nValidating ALL stages for {len(pdbs)} PDBs")
        print("This will create 8 result files, one per stage\n")
        
        all_results = {}
        for stage_id, stage_info in STAGES.items():
            results = validate_stage(stage_id, stage_info, pdbs, legacy_dir, modern_dir, 
                                    output_dir, args.limit)
            all_results[stage_id] = results["summary"]
        
        # Print overall summary
        print(f"\n{'='*60}")
        print("OVERALL SUMMARY")
        print(f"{'='*60}")
        for stage_id, summary in all_results.items():
            stage_name = STAGES[stage_id]["name"]
            print(f"{stage_name:20s}: {summary['passed_count']:4d} passed, "
                  f"{summary['failed_count']:4d} failed, "
                  f"{summary['pass_rate']:6.2f}% pass rate")


if __name__ == "__main__":
    main()

