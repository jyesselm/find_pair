#!/usr/bin/env python3
"""
Compare generated JSON files with legacy JSON files.

This script compares JSON files in data/json/ with corresponding files
in data/json_legacy/ to verify that the generated files match the legacy
format and values.

Usage:
    python scripts/compare_json_files.py [--tolerance TOL] [--threads N] [--verbose]
"""

import json
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from threading import Lock
import argparse
from collections import defaultdict
import multiprocessing


class JsonComparator:
    """Compare JSON files with tolerance for floating point values."""
    
    def __init__(self, tolerance: float = 0.001):
        self.tolerance = tolerance
        self.print_lock = Lock()
    
    def approximately_equal(self, a: float, b: float) -> bool:
        """Check if two floats are approximately equal within tolerance."""
        return abs(a - b) < self.tolerance
    
    def compare_numbers(self, val1: Any, val2: Any, path: str) -> List[str]:
        """Compare two numeric values."""
        errors = []
        try:
            num1 = float(val1)
            num2 = float(val2)
            if not self.approximately_equal(num1, num2):
                errors.append(f"{path}: {num1} != {num2} (diff: {abs(num1 - num2)})")
        except (ValueError, TypeError):
            if val1 != val2:
                errors.append(f"{path}: {val1} != {val2}")
        return errors
    
    def compare_arrays(self, arr1: List, arr2: List, path: str) -> List[str]:
        """Compare two arrays."""
        errors = []
        if len(arr1) != len(arr2):
            errors.append(f"{path}: Array length mismatch ({len(arr1)} != {len(arr2)})")
            return errors
        
        for i, (item1, item2) in enumerate(zip(arr1, arr2)):
            errors.extend(self.compare_value(item1, item2, f"{path}[{i}]"))
        
        return errors
    
    def compare_objects(self, obj1: Dict, obj2: Dict, path: str) -> List[str]:
        """Compare two JSON objects."""
        errors = []
        
        # Get all keys from both objects
        all_keys = set(obj1.keys()) | set(obj2.keys())
        
        for key in all_keys:
            key_path = f"{path}.{key}" if path else key
            
            if key not in obj1:
                errors.append(f"{key_path}: Missing in generated JSON")
                continue
            
            if key not in obj2:
                errors.append(f"{key_path}: Missing in legacy JSON")
                continue
            
            errors.extend(self.compare_value(obj1[key], obj2[key], key_path))
        
        return errors
    
    def compare_value(self, val1: Any, val2: Any, path: str = "") -> List[str]:
        """Compare two JSON values recursively."""
        errors = []
        
        # Type mismatch
        if type(val1) != type(val2):
            # Allow int/float comparison
            if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                return self.compare_numbers(val1, val2, path)
            errors.append(f"{path}: Type mismatch ({type(val1).__name__} != {type(val2).__name__})")
            return errors
        
        # Compare based on type
        if isinstance(val1, dict):
            errors.extend(self.compare_objects(val1, val2, path))
        elif isinstance(val1, list):
            errors.extend(self.compare_arrays(val1, val2, path))
        elif isinstance(val1, (int, float)):
            errors.extend(self.compare_numbers(val1, val2, path))
        elif val1 != val2:
            errors.append(f"{path}: {val1} != {val2}")
        
        return errors
    
    def compare_calculations(self, gen_calcs: List[Dict], leg_calcs: List[Dict]) -> List[str]:
        """Compare calculations arrays, matching by type and index."""
        errors = []
        
        # Group by type
        gen_by_type = defaultdict(list)
        leg_by_type = defaultdict(list)
        
        for calc in gen_calcs:
            if "type" in calc:
                gen_by_type[calc["type"]].append(calc)
        
        for calc in leg_calcs:
            if "type" in calc:
                leg_by_type[calc["type"]].append(calc)
        
        # Compare each type
        all_types = set(gen_by_type.keys()) | set(leg_by_type.keys())
        
        for calc_type in all_types:
            gen_list = gen_by_type[calc_type]
            leg_list = leg_by_type[calc_type]
            
            if len(gen_list) != len(leg_list):
                errors.append(
                    f"calculations[{calc_type}]: Count mismatch "
                    f"({len(gen_list)} != {len(leg_list)})"
                )
                continue
            
            # Compare each record of this type
            for i, (gen_calc, leg_calc) in enumerate(zip(gen_list, leg_list)):
                errors.extend(
                    self.compare_value(
                        gen_calc, 
                        leg_calc, 
                        f"calculations[{calc_type}][{i}]"
                    )
                )
        
        return errors
    
    def compare_files(self, gen_file: Path, leg_file: Path) -> Tuple[bool, List[str]]:
        """Compare two JSON files."""
        errors = []
        
        try:
            # Load JSON files
            with open(gen_file, 'r') as f:
                gen_json = json.load(f)
            
            with open(leg_file, 'r') as f:
                leg_json = json.load(f)
            
            # Compare top-level fields
            if "pdb_name" in gen_json and "pdb_name" in leg_json:
                if gen_json["pdb_name"] != leg_json["pdb_name"]:
                    errors.append("pdb_name mismatch")
            
            # Compare calculations arrays
            if "calculations" in gen_json and "calculations" in leg_json:
                errors.extend(
                    self.compare_calculations(
                        gen_json["calculations"],
                        leg_json["calculations"]
                    )
                )
            elif "calculations" in gen_json:
                errors.append("calculations missing in legacy JSON")
            elif "calculations" in leg_json:
                errors.append("calculations missing in generated JSON")
            
            # Compare other top-level fields
            for key in set(gen_json.keys()) | set(leg_json.keys()):
                if key == "calculations":
                    continue  # Already compared
                if key not in gen_json:
                    errors.append(f"{key}: Missing in generated JSON")
                elif key not in leg_json:
                    errors.append(f"{key}: Missing in legacy JSON")
                else:
                    errors.extend(self.compare_value(gen_json[key], leg_json[key], key))
            
            return len(errors) == 0, errors
            
        except json.JSONDecodeError as e:
            return False, [f"JSON decode error: {e}"]
        except Exception as e:
            return False, [f"Error comparing files: {e}"]


def find_json_pairs(gen_dir: Path, leg_dir: Path) -> List[Tuple[Path, Path, str]]:
    """Find pairs of generated and legacy JSON files."""
    pairs = []
    
    if not gen_dir.exists():
        print(f"Warning: Generated JSON directory does not exist: {gen_dir}")
        return pairs
    
    if not leg_dir.exists():
        print(f"Warning: Legacy JSON directory does not exist: {leg_dir}")
        return pairs
    
    # Find all generated JSON files
    for gen_file in gen_dir.glob("*.json"):
        # Skip globals files
        if "_globals" in gen_file.name:
            continue
        
        pdb_name = gen_file.stem
        leg_file = leg_dir / f"{pdb_name}.json"
        
        if leg_file.exists():
            pairs.append((gen_file, leg_file, pdb_name))
    
    return pairs


def compare_single_pair_worker(
    gen_file_str: str,
    leg_file_str: str,
    pdb_name: str,
    tolerance: float
) -> Tuple[str, bool, List[str]]:
    """Worker function for comparing a single pair of JSON files.
    
    This function is used by multiprocessing.
    It creates its own comparator to avoid sharing state.
    Accepts string paths (not Path objects) for pickling.
    """
    gen_file = Path(gen_file_str)
    leg_file = Path(leg_file_str)
    comparator = JsonComparator(tolerance=tolerance)
    success, errors = comparator.compare_files(gen_file, leg_file)
    return pdb_name, success, errors


def compare_single_pair(
    gen_file: Path,
    leg_file: Path,
    pdb_name: str,
    comparator: JsonComparator,
    verbose: bool = False
) -> Tuple[str, bool, List[str]]:
    """Compare a single pair of JSON files (threading version with shared comparator)."""
    success, errors = comparator.compare_files(gen_file, leg_file)
    
    if verbose or not success:
        with comparator.print_lock:
            status = "✓" if success else "✗"
            print(f"{status} {pdb_name}: ", end="")
            if success:
                print("MATCH")
            else:
                print(f"FAILED ({len(errors)} errors)")
                if verbose and errors:
                    for error in errors[:10]:  # Print first 10 errors
                        print(f"    {error}")
                    if len(errors) > 10:
                        print(f"    ... and {len(errors) - 10} more errors")
    
    return pdb_name, success, errors


def main():
    parser = argparse.ArgumentParser(
        description="Compare generated JSON files with legacy JSON files"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.001,
        help="Tolerance for floating point comparisons (default: 0.001)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads/processes (default: CPU count)"
    )
    parser.add_argument(
        "--use-processes",
        action="store_true",
        help="Use multiprocessing (default, better for CPU-bound JSON parsing)"
    )
    parser.add_argument(
        "--use-threads",
        action="store_true",
        help="Use threading instead of multiprocessing (better for I/O-bound tasks)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed error messages"
    )
    parser.add_argument(
        "--gen-dir",
        type=Path,
        default=Path("data/json"),
        help="Directory containing generated JSON files (default: data/json)"
    )
    parser.add_argument(
        "--leg-dir",
        type=Path,
        default=Path("data/json_legacy"),
        help="Directory containing legacy JSON files (default: data/json_legacy)"
    )
    
    args = parser.parse_args()
    
    # Find JSON file pairs
    pairs = find_json_pairs(args.gen_dir, args.leg_dir)
    
    if not pairs:
        print("No JSON file pairs found.")
        print(f"  Generated JSON dir: {args.gen_dir}")
        print(f"  Legacy JSON dir: {args.leg_dir}")
        sys.exit(1)
    
    print(f"Found {len(pairs)} JSON file pairs to compare")
    print(f"Tolerance: {args.tolerance}")
    
    # Determine execution method
    use_processes = args.use_processes or (not args.use_threads)
    num_workers = args.threads or multiprocessing.cpu_count()
    
    print(f"Workers: {num_workers} ({'processes' if use_processes else 'threads'})")
    print()
    
    # Compare files in parallel
    results = []
    
    if use_processes:
        # Use multiprocessing (better for CPU-bound JSON parsing)
        # Convert Path objects to strings for pickling
        tasks = [
            (str(gen_file), str(leg_file), pdb_name, args.tolerance)
            for gen_file, leg_file, pdb_name in pairs
        ]
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {
                executor.submit(compare_single_pair_worker, *task): task[2]  # task[2] is pdb_name
                for task in tasks
            }
            
            completed = 0
            for future in as_completed(futures):
                pdb_name, success, errors = future.result()
                results.append((pdb_name, success, errors))
                completed += 1
                
                # Print progress
                if args.verbose or not success:
                    status = "✓" if success else "✗"
                    print(f"{status} {pdb_name}: ", end="")
                    if success:
                        print("MATCH")
                    else:
                        print(f"FAILED ({len(errors)} errors)")
                        if args.verbose and errors:
                            for error in errors[:10]:
                                print(f"    {error}")
                            if len(errors) > 10:
                                print(f"    ... and {len(errors) - 10} more errors")
                elif completed % 10 == 0:
                    # Print progress every 10 files
                    print(f"Progress: {completed}/{len(pairs)} files compared...")
    else:
        # Use threading (better for I/O-bound tasks, but JSON parsing is CPU-bound)
        comparator = JsonComparator(tolerance=args.tolerance)
        
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = {
                executor.submit(
                    compare_single_pair,
                    gen_file,
                    leg_file,
                    pdb_name,
                    comparator,
                    args.verbose
                ): (gen_file, leg_file, pdb_name)
                for gen_file, leg_file, pdb_name in pairs
            }
            
            completed = 0
            for future in as_completed(futures):
                pdb_name, success, errors = future.result()
                results.append((pdb_name, success, errors))
                completed += 1
                
                # Print progress for non-verbose mode
                if not args.verbose and completed % 10 == 0:
                    print(f"Progress: {completed}/{len(pairs)} files compared...")
    
    # Print summary
    print()
    print("=" * 60)
    print("Comparison Summary")
    print("=" * 60)
    
    total = len(results)
    passed = sum(1 for _, success, _ in results if success)
    failed = total - passed
    
    print(f"Total files compared: {total}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    
    if failed > 0:
        print()
        print("Failed files:")
        for pdb_name, success, errors in results:
            if not success:
                print(f"  {pdb_name}: {len(errors)} errors")
                if args.verbose:
                    for error in errors[:5]:
                        print(f"    - {error}")
                    if len(errors) > 5:
                        print(f"    ... and {len(errors) - 5} more")
    
    # Exit with error code if any failures
    sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()

