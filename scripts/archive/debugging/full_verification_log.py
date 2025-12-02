#!/usr/bin/env python3
"""
Full verification with detailed logging of all differences.
"""

import os
import sys
import json
import csv
import subprocess
import numpy as np
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from dataclasses import dataclass, field
from typing import List, Set, Tuple, Dict, Optional

KNOWN_LEGACY_ISSUES = {'1EFW', '1QRU', '1TN1', '1TN2', '1TTT'}


@dataclass
class PairDiff:
    pdb_id: str
    pair_i: int
    pair_j: int
    diff_type: str  # 'missing' or 'extra'
    bp_type: str = ""


@dataclass 
class PDBResult:
    pdb_id: str
    status: str
    legacy_frames: int = 0
    modern_frames: int = 0
    common_frames: int = 0
    max_origin_diff: float = 0.0
    max_orient_diff: float = 0.0
    legacy_pairs: int = 0
    modern_pairs: int = 0
    missing_pairs: int = 0
    extra_pairs: int = 0
    pair_diffs: List[PairDiff] = field(default_factory=list)
    error_message: str = ""


def generate_modern_json(pdb_id: str, base_dir: str) -> Tuple[bool, str]:
    """Generate modern JSON for a PDB."""
    pdb_path = os.path.join(base_dir, f'data/pdb/{pdb_id}.pdb')
    output_dir = os.path.join(base_dir, 'data/json')
    
    if not os.path.exists(pdb_path):
        return False, f"PDB file not found"
    
    cmd = [
        os.path.join(base_dir, 'build/generate_modern_json'),
        pdb_path, output_dir, '--fix-indices'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            return False, result.stderr[:100]
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)[:100]


def compare_frames(pdb_id: str, base_dir: str) -> Optional[Tuple[int, int, int, float, float]]:
    """Compare frame origins and orientations."""
    legacy_path = os.path.join(base_dir, f'data/json_legacy/base_frame_calc/{pdb_id}.json')
    modern_path = os.path.join(base_dir, f'data/json/base_frame_calc/{pdb_id}.json')
    
    if not os.path.exists(legacy_path) or not os.path.exists(modern_path):
        return None
    
    with open(legacy_path) as f:
        legacy = json.load(f)
    with open(modern_path) as f:
        modern = json.load(f)
    
    legacy_by_idx = {r['residue_idx']: r for r in legacy}
    modern_by_idx = {r['residue_idx']: r for r in modern}
    common = set(legacy_by_idx.keys()) & set(modern_by_idx.keys())
    
    max_origin_diff = 0.0
    max_orient_diff = 0.0
    
    for idx in common:
        l, m = legacy_by_idx[idx], modern_by_idx[idx]
        if 'origin' in l and 'origin' in m:
            diff = np.sqrt(sum((l['origin'][i] - m['origin'][i])**2 for i in range(3)))
            max_origin_diff = max(max_origin_diff, diff)
        if 'ref_frame' in l and 'ref_frame' in m:
            for row in range(3):
                diff = np.sqrt(sum((l['ref_frame'][row][i] - m['ref_frame'][row][i])**2 for i in range(3)))
                max_orient_diff = max(max_orient_diff, diff)
    
    return (len(legacy), len(modern), len(common), max_origin_diff, max_orient_diff)


def compare_pairs(pdb_id: str, base_dir: str) -> Optional[Tuple[int, int, List[PairDiff]]]:
    """Compare pair selection and return detailed diffs."""
    legacy_path = os.path.join(base_dir, f'data/json_legacy/find_bestpair_selection/{pdb_id}.json')
    modern_path = os.path.join(base_dir, f'data/json/find_bestpair_selection/{pdb_id}.json')
    
    if not os.path.exists(legacy_path) or not os.path.exists(modern_path):
        return None
    
    def load_pairs(filepath):
        with open(filepath) as f:
            data = json.load(f)
        pairs = set()
        for rec in data:
            if 'pairs' in rec:
                for pair in rec['pairs']:
                    if len(pair) >= 2:
                        pairs.add((min(pair[0], pair[1]), max(pair[0], pair[1])))
        return pairs
    
    legacy = load_pairs(legacy_path)
    modern = load_pairs(modern_path)
    
    diffs = []
    for pair in sorted(legacy - modern):
        diffs.append(PairDiff(pdb_id, pair[0], pair[1], 'missing'))
    for pair in sorted(modern - legacy):
        diffs.append(PairDiff(pdb_id, pair[0], pair[1], 'extra'))
    
    return (len(legacy), len(modern), diffs)


def process_pdb(pdb_id: str, base_dir: str, generate: bool = True) -> PDBResult:
    """Process a single PDB."""
    result = PDBResult(pdb_id=pdb_id, status='error')
    
    if pdb_id in KNOWN_LEGACY_ISSUES:
        result.status = 'skipped'
        result.error_message = 'Known legacy data issue'
        return result
    
    try:
        if generate:
            modern_path = os.path.join(base_dir, f'data/json/base_frame_calc/{pdb_id}.json')
            if not os.path.exists(modern_path):
                success, error = generate_modern_json(pdb_id, base_dir)
                if not success:
                    result.status = 'error'
                    result.error_message = f'Generation: {error}'
                    return result
        
        frame_result = compare_frames(pdb_id, base_dir)
        if frame_result:
            result.legacy_frames, result.modern_frames, result.common_frames, \
                result.max_origin_diff, result.max_orient_diff = frame_result
        
        pair_result = compare_pairs(pdb_id, base_dir)
        if pair_result:
            result.legacy_pairs, result.modern_pairs, result.pair_diffs = pair_result
            result.missing_pairs = sum(1 for d in result.pair_diffs if d.diff_type == 'missing')
            result.extra_pairs = sum(1 for d in result.pair_diffs if d.diff_type == 'extra')
        
        frames_ok = result.max_origin_diff < 0.0001 and result.max_orient_diff < 0.0001
        pairs_ok = result.missing_pairs == 0 and result.extra_pairs == 0
        
        if frames_ok and pairs_ok:
            result.status = 'perfect'
        elif frame_result is None and pair_result is None:
            result.status = 'error'
            result.error_message = 'No data'
        else:
            result.status = 'difference'
            
    except json.JSONDecodeError as e:
        result.status = 'error'
        result.error_message = f'JSON parse error: {str(e)[:50]}'
    except Exception as e:
        result.status = 'error'
        result.error_message = str(e)[:100]
    
    return result


class BatchProcessor:
    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        self.lock = Lock()
        self.processed = 0
        self.total = 0
        self.all_pair_diffs: List[PairDiff] = []
    
    def process_batch(self, pdbs: List[str], generate: bool = True, threads: int = 4) -> List[PDBResult]:
        self.total = len(pdbs)
        self.processed = 0
        results = []
        
        def process_one(pdb_id):
            result = process_pdb(pdb_id, self.base_dir, generate)
            with self.lock:
                self.processed += 1
                self.all_pair_diffs.extend(result.pair_diffs)
                if self.processed % 50 == 0 or self.processed == self.total:
                    print(f"  Progress: {self.processed}/{self.total} ({100*self.processed/self.total:.1f}%)")
            return result
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(process_one, pdb): pdb for pdb in pdbs}
            for future in as_completed(futures):
                results.append(future.result())
        
        return results


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', '-l', type=int, default=100)
    parser.add_argument('--offset', '-o', type=int, default=0)
    parser.add_argument('--threads', '-t', type=int, default=4)
    parser.add_argument('--no-generate', action='store_true')
    args = parser.parse_args()
    
    base_dir = '/Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2'
    
    # Load already processed
    results_file = os.path.join(base_dir, 'data/batch_results.csv')
    processed = set()
    if os.path.exists(results_file):
        with open(results_file) as f:
            processed = set(r['pdb_id'] for r in csv.DictReader(f))
    
    # Find remaining PDBs
    legacy_frames = set(f.replace('.json','') for f in os.listdir(f'{base_dir}/data/json_legacy/base_frame_calc') if f.endswith('.json'))
    pdb_files = set(f.replace('.pdb','') for f in os.listdir(f'{base_dir}/data/pdb') if f.endswith('.pdb'))
    remaining = sorted((legacy_frames & pdb_files) - processed)
    
    print("=" * 70)
    print("FULL VERIFICATION WITH DETAILED LOGGING")
    print("=" * 70)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Already processed: {len(processed)}")
    print(f"Remaining: {len(remaining)}")
    print(f"Processing: {args.offset} to {args.offset + args.limit}")
    print("-" * 70)
    
    batch = remaining[args.offset:args.offset + args.limit]
    if not batch:
        print("No PDBs to process!")
        return
    
    print(f"\nProcessing {len(batch)} PDBs...")
    processor = BatchProcessor(base_dir)
    results = processor.process_batch(batch, generate=not args.no_generate, threads=args.threads)
    
    # Categorize
    perfect = [r for r in results if r.status == 'perfect']
    differences = [r for r in results if r.status == 'difference']
    errors = [r for r in results if r.status == 'error']
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Perfect: {len(perfect)}, Differences: {len(differences)}, Errors: {len(errors)}")
    
    # Save main results (append)
    with open(results_file, 'a', newline='') as f:
        writer = csv.writer(f)
        for r in sorted(results, key=lambda x: x.pdb_id):
            writer.writerow([r.pdb_id, r.status, r.legacy_frames, r.modern_frames, r.common_frames,
                            f'{r.max_origin_diff:.10f}', f'{r.max_orient_diff:.10f}',
                            r.legacy_pairs, r.modern_pairs, r.missing_pairs, r.extra_pairs,
                            r.error_message])
    
    # Save detailed pair differences (append)
    diff_file = os.path.join(base_dir, 'data/all_pair_differences.csv')
    file_exists = os.path.exists(diff_file)
    with open(diff_file, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['pdb_id', 'pair_i', 'pair_j', 'diff_type', 'bp_type'])
        for d in processor.all_pair_diffs:
            writer.writerow([d.pdb_id, d.pair_i, d.pair_j, d.diff_type, d.bp_type])
    
    print(f"\nResults appended to: {results_file}")
    print(f"Pair diffs appended to: {diff_file}")
    
    # Show new differences
    if differences:
        print("\n" + "-" * 70)
        print(f"NEW PDBs WITH DIFFERENCES ({len(differences)})")
        print("-" * 70)
        for r in differences[:15]:
            print(f"  {r.pdb_id}: missing={r.missing_pairs}, extra={r.extra_pairs}")
            for d in r.pair_diffs[:3]:
                print(f"    {d.diff_type}: ({d.pair_i}, {d.pair_j})")


if __name__ == '__main__':
    main()

