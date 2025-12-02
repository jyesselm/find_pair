#!/usr/bin/env python3
"""
Batch generate modern JSON and verify against legacy.
Processes PDBs in small batches to save space.
"""

import os
import sys
import json
import subprocess
import csv
import numpy as np
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from dataclasses import dataclass
from typing import List, Set, Tuple, Optional

# Known legacy data issues - skip these
KNOWN_LEGACY_ISSUES = {'1EFW', '1QRU', '1TN1', '1TN2', '1TTT'}


@dataclass
class PDBResult:
    pdb_id: str
    status: str  # 'perfect', 'difference', 'error', 'skipped'
    legacy_frames: int = 0
    modern_frames: int = 0
    common_frames: int = 0
    max_origin_diff: float = 0.0
    max_orient_diff: float = 0.0
    legacy_pairs: int = 0
    modern_pairs: int = 0
    missing_pairs: int = 0
    extra_pairs: int = 0
    error_message: str = ""


def generate_modern_json(pdb_id: str, base_dir: str) -> Tuple[bool, str]:
    """Generate modern JSON for a PDB."""
    pdb_path = os.path.join(base_dir, f'data/pdb/{pdb_id}.pdb')
    output_dir = os.path.join(base_dir, 'data/json')
    
    if not os.path.exists(pdb_path):
        return False, f"PDB file not found: {pdb_path}"
    
    cmd = [
        os.path.join(base_dir, 'build/generate_modern_json'),
        pdb_path,
        output_dir,
        '--fix-indices'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            return False, result.stderr[:200]
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)


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
        l = legacy_by_idx[idx]
        m = modern_by_idx[idx]
        
        if 'origin' in l and 'origin' in m:
            diff = np.sqrt(sum((l['origin'][i] - m['origin'][i])**2 for i in range(3)))
            max_origin_diff = max(max_origin_diff, diff)
        
        if 'ref_frame' in l and 'ref_frame' in m:
            for row in range(3):
                diff = np.sqrt(sum((l['ref_frame'][row][i] - m['ref_frame'][row][i])**2 for i in range(3)))
                max_orient_diff = max(max_orient_diff, diff)
    
    return (len(legacy), len(modern), len(common), max_origin_diff, max_orient_diff)


def compare_pairs(pdb_id: str, base_dir: str) -> Optional[Tuple[int, int, int, int]]:
    """Compare pair selection."""
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
    
    missing = len(legacy - modern)
    extra = len(modern - legacy)
    
    return (len(legacy), len(modern), missing, extra)


def process_pdb(pdb_id: str, base_dir: str, generate: bool = True) -> PDBResult:
    """Process a single PDB: generate (if needed) and verify."""
    result = PDBResult(pdb_id=pdb_id, status='error')
    
    if pdb_id in KNOWN_LEGACY_ISSUES:
        result.status = 'skipped'
        result.error_message = 'Known legacy data issue'
        return result
    
    try:
        # Generate modern JSON if requested
        if generate:
            modern_frame_path = os.path.join(base_dir, f'data/json/base_frame_calc/{pdb_id}.json')
            if not os.path.exists(modern_frame_path):
                success, error = generate_modern_json(pdb_id, base_dir)
                if not success:
                    result.status = 'error'
                    result.error_message = f'Generation failed: {error}'
                    return result
        
        # Compare frames
        frame_result = compare_frames(pdb_id, base_dir)
        if frame_result:
            result.legacy_frames, result.modern_frames, result.common_frames, \
                result.max_origin_diff, result.max_orient_diff = frame_result
        
        # Compare pairs
        pair_result = compare_pairs(pdb_id, base_dir)
        if pair_result:
            result.legacy_pairs, result.modern_pairs, result.missing_pairs, result.extra_pairs = pair_result
        
        # Determine status
        frames_ok = result.max_origin_diff < 0.0001 and result.max_orient_diff < 0.0001
        pairs_ok = result.missing_pairs == 0 and result.extra_pairs == 0
        
        if frames_ok and pairs_ok:
            result.status = 'perfect'
        elif frame_result is None and pair_result is None:
            result.status = 'error'
            result.error_message = 'No data to compare'
        else:
            result.status = 'difference'
            
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
    
    def process_batch(self, pdbs: List[str], generate: bool = True, threads: int = 4) -> List[PDBResult]:
        """Process a batch of PDBs."""
        self.total = len(pdbs)
        self.processed = 0
        results = []
        
        def process_one(pdb_id):
            result = process_pdb(pdb_id, self.base_dir, generate)
            with self.lock:
                self.processed += 1
                if self.processed % 10 == 0 or self.processed == self.total:
                    print(f"  Progress: {self.processed}/{self.total}")
            return result
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(process_one, pdb): pdb for pdb in pdbs}
            for future in as_completed(futures):
                results.append(future.result())
        
        return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Batch generate and verify PDBs')
    parser.add_argument('--limit', '-l', type=int, default=50, help='Number of PDBs to process')
    parser.add_argument('--offset', '-o', type=int, default=0, help='Start offset')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads')
    parser.add_argument('--no-generate', action='store_true', help='Skip generation, only verify')
    parser.add_argument('--output', default='data/batch_results.csv', help='Output CSV file')
    args = parser.parse_args()
    
    base_dir = '/Users/jyesselman2/Library/CloudStorage/Dropbox/2_code/cpp/find_pair_2'
    
    # Find PDBs that need processing (have legacy but no modern)
    legacy_frames = set(f.replace('.json','') for f in os.listdir(f'{base_dir}/data/json_legacy/base_frame_calc') if f.endswith('.json'))
    modern_frames = set(f.replace('.json','') for f in os.listdir(f'{base_dir}/data/json/base_frame_calc') if f.endswith('.json'))
    
    # PDBs with legacy frames but no modern frames
    need_processing = sorted(legacy_frames - modern_frames)
    
    # Also check which have PDB files
    pdb_dir = f'{base_dir}/data/pdb'
    available_pdbs = set(f.replace('.pdb','') for f in os.listdir(pdb_dir) if f.endswith('.pdb'))
    need_processing = [p for p in need_processing if p in available_pdbs]
    
    print("=" * 70)
    print("BATCH GENERATE AND VERIFY")
    print("=" * 70)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"PDBs needing processing: {len(need_processing)}")
    print(f"Processing: {args.offset} to {args.offset + args.limit}")
    print("-" * 70)
    
    # Select batch
    batch = need_processing[args.offset:args.offset + args.limit]
    
    if not batch:
        print("No PDBs to process!")
        return
    
    print(f"\nProcessing {len(batch)} PDBs...")
    
    processor = BatchProcessor(base_dir)
    results = processor.process_batch(batch, generate=not args.no_generate, threads=args.threads)
    
    # Categorize results
    perfect = [r for r in results if r.status == 'perfect']
    differences = [r for r in results if r.status == 'difference']
    errors = [r for r in results if r.status == 'error']
    skipped = [r for r in results if r.status == 'skipped']
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total processed: {len(results)}")
    print(f"  ✅ Perfect:     {len(perfect)}")
    print(f"  ⚠️ Differences: {len(differences)}")
    print(f"  ❌ Errors:      {len(errors)}")
    print(f"  ⏭️ Skipped:     {len(skipped)}")
    
    # Save to CSV (append mode)
    csv_path = os.path.join(base_dir, args.output)
    file_exists = os.path.exists(csv_path)
    
    with open(csv_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(['pdb_id', 'status', 'legacy_frames', 'modern_frames', 'common_frames',
                            'max_origin_diff', 'max_orient_diff', 'legacy_pairs', 'modern_pairs',
                            'missing_pairs', 'extra_pairs', 'error'])
        for r in sorted(results, key=lambda x: x.pdb_id):
            writer.writerow([r.pdb_id, r.status, r.legacy_frames, r.modern_frames, r.common_frames,
                            f'{r.max_origin_diff:.10f}', f'{r.max_orient_diff:.10f}',
                            r.legacy_pairs, r.modern_pairs, r.missing_pairs, r.extra_pairs,
                            r.error_message])
    
    print(f"\nResults appended to: {csv_path}")
    
    # Show differences
    if differences:
        print("\n" + "-" * 70)
        print("PDBs WITH DIFFERENCES")
        print("-" * 70)
        for r in differences[:10]:
            print(f"  {r.pdb_id}: origin={r.max_origin_diff:.6f}, orient={r.max_orient_diff:.6f}, "
                  f"missing={r.missing_pairs}, extra={r.extra_pairs}")
    
    if errors:
        print("\n" + "-" * 70)
        print("ERRORS")
        print("-" * 70)
        for r in errors[:10]:
            print(f"  {r.pdb_id}: {r.error_message}")


if __name__ == '__main__':
    main()


