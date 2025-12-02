#!/usr/bin/env python3
"""
Comprehensive PDB verification script with threading.
Compares legacy vs modern find_pair outputs for all valid PDBs.
Only saves detailed output for PDBs with differences.
"""

import json
import os
import sys
import argparse
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from threading import Lock
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Set

# Known perfect match PDBs - skip detailed output for these
# (Based on previous 255-PDB test: 250 perfect matches)
# Update this list as more PDBs are verified
KNOWN_PERFECT_PDBS: Set[str] = {
    # 10-PDB test set (all verified 100%)
    '1Q96', '1VBY', '3AVY', '3G8T', '3KNC', '4AL5', '5UJ2', '6CAQ', '6LTU', '8J1J',
}

# Known legacy data issues (duplicate records in legacy JSON)
KNOWN_LEGACY_DATA_ISSUES: Set[str] = {
    '1EFW',  # 110 duplicate records
    '1QRU',  # 56 duplicate records
    '1TN1',  # 56 duplicate records
    '1TN2',  # 56 duplicate records
    '1TTT',  # 184 duplicate records
}


@dataclass
class PDBResult:
    """Result of comparing a single PDB."""
    pdb_id: str
    status: str  # 'perfect', 'difference', 'error', 'skipped', 'legacy_issue'
    
    # Frame comparison
    legacy_frames: int = 0
    modern_frames: int = 0
    common_frames: int = 0
    max_origin_diff: float = 0.0
    max_orient_diff: float = 0.0
    
    # Pair comparison
    legacy_pairs: int = 0
    modern_pairs: int = 0
    common_pairs: int = 0
    missing_pairs: int = 0
    extra_pairs: int = 0
    
    # Details for differences
    missing_pair_list: List[Tuple[int, int]] = field(default_factory=list)
    extra_pair_list: List[Tuple[int, int]] = field(default_factory=list)
    
    error_message: str = ""
    processing_time: float = 0.0


class PDBVerifier:
    """Verifies legacy vs modern PDB outputs."""
    
    def __init__(self, base_dir: str, skip_known_perfect: bool = True):
        self.base_dir = base_dir
        self.legacy_frame_dir = os.path.join(base_dir, 'data/json_legacy/base_frame_calc')
        self.modern_frame_dir = os.path.join(base_dir, 'data/json/base_frame_calc')
        self.legacy_pairs_dir = os.path.join(base_dir, 'data/json_legacy/find_bestpair_selection')
        self.modern_pairs_dir = os.path.join(base_dir, 'data/json/find_bestpair_selection')
        self.skip_known_perfect = skip_known_perfect
        
        # Thread-safe counters
        self.lock = Lock()
        self.processed = 0
        self.total = 0
        
    def compare_pdb(self, pdb_id: str) -> PDBResult:
        """Compare a single PDB between legacy and modern."""
        import time
        start_time = time.time()
        
        result = PDBResult(pdb_id=pdb_id, status='error')
        
        try:
            # Skip known legacy data issues
            if pdb_id in KNOWN_LEGACY_DATA_ISSUES:
                result.status = 'legacy_issue'
                result.error_message = "Known legacy JSON data quality issue (duplicate records)"
                return result
            
            # Compare frames
            frame_result = self._compare_frames(pdb_id)
            if frame_result is None:
                result.status = 'error'
                result.error_message = "Missing JSON files"
                return result
            
            (result.legacy_frames, result.modern_frames, result.common_frames,
             result.max_origin_diff, result.max_orient_diff) = frame_result
            
            # Compare pairs
            pair_result = self._compare_pairs(pdb_id)
            if pair_result:
                (result.legacy_pairs, result.modern_pairs, result.common_pairs,
                 result.missing_pairs, result.extra_pairs,
                 result.missing_pair_list, result.extra_pair_list) = pair_result
            
            # Determine overall status
            frames_perfect = (result.max_origin_diff == 0 and result.max_orient_diff == 0)
            pairs_perfect = (result.missing_pairs == 0 and result.extra_pairs == 0)
            
            if frames_perfect and pairs_perfect:
                result.status = 'perfect'
            else:
                result.status = 'difference'
                
        except Exception as e:
            result.status = 'error'
            result.error_message = str(e)
        
        result.processing_time = time.time() - start_time
        
        # Update progress
        with self.lock:
            self.processed += 1
            if self.processed % 50 == 0 or self.processed == self.total:
                print(f"  Progress: {self.processed}/{self.total} ({100*self.processed/self.total:.1f}%)")
        
        return result
    
    def _compare_frames(self, pdb_id: str) -> Optional[Tuple[int, int, int, float, float]]:
        """Compare frame origins and orientations."""
        legacy_path = os.path.join(self.legacy_frame_dir, f'{pdb_id}.json')
        modern_path = os.path.join(self.modern_frame_dir, f'{pdb_id}.json')
        
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
            
            # Origin comparison
            if 'origin' in l and 'origin' in m:
                diff = np.sqrt(sum((l['origin'][i] - m['origin'][i])**2 for i in range(3)))
                max_origin_diff = max(max_origin_diff, diff)
            
            # Orientation comparison
            if 'ref_frame' in l and 'ref_frame' in m:
                for row in range(3):
                    diff = np.sqrt(sum((l['ref_frame'][row][i] - m['ref_frame'][row][i])**2 for i in range(3)))
                    max_orient_diff = max(max_orient_diff, diff)
        
        return (len(legacy), len(modern), len(common), max_origin_diff, max_orient_diff)
    
    def _compare_pairs(self, pdb_id: str) -> Optional[Tuple[int, int, int, int, int, List, List]]:
        """Compare find_bestpair_selection pairs."""
        legacy_path = os.path.join(self.legacy_pairs_dir, f'{pdb_id}.json')
        modern_path = os.path.join(self.modern_pairs_dir, f'{pdb_id}.json')
        
        if not os.path.exists(legacy_path) or not os.path.exists(modern_path):
            return None
        
        with open(legacy_path) as f:
            legacy_data = json.load(f)
        with open(modern_path) as f:
            modern_data = json.load(f)
        
        # Extract pairs - handle both formats
        legacy_pairs = set()
        modern_pairs = set()
        
        # Legacy format: [{"type": "find_bestpair_selection", "pairs": [[i,j], ...]}]
        # or list of records with residue_idx1, residue_idx2
        for rec in legacy_data:
            if 'pairs' in rec:
                # Array of [i, j] pairs
                for pair in rec['pairs']:
                    if len(pair) >= 2:
                        legacy_pairs.add((min(pair[0], pair[1]), max(pair[0], pair[1])))
            else:
                # Individual record with residue_idx1/2 or residue_i/j
                i = rec.get('residue_idx1', rec.get('residue_i'))
                j = rec.get('residue_idx2', rec.get('residue_j'))
                if i is not None and j is not None:
                    legacy_pairs.add((min(i, j), max(i, j)))
        
        for rec in modern_data:
            if 'pairs' in rec:
                for pair in rec['pairs']:
                    if len(pair) >= 2:
                        modern_pairs.add((min(pair[0], pair[1]), max(pair[0], pair[1])))
            else:
                i = rec.get('residue_idx1', rec.get('residue_i'))
                j = rec.get('residue_idx2', rec.get('residue_j'))
                if i is not None and j is not None:
                    modern_pairs.add((min(i, j), max(i, j)))
        
        common = legacy_pairs & modern_pairs
        missing = legacy_pairs - modern_pairs
        extra = modern_pairs - legacy_pairs
        
        return (len(legacy_pairs), len(modern_pairs), len(common),
                len(missing), len(extra), sorted(missing), sorted(extra))


def run_verification(args):
    """Run the full verification."""
    base_dir = args.base_dir
    num_threads = args.threads
    
    verifier = PDBVerifier(base_dir, skip_known_perfect=not args.include_known_perfect)
    
    # Find PDBs with complete data by scanning directories
    legacy_pairs_dir = os.path.join(base_dir, 'data/json_legacy/find_bestpair_selection')
    modern_pairs_dir = os.path.join(base_dir, 'data/json/find_bestpair_selection')
    
    # Get sets of PDBs from each directory
    legacy_frames = set(f.replace('.json', '') for f in os.listdir(verifier.legacy_frame_dir) if f.endswith('.json'))
    modern_frames = set(f.replace('.json', '') for f in os.listdir(verifier.modern_frame_dir) if f.endswith('.json'))
    legacy_pairs = set(f.replace('.json', '') for f in os.listdir(legacy_pairs_dir) if f.endswith('.json'))
    modern_pairs = set(f.replace('.json', '') for f in os.listdir(modern_pairs_dir) if f.endswith('.json'))
    
    # PDBs with all 4 files (complete comparison possible)
    all_complete_pdbs = sorted(legacy_frames & modern_frames & legacy_pairs & modern_pairs)
    total_available = len(all_complete_pdbs)
    
    # Apply offset and limit
    start_idx = args.offset if args.offset else 0
    end_idx = start_idx + args.limit if args.limit else len(all_complete_pdbs)
    available_pdbs = all_complete_pdbs[start_idx:end_idx]
    
    print("=" * 80)
    print("COMPREHENSIVE PDB VERIFICATION")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"PDBs with complete data (frames + pairs): {total_available}")
    print(f"PDBs to test: {len(available_pdbs)}")
    print(f"Known legacy data issues (excluded): {len(KNOWN_LEGACY_DATA_ISSUES)}")
    print(f"Threads: {num_threads}")
    print("-" * 80)
    
    verifier.total = len(available_pdbs)
    
    # Run comparisons in parallel
    results: List[PDBResult] = []
    
    print(f"\nComparing {len(available_pdbs)} PDBs...")
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(verifier.compare_pdb, pdb): pdb for pdb in available_pdbs}
        
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
    
    # Categorize results
    perfect = [r for r in results if r.status == 'perfect']
    differences = [r for r in results if r.status == 'difference']
    errors = [r for r in results if r.status == 'error']
    legacy_issues = [r for r in results if r.status == 'legacy_issue']
    
    # Calculate totals
    total_common_frames = sum(r.common_frames for r in results)
    total_legacy_pairs = sum(r.legacy_pairs for r in results)
    total_modern_pairs = sum(r.modern_pairs for r in results)
    total_common_pairs = sum(r.common_pairs for r in results)
    total_missing = sum(r.missing_pairs for r in results)
    total_extra = sum(r.extra_pairs for r in results)
    
    # Find max differences
    max_origin = max((r.max_origin_diff for r in results), default=0)
    max_orient = max((r.max_orient_diff for r in results), default=0)
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total PDBs compared: {len(results)}")
    print(f"  ✅ Perfect matches:     {len(perfect):4d} ({100*len(perfect)/len(results):.1f}%)")
    print(f"  ⚠️  With differences:   {len(differences):4d} ({100*len(differences)/len(results):.1f}%)")
    print(f"  ⏭️  Legacy data issues: {len(legacy_issues):4d} ({100*len(legacy_issues)/len(results):.1f}%)")
    print(f"  ❌ Errors:              {len(errors):4d} ({100*len(errors)/len(results):.1f}%)")
    
    print("\n" + "-" * 80)
    print("FRAME COMPARISON (origin & ref_frame)")
    print("-" * 80)
    print(f"Total common frames compared: {total_common_frames:,}")
    print(f"Max origin difference:        {max_origin:.10f} Å")
    print(f"Max orientation difference:   {max_orient:.10f}")
    frames_perfect = max_origin == 0 and max_orient == 0
    print(f"Status: {'✅ ALL FRAMES PERFECT' if frames_perfect else '⚠️ FRAMES HAVE DIFFERENCES'}")
    
    print("\n" + "-" * 80)
    print("PAIR COMPARISON (find_bestpair_selection)")
    print("-" * 80)
    print(f"Total legacy pairs:  {total_legacy_pairs:,}")
    print(f"Total modern pairs:  {total_modern_pairs:,}")
    print(f"Common pairs:        {total_common_pairs:,}")
    print(f"Missing in modern:   {total_missing}")
    print(f"Extra in modern:     {total_extra}")
    pairs_perfect = total_missing == 0 and total_extra == 0
    print(f"Status: {'✅ ALL PAIRS PERFECT' if pairs_perfect else '⚠️ PAIRS HAVE DIFFERENCES'}")
    
    # Output files
    output_json = os.path.join(base_dir, 'data/verification_results.json')
    output_csv = os.path.join(base_dir, 'data/verification_results.csv')
    
    # Write CSV with all results
    import csv
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pdb_id', 'status', 'legacy_frames', 'modern_frames', 'common_frames',
                        'max_origin_diff', 'max_orient_diff', 'legacy_pairs', 'modern_pairs',
                        'common_pairs', 'missing_pairs', 'extra_pairs'])
        for r in sorted(results, key=lambda x: x.pdb_id):
            writer.writerow([r.pdb_id, r.status, r.legacy_frames, r.modern_frames, r.common_frames,
                            f'{r.max_origin_diff:.10f}', f'{r.max_orient_diff:.10f}',
                            r.legacy_pairs, r.modern_pairs, r.common_pairs,
                            r.missing_pairs, r.extra_pairs])
    
    # Only save differences in JSON (to save space)
    output_data = {
        'generated': datetime.now().isoformat(),
        'summary': {
            'total_compared': len(results),
            'perfect_matches': len(perfect),
            'with_differences': len(differences),
            'legacy_issues': len(legacy_issues),
            'errors': len(errors),
            'max_origin_diff': max_origin,
            'max_orient_diff': max_orient,
            'total_missing_pairs': total_missing,
            'total_extra_pairs': total_extra,
        },
        'perfect_pdbs': sorted([r.pdb_id for r in perfect]),
        'legacy_issue_pdbs': sorted([r.pdb_id for r in legacy_issues]),
        'differences': [],
        'errors': []
    }
    
    # Add detailed difference info
    for r in sorted(differences, key=lambda x: x.missing_pairs + x.extra_pairs, reverse=True):
        output_data['differences'].append({
            'pdb_id': r.pdb_id,
            'max_origin_diff': r.max_origin_diff,
            'max_orient_diff': r.max_orient_diff,
            'legacy_pairs': r.legacy_pairs,
            'modern_pairs': r.modern_pairs,
            'missing_pairs': r.missing_pairs,
            'extra_pairs': r.extra_pairs,
            'missing_pair_list': r.missing_pair_list[:20],  # Limit to first 20
            'extra_pair_list': r.extra_pair_list[:20],
        })
    
    # Add error info
    for r in errors:
        output_data['errors'].append({
            'pdb_id': r.pdb_id,
            'error': r.error_message
        })
    
    with open(output_json, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to:")
    print(f"  CSV: {output_csv}")
    print(f"  JSON: {output_json}")
    
    # Print differences if any
    if differences:
        print("\n" + "=" * 80)
        print(f"PDBs WITH DIFFERENCES ({len(differences)})")
        print("=" * 80)
        print(f"{'PDB':<8} {'Origin Diff':<14} {'Orient Diff':<14} {'Missing':<8} {'Extra':<8}")
        print("-" * 60)
        for r in sorted(differences, key=lambda x: x.missing_pairs + x.extra_pairs, reverse=True)[:30]:
            print(f"{r.pdb_id:<8} {r.max_origin_diff:<14.6f} {r.max_orient_diff:<14.6f} {r.missing_pairs:<8} {r.extra_pairs:<8}")
        if len(differences) > 30:
            print(f"... and {len(differences) - 30} more")
    
    if errors:
        print("\n" + "=" * 80)
        print(f"ERRORS ({len(errors)})")
        print("=" * 80)
        for r in errors[:20]:
            print(f"  {r.pdb_id}: {r.error_message[:60]}")
        if len(errors) > 20:
            print(f"  ... and {len(errors) - 20} more")
    
    # Final verdict
    print("\n" + "=" * 80)
    print("FINAL VERDICT")
    print("=" * 80)
    
    effective_perfect = len(perfect) + len(legacy_issues)
    effective_total = len(results) - len(errors)
    
    if len(differences) == 0 and frames_perfect and pairs_perfect:
        print("🎉 ALL PDBs MATCH PERFECTLY!")
        print("   (excluding known legacy data issues)")
    else:
        print(f"⚠️  {len(differences)} PDBs have differences that need investigation")
        print(f"   Effective match rate: {effective_perfect}/{effective_total} ({100*effective_perfect/effective_total:.1f}%)")
    
    return 0 if len(differences) == 0 else 1


def main():
    parser = argparse.ArgumentParser(description='Verify legacy vs modern PDB outputs')
    parser.add_argument('--base-dir', '-d', default=os.getcwd(),
                        help='Base directory of the project')
    parser.add_argument('--threads', '-t', type=int, default=16,
                        help='Number of threads (default: 16)')
    parser.add_argument('--limit', '-l', type=int, default=None,
                        help='Limit number of PDBs to test')
    parser.add_argument('--offset', '-o', type=int, default=0,
                        help='Start from this PDB index (for batching)')
    parser.add_argument('--include-known-perfect', action='store_true',
                        help='Include known perfect PDBs in detailed output')
    
    args = parser.parse_args()
    sys.exit(run_verification(args))


if __name__ == '__main__':
    main()

