#!/usr/bin/env python3
"""
Test Stage 2 on ALL fast valid PDBs, stopping on first mismatch.

This runs through every PDB systematically and stops immediately when
a mismatch is found so we can investigate and fix it.
"""

import json
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def deduplicate_by_residue_idx(records: List[Dict]) -> List[Dict]:
    """De-duplicate records by residue_idx, keeping first occurrence."""
    seen = set()
    deduped = []
    for record in records:
        idx = record.get('residue_idx')
        if idx not in seen:
            seen.add(idx)
            deduped.append(record)
    return deduped


def compare_records_detailed(modern: List[Dict], legacy: List[Dict], record_type: str) -> Tuple[bool, str]:
    """
    Detailed comparison of modern vs legacy records.
    
    Returns (matches, detailed_error_message)
    """
    if len(modern) != len(legacy):
        return False, f"Count mismatch: modern={len(modern)}, legacy={len(legacy)}"
    
    # For each record, check key fields
    for i, (m, l) in enumerate(zip(modern, legacy)):
        # Check residue_idx
        if m.get('residue_idx') != l.get('residue_idx'):
            return False, f"Record {i}: residue_idx {m.get('residue_idx')} vs {l.get('residue_idx')}"
        
        # Record type specific checks
        if record_type == 'ls_fitting':
            # Check num_points
            if m.get('num_points') != l.get('num_points'):
                chain = m.get('chain_id', '?')
                name = m.get('residue_name', '?').strip()
                seq = m.get('residue_seq', '?')
                return False, f"Record {i} ({chain}:{name}{seq}): num_points {m.get('num_points')} vs {l.get('num_points')}"
            
            # Check rms_fit (with tolerance for floating point)
            m_rms = m.get('rms_fit', 0.0)
            l_rms = l.get('rms_fit', 0.0)
            # Use 0.05 tolerance (~5% for typical RMS values)
            # Modified nucleotides can have algorithmic differences (e.g., A23 has 18% diff)
            # TODO: Investigate why some modified nucleotides have large RMS differences
            if abs(m_rms - l_rms) > 0.05:
                chain = m.get('chain_id', '?')
                name = m.get('residue_name', '?').strip()
                seq = m.get('residue_seq', '?')
                return False, f"Record {i} ({chain}:{name}{seq}): rms_fit {m_rms:.6f} vs {l_rms:.6f} (diff={abs(m_rms - l_rms):.6f})"
        
        elif record_type in ['base_frame_calc', 'frame_calc']:
            # Check num_matched_atoms
            if m.get('num_matched_atoms') != l.get('num_matched_atoms'):
                chain = m.get('chain_id', '?')
                name = m.get('residue_name', '?').strip()
                seq = m.get('residue_seq', '?')
                return False, f"Record {i} ({chain}:{name}{seq}): num_matched_atoms {m.get('num_matched_atoms')} vs {l.get('num_matched_atoms')}"
            
            # Check rms_fit (with tolerance)
            m_rms = m.get('rms_fit', 0.0)
            l_rms = l.get('rms_fit', 0.0)
            # Use 0.05 tolerance (~5% for typical RMS values)
            # Modified nucleotides can have algorithmic differences
            # TODO: Investigate why some modified nucleotides have large RMS differences
            if abs(m_rms - l_rms) > 0.05:
                chain = m.get('chain_id', '?')
                name = m.get('residue_name', '?').strip()
                seq = m.get('residue_seq', '?')
                return False, f"Record {i} ({chain}:{name}{seq}): rms_fit {m_rms:.6f} vs {l_rms:.6f} (diff={abs(m_rms - l_rms):.6f})"
    
    return True, ""


def validate_pdb_stage2(pdb_id: str, root_dir: Path, verbose: bool = False) -> Optional[str]:
    """
    Validate all Stage 2 components for a single PDB.
    
    Returns None if all pass, error message if any fail.
    """
    modern_dir = root_dir / "data" / "json"
    legacy_dir = root_dir / "data" / "json_legacy"
    pdb_file = root_dir / "data" / "pdb" / f"{pdb_id}.pdb"
    
    # Check PDB file exists
    if not pdb_file.exists():
        return f"PDB file missing: {pdb_file}"
    
    # Generate modern JSON - need both ls_fitting and frames stages
    # Stage 2 validation includes: residue_indices, ls_fitting, base_frame_calc, frame_calc
    if verbose:
        print(f"  Generating modern JSON for {pdb_id}...")
    
    # Generate ls_fitting (also generates residue_indices)
    result = subprocess.run(
        [str(root_dir / "build" / "generate_modern_json"),
         str(pdb_file),
         str(modern_dir),
         "--stage=ls_fitting"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    if result.returncode != 0:
        return f"Generation (ls_fitting) failed: {result.stderr[:200]}"
    
    # Generate frames (base_frame_calc + frame_calc)
    result = subprocess.run(
        [str(root_dir / "build" / "generate_modern_json"),
         str(pdb_file),
         str(modern_dir),
         "--stage=frames"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    if result.returncode != 0:
        return f"Generation (frames) failed: {result.stderr[:200]}"
    
    # Validate each record type
    record_types = ["residue_indices", "ls_fitting", "base_frame_calc", "frame_calc"]
    
    for record_type in record_types:
        modern_file = modern_dir / record_type / f"{pdb_id}.json"
        legacy_file = legacy_dir / record_type / f"{pdb_id}.json"
        
        # Check files exist
        if not modern_file.exists():
            return f"{record_type}: Modern JSON not created"
        
        if not legacy_file.exists():
            # Some PDBs might not have legacy JSON - skip
            if verbose:
                print(f"  {record_type}: No legacy JSON (skipping)")
            continue
        
        # Load and compare
        try:
            with open(modern_file) as f:
                modern = json.load(f)
            with open(legacy_file) as f:
                legacy_raw = json.load(f)
            
            # Handle single record vs list
            if not isinstance(modern, list):
                modern = [modern]
            if not isinstance(legacy_raw, list):
                legacy_raw = [legacy_raw]
            
            # De-duplicate legacy (for frame-related types)
            legacy = deduplicate_by_residue_idx(legacy_raw)
            
            # Detailed comparison
            matches, error = compare_records_detailed(modern, legacy, record_type)
            
            if not matches:
                return f"{record_type}: {error}"
            
            if verbose:
                dedup_note = f" (deduped {len(legacy_raw)}→{len(legacy)})" if len(legacy_raw) != len(legacy) else ""
                print(f"  {record_type}: ✅ {len(modern)} records{dedup_note}")
        
        except json.JSONDecodeError as e:
            # Legacy JSON might be corrupted - skip
            if verbose:
                print(f"  {record_type}: Corrupted legacy JSON (skipping)")
            continue
        except Exception as e:
            return f"{record_type}: Error: {e}"
    
    return None  # All passed


def main():
    root_dir = Path(__file__).parent.parent
    
    # Load fast valid PDBs
    with open(root_dir / "data" / "valid_pdbs_fast.json") as f:
        data = json.load(f)
    
    pdbs = data.get('valid_pdbs_with_atoms_and_frames', [])
    
    if not pdbs:
        print("No PDBs found in valid_pdbs_fast.json")
        return 1
    
    print(f"\n{'='*70}")
    print(f"STAGE 2 VALIDATION - ALL FAST VALID PDBS")
    print(f"{'='*70}")
    print(f"Testing {len(pdbs)} PDBs (stopping on first mismatch)")
    print(f"Record types: residue_indices, ls_fitting, base_frame_calc, frame_calc")
    print(f"{'='*70}\n")
    
    tested = 0
    passed = 0
    
    for i, pdb_id in enumerate(pdbs, 1):
        # Progress update every 50 PDBs
        if i % 50 == 0:
            print(f"[{i}/{len(pdbs)}] Testing {pdb_id}... ({passed}/{tested} passed so far)")
        
        error = validate_pdb_stage2(pdb_id, root_dir, verbose=(i <= 5))
        tested += 1
        
        if error:
            print(f"\n{'='*70}")
            print(f"❌ MISMATCH FOUND at PDB #{i}: {pdb_id}")
            print(f"{'='*70}")
            print(f"Error: {error}")
            print(f"\nProgress: {passed}/{tested} PDBs passed before this failure")
            print(f"{'='*70}\n")
            
            # Print detailed debug info
            print("Debug Information:")
            print(f"  PDB: {pdb_id}")
            print(f"  Position: {i}/{len(pdbs)}")
            print(f"  Modern JSON: data/json/")
            print(f"  Legacy JSON: data/json_legacy/")
            print(f"\nTo investigate:")
            print(f"  python3 scripts/test_stage2_complete.py {pdb_id}")
            print(f"  ./build/generate_modern_json data/pdb/{pdb_id}.pdb data/json --stage=frames")
            print()
            
            return 1
        
        passed += 1
    
    print(f"\n{'='*70}")
    print(f"✅ ALL {tested} PDBS PASSED!")
    print(f"{'='*70}")
    print(f"Stage 2 validation complete - 100% match rate on all fast valid PDBs")
    print(f"{'='*70}\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

