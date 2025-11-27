#!/usr/bin/env python3
"""
Analyze mismatched pairs with detailed quality score and H-bond analysis.

This script extends the existing comparison infrastructure to provide
detailed analysis of mismatched pairs, including quality scores, H-bonds,
and root cause identification.
"""

import json
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import comparison functions directly to avoid dependency issues
try:
    from x3dna_json_compare.find_bestpair_comparison import compare_find_bestpair_selections
except ImportError:
    # Fallback implementation
    from typing import Set
    from dataclasses import dataclass, field
    
    @dataclass
    class FindBestpairComparison:
        missing_in_modern: List[tuple] = field(default_factory=list)
        extra_in_modern: List[tuple] = field(default_factory=list)
        total_legacy: int = 0
        total_modern: int = 0
        common_count: int = 0
    
    def compare_find_bestpair_selections(legacy_records, modern_records):
        result = FindBestpairComparison()
        legacy_pairs = set()
        modern_pairs = set()
        
        for rec in legacy_records:
            if rec.get('type') != 'find_bestpair_selection':
                continue
            for pair in rec.get('pairs', []):
                if isinstance(pair, list) and len(pair) >= 2:
                    legacy_pairs.add(normalize_pair(pair[0], pair[1]))
        
        for rec in modern_records:
            if rec.get('type') != 'find_bestpair_selection':
                continue
            for pair in rec.get('pairs', []):
                if isinstance(pair, list) and len(pair) >= 2:
                    modern_pairs.add(normalize_pair(pair[0], pair[1]))
        
        result.total_legacy = len(legacy_pairs)
        result.total_modern = len(modern_pairs)
        result.common_count = len(legacy_pairs & modern_pairs)
        result.missing_in_modern = sorted(legacy_pairs - modern_pairs)
        result.extra_in_modern = sorted(modern_pairs - legacy_pairs)
        return result


@dataclass
class PairAnalysis:
    """Analysis of a single mismatched pair."""
    residue1: int
    residue2: int
    is_missing: bool  # True if missing in modern, False if extra in modern
    
    # Validation status
    exists_in_validation_legacy: bool = False
    exists_in_validation_modern: bool = False
    is_valid_legacy: bool = False
    is_valid_modern: bool = False
    
    # Quality scores
    base_score_legacy: float = 0.0
    base_score_modern: float = 0.0
    adjust_pairQuality_legacy: float = 0.0
    adjust_pairQuality_modern: float = 0.0
    bp_type_id_legacy: int = 0
    bp_type_id_modern: int = 0
    final_score_legacy: float = 0.0
    final_score_modern: float = 0.0
    
    # H-bond info
    num_good_hb_legacy: int = 0
    num_good_hb_modern: int = 0
    num_total_hb_legacy: int = 0
    num_total_hb_modern: int = 0
    
    # Root cause analysis
    root_cause: str = ""
    priority: str = "MEDIUM"


def normalize_pair(i: int, j: int) -> Tuple[int, int]:
    """Normalize pair to (min, max)."""
    return (min(i, j), max(i, j))


def find_json_file(pdb_id: str, record_type: str, is_legacy: bool) -> Optional[Path]:
    """Find JSON file for a given PDB ID and record type."""
    base_dir = project_root / ("data/json_legacy" if is_legacy else "data/json")
    
    # Try segmented directory structure first
    segmented_file = base_dir / record_type / f"{pdb_id}.json"
    if segmented_file.exists():
        return segmented_file
    
    # Fall back to old format
    suffix_file = base_dir / f"{pdb_id}_{record_type}.json"
    if suffix_file.exists():
        return suffix_file
    
    return None


def load_json_records(pdb_id: str, record_type: str, is_legacy: bool) -> List[Dict]:
    """Load JSON records of a specific type."""
    json_file = find_json_file(pdb_id, record_type, is_legacy)
    if not json_file or not json_file.exists():
        return []
    
    with open(json_file) as f:
        data = json.load(f)
    
    if isinstance(data, list):
        return data
    elif isinstance(data, dict) and 'calculations' in data:
        calc = data['calculations']
        if isinstance(calc, list):
            return [r for r in calc if r.get('type') == record_type]
        elif isinstance(calc, dict) and record_type in calc:
            return calc[record_type]
    
    return []


def extract_quality_info(records: List[Dict], idx1: int, idx2: int) -> Dict:
    """Extract quality score information from validation records."""
    info = {
        'exists': False,
        'is_valid': False,
        'base_score': 0.0,
        'bp_type_id': 0,
        'dorg': 0.0,
        'd_v': 0.0,
        'plane_angle': 0.0,
        'dNN': 0.0,
    }
    
    pair_key = normalize_pair(idx1, idx2)
    
    for rec in records:
        if rec.get('type') != 'pair_validation':
            continue
        
        base_i = rec.get('base_i') or rec.get('residue1_idx')
        base_j = rec.get('base_j') or rec.get('residue2_idx')
        
        if base_i and base_j and normalize_pair(base_i, base_j) == pair_key:
            info['exists'] = True
            
            if rec.get('is_valid'):
                info['is_valid'] = bool(rec.get('is_valid'))
            
            calc = rec.get('calculated_values', {})
            if calc:
                info['base_score'] = calc.get('quality_score', 0.0)
                info['dorg'] = calc.get('dorg', 0.0)
                info['d_v'] = calc.get('d_v', 0.0)
                info['plane_angle'] = calc.get('plane_angle', 0.0)
                info['dNN'] = calc.get('dNN', 0.0)
            
            info['bp_type_id'] = rec.get('bp_type_id', 0)
            break
    
    return info


def extract_hbond_info(records: List[Dict], idx1: int, idx2: int) -> Dict:
    """Extract H-bond information."""
    info = {
        'num_total': 0,
        'num_good': 0,
    }
    
    pair_key = normalize_pair(idx1, idx2)
    
    for rec in records:
        if rec.get('type') != 'hbond_list':
            continue
        
        base_i = rec.get('base_i') or rec.get('residue1_idx')
        base_j = rec.get('base_j') or rec.get('residue2_idx')
        
        if base_i and base_j and normalize_pair(base_i, base_j) == pair_key:
            hbonds = rec.get('hbonds', [])
            info['num_total'] = len(hbonds)
            
            # Count good H-bonds (type='-' and distance in [2.5, 3.5])
            for hb in hbonds:
                hb_type = hb.get('type', ' ')
                dist = hb.get('distance', 0.0)
                if hb_type == '-' and 2.5 <= dist <= 3.5:
                    info['num_good'] += 1
            
            break
    
    return info


def calculate_adjust_pairQuality(num_good_hb: int) -> float:
    """Calculate adjust_pairQuality from number of good H-bonds."""
    return -3.0 * num_good_hb


def calculate_final_score(base_score: float, adjust: float, bp_type_id: int) -> float:
    """Calculate final quality score."""
    final = base_score + adjust
    if bp_type_id == 2:
        final -= 2.0
    return final


def analyze_pair(
    pdb_id: str,
    pair: Tuple[int, int],
    is_missing: bool,
    project_root: Path
) -> PairAnalysis:
    """Analyze a single mismatched pair."""
    idx1, idx2 = pair
    
    # Load validation and H-bond records
    legacy_validation = load_json_records(pdb_id, 'pair_validation', True)
    modern_validation = load_json_records(pdb_id, 'pair_validation', False)
    legacy_hbond = load_json_records(pdb_id, 'hbond_list', True)
    modern_hbond = load_json_records(pdb_id, 'hbond_list', False)
    
    # Extract information
    legacy_quality = extract_quality_info(legacy_validation, idx1, idx2)
    modern_quality = extract_quality_info(modern_validation, idx1, idx2)
    legacy_hb = extract_hbond_info(legacy_hbond, idx1, idx2)
    modern_hb = extract_hbond_info(modern_hbond, idx1, idx2)
    
    # Calculate adjust_pairQuality and final scores
    legacy_adjust = calculate_adjust_pairQuality(legacy_hb['num_good'])
    modern_adjust = calculate_adjust_pairQuality(modern_hb['num_good'])
    
    legacy_final = calculate_final_score(
        legacy_quality['base_score'],
        legacy_adjust,
        legacy_quality['bp_type_id']
    )
    modern_final = calculate_final_score(
        modern_quality['base_score'],
        modern_adjust,
        modern_quality['bp_type_id']
    )
    
    # Create analysis
    analysis = PairAnalysis(
        residue1=idx1,
        residue2=idx2,
        is_missing=is_missing,
        exists_in_validation_legacy=legacy_quality['exists'],
        exists_in_validation_modern=modern_quality['exists'],
        is_valid_legacy=legacy_quality['is_valid'],
        is_valid_modern=modern_quality['is_valid'],
        base_score_legacy=legacy_quality['base_score'],
        base_score_modern=modern_quality['base_score'],
        adjust_pairQuality_legacy=legacy_adjust,
        adjust_pairQuality_modern=modern_adjust,
        bp_type_id_legacy=legacy_quality['bp_type_id'],
        bp_type_id_modern=modern_quality['bp_type_id'],
        final_score_legacy=legacy_final,
        final_score_modern=modern_final,
        num_good_hb_legacy=legacy_hb['num_good'],
        num_good_hb_modern=modern_hb['num_good'],
        num_total_hb_legacy=legacy_hb['num_total'],
        num_total_hb_modern=modern_hb['num_total'],
    )
    
    # Determine root cause and priority
    if not analysis.exists_in_validation_modern:
        analysis.root_cause = "Pair not found in modern validation records"
        analysis.priority = "HIGH"
    elif not analysis.is_valid_modern:
        analysis.root_cause = "Pair exists but is invalid in modern (is_valid=0)"
        analysis.priority = "HIGH"
    else:
        base_diff = abs(analysis.base_score_legacy - analysis.base_score_modern)
        adjust_diff = abs(analysis.adjust_pairQuality_legacy - analysis.adjust_pairQuality_modern)
        final_diff = abs(analysis.final_score_legacy - analysis.final_score_modern)
        
        if base_diff > 0.001:
            analysis.root_cause = f"Base quality score difference: {base_diff:.6f}"
            analysis.priority = "HIGH"
        elif adjust_diff > 0.001:
            analysis.root_cause = f"adjust_pairQuality difference: {adjust_diff:.6f}"
            if analysis.num_good_hb_legacy != analysis.num_good_hb_modern:
                analysis.root_cause += (
                    f" (H-bond count difference: legacy={analysis.num_good_hb_legacy}, "
                    f"modern={analysis.num_good_hb_modern})"
                )
            analysis.priority = "HIGH"
        elif analysis.bp_type_id_legacy != analysis.bp_type_id_modern:
            analysis.root_cause = (
                f"bp_type_id difference: legacy={analysis.bp_type_id_legacy}, "
                f"modern={analysis.bp_type_id_modern}"
            )
            analysis.priority = "MEDIUM"
        elif final_diff > 0.001:
            analysis.root_cause = f"Final quality score difference: {final_diff:.6f}"
            analysis.priority = "HIGH"
        else:
            analysis.root_cause = "Quality scores match - likely tie-breaking or iteration order issue"
            analysis.priority = "LOW"
    
    return analysis


def generate_report(
    pdb_id: str,
    missing_analyses: List[PairAnalysis],
    extra_analyses: List[PairAnalysis],
    output_file: Optional[Path] = None
):
    """Generate markdown report."""
    from datetime import datetime
    
    lines = [
        f"# Mismatched Pairs Analysis: {pdb_id}",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Summary",
        "",
        f"- **Missing Pairs (in legacy but not modern)**: {len(missing_analyses)}",
        f"- **Extra Pairs (in modern but not legacy)**: {len(extra_analyses)}",
        "",
    ]
    
    if missing_analyses:
        lines.extend([
            f"## Missing Pairs ({len(missing_analyses)})",
            "",
        ])
        
        for i, analysis in enumerate(missing_analyses, 1):
            lines.extend([
                f"### Pair ({analysis.residue1}, {analysis.residue2})",
                "",
                f"- **Status**: ",
            ])
            
            if analysis.exists_in_validation_modern:
                lines[-1] += f"✅ EXISTS in validation (is_valid={1 if analysis.is_valid_modern else 0})"
            else:
                lines[-1] += "❌ NOT FOUND in validation"
            
            lines.extend([
                f"- **Priority**: {analysis.priority}",
                f"- **Root Cause**: {analysis.root_cause}",
                "",
            ])
            
            if analysis.exists_in_validation_modern:
                lines.extend([
                    "**Quality Scores**:",
                    f"- Base Score: Legacy={analysis.base_score_legacy:.6f}, "
                    f"Modern={analysis.base_score_modern:.6f} "
                    f"(diff={abs(analysis.base_score_legacy - analysis.base_score_modern):.6f})",
                    f"- adjust_pairQuality: Legacy={analysis.adjust_pairQuality_legacy:.6f}, "
                    f"Modern={analysis.adjust_pairQuality_modern:.6f} "
                    f"(diff={abs(analysis.adjust_pairQuality_legacy - analysis.adjust_pairQuality_modern):.6f})",
                    f"- bp_type_id: Legacy={analysis.bp_type_id_legacy}, Modern={analysis.bp_type_id_modern}",
                    f"- Final Score: Legacy={analysis.final_score_legacy:.6f}, "
                    f"Modern={analysis.final_score_modern:.6f} "
                    f"(diff={abs(analysis.final_score_legacy - analysis.final_score_modern):.6f})",
                    "",
                    "**H-bonds**:",
                    f"- Good H-bonds: Legacy={analysis.num_good_hb_legacy}, Modern={analysis.num_good_hb_modern}",
                    f"- Total H-bonds: Legacy={analysis.num_total_hb_legacy}, Modern={analysis.num_total_hb_modern}",
                    "",
                ])
                
                if (abs(analysis.base_score_legacy - analysis.base_score_modern) > 0.001 or
                    abs(analysis.adjust_pairQuality_legacy - analysis.adjust_pairQuality_modern) > 0.001 or
                    analysis.bp_type_id_legacy != analysis.bp_type_id_modern):
                    lines.append("⚠️ **QUALITY SCORE DIFFERENCE DETECTED**")
                    lines.append("")
    
    if extra_analyses:
        lines.extend([
            "",
            f"## Extra Pairs ({len(extra_analyses)})",
            "",
        ])
        
        for i, analysis in enumerate(extra_analyses, 1):
            lines.extend([
                f"### Pair ({analysis.residue1}, {analysis.residue2})",
                "",
                f"- **Status**: ",
            ])
            
            if analysis.exists_in_validation_legacy:
                lines[-1] += f"✅ EXISTS in legacy validation (is_valid={1 if analysis.is_valid_legacy else 0})"
            else:
                lines[-1] += "❌ NOT FOUND in legacy validation"
            
            lines.extend([
                f"- **Priority**: {analysis.priority}",
                f"- **Root Cause**: {analysis.root_cause}",
                "",
            ])
    
    report_text = "\n".join(lines)
    
    if output_file:
        output_file.write_text(report_text)
        print(f"Report written to: {output_file}")
    else:
        print(report_text)


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/analyze_mismatched_pairs.py <pdb_id> [--output <file>]")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    output_file = None
    
    if '--output' in sys.argv:
        idx = sys.argv.index('--output')
        if idx + 1 < len(sys.argv):
            output_file = Path(sys.argv[idx + 1])
    
    # Check if modern JSON exists, generate if missing
    modern_json_file = find_json_file(pdb_id, 'find_bestpair_selection', False)
    pdb_file = project_root / "data/pdb" / f"{pdb_id}.pdb"
    
    if not modern_json_file or not modern_json_file.exists():
        if pdb_file.exists():
            print(f"Modern JSON not found for {pdb_id}. Generating...")
            try:
                result = subprocess.run(
                    ['build/generate_modern_json', str(pdb_file), str(project_root / 'data/json')],
                    capture_output=True,
                    text=True,
                    timeout=600
                )
                if result.returncode != 0:
                    print(f"Warning: Failed to generate modern JSON: {result.stderr[:200]}")
            except subprocess.TimeoutExpired:
                print(f"Warning: Timeout generating modern JSON for {pdb_id}")
            except Exception as e:
                print(f"Warning: Error generating modern JSON: {e}")
        else:
            print(f"Error: PDB file not found: {pdb_file}")
            sys.exit(1)
    
    # Use existing comparison infrastructure
    legacy_records = load_json_records(pdb_id, 'find_bestpair_selection', True)
    modern_records = load_json_records(pdb_id, 'find_bestpair_selection', False)
    
    if not legacy_records and not modern_records:
        print(f"Error: No find_bestpair_selection records found for {pdb_id}")
        sys.exit(1)
    
    fbc = compare_find_bestpair_selections(legacy_records, modern_records)
    
    if not fbc.missing_in_modern and not fbc.extra_in_modern:
        print(f"✅ Perfect match! No mismatched pairs for {pdb_id}")
        sys.exit(0)
    
    print(f"Found {len(fbc.missing_in_modern)} missing pairs and {len(fbc.extra_in_modern)} extra pairs")
    
    # Analyze missing pairs
    missing_analyses = []
    for pair in fbc.missing_in_modern:
        analysis = analyze_pair(pdb_id, pair, True, project_root)
        missing_analyses.append(analysis)
    
    # Analyze extra pairs
    extra_analyses = []
    for pair in fbc.extra_in_modern:
        analysis = analyze_pair(pdb_id, pair, False, project_root)
        extra_analyses.append(analysis)
    
    # Sort by priority
    missing_analyses.sort(key=lambda x: {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2}[x.priority])
    extra_analyses.sort(key=lambda x: {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2}[x.priority])
    
    # Generate report
    if output_file is None:
        output_file = Path(f"{pdb_id}_mismatched_pairs_analysis.md")
    
    generate_report(pdb_id, missing_analyses, extra_analyses, output_file)


if __name__ == '__main__':
    main()

