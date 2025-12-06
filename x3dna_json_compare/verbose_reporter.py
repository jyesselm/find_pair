"""
Verbose comparison reporter for detailed field-by-field analysis.

Generates human-readable detailed comparison reports showing:
- Field-by-field comparisons
- Diff highlighting
- JSON source paths
- Related record lookups
- Calculation provenance (when available)
"""

from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime
import json


@dataclass
class FieldComparison:
    """Result of comparing a single field."""
    field_name: str
    legacy_value: Any
    modern_value: Any
    matches: bool
    diff: Optional[float] = None
    tolerance: Optional[float] = None
    notes: List[str] = None
    
    def __post_init__(self):
        if self.notes is None:
            self.notes = []


@dataclass
class RecordComparison:
    """Result of comparing a single record (e.g., one base pair)."""
    record_key: Any  # e.g., (base_i, base_j) or residue_key
    record_type: str  # e.g., "distance_checks", "base_pair"
    matches: bool
    field_comparisons: List[FieldComparison]
    legacy_source: Optional[str] = None
    modern_source: Optional[str] = None
    related_records: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.related_records is None:
            self.related_records = {}


class VerboseReporter:
    """Generate detailed verbose comparison reports."""
    
    def __init__(
        self,
        tolerance: float = 1e-6,
        show_provenance: bool = False,
        show_related: bool = True,
        max_mismatches_per_stage: int = 20,
        diff_only: bool = False
    ):
        """
        Initialize verbose reporter.
        
        Args:
            tolerance: Default numerical tolerance
            show_provenance: Show calculation source info (if available)
            show_related: Show related records for context
            max_mismatches_per_stage: Maximum mismatches to show per stage
            diff_only: Only show records with differences
        """
        self.tolerance = tolerance
        self.show_provenance = show_provenance
        self.show_related = show_related
        self.max_mismatches_per_stage = max_mismatches_per_stage
        self.diff_only = diff_only
        self.sections: List[str] = []
    
    def add_header(self, pdb_id: str, stages: List[str] = None):
        """Add report header."""
        lines = []
        lines.append("=" * 80)
        lines.append(f"VERBOSE COMPARISON: {pdb_id}")
        lines.append("=" * 80)
        lines.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Tolerance: {self.tolerance}")
        if stages:
            lines.append(f"Stages: {', '.join(stages)}")
        lines.append(f"Mode: {'Differences only' if self.diff_only else 'All records'}")
        lines.append("")
        
        self.sections.append("\n".join(lines))
    
    def add_stage_header(self, stage_name: str, stage_number: int = None):
        """Add stage section header."""
        lines = []
        lines.append("-" * 80)
        if stage_number is not None:
            lines.append(f"STAGE {stage_number}: {stage_name}")
        else:
            lines.append(f"STAGE: {stage_name}")
        lines.append("-" * 80)
        
        self.sections.append("\n".join(lines))
    
    def add_stage_summary(
        self,
        total_legacy: int,
        total_modern: int,
        common: int,
        missing_in_modern: int,
        extra_in_modern: int,
        mismatches: int
    ):
        """Add stage-level summary statistics."""
        lines = []
        lines.append(f"Total legacy records: {total_legacy}")
        lines.append(f"Total modern records: {total_modern}")
        lines.append(f"Common records: {common}")
        
        if missing_in_modern > 0:
            lines.append(f"⚠️  Missing in modern: {missing_in_modern}")
        if extra_in_modern > 0:
            lines.append(f"⚠️  Extra in modern: {extra_in_modern}")
        if mismatches > 0:
            lines.append(f"❌ Mismatched records: {mismatches}")
        
        if missing_in_modern == 0 and extra_in_modern == 0 and mismatches == 0:
            lines.append("✅ All records match perfectly")
        
        lines.append("")
        
        self.sections.append("\n".join(lines))
    
    def add_record_comparison(self, comparison: RecordComparison):
        """Add detailed record comparison."""
        if self.diff_only and comparison.matches:
            return  # Skip matching records in diff-only mode
        
        lines = []
        
        # Record header
        status = "✅ MATCH" if comparison.matches else "❌ MISMATCH"
        lines.append(f"{status} {self._format_record_key(comparison.record_key, comparison.record_type)}")
        
        # Field comparisons
        if comparison.field_comparisons:
            lines.append("")
            lines.append("  Fields:")
            
            for field_comp in comparison.field_comparisons:
                field_line = self._format_field_comparison(field_comp)
                lines.append(f"    {field_line}")
                
                # Add notes if any
                for note in field_comp.notes:
                    lines.append(f"      ℹ️  {note}")
        
        # Related records (if enabled and available)
        if self.show_related and comparison.related_records:
            lines.append("")
            lines.append("  Related records:")
            for record_type, info in comparison.related_records.items():
                lines.append(f"    └─ {record_type}: {self._format_related_info(info)}")
        
        lines.append("")
        
        self.sections.append("\n".join(lines))
    
    def add_summary(
        self,
        stages_compared: int,
        perfect_matches: int,
        stages_with_diffs: int,
        diff_details: List[str]
    ):
        """Add final summary."""
        lines = []
        lines.append("-" * 80)
        lines.append("SUMMARY")
        lines.append("-" * 80)
        lines.append(f"Stages compared: {stages_compared}")
        lines.append(f"Perfect matches: {perfect_matches}")
        lines.append(f"Stages with differences: {stages_with_diffs}")
        
        if diff_details:
            lines.append("")
            lines.append("Differences found:")
            for detail in diff_details:
                lines.append(f"  - {detail}")
        
        lines.append("")
        
        if stages_with_diffs == 0:
            lines.append("Overall: ✅ ALL STAGES MATCH PERFECTLY")
        else:
            lines.append("Overall: ⚠️  DIFFERENCES FOUND")
        
        lines.append("")
        
        self.sections.append("\n".join(lines))
    
    def generate_report(self) -> str:
        """Generate final verbose report."""
        return "\n".join(self.sections)
    
    def _format_record_key(self, key: Any, record_type: str) -> str:
        """Format record key for display."""
        if isinstance(key, tuple):
            if len(key) == 2 and all(isinstance(x, int) for x in key):
                # Pair key: (base_i, base_j)
                return f"(base_i={key[0]}, base_j={key[1]})"
            elif len(key) == 3:
                # Residue key: (chain_id, residue_seq, insertion)
                chain_id, residue_seq, insertion = key
                ins_str = f":{insertion}" if insertion and insertion != ' ' else ""
                return f"(chain {chain_id}, seq {residue_seq}{ins_str})"
        
        return str(key)
    
    def _format_field_comparison(self, field: FieldComparison) -> str:
        """Format a single field comparison."""
        name_width = 20
        value_width = 15
        
        # Format field name
        name_str = f"{field.field_name}:".ljust(name_width)
        
        # Format values
        legacy_str = self._format_value(field.legacy_value, value_width)
        modern_str = self._format_value(field.modern_value, value_width)
        
        # Status indicator
        if field.matches:
            status = "✓"
            comparison = "=="
        else:
            status = "✗"
            comparison = "vs"
        
        result = f"{name_str} {legacy_str} {comparison} {modern_str} {status}"
        
        # Add diff info for mismatches
        if not field.matches and field.diff is not None:
            if field.tolerance is not None:
                result += f" (diff: {field.diff:.6e}, tolerance: {field.tolerance:.6e})"
            else:
                result += f" (diff: {field.diff:.6e})"
        
        return result
    
    def _format_value(self, value: Any, width: int = 15) -> str:
        """Format a value for display."""
        if value is None:
            return "None".ljust(width)
        elif isinstance(value, float):
            # Format float with appropriate precision
            if abs(value) < 1e-10:
                return "0.0".ljust(width)
            elif abs(value) < 0.001 or abs(value) > 1e6:
                return f"{value:.6e}".ljust(width)
            else:
                return f"{value:.6f}".ljust(width)
        elif isinstance(value, (list, tuple)) and len(value) <= 3:
            # Format short arrays (e.g., xyz coordinates)
            formatted = "[" + ", ".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in value) + "]"
            return formatted.ljust(width) if len(formatted) <= width else formatted
        elif isinstance(value, str):
            # Truncate long strings
            if len(value) > width:
                return value[:width-3] + "..."
            return value.ljust(width)
        else:
            # Default string representation
            s = str(value)
            if len(s) > width:
                return s[:width-3] + "..."
            return s.ljust(width)
    
    def _format_related_info(self, info: Any) -> str:
        """Format related record information."""
        if isinstance(info, dict):
            # Show key fields
            parts = []
            for key, value in list(info.items())[:3]:  # Limit to 3 fields
                parts.append(f"{key}={self._format_value(value, 10).strip()}")
            return ", ".join(parts)
        return str(info)


def compare_values_verbose(
    field_name: str,
    legacy_val: Any,
    modern_val: Any,
    tolerance: float = 1e-6
) -> FieldComparison:
    """
    Compare two values and return detailed FieldComparison.
    
    Args:
        field_name: Name of the field being compared
        legacy_val: Legacy value
        modern_val: Modern value
        tolerance: Numerical tolerance for floats
    
    Returns:
        FieldComparison with detailed results
    """
    notes = []
    
    # Handle None values
    if legacy_val is None and modern_val is None:
        return FieldComparison(
            field_name=field_name,
            legacy_value=legacy_val,
            modern_value=modern_val,
            matches=True
        )
    
    if legacy_val is None or modern_val is None:
        notes.append(f"One value is None")
        return FieldComparison(
            field_name=field_name,
            legacy_value=legacy_val,
            modern_value=modern_val,
            matches=False,
            notes=notes
        )
    
    # Numerical comparison
    if isinstance(legacy_val, (int, float)) and isinstance(modern_val, (int, float)):
        diff = abs(float(legacy_val) - float(modern_val))
        matches = diff <= tolerance
        
        if not matches and diff > 0:
            notes.append(f"Exceeds tolerance by {diff - tolerance:.6e}")
        
        return FieldComparison(
            field_name=field_name,
            legacy_value=legacy_val,
            modern_value=modern_val,
            matches=matches,
            diff=diff,
            tolerance=tolerance,
            notes=notes
        )
    
    # String comparison
    if isinstance(legacy_val, str) and isinstance(modern_val, str):
        matches = legacy_val.strip() == modern_val.strip()
        return FieldComparison(
            field_name=field_name,
            legacy_value=legacy_val,
            modern_value=modern_val,
            matches=matches,
            notes=notes
        )
    
    # List/array comparison
    if isinstance(legacy_val, (list, tuple)) and isinstance(modern_val, (list, tuple)):
        if len(legacy_val) != len(modern_val):
            notes.append(f"Length mismatch: {len(legacy_val)} vs {len(modern_val)}")
            return FieldComparison(
                field_name=field_name,
                legacy_value=legacy_val,
                modern_value=modern_val,
                matches=False,
                notes=notes
            )
        
        # Element-wise comparison
        all_match = True
        max_diff = 0.0
        
        for i, (l, m) in enumerate(zip(legacy_val, modern_val)):
            if isinstance(l, (int, float)) and isinstance(m, (int, float)):
                diff = abs(float(l) - float(m))
                max_diff = max(max_diff, diff)
                if diff > tolerance:
                    all_match = False
                    notes.append(f"Element {i}: diff={diff:.6e} > tolerance={tolerance:.6e}")
            elif l != m:
                all_match = False
                notes.append(f"Element {i}: {l} != {m}")
        
        return FieldComparison(
            field_name=field_name,
            legacy_value=legacy_val,
            modern_value=modern_val,
            matches=all_match,
            diff=max_diff if max_diff > 0 else None,
            tolerance=tolerance,
            notes=notes
        )
    
    # Direct equality for other types
    matches = legacy_val == modern_val
    if not matches:
        notes.append(f"Type mismatch or unequal: {type(legacy_val).__name__} vs {type(modern_val).__name__}")
    
    return FieldComparison(
        field_name=field_name,
        legacy_value=legacy_val,
        modern_value=modern_val,
        matches=matches,
        notes=notes
    )


def create_record_comparison_from_dicts(
    record_key: Any,
    record_type: str,
    legacy_record: Dict,
    modern_record: Dict,
    fields_to_compare: List[str],
    tolerance: float = 1e-6,
    legacy_source: Optional[str] = None,
    modern_source: Optional[str] = None,
    related_records: Optional[Dict[str, Any]] = None
) -> RecordComparison:
    """
    Create a RecordComparison by comparing specified fields in two record dicts.
    
    Args:
        record_key: Key identifying the record (e.g., pair tuple)
        record_type: Type of record (e.g., "distance_checks")
        legacy_record: Legacy record dictionary
        modern_record: Modern record dictionary
        fields_to_compare: List of field names to compare
        tolerance: Numerical tolerance
        legacy_source: Source path for legacy record
        modern_source: Source path for modern record
        related_records: Related records for context
    
    Returns:
        RecordComparison with field-by-field results
    """
    field_comparisons = []
    all_match = True
    
    for field_name in fields_to_compare:
        legacy_val = legacy_record.get(field_name)
        modern_val = modern_record.get(field_name)
        
        field_comp = compare_values_verbose(field_name, legacy_val, modern_val, tolerance)
        field_comparisons.append(field_comp)
        
        if not field_comp.matches:
            all_match = False
    
    return RecordComparison(
        record_key=record_key,
        record_type=record_type,
        matches=all_match,
        field_comparisons=field_comparisons,
        legacy_source=legacy_source,
        modern_source=modern_source,
        related_records=related_records
    )


def generate_full_verbose_report(pdb_id: str, result, tolerance: float = 1e-6, 
                                  diff_only: bool = False) -> str:
    """Generate a complete verbose report for a comparison result.
    
    Args:
        pdb_id: PDB identifier
        result: ComparisonResult object
        tolerance: Numerical comparison tolerance
        diff_only: If True, only show mismatches. If False, show all records.
    
    Returns:
        Formatted verbose report string
    """
    reporter = VerboseReporter(
        tolerance=tolerance,
        show_provenance=False,
        show_related=True,
        diff_only=diff_only
    )
    
    # Build stages list
    stages = []
    if result.atom_comparison:
        stages.append("atoms")
    if result.frame_comparison:
        stages.append("frames")
    if result.step_comparison:
        stages.append("steps")
    if result.helical_comparison:
        stages.append("helical")
    if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
        stages.append("distance_checks")
    if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
        stages.append("hbond_list")
    if hasattr(result, 'pair_validation_comparison') and result.pair_validation_comparison:
        stages.append("pair_validation")
    if hasattr(result, 'find_bestpair_comparison') and result.find_bestpair_comparison:
        stages.append("find_bestpair_selection")
    if hasattr(result, 'base_pair_comparison') and result.base_pair_comparison:
        stages.append("base_pair")
    
    reporter.add_header(pdb_id, stages)
    
    # Add atoms - load JSON to show ALL records
    if result.atom_comparison:
        import json
        from pathlib import Path
        
        ac = result.atom_comparison
        reporter.add_stage_header("atoms", 1)
        reporter.add_stage_summary(
            ac.total_legacy, ac.total_modern,
            ac.common_count,
            len(ac.missing_in_modern), len(ac.extra_in_modern),
            len(ac.mismatched_fields)
        )
        
        # Load all atom records to show every record
        if not diff_only:
            legacy_file = Path(f"data/json_legacy/pdb_atoms/{pdb_id}.json")
            modern_file = Path(f"data/json/pdb_atoms/{pdb_id}.json")
            
            if legacy_file.exists() and modern_file.exists():
                try:
                    with open(legacy_file) as f:
                        legacy_data = json.load(f)
                    with open(modern_file) as f:
                        modern_data = json.load(f)
                    
                    # Extract atoms array from the wrapper structure
                    # Legacy: [{num_atoms: ..., atoms: [...]}]
                    # Modern: {atoms: [...]}
                    if isinstance(legacy_data, list) and len(legacy_data) > 0:
                        legacy_atoms = legacy_data[0].get('atoms', [])
                    elif isinstance(legacy_data, dict):
                        legacy_atoms = legacy_data.get('atoms', [])
                    else:
                        legacy_atoms = []
                    
                    if isinstance(modern_data, list) and len(modern_data) > 0:
                        modern_atoms = modern_data[0].get('atoms', [])
                    elif isinstance(modern_data, dict):
                        modern_atoms = modern_data.get('atoms', [])
                    else:
                        modern_atoms = []
                    
                    # Create dict for quick modern lookup
                    # Atoms are identified by atom_idx
                    modern_dict = {}
                    for rec in modern_atoms:
                        atom_idx = rec.get('atom_idx')
                        if atom_idx:
                            modern_dict[atom_idx] = rec
                    
                    # Compare all legacy records
                    shown_count = 0
                    max_to_show = 1000
                    
                    for legacy_rec in legacy_atoms:
                        if shown_count >= max_to_show:
                            break
                        
                        atom_idx = legacy_rec.get('atom_idx')
                        
                        # Check if modern has this atom
                        if atom_idx not in modern_dict:
                            # Add a warning about missing atom
                            reporter.sections.append(
                                f"⚠️  WARNING: atom_{atom_idx} - "
                                f"Legacy atom (atom_idx={atom_idx}) not found in modern output!\n"
                            )
                            shown_count += 1
                            continue
                        
                        modern_rec = modern_dict.get(atom_idx)
                        
                        # Fields to compare for atoms
                        fields = ['atom_name', 'residue_name', 'chain_id', 'residue_seq', 'xyz']
                        
                        rec_comp = create_record_comparison_from_dicts(
                            record_key=f"atom_{atom_idx}",
                            record_type="pdb_atoms",
                            legacy_record=legacy_rec,
                            modern_record=modern_rec,
                            fields_to_compare=fields,
                            tolerance=tolerance,
                            legacy_source=str(legacy_file),
                            modern_source=str(modern_file)
                        )
                        reporter.add_record_comparison(rec_comp)
                        shown_count += 1
                except Exception as e:
                    # If loading fails, just continue without atom details
                    pass
    
    # Add distance checks - load JSON to show ALL records
    if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
        import json
        from pathlib import Path
        
        dcc = result.distance_checks_comparison
        reporter.add_stage_header("distance_checks", 3)
        reporter.add_stage_summary(
            dcc.total_legacy, dcc.total_modern,
            dcc.total_legacy - len(dcc.missing_in_modern),
            len(dcc.missing_in_modern), len(dcc.extra_in_modern),
            len(dcc.mismatched_checks)
        )
        
        # Load all distance check records to show every record
        if not diff_only:
            legacy_file = Path(f"data/json_legacy/distance_checks/{pdb_id}.json")
            modern_file = Path(f"data/json/distance_checks/{pdb_id}.json")
            
            if legacy_file.exists() and modern_file.exists():
                try:
                    with open(legacy_file) as f:
                        legacy_data = json.load(f)
                    with open(modern_file) as f:
                        modern_data = json.load(f)
                    
                    # Create dict for quick modern lookup
                    modern_dict = {}
                    for rec in modern_data:
                        key = (rec.get('base_i'), rec.get('base_j'))
                        modern_dict[key] = rec
                    
                    # Compare all legacy records
                    shown_count = 0
                    max_to_show = 1000
                    
                    for legacy_rec in legacy_data:
                        if shown_count >= max_to_show:
                            break
                        
                        pair_key = (legacy_rec.get('base_i'), legacy_rec.get('base_j'))
                        
                        # Check if modern has this pair
                        if pair_key not in modern_dict:
                            reporter.sections.append(
                                f"⚠️  WARNING: (base_i={pair_key[0]}, base_j={pair_key[1]}) - "
                                f"Legacy pair not found in modern output!\n"
                            )
                            shown_count += 1
                            continue
                        
                        modern_rec = modern_dict.get(pair_key)
                        
                        rec_comp = create_record_comparison_from_dicts(
                            record_key=pair_key,
                            record_type="distance_checks",
                            legacy_record=legacy_rec,
                            modern_record=modern_rec,
                            fields_to_compare=['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area'],
                            tolerance=tolerance,
                            legacy_source=str(legacy_file),
                            modern_source=str(modern_file)
                        )
                        reporter.add_record_comparison(rec_comp)
                        shown_count += 1
                except Exception:
                    pass
        
        # Fall back to mismatch-only if needed
        if diff_only or len([s for s in reporter.sections if 'MISMATCH' in s or 'MATCH' in s]) == 0:
            for mismatch in dcc.mismatched_checks[:reporter.max_mismatches_per_stage]:
                pair_key = (mismatch.get('base_i'), mismatch.get('base_j'))
                rec_comp = create_record_comparison_from_dicts(
                    record_key=pair_key,
                    record_type="distance_checks",
                    legacy_record=mismatch.get('legacy_record', {}),
                    modern_record=mismatch.get('modern_record', {}),
                    fields_to_compare=['dorg', 'dNN', 'plane_angle', 'd_v', 'overlap_area'],
                    tolerance=tolerance,
                    legacy_source=f"data/json_legacy/distance_checks/{pdb_id}.json",
                    modern_source=f"data/json/distance_checks/{pdb_id}.json"
                )
                reporter.add_record_comparison(rec_comp)
    
    # Add hbond list - load JSON to show ALL records
    if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
        import json
        from pathlib import Path
        
        hlc = result.hbond_list_comparison
        reporter.add_stage_header("hbond_list", 4)
        reporter.add_stage_summary(
            hlc.total_legacy, hlc.total_modern,
            hlc.total_legacy - len(hlc.missing_in_modern),
            len(hlc.missing_in_modern), len(hlc.extra_in_modern),
            len(hlc.mismatched_pairs)
        )
        
        # Load all hbond records to show every record
        if not diff_only:
            legacy_file = Path(f"data/json_legacy/hbond_list/{pdb_id}.json")
            modern_file = Path(f"data/json/hbond_list/{pdb_id}.json")
            
            if legacy_file.exists() and modern_file.exists():
                try:
                    with open(legacy_file) as f:
                        legacy_data = json.load(f)
                    with open(modern_file) as f:
                        modern_data = json.load(f)
                    
                    # Create dict for quick modern lookup
                    modern_dict = {}
                    for rec in modern_data:
                        key = (rec.get('base_i'), rec.get('base_j'))
                        modern_dict[key] = rec
                    
                    # Compare all legacy records
                    shown_count = 0
                    max_to_show = 1000
                    
                    for legacy_rec in legacy_data:
                        if shown_count >= max_to_show:
                            break
                        
                        pair_key = (legacy_rec.get('base_i'), legacy_rec.get('base_j'))
                        
                        # Check if modern has this pair
                        if pair_key not in modern_dict:
                            reporter.sections.append(
                                f"⚠️  WARNING: (base_i={pair_key[0]}, base_j={pair_key[1]}) - "
                                f"Legacy hbond pair not found in modern output!\n"
                            )
                            shown_count += 1
                            continue
                        
                        modern_rec = modern_dict.get(pair_key)
                        
                        rec_comp = create_record_comparison_from_dicts(
                            record_key=pair_key,
                            record_type="hbond_list",
                            legacy_record=legacy_rec,
                            modern_record=modern_rec,
                            fields_to_compare=['num_hbonds'],
                            tolerance=tolerance,
                            legacy_source=str(legacy_file),
                            modern_source=str(modern_file)
                        )
                        reporter.add_record_comparison(rec_comp)
                        shown_count += 1
                except Exception:
                    pass
        
        # Fall back to mismatch-only if needed
        if diff_only or len([s for s in reporter.sections if 'MISMATCH' in s or 'MATCH' in s]) == 0:
            for mismatch in hlc.mismatched_pairs[:reporter.max_mismatches_per_stage]:
                pair_key = (mismatch.get('base_i'), mismatch.get('base_j'))
                rec_comp = create_record_comparison_from_dicts(
                    record_key=pair_key,
                    record_type="hbond_list",
                    legacy_record=mismatch.get('legacy_record', {}),
                    modern_record=mismatch.get('modern_record', {}),
                    fields_to_compare=['num_hbonds'],
                    tolerance=tolerance,
                    legacy_source=f"data/json_legacy/hbond_list/{pdb_id}.json",
                    modern_source=f"data/json/hbond_list/{pdb_id}.json"
                )
                reporter.add_record_comparison(rec_comp)
    
    # Add frames - load JSON files to show ALL records
    if result.frame_comparison:
        import json
        from pathlib import Path
        
        fc = result.frame_comparison
        reporter.add_stage_header("frames", 2)
        reporter.add_stage_summary(
            fc.total_legacy, fc.total_modern,
            fc.total_legacy - len(fc.missing_residues),
            len(fc.missing_residues),
            0,
            len(fc.mismatched_calculations)
        )
        
        # Load all frame records from JSON files to show every record
        if not diff_only:
            # Try to load base_frame_calc files first
            for record_type in ['base_frame_calc', 'ls_fitting', 'frame_calc']:
                legacy_file = Path(f"data/json_legacy/{record_type}/{pdb_id}.json")
                modern_file = Path(f"data/json/{record_type}/{pdb_id}.json")
                
                if legacy_file.exists() and modern_file.exists():
                    try:
                        with open(legacy_file) as f:
                            legacy_data = json.load(f)
                        with open(modern_file) as f:
                            modern_data = json.load(f)
                        
                        # Determine fields to compare based on record type
                        if record_type == 'base_frame_calc':
                            fields = ['rms_fit', 'num_matched_atoms', 'base_type']
                        elif record_type == 'ls_fitting':
                            fields = ['rms_fit', 'num_points']
                        else:
                            fields = ['rms_fit', 'num_matched_atoms']
                        
                        # Create a dict for quick modern lookup
                        modern_dict = {}
                        for rec in modern_data:
                            key = (rec.get('chain_id'), rec.get('residue_seq'), rec.get('insertion_code', ' '))
                            modern_dict[key] = rec
                        
                        # Compare all legacy records
                        shown_count = 0
                        max_to_show = 1000  # Safety limit
                        
                        for legacy_rec in legacy_data:
                            if shown_count >= max_to_show:
                                break
                            
                            residue_key = (legacy_rec.get('chain_id'), 
                                         legacy_rec.get('residue_seq'), 
                                         legacy_rec.get('insertion_code', ' '))
                            
                            # Check if modern has this residue
                            if residue_key not in modern_dict:
                                reporter.sections.append(
                                    f"⚠️  WARNING: (chain {residue_key[0]}, seq {residue_key[1]}) - "
                                    f"Legacy residue not found in modern output!\n"
                                )
                                shown_count += 1
                                continue
                            
                            modern_rec = modern_dict.get(residue_key)
                            
                            rec_comp = create_record_comparison_from_dicts(
                                record_key=residue_key,
                                record_type=record_type,
                                legacy_record=legacy_rec,
                                modern_record=modern_rec,
                                fields_to_compare=fields,
                                tolerance=tolerance,
                                legacy_source=str(legacy_file),
                                modern_source=str(modern_file)
                            )
                            reporter.add_record_comparison(rec_comp)
                            shown_count += 1
                        
                        # Only show one record type
                        break
                    except Exception as e:
                        # Fall back to mismatch-only if loading fails
                        pass
        
        # If diff_only or file loading failed, show just mismatches
        if diff_only or not reporter.sections or len(reporter.sections) < 3:
            for mismatch in fc.mismatched_calculations[:reporter.max_mismatches_per_stage]:
                residue_key = mismatch.residue_key
                record_type = "unknown"
                mismatches_dict = mismatch.mismatches
                if 'base_frame_calc' in mismatches_dict:
                    record_type = "base_frame_calc"
                    fields = ['rms_fit', 'num_matched_atoms', 'base_type']
                elif 'ls_fitting' in mismatches_dict:
                    record_type = "ls_fitting"
                    fields = ['rms_fit', 'num_points']
                elif 'frame_calc' in mismatches_dict:
                    record_type = "frame_calc"
                    fields = ['rms_fit', 'num_matched_atoms']
                else:
                    fields = ['rms_fit']
                
                rec_comp = create_record_comparison_from_dicts(
                    record_key=residue_key,
                    record_type=record_type,
                    legacy_record=mismatch.legacy_record,
                    modern_record=mismatch.modern_record,
                    fields_to_compare=fields,
                    tolerance=tolerance,
                    legacy_source=f"data/json_legacy/{record_type}/{pdb_id}.json",
                    modern_source=f"data/json/{record_type}/{pdb_id}.json"
                )
                reporter.add_record_comparison(rec_comp)
    
    # Add steps (bpstep_params)
    if result.step_comparison:
        import json
        from pathlib import Path
        
        sc = result.step_comparison
        reporter.add_stage_header("steps", 5)
        reporter.add_stage_summary(
            sc.total_legacy, sc.total_modern,
            sc.total_legacy - len(sc.missing_steps),
            len(sc.missing_steps), 0,
            len(sc.mismatched_steps)
        )
        
        # Load all step records to show every record
        if not diff_only:
            legacy_file = Path(f"data/json_legacy/bpstep_params/{pdb_id}.json")
            modern_file = Path(f"data/json/bpstep_params/{pdb_id}.json")
            
            if legacy_file.exists() and modern_file.exists():
                try:
                    with open(legacy_file) as f:
                        legacy_data = json.load(f)
                    with open(modern_file) as f:
                        modern_data = json.load(f)
                    
                    # Create dict for quick modern lookup
                    modern_dict = {}
                    for rec in modern_data:
                        step_id = rec.get('step_id')
                        if step_id:
                            modern_dict[step_id] = rec
                    
                    # Compare all legacy records
                    shown_count = 0
                    max_to_show = 1000
                    
                    # Common step parameter fields
                    fields = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
                    
                    for legacy_rec in legacy_data:
                        if shown_count >= max_to_show:
                            break
                        
                        step_id = legacy_rec.get('step_id')
                        
                        # Check if modern has this step
                        if step_id not in modern_dict:
                            reporter.sections.append(
                                f"⚠️  WARNING: step {step_id} - "
                                f"Legacy step not found in modern output!\n"
                            )
                            shown_count += 1
                            continue
                        
                        modern_rec = modern_dict.get(step_id)
                        
                        rec_comp = create_record_comparison_from_dicts(
                            record_key=step_id,
                            record_type="bpstep_params",
                            legacy_record=legacy_rec,
                            modern_record=modern_rec,
                            fields_to_compare=fields,
                            tolerance=tolerance,
                            legacy_source=str(legacy_file),
                            modern_source=str(modern_file)
                        )
                        reporter.add_record_comparison(rec_comp)
                        shown_count += 1
                except Exception:
                    pass
    
    # Add helical params
    if result.helical_comparison:
        import json
        from pathlib import Path
        
        hc = result.helical_comparison
        reporter.add_stage_header("helical", 6)
        reporter.add_stage_summary(
            hc.total_legacy, hc.total_modern,
            hc.total_legacy - len(hc.missing_steps),
            len(hc.missing_steps), 0,
            len(hc.mismatched_steps)
        )
        
        # Load all helical records to show every record
        if not diff_only:
            legacy_file = Path(f"data/json_legacy/helical_params/{pdb_id}.json")
            modern_file = Path(f"data/json/helical_params/{pdb_id}.json")
            
            if legacy_file.exists() and modern_file.exists():
                try:
                    with open(legacy_file) as f:
                        legacy_data = json.load(f)
                    with open(modern_file) as f:
                        modern_data = json.load(f)
                    
                    # Create dict for quick modern lookup
                    modern_dict = {}
                    for rec in modern_data:
                        step_id = rec.get('step_id')
                        if step_id:
                            modern_dict[step_id] = rec
                    
                    # Compare all legacy records
                    shown_count = 0
                    max_to_show = 1000
                    
                    # Common helical parameter fields
                    fields = ['x_disp', 'y_disp', 'h_rise', 'inclination', 'tip', 'h_twist']
                    
                    for legacy_rec in legacy_data:
                        if shown_count >= max_to_show:
                            break
                        
                        step_id = legacy_rec.get('step_id')
                        
                        # Check if modern has this step
                        if step_id not in modern_dict:
                            reporter.sections.append(
                                f"⚠️  WARNING: helical step {step_id} - "
                                f"Legacy helical step not found in modern output!\n"
                            )
                            shown_count += 1
                            continue
                        
                        modern_rec = modern_dict.get(step_id)
                        
                        rec_comp = create_record_comparison_from_dicts(
                            record_key=step_id,
                            record_type="helical_params",
                            legacy_record=legacy_rec,
                            modern_record=modern_rec,
                            fields_to_compare=fields,
                            tolerance=tolerance,
                            legacy_source=str(legacy_file),
                            modern_source=str(modern_file)
                        )
                        reporter.add_record_comparison(rec_comp)
                        shown_count += 1
                except Exception:
                    pass
    
    # Add summary
    stages_compared = len(stages)
    stages_with_diffs = 0
    diff_details = []
    
    if result.atom_comparison and len(result.atom_comparison.mismatched_fields) > 0:
        stages_with_diffs += 1
        diff_details.append(f"atoms: {len(result.atom_comparison.mismatched_fields)} mismatches")
    
    if result.frame_comparison and len(result.frame_comparison.mismatched_calculations) > 0:
        stages_with_diffs += 1
        diff_details.append(f"frames: {len(result.frame_comparison.mismatched_calculations)} mismatches")
    
    if result.step_comparison and len(result.step_comparison.mismatched_steps) > 0:
        stages_with_diffs += 1
        diff_details.append(f"steps: {len(result.step_comparison.mismatched_steps)} mismatches")
    
    if result.helical_comparison and len(result.helical_comparison.mismatched_steps) > 0:
        stages_with_diffs += 1
        diff_details.append(f"helical: {len(result.helical_comparison.mismatched_steps)} mismatches")
    
    if hasattr(result, 'distance_checks_comparison') and result.distance_checks_comparison:
        if len(result.distance_checks_comparison.mismatched_checks) > 0:
            stages_with_diffs += 1
            diff_details.append(f"distance_checks: {len(result.distance_checks_comparison.mismatched_checks)} mismatches")
    
    if hasattr(result, 'hbond_list_comparison') and result.hbond_list_comparison:
        if len(result.hbond_list_comparison.mismatched_pairs) > 0:
            stages_with_diffs += 1
            diff_details.append(f"hbond_list: {len(result.hbond_list_comparison.mismatched_pairs)} mismatches")
    
    reporter.add_summary(
        stages_compared=stages_compared,
        perfect_matches=stages_compared - stages_with_diffs,
        stages_with_diffs=stages_with_diffs,
        diff_details=diff_details
    )
    
    return reporter.generate_report()

