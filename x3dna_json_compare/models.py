"""
Data models for comparison results.

Provides structured data classes for representing comparison results
between legacy and modern JSON outputs.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, TYPE_CHECKING
from pathlib import Path

if TYPE_CHECKING:
    from .pair_comparison import PairValidationComparison, DistanceChecksComparison
    from .base_pair_comparison import BasePairComparison
    from .find_bestpair_comparison import FindBestpairComparison
    from .hbond_comparison import HBondListComparison


@dataclass
class AtomInfo:
    """Information about a single atom."""
    atom_name: str
    chain_id: str
    residue_seq: int
    insertion: str
    residue_name: str
    line_number: Optional[int] = None
    pdb_line: Optional[str] = None
    legacy_atom_idx: Optional[int] = None  # Legacy atom index for direct comparison
    atom_data: Dict[str, Any] = field(default_factory=dict)


@dataclass
class FieldMismatch:
    """Details of a mismatched field between legacy and modern."""
    field_name: str
    legacy_value: Any
    modern_value: Any
    atom_key: Tuple[str, int, str, str]  # (chain_id, residue_seq, insertion, atom_name)
    legacy_line_number: Optional[int] = None
    modern_line_number: Optional[int] = None
    legacy_pdb_line: Optional[str] = None
    modern_pdb_line: Optional[str] = None


@dataclass
class AtomComparison:
    """Result of atom-level comparison."""
    missing_in_modern: List[AtomInfo] = field(default_factory=list)
    extra_in_modern: List[AtomInfo] = field(default_factory=list)
    mismatched_fields: List[FieldMismatch] = field(default_factory=list)
    count_difference: bool = False
    total_legacy: int = 0
    total_modern: int = 0
    common_count: int = 0


@dataclass
class FrameRecord:
    """Frame calculation record for a residue."""
    residue_name: str
    chain_id: str
    residue_seq: int
    insertion: str
    base_type: str
    rms_fit: float
    num_matched_atoms: int
    matched_atoms: List[str] = field(default_factory=list)
    template_file: Optional[str] = None
    residue_idx: Optional[int] = None
    legacy_residue_idx: Optional[int] = None  # Legacy residue index for direct comparison
    rotation_matrix: Optional[List[List[float]]] = None
    translation: Optional[List[float]] = None


@dataclass
class AtomLineInfo:
    """PDB line information for an atom."""
    atom_name: str
    line_number: int
    pdb_line: str
    legacy_atom_idx: Optional[int] = None  # Legacy atom index for direct comparison


@dataclass
class FrameMismatch:
    """Details of a mismatched frame calculation."""
    residue_key: Tuple[str, int, str]  # (chain_id, residue_seq, insertion)
    legacy_record: Dict[str, Any]
    modern_record: Dict[str, Any]
    mismatches: Dict[str, Any] = field(default_factory=dict)
    atom_pdb_lines: List[AtomLineInfo] = field(default_factory=list)
    legacy_matched_atoms: List[str] = field(default_factory=list)
    modern_matched_atoms: List[str] = field(default_factory=list)


@dataclass
class FrameComparison:
    """Result of frame calculation comparison."""
    missing_residues: List[FrameRecord] = field(default_factory=list)
    mismatched_calculations: List[FrameMismatch] = field(default_factory=list)
    total_legacy: int = 0
    total_modern: int = 0


@dataclass
class ComparisonResult:
    """Complete comparison result for a PDB."""
    pdb_id: str
    status: str  # 'match', 'diff', 'error'
    atom_comparison: Optional[AtomComparison] = None
    frame_comparison: Optional[FrameComparison] = None
    step_comparison: Optional['StepComparison'] = None  # Forward reference - bpstep_params comparison (for backward compatibility)
    helical_comparison: Optional['StepComparison'] = None  # helical_params comparison
    find_bestpair_comparison: Optional[Any] = None  # FindBestpairComparison (actual selected pairs)
    base_pair_comparison: Optional[Any] = None  # BasePairComparison from base_pair_comparison (pairs that pass all checks)
    pair_validation_comparison: Optional[Any] = None  # PairValidationComparison from pair_comparison
    distance_checks_comparison: Optional[Any] = None  # DistanceChecksComparison from pair_comparison
    hbond_list_comparison: Optional[Any] = None  # HBondListComparison from hbond_comparison
    residue_indices_comparison: Optional[Any] = None  # ResidueIndicesComparison from residue_indices_comparison
    errors: List[str] = field(default_factory=list)
    pdb_file_path: Optional[str] = None
    
    def has_differences(self) -> bool:
        """Check if this result contains any differences."""
        if self.status == 'error':
            return False
        
        if self.atom_comparison:
            if (self.atom_comparison.missing_in_modern or
                self.atom_comparison.extra_in_modern or
                self.atom_comparison.mismatched_fields or
                self.atom_comparison.count_difference):
                return True
        
        if self.frame_comparison:
            if (self.frame_comparison.missing_residues or
                self.frame_comparison.mismatched_calculations):
                return True
        
        if self.step_comparison:
            if (self.step_comparison.missing_steps or
                self.step_comparison.mismatched_steps):
                return True
        
        if self.helical_comparison:
            if (self.helical_comparison.missing_steps or
                self.helical_comparison.mismatched_steps):
                return True
        
        if self.find_bestpair_comparison:
            if (self.find_bestpair_comparison.missing_in_modern or
                self.find_bestpair_comparison.extra_in_modern):
                return True
        if self.base_pair_comparison:
            if (self.base_pair_comparison.missing_in_modern or
                self.base_pair_comparison.extra_in_modern or
                self.base_pair_comparison.mismatched_pairs):
                return True
        
        if self.pair_validation_comparison:
            if (self.pair_validation_comparison.missing_in_modern or
                self.pair_validation_comparison.extra_in_modern or
                self.pair_validation_comparison.mismatched_validations):
                return True
        
        if self.distance_checks_comparison:
            if (self.distance_checks_comparison.missing_in_modern or
                self.distance_checks_comparison.extra_in_modern or
                self.distance_checks_comparison.mismatched_checks):
                return True
        
        if self.hbond_list_comparison:
            if (self.hbond_list_comparison.missing_in_modern or
                self.hbond_list_comparison.extra_in_modern or
                self.hbond_list_comparison.mismatched_pairs):
                return True
        
        if self.residue_indices_comparison:
            if (self.residue_indices_comparison.missing_in_modern or
                self.residue_indices_comparison.extra_in_modern or
                self.residue_indices_comparison.mismatched_entries or
                not self.residue_indices_comparison.num_residue_match):
                return True
        
        return False
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {
            'pdb_id': self.pdb_id,
            'status': self.status,
            'timestamp': self.timestamp,
            'errors': self.errors,
        }
        
        if self.cache_key:
            result['cache_key'] = self.cache_key
        
        if self.pdb_file_path:
            result['pdb_file_path'] = self.pdb_file_path
        
        if self.atom_comparison:
            result['atom_differences'] = {
                'missing_in_modern': len(self.atom_comparison.missing_in_modern),
                'extra_in_modern': len(self.atom_comparison.extra_in_modern),
                'mismatched_fields': len(self.atom_comparison.mismatched_fields),
                'count_difference': self.atom_comparison.count_difference,
                'total_legacy': self.atom_comparison.total_legacy,
                'total_modern': self.atom_comparison.total_modern,
                'common_count': self.atom_comparison.common_count,
            }
        
        if self.frame_comparison:
            result['frame_differences'] = {
                'missing_residues': len(self.frame_comparison.missing_residues),
                'mismatched_calculations': len(self.frame_comparison.mismatched_calculations),
                'total_legacy': self.frame_comparison.total_legacy,
                'total_modern': self.frame_comparison.total_modern,
            }
        
        if self.step_comparison:
            result['step_differences'] = {
                'missing_steps': len(self.step_comparison.missing_steps),
                'mismatched_steps': len(self.step_comparison.mismatched_steps),
                'total_legacy': self.step_comparison.total_legacy,
                'total_modern': self.step_comparison.total_modern,
                'parameter_type': self.step_comparison.parameter_type,
            }
        
        if self.helical_comparison:
            result['helical_differences'] = {
                'missing_steps': len(self.helical_comparison.missing_steps),
                'mismatched_steps': len(self.helical_comparison.mismatched_steps),
                'total_legacy': self.helical_comparison.total_legacy,
                'total_modern': self.helical_comparison.total_modern,
                'parameter_type': self.helical_comparison.parameter_type,
            }
        
        if self.find_bestpair_comparison:
            result['find_bestpair_differences'] = {
                'missing_in_modern': len(self.find_bestpair_comparison.missing_in_modern),
                'extra_in_modern': len(self.find_bestpair_comparison.extra_in_modern),
                'total_legacy': self.find_bestpair_comparison.total_legacy,
                'total_modern': self.find_bestpair_comparison.total_modern,
                'common_count': self.find_bestpair_comparison.common_count,
            }
        
        if self.base_pair_comparison:
            result['base_pair_differences'] = {
                'missing_in_modern': len(self.base_pair_comparison.missing_in_modern),
                'extra_in_modern': len(self.base_pair_comparison.extra_in_modern),
                'mismatched_pairs': len(self.base_pair_comparison.mismatched_pairs),
                'total_legacy': self.base_pair_comparison.total_legacy,
                'total_modern': self.base_pair_comparison.total_modern,
                'common_count': self.base_pair_comparison.common_count,
            }
        
        if self.pair_validation_comparison:
            result['pair_validation_differences'] = {
                'missing_in_modern': len(self.pair_validation_comparison.missing_in_modern),
                'extra_in_modern': len(self.pair_validation_comparison.extra_in_modern),
                'mismatched_validations': len(self.pair_validation_comparison.mismatched_validations),
                'total_legacy': self.pair_validation_comparison.total_legacy,
                'total_modern': self.pair_validation_comparison.total_modern,
                'common_count': self.pair_validation_comparison.common_count,
            }
        
        if self.distance_checks_comparison:
            result['distance_checks_differences'] = {
                'missing_in_modern': len(self.distance_checks_comparison.missing_in_modern),
                'extra_in_modern': len(self.distance_checks_comparison.extra_in_modern),
                'mismatched_checks': len(self.distance_checks_comparison.mismatched_checks),
                'total_legacy': self.distance_checks_comparison.total_legacy,
                'total_modern': self.distance_checks_comparison.total_modern,
                'common_count': self.distance_checks_comparison.common_count,
            }
        
        if self.hbond_list_comparison:
            result['hbond_list_differences'] = {
                'missing_in_modern': len(self.hbond_list_comparison.missing_in_modern),
                'extra_in_modern': len(self.hbond_list_comparison.extra_in_modern),
                'mismatched_pairs': len(self.hbond_list_comparison.mismatched_pairs),
                'total_legacy': self.hbond_list_comparison.total_legacy,
                'total_modern': self.hbond_list_comparison.total_modern,
                'common_count': self.hbond_list_comparison.common_count,
            }
        
        if self.residue_indices_comparison:
            result['residue_indices_differences'] = {
                'missing_in_modern': 1 if self.residue_indices_comparison.missing_in_modern else 0,
                'extra_in_modern': 1 if self.residue_indices_comparison.extra_in_modern else 0,
                'num_residue_match': self.residue_indices_comparison.num_residue_match,
                'mismatched_entries': len(self.residue_indices_comparison.mismatched_entries),
                'legacy_num_residue': self.residue_indices_comparison.legacy_num_residue,
                'modern_num_residue': self.residue_indices_comparison.modern_num_residue,
            }
        
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ComparisonResult':
        """Create from dictionary (for cache deserialization)."""
        # Reconstruct minimal comparison objects from summary data
        # so has_differences() works correctly
        atom_comparison = None
        if 'atom_differences' in data:
            atom_diff = data['atom_differences']
            atom_comparison = AtomComparison(
                missing_in_modern=[],  # Empty list but count preserved
                extra_in_modern=[],    # Empty list but count preserved
                mismatched_fields=[],  # Empty list but count preserved
                count_difference=atom_diff.get('count_difference', False),
                total_legacy=atom_diff.get('total_legacy', 0),
                total_modern=atom_diff.get('total_modern', 0),
                common_count=atom_diff.get('common_count', 0),
            )
            # Set lists to non-empty if there were differences
            # This ensures has_differences() returns True when needed
            # We only need one placeholder per type since has_differences() just checks if list is non-empty
            if atom_diff.get('missing_in_modern', 0) > 0:
                atom_comparison.missing_in_modern = [
                    AtomInfo('?', '', 0, ' ', '?')  # Minimal placeholder
                ]
            if atom_diff.get('extra_in_modern', 0) > 0:
                atom_comparison.extra_in_modern = [
                    AtomInfo('?', '', 0, ' ', '?')  # Minimal placeholder
                ]
            if atom_diff.get('mismatched_fields', 0) > 0:
                atom_comparison.mismatched_fields = [
                    FieldMismatch('?', None, None, ('', 0, ' ', '?'))  # Minimal placeholder
                ]
        
        frame_comparison = None
        if 'frame_differences' in data:
            frame_diff = data['frame_differences']
            frame_comparison = FrameComparison(
                missing_residues=[],  # Empty list but count preserved
                mismatched_calculations=[],  # Empty list but count preserved
                total_legacy=frame_diff.get('total_legacy', 0),
                total_modern=frame_diff.get('total_modern', 0),
            )
            # Set lists to non-empty if there were differences
            # This ensures has_differences() returns True when needed
            # We only need one placeholder per type since has_differences() just checks if list is non-empty
            if frame_diff.get('missing_residues', 0) > 0:
                frame_comparison.missing_residues = [
                    FrameRecord('?', '', 0, ' ', '?', 0.0, 0)  # Minimal placeholder
                ]
            if frame_diff.get('mismatched_calculations', 0) > 0:
                frame_comparison.mismatched_calculations = [
                    FrameMismatch(('', 0, ' '), {}, {}, {})  # Minimal placeholder
                ]
        
        step_comparison = None
        if 'step_differences' in data:
            step_diff = data['step_differences']
            # Import StepComparison here to avoid circular dependency
            from .step_comparison import StepComparison
            step_comparison = StepComparison(
                missing_steps=[],  # Empty list but count preserved
                mismatched_steps=[],  # Empty list but count preserved
                total_legacy=step_diff.get('total_legacy', 0),
                total_modern=step_diff.get('total_modern', 0),
                parameter_type=step_diff.get('parameter_type', 'bpstep_params'),
            )
            # Set lists to non-empty if there were differences
            if step_diff.get('missing_steps', 0) > 0:
                step_comparison.missing_steps = [{}]  # Minimal placeholder
            if step_diff.get('mismatched_steps', 0) > 0:
                step_comparison.mismatched_steps = [{}]  # Minimal placeholder
        
        helical_comparison = None
        if 'helical_differences' in data:
            helical_diff = data['helical_differences']
            # Import StepComparison here to avoid circular dependency
            from .step_comparison import StepComparison
            helical_comparison = StepComparison(
                missing_steps=[],  # Empty list but count preserved
                mismatched_steps=[],  # Empty list but count preserved
                total_legacy=helical_diff.get('total_legacy', 0),
                total_modern=helical_diff.get('total_modern', 0),
                parameter_type=helical_diff.get('parameter_type', 'helical_params'),
            )
            # Set lists to non-empty if there were differences
            if helical_diff.get('missing_steps', 0) > 0:
                helical_comparison.missing_steps = [{}]  # Minimal placeholder
            if helical_diff.get('mismatched_steps', 0) > 0:
                helical_comparison.mismatched_steps = [{}]  # Minimal placeholder
        
        # Determine correct status from differences (recalculate from summary data)
        # This fixes cases where cached status might be incorrect
        status = data.get('status', 'error')
        if data.get('errors'):
            status = 'error'
        else:
            # Check if there are any differences based on summary data
            has_diffs = False
            if atom_comparison:
                has_diffs = (atom_comparison.missing_in_modern or
                           atom_comparison.extra_in_modern or
                           atom_comparison.mismatched_fields or
                           atom_comparison.count_difference)
            if not has_diffs and frame_comparison:
                has_diffs = (frame_comparison.missing_residues or
                           frame_comparison.mismatched_calculations)
            
            if has_diffs:
                status = 'diff'
            else:
                status = 'match'
        
        return cls(
            pdb_id=data['pdb_id'],
            status=status,  # Use recalculated status
            timestamp=data.get('timestamp', 0.0),
            errors=data.get('errors', []),
            cache_key=data.get('cache_key'),
            pdb_file_path=data.get('pdb_file_path'),
            atom_comparison=atom_comparison,
            frame_comparison=frame_comparison,
            step_comparison=step_comparison,
            helical_comparison=helical_comparison,
        )

