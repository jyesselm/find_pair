"""
Data models for comparison results.

Provides structured data classes for representing comparison results
between legacy and modern JSON outputs.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path


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
    rotation_matrix: Optional[List[List[float]]] = None
    translation: Optional[List[float]] = None


@dataclass
class AtomLineInfo:
    """PDB line information for an atom."""
    atom_name: str
    line_number: int
    pdb_line: str


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
    errors: List[str] = field(default_factory=list)
    cache_key: Optional[str] = None
    timestamp: float = 0.0
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
        
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ComparisonResult':
        """Create from dictionary (for cache deserialization)."""
        # This is a simplified version - full deserialization would restore all fields
        return cls(
            pdb_id=data['pdb_id'],
            status=data['status'],
            timestamp=data.get('timestamp', 0.0),
            errors=data.get('errors', []),
            cache_key=data.get('cache_key'),
            pdb_file_path=data.get('pdb_file_path'),
        )

