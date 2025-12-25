"""cWW miss annotator tool for analyzing missed Watson-Crick pairs."""

from .diagnostics import (
    ExpectedHBond,
    HBondDiagnostics,
    GeometricDiagnostics,
    MissAnnotation,
    PDBReport,
)
from .hbond_analyzer import (
    HBondAnalyzer,
    parse_dssr_hbonds,
    CWW_EXPECTED_PATTERNS,
)
from .loaders import (
    DSSRPair,
    SlotHBond,
    dssr_to_slot_id,
    load_dssr_pairs,
    load_slot_hbonds,
    slot_to_dssr_id,
)
from .report import (
    AggregateStats,
    AggregateReporter,
    save_pdb_report,
    save_aggregate_report,
    load_pdb_report,
    load_pdb_reports,
)
from .geometric_analyzer import GeometricAnalyzer
from .annotator import MissAnnotator

__all__ = [
    # Diagnostics dataclasses
    "ExpectedHBond",
    "HBondDiagnostics",
    "GeometricDiagnostics",
    "MissAnnotation",
    "PDBReport",
    # H-bond analyzer
    "HBondAnalyzer",
    "parse_dssr_hbonds",
    "CWW_EXPECTED_PATTERNS",
    # Geometric analyzer
    "GeometricAnalyzer",
    # Main annotator
    "MissAnnotator",
    # Loaders
    "DSSRPair",
    "SlotHBond",
    "dssr_to_slot_id",
    "slot_to_dssr_id",
    "load_dssr_pairs",
    "load_slot_hbonds",
    # Report generation
    "AggregateStats",
    "AggregateReporter",
    "save_pdb_report",
    "save_aggregate_report",
    "load_pdb_report",
    "load_pdb_reports",
]
