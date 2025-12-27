"""Analysis and diagnostic tools for pair classification.

This package provides tools for analyzing and debugging pair classification:
- MissAnnotator: Annotates false negatives with detailed diagnostics
- Diagnostic dataclasses: HBondDiagnostics, GeometricDiagnostics, etc.
- Reporting: AggregateStats, PDBReport, report generation
"""

from prototypes.pair_identification.analysis.diagnostics import (
    ExpectedHBond,
    HBondDiagnostics,
    GeometricDiagnostics,
    MissAnnotation,
    PDBReport,
)
from prototypes.pair_identification.analysis.annotator import (
    MissAnnotator,
    DSSRPair,
    SlotHBond,
    load_dssr_pairs,
    load_slot_hbonds,
)
from prototypes.pair_identification.analysis.report import (
    AggregateStats,
    AggregateReporter,
    save_pdb_report,
    save_aggregate_report,
    load_pdb_report,
    load_pdb_reports,
)

__all__ = [
    # Diagnostics
    "ExpectedHBond",
    "HBondDiagnostics",
    "GeometricDiagnostics",
    "MissAnnotation",
    "PDBReport",
    # Annotator
    "MissAnnotator",
    "DSSRPair",
    "SlotHBond",
    "load_dssr_pairs",
    "load_slot_hbonds",
    # Reporting
    "AggregateStats",
    "AggregateReporter",
    "save_pdb_report",
    "save_aggregate_report",
    "load_pdb_report",
    "load_pdb_reports",
]
