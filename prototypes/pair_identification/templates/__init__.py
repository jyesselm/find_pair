"""Template-based pair classification.

This package provides template loading, alignment, and classification
for Leontis-Westhof base pair identification.

Modules:
    loader: Template file discovery and loading
    aligner: Template alignment using Kabsch algorithm
    classifier: Combined RMSD + H-bond LW classification
    stats: Statistical geometry for validation

Usage:
    from templates import TemplateLoader, TemplateAligner, LWClassifier
    from templates.stats import TemplateStatsLoader

    # Load templates and classify
    classifier = LWClassifier()
    result = classifier.classify(res1, res2, hbonds)
    print(f"Best LW: {result.best_lw} (RMSD={result.best_rmsd:.3f})")
"""

from prototypes.pair_identification.templates.loader import TemplateLoader
from prototypes.pair_identification.templates.aligner import (
    TemplateAligner,
    AlignmentResult,
    ClassificationResult,
)
from prototypes.pair_identification.templates.classifier import (
    LWClassifier,
    LWClassificationResult,
)
from prototypes.pair_identification.templates.stats import (
    TemplateStatsLoader,
    GeometryStats,
)

__all__ = [
    "TemplateLoader",
    "TemplateAligner",
    "AlignmentResult",
    "ClassificationResult",
    "LWClassifier",
    "LWClassificationResult",
    "TemplateStatsLoader",
    "GeometryStats",
]
