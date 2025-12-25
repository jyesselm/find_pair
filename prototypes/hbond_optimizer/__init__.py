"""
H-Bond Optimizer Prototype Package.

Provides slot-based H-bond detection with bifurcation support.
Achieves ~96% recall, ~94% precision vs DSSR on 100 PDB test set.

Usage:
    from prototypes.hbond_optimizer import HBondOptimizer, Residue
    from prototypes.hbond_optimizer import parse_pdb_residues, run_dssr, parse_dssr_output
"""

from .optimizer import (
    HBondOptimizer,
    Residue,
    HBond,
    HBondCandidate,
    format_hbond,
)

from .geometry import (
    HSlot,
    LPSlot,
    DONOR_CAPACITY,
    ACCEPTOR_CAPACITY,
    predict_h_slots,
    predict_lp_slots,
    score_hbond_alignment,
    compute_base_normal,
    normalize,
    angle_between,
    get_donor_capacity,
    get_acceptor_capacity,
)

from .compare_with_dssr import (
    parse_pdb_residues,
    run_dssr,
    parse_dssr_output,
    DSSRHBond,
)

from .modified_registry import (
    ModifiedResidueRegistry,
    get_registry,
    get_parent_base,
)

__all__ = [
    # Core optimizer
    'HBondOptimizer',
    'Residue',
    'HBond',
    'HBondCandidate',
    'format_hbond',
    # Geometry
    'HSlot',
    'LPSlot',
    'DONOR_CAPACITY',
    'ACCEPTOR_CAPACITY',
    'predict_h_slots',
    'predict_lp_slots',
    'score_hbond_alignment',
    'compute_base_normal',
    'normalize',
    'angle_between',
    'get_donor_capacity',
    'get_acceptor_capacity',
    # DSSR comparison
    'parse_pdb_residues',
    'run_dssr',
    'parse_dssr_output',
    'DSSRHBond',
]
