"""Core data structures for pair identification."""

from .residue import (
    Atom,
    Residue,
    is_purine,
    is_pyrimidine,
    normalize_base_type,
)

from .constants import (
    PURINE_RING_ATOMS,
    PYRIMIDINE_RING_ATOMS,
    GLYCOSIDIC_N,
    BASE_ATOMS,
    PURINES,
    PYRIMIDINES,
    DNA_BASES,
    RNA_BASES,
    MODIFIED_BASES,
    CANONICAL_WC_SEQUENCES,
    WOBBLE_SEQUENCES,
    CANONICAL_SEQUENCES,
    MAX_DORG,
    MAX_D_V,
    MAX_PLANE_ANGLE,
    MIN_DNN,
    OVERLAP_THRESHOLD,
    D_V_WEIGHT,
    PLANE_ANGLE_DIVISOR,
    MODIFIED_BASE_MAP,
    DNA_PREFIX_MAP,
)

from .alignment import (
    kabsch_align,
    compute_rmsd,
    align_and_compute_rmsd,
    align_atom_dicts,
    transform_points,
)

from .identifiers import (
    dssr_to_res_id,
    res_id_to_dssr,
    parse_res_id,
    parse_dssr_id,
    make_res_id,
    normalize_chain,
    extract_sequence,
    is_standard_wc_sequence,
    is_canonical_pair,
)

from .pdb_parser import (
    parse_pdb,
    parse_template_pdb,
    write_pair_pdb,
    extract_pair_from_pdb,
)

__all__ = [
    # Data structures
    "Atom",
    "Residue",
    "is_purine",
    "is_pyrimidine",
    "normalize_base_type",
    # Constants
    "PURINE_RING_ATOMS",
    "PYRIMIDINE_RING_ATOMS",
    "GLYCOSIDIC_N",
    "BASE_ATOMS",
    "PURINES",
    "PYRIMIDINES",
    "DNA_BASES",
    "RNA_BASES",
    "MODIFIED_BASES",
    "CANONICAL_WC_SEQUENCES",
    "WOBBLE_SEQUENCES",
    "CANONICAL_SEQUENCES",
    "MAX_DORG",
    "MAX_D_V",
    "MAX_PLANE_ANGLE",
    "MIN_DNN",
    "OVERLAP_THRESHOLD",
    "D_V_WEIGHT",
    "PLANE_ANGLE_DIVISOR",
    "MODIFIED_BASE_MAP",
    "DNA_PREFIX_MAP",
    # Alignment
    "kabsch_align",
    "compute_rmsd",
    "align_and_compute_rmsd",
    "align_atom_dicts",
    "transform_points",
    # Identifiers
    "dssr_to_res_id",
    "res_id_to_dssr",
    "parse_res_id",
    "parse_dssr_id",
    "make_res_id",
    "normalize_chain",
    "extract_sequence",
    "is_standard_wc_sequence",
    "is_canonical_pair",
    # PDB parsing
    "parse_pdb",
    "parse_template_pdb",
    "write_pair_pdb",
    "extract_pair_from_pdb",
]
