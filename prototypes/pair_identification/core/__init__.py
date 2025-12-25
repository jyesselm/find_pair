"""Core data structures for pair identification."""

from .residue import (
    Atom,
    Residue,
    is_purine,
    is_pyrimidine,
    normalize_base_type,
    PURINE_RING_ATOMS,
    PYRIMIDINE_RING_ATOMS,
    GLYCOSIDIC_N,
    BASE_ATOMS,
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
    "Atom",
    "Residue",
    "is_purine",
    "is_pyrimidine",
    "normalize_base_type",
    "PURINE_RING_ATOMS",
    "PYRIMIDINE_RING_ATOMS",
    "GLYCOSIDIC_N",
    "BASE_ATOMS",
    "kabsch_align",
    "compute_rmsd",
    "align_and_compute_rmsd",
    "align_atom_dicts",
    "transform_points",
    "dssr_to_res_id",
    "res_id_to_dssr",
    "parse_res_id",
    "parse_dssr_id",
    "make_res_id",
    "normalize_chain",
    "extract_sequence",
    "is_standard_wc_sequence",
    "is_canonical_pair",
    "parse_pdb",
    "parse_template_pdb",
    "write_pair_pdb",
    "extract_pair_from_pdb",
]
