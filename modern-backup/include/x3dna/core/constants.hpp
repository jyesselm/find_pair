/**
 * @file constants.hpp
 * @brief Named constants for the x3dna library
 *
 * Centralizes magic numbers and configuration values to improve
 * code readability and maintainability.
 */

#pragma once

#include <cstddef>
#include <cctype>
#include <string>
#include <string_view>
#include <vector>
#include <x3dna/core/residue_type.hpp>
#include <x3dna/core/molecule_type.hpp>

namespace x3dna {
namespace constants {

/**
 * @brief PDB file format column definitions
 *
 * Based on the official PDB file format specification.
 * Column indices are 0-based for use with std::string::substr().
 * See: https://www.wwpdb.org/documentation/file-format
 */
namespace pdb_columns {
// ATOM/HETATM record columns (0-based indices for substr)
constexpr size_t RECORD_TYPE_START = 0;
constexpr size_t RECORD_TYPE_LEN = 6;

constexpr size_t ATOM_SERIAL_START = 6;
constexpr size_t ATOM_SERIAL_LEN = 5;

constexpr size_t ATOM_NAME_START = 12;
constexpr size_t ATOM_NAME_LEN = 4;

constexpr size_t ALT_LOC = 16;

constexpr size_t RESIDUE_NAME_START = 17;
constexpr size_t RESIDUE_NAME_LEN = 3;

constexpr size_t CHAIN_ID = 21;

constexpr size_t RESIDUE_SEQ_START = 22;
constexpr size_t RESIDUE_SEQ_LEN = 4;

constexpr size_t INSERTION_CODE = 26;

constexpr size_t X_COORD_START = 30;
constexpr size_t X_COORD_LEN = 8;

constexpr size_t Y_COORD_START = 38;
constexpr size_t Y_COORD_LEN = 8;

constexpr size_t Z_COORD_START = 46;
constexpr size_t Z_COORD_LEN = 8;

constexpr size_t OCCUPANCY_START = 54;
constexpr size_t OCCUPANCY_LEN = 6;

constexpr size_t B_FACTOR_START = 60;
constexpr size_t B_FACTOR_LEN = 6;

constexpr size_t ELEMENT_START = 76;
constexpr size_t ELEMENT_LEN = 2;

// Minimum line lengths for different record types
constexpr size_t MIN_ATOM_LINE = 52;      // Minimum for coordinates
constexpr size_t MIN_FULL_ATOM_LINE = 78; // Full ATOM record

// MODEL record
constexpr size_t MODEL_NUM_START = 6;
constexpr size_t MODEL_NUM_LEN = 4;

// HEADER record
constexpr size_t HEADER_PDB_ID_START = 62;
constexpr size_t HEADER_PDB_ID_LEN = 4;
} // namespace pdb_columns

// Note: Geometric thresholds now centralized in config/parameters_generated.hpp
// Use config::params::COVALENT_BOND_MAX, etc.

/**
 * @brief Nucleic acid base constants
 */
namespace nucleotides {
// Ring atom counts
constexpr size_t PURINE_RING_ATOM_COUNT = 9;     // N1, C2, N3, C4, C5, C6, N7, C8, N9
constexpr size_t PYRIMIDINE_RING_ATOM_COUNT = 6; // N1, C2, N3, C4, C5, C6
constexpr size_t MIN_ATOMS_FOR_FIT = 3;          // Minimum atoms for least-squares fitting

// Ring atom names (trimmed, without padding)
// Purine ring atoms: fused 6+5 ring system (A, G, I)
inline const std::vector<std::string>& purine_ring_atoms() {
    static const std::vector<std::string> ATOMS = {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"};
    return ATOMS;
}

// Pyrimidine ring atoms: single 6-membered ring (C, U, T, P)
inline const std::vector<std::string>& pyrimidine_ring_atoms() {
    static const std::vector<std::string> ATOMS = {"N1", "C2", "N3", "C4", "C5", "C6"};
    return ATOMS;
}

/**
 * @brief Check if an atom name is a ring atom
 * @param atom_name Atom name (leading/trailing spaces are trimmed)
 * @return true if this is a base ring atom
 */
inline bool is_ring_atom(std::string_view atom_name) {
    // Trim whitespace
    while (!atom_name.empty() && std::isspace(static_cast<unsigned char>(atom_name.front()))) {
        atom_name.remove_prefix(1);
    }
    while (!atom_name.empty() && std::isspace(static_cast<unsigned char>(atom_name.back()))) {
        atom_name.remove_suffix(1);
    }

    // Check against purine atoms (superset includes all pyrimidine atoms)
    for (const auto& ring_atom : purine_ring_atoms()) {
        if (atom_name == ring_atom) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Check if a residue type is a purine
 * @param type ResidueType enum value
 * @return true for ADENINE, GUANINE, INOSINE; false otherwise
 */
inline bool is_purine(core::ResidueType type) {
    return type == core::ResidueType::ADENINE || type == core::ResidueType::GUANINE ||
           type == core::ResidueType::INOSINE;
}

/**
 * @brief Get ring atom names for a residue type
 * @param type ResidueType enum value
 * @return Ring atom names for that type (purine or pyrimidine set)
 */
inline const std::vector<std::string>& ring_atoms_for_type(core::ResidueType type) {
    if (is_purine(type)) {
        return purine_ring_atoms();
    }
    return pyrimidine_ring_atoms();
}

/**
 * @brief Check if a BaseType is a purine
 * @param type BaseType enum value
 * @return true for ADENINE, GUANINE, INOSINE; false otherwise
 */
inline bool is_purine(core::BaseType type) {
    return type == core::BaseType::ADENINE || type == core::BaseType::GUANINE || type == core::BaseType::INOSINE;
}

/**
 * @brief Check if a BaseType is a pyrimidine
 * @param type BaseType enum value
 * @return true for CYTOSINE, THYMINE, URACIL, PSEUDOURIDINE; false otherwise
 */
inline bool is_pyrimidine(core::BaseType type) {
    return type == core::BaseType::CYTOSINE || type == core::BaseType::THYMINE || type == core::BaseType::URACIL ||
           type == core::BaseType::PSEUDOURIDINE;
}

/**
 * @brief Check if a ResidueType is a standard nucleotide base
 * @param type ResidueType enum value
 * @return true for ADENINE, GUANINE, CYTOSINE, THYMINE, URACIL; false otherwise
 */
inline bool is_standard_base(core::ResidueType type) {
    return type == core::ResidueType::ADENINE || type == core::ResidueType::GUANINE ||
           type == core::ResidueType::CYTOSINE || type == core::ResidueType::THYMINE ||
           type == core::ResidueType::URACIL;
}

/**
 * @brief Check if a ResidueType is a special nucleotide (not one of the standard 5)
 * @param type ResidueType enum value
 * @return true for INOSINE, PSEUDOURIDINE, NONCANONICAL_RNA; false otherwise
 */
inline bool is_special_base(core::ResidueType type) {
    return type == core::ResidueType::INOSINE || type == core::ResidueType::PSEUDOURIDINE ||
           type == core::ResidueType::NONCANONICAL_RNA;
}

/**
 * @brief Get ring atom names for a BaseType
 * @param type BaseType enum value
 * @return Ring atom names for that type (purine or pyrimidine set)
 */
inline const std::vector<std::string>& ring_atoms_for_type(core::BaseType type) {
    if (is_purine(type)) {
        return purine_ring_atoms();
    }
    return pyrimidine_ring_atoms();
}
} // namespace nucleotides

/**
 * @brief Output formatting constants
 */
namespace formatting {
constexpr int COORDINATE_PRECISION = 3; // Decimal places for coordinates
constexpr int ANGLE_PRECISION = 2;      // Decimal places for angles
constexpr int PARAMETER_PRECISION = 2;  // Decimal places for step parameters
} // namespace formatting

/**
 * @brief Ring geometry data for nucleotide bases
 *
 * Standard coordinates and atom names for RMSD-based nucleotide detection.
 * Order matches legacy RA_LIST: C4, N3, C2, N1, C6, C5, N7, C8, N9
 */
namespace ring_data {
// Standard ring atom names (matches legacy RA_LIST order)
constexpr std::array<const char*, 9> RING_ATOM_NAMES = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ",
                                                        " C5 ", " N7 ", " C8 ", " N9 "};

// Standard ring geometry for RMSD fitting (matches RA_LIST order)
// From legacy xyz_ring array - standard purine/pyrimidine ring system
constexpr std::array<std::array<double, 3>, 9> STANDARD_RING_GEOMETRY = {{
    {{-1.265, 3.177, 0.000}}, // C4
    {{-2.342, 2.364, 0.001}}, // N3
    {{-1.999, 1.087, 0.000}}, // C2
    {{-0.700, 0.641, 0.000}}, // N1
    {{0.424, 1.460, 0.000}},  // C6
    {{0.071, 2.833, 0.000}},  // C5
    {{0.870, 3.969, 0.000}},  // N7 (purine only)
    {{0.023, 4.962, 0.000}},  // C8 (purine only)
    {{-1.289, 4.551, 0.000}}, // N9 (purine only)
}};

// Common ring atoms (pyrimidine + purine share these 6 atoms)
constexpr std::array<const char*, 6> COMMON_RING_ATOMS = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "};

// Additional purine-only ring atoms (3-membered imidazole ring)
constexpr std::array<const char*, 3> PURINE_RING_ATOMS = {" N7 ", " C8 ", " N9 "};
} // namespace ring_data

/**
 * @brief Nucleotide residue name lists
 *
 * Standard nucleotide names and molecules to exclude from analysis.
 */
namespace nucleotide_lists {
// Standard nucleotide residue names (matches legacy NT_LIST)
// DNA: DA, DC, DG, DT, DU, A, C, G, T, U
// RNA: A, C, G, U, PSU (pseudouridine), I (inosine), P5P
// Nucleotide derivatives: ADP, GDP, CDP, UDP, TDP, PU, DI
constexpr std::array<const char*, 20> NT_LIST = {"A",   "C",   "G",   "T",   "U",   "PSU", "P5P", "PU", "I",  "DI",
                                                 "ADP", "GDP", "CDP", "UDP", "TDP", "DA",  "DC",  "DG", "DT", "DU"};

// Non-nucleotide molecules to exclude (buffers, salts, etc.)
// These are common crystallization additives that should not be analyzed
constexpr std::array<const char*, 11> EXCLUDED_MOLECULES = {"MES", "HEPES", "TRIS", "EDO", "GOL", "SO4",
                                                            "PO4", "ACT",   "FMT",  "EFZ", "LYA"};
} // namespace nucleotide_lists

/**
 * @brief Hydrogen bond classification data
 *
 * Constants for H-bond quality scoring and Watson-Crick pair detection.
 */
namespace hbond_data {
// Watson-Crick pair types (for quality scoring)
// Includes standard WC pairs and inosine-cytosine pairs
constexpr std::array<const char*, 9> WC_PAIR_LIST = {"XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"};

// Backbone atoms excluded from H-bond quality checks
// These atoms are part of sugar-phosphate backbone, not base pairing
// Includes phosphate oxygens (O1P, O2P), sugar oxygens (O3', O4', O5'), and N7
constexpr std::array<const char*, 6> HBOND_EXCLUDED_ATOMS = {" O1P", " O2P", " O3'", " O4'", " O5'", " N7 "};
} // namespace hbond_data

} // namespace constants
} // namespace x3dna
