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
}  // namespace pdb_columns

/**
 * @brief Geometric thresholds for structure analysis
 */
namespace geometry {
    // Distance thresholds (Angstroms)
    constexpr double COVALENT_BOND_MAX = 2.0;          // Maximum covalent bond distance
    constexpr double HYDROGEN_BOND_MIN = 1.8;          // Minimum H-bond distance
    constexpr double HYDROGEN_BOND_MAX = 3.5;          // Maximum H-bond distance
    constexpr double SAME_ATOM_THRESHOLD = 0.1;        // Below this, atoms are the same

    // Angle thresholds (degrees)
    constexpr double MAX_PLANE_ANGLE = 65.0;           // Maximum base plane angle

    // Numeric constants
    constexpr double LARGE_NUMBER = 1.0e+18;           // Very large number for initialization
    constexpr double SMALL_NUMBER = 1.0e-10;           // Very small number for comparisons
}  // namespace geometry

/**
 * @brief Nucleic acid base constants
 */
namespace nucleotides {
    // Ring atom counts
    constexpr size_t PURINE_RING_ATOM_COUNT = 9;       // N1, C2, N3, C4, C5, C6, N7, C8, N9
    constexpr size_t PYRIMIDINE_RING_ATOM_COUNT = 6;   // N1, C2, N3, C4, C5, C6
    constexpr size_t MIN_ATOMS_FOR_FIT = 3;            // Minimum atoms for least-squares fitting

    // Ring atom names (trimmed, without padding)
    // Purine ring atoms: fused 6+5 ring system (A, G, I)
    inline const std::vector<std::string>& purine_ring_atoms() {
        static const std::vector<std::string> atoms = {
            "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"
        };
        return atoms;
    }

    // Pyrimidine ring atoms: single 6-membered ring (C, U, T, P)
    inline const std::vector<std::string>& pyrimidine_ring_atoms() {
        static const std::vector<std::string> atoms = {
            "N1", "C2", "N3", "C4", "C5", "C6"
        };
        return atoms;
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
        return type == core::ResidueType::ADENINE ||
               type == core::ResidueType::GUANINE ||
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
}  // namespace nucleotides

/**
 * @brief Output formatting constants
 */
namespace formatting {
    constexpr int COORDINATE_PRECISION = 3;            // Decimal places for coordinates
    constexpr int ANGLE_PRECISION = 2;                 // Decimal places for angles
    constexpr int PARAMETER_PRECISION = 2;             // Decimal places for step parameters
}  // namespace formatting

}  // namespace constants
}  // namespace x3dna
