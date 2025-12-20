/**
 * @file molecule_type.hpp
 * @brief Top-level molecule type classification
 */

#pragma once

namespace x3dna {
namespace core {
namespace typing {

/**
 * @enum MoleculeType
 * @brief Top-level classification of molecular entities
 */
enum class MoleculeType {
    UNKNOWN,
    NUCLEIC_ACID,  ///< RNA or DNA nucleotides
    PROTEIN,       ///< Amino acids
    LIPID,         ///< Lipid molecules
    WATER,         ///< Water molecules (HOH, WAT)
    ION,           ///< Metal ions and small charged species
    LIGAND         ///< Other small molecules, drugs, cofactors
};

/**
 * @brief Convert MoleculeType to string
 */
[[nodiscard]] inline const char* to_string(MoleculeType type) {
    switch (type) {
        case MoleculeType::UNKNOWN: return "UNKNOWN";
        case MoleculeType::NUCLEIC_ACID: return "NUCLEIC_ACID";
        case MoleculeType::PROTEIN: return "PROTEIN";
        case MoleculeType::LIPID: return "LIPID";
        case MoleculeType::WATER: return "WATER";
        case MoleculeType::ION: return "ION";
        case MoleculeType::LIGAND: return "LIGAND";
    }
    return "UNKNOWN";
}

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::MoleculeType;

} // namespace core
} // namespace x3dna
