/**
 * @file atom_type.hpp
 * @brief Atom-level type classifications
 */

#pragma once

namespace x3dna {
namespace core {
namespace typing {

/**
 * @enum ElementType
 * @brief Classification of chemical elements
 */
enum class ElementType {
    UNKNOWN,
    // Common organic elements
    CARBON,
    NITROGEN,
    OXYGEN,
    HYDROGEN,
    PHOSPHORUS,
    SULFUR,
    // Halogens
    FLUORINE,
    CHLORINE,
    BROMINE,
    IODINE,
    // Common metal ions
    SODIUM,
    POTASSIUM,
    MAGNESIUM,
    CALCIUM,
    ZINC,
    IRON,
    MANGANESE,
    COPPER,
    COBALT,
    NICKEL,
    // Less common elements in biomolecules
    SELENIUM,
    CADMIUM,
    BARIUM,
    STRONTIUM,
    RUBIDIUM,
    CESIUM,
    LITHIUM
};

/**
 * @enum AtomLocation
 * @brief Classification of atom location within a residue
 */
enum class AtomLocation {
    UNKNOWN,
    // Nucleotide locations
    BACKBONE,      ///< Phosphate backbone: P, OP1, OP2, O5', O3'
    SUGAR,         ///< Ribose sugar: C1'-C5', O4', O2'
    NUCLEOBASE,    ///< Base ring atoms + exocyclic groups
    // Protein locations
    MAINCHAIN,     ///< Protein backbone: N, CA, C, O
    SIDECHAIN,     ///< Amino acid side chain atoms
    // Other
    SOLVENT        ///< Water/ion atoms
};

/**
 * @enum HBondRole
 * @brief Classification of hydrogen bonding capability
 */
enum class HBondRole {
    UNKNOWN,
    DONOR,         ///< Can only donate H-bond (e.g., N-H)
    ACCEPTOR,      ///< Can only accept H-bond (e.g., C=O)
    BOTH,          ///< Can donate and accept (e.g., O-H)
    NONE           ///< Cannot participate in H-bonding
};

/**
 * @brief Convert ElementType to string
 */
[[nodiscard]] inline const char* to_string(ElementType type) {
    switch (type) {
        case ElementType::UNKNOWN: return "UNKNOWN";
        case ElementType::CARBON: return "C";
        case ElementType::NITROGEN: return "N";
        case ElementType::OXYGEN: return "O";
        case ElementType::HYDROGEN: return "H";
        case ElementType::PHOSPHORUS: return "P";
        case ElementType::SULFUR: return "S";
        case ElementType::FLUORINE: return "F";
        case ElementType::CHLORINE: return "CL";
        case ElementType::BROMINE: return "BR";
        case ElementType::IODINE: return "I";
        case ElementType::SODIUM: return "NA";
        case ElementType::POTASSIUM: return "K";
        case ElementType::MAGNESIUM: return "MG";
        case ElementType::CALCIUM: return "CA";
        case ElementType::ZINC: return "ZN";
        case ElementType::IRON: return "FE";
        case ElementType::MANGANESE: return "MN";
        case ElementType::COPPER: return "CU";
        case ElementType::COBALT: return "CO";
        case ElementType::NICKEL: return "NI";
        case ElementType::SELENIUM: return "SE";
        case ElementType::CADMIUM: return "CD";
        case ElementType::BARIUM: return "BA";
        case ElementType::STRONTIUM: return "SR";
        case ElementType::RUBIDIUM: return "RB";
        case ElementType::CESIUM: return "CS";
        case ElementType::LITHIUM: return "LI";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert AtomLocation to string
 */
[[nodiscard]] inline const char* to_string(AtomLocation loc) {
    switch (loc) {
        case AtomLocation::UNKNOWN: return "UNKNOWN";
        case AtomLocation::BACKBONE: return "BACKBONE";
        case AtomLocation::SUGAR: return "SUGAR";
        case AtomLocation::NUCLEOBASE: return "NUCLEOBASE";
        case AtomLocation::MAINCHAIN: return "MAINCHAIN";
        case AtomLocation::SIDECHAIN: return "SIDECHAIN";
        case AtomLocation::SOLVENT: return "SOLVENT";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert HBondRole to string
 */
[[nodiscard]] inline const char* to_string(HBondRole role) {
    switch (role) {
        case HBondRole::UNKNOWN: return "UNKNOWN";
        case HBondRole::DONOR: return "DONOR";
        case HBondRole::ACCEPTOR: return "ACCEPTOR";
        case HBondRole::BOTH: return "BOTH";
        case HBondRole::NONE: return "NONE";
    }
    return "UNKNOWN";
}

/**
 * @brief Get legacy element index for backwards compatibility
 *
 * Element indices (legacy asym_idx compatible):
 * - 0: Unknown
 * - 1: Carbon (C)
 * - 2: Oxygen (O)
 * - 3: Hydrogen (H)
 * - 4: Nitrogen (N)
 * - 5: Sulfur (S)
 * - 6: Phosphorus (P)
 */
[[nodiscard]] inline int get_legacy_element_index(ElementType type) {
    switch (type) {
        case ElementType::CARBON: return 1;
        case ElementType::OXYGEN: return 2;
        case ElementType::HYDROGEN: return 3;
        case ElementType::NITROGEN: return 4;
        case ElementType::SULFUR: return 5;
        case ElementType::PHOSPHORUS: return 6;
        default: return 0;
    }
}

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::ElementType;
using typing::AtomLocation;
using typing::HBondRole;

} // namespace core
} // namespace x3dna
