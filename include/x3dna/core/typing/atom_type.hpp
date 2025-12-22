/**
 * @file atom_type.hpp
 * @brief Atom-level type classifications
 */

#pragma once

#include <cstddef>
#include <cstdint>

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
    BACKBONE,   ///< Phosphate backbone: P, OP1, OP2, O5', O3'
    SUGAR,      ///< Ribose sugar: C1'-C5', O4', O2'
    NUCLEOBASE, ///< Base ring atoms + exocyclic groups
    // Protein locations
    MAINCHAIN, ///< Protein backbone: N, CA, C, O
    SIDECHAIN, ///< Amino acid side chain atoms
    // Other
    SOLVENT ///< Water/ion atoms
};

/**
 * @enum AtomType
 * @brief Standard atom types for fast integer comparison
 *
 * Covers all common atoms in nucleotides and amino acids.
 * Using enum comparison is ~100x faster than string comparison.
 * UNKNOWN is used for non-standard/modified atoms.
 */
enum class AtomType : uint8_t {
    UNKNOWN = 0,

    // === Nucleotide ring atoms (indices 1-9 for array lookup) ===
    C4 = 1, ///< Pyrimidine/purine ring
    N3 = 2, ///< Pyrimidine/purine ring
    C2 = 3, ///< Pyrimidine/purine ring
    N1 = 4, ///< Pyrimidine/purine ring
    C6 = 5, ///< Pyrimidine/purine ring
    C5 = 6, ///< Pyrimidine/purine ring
    N7 = 7, ///< Purine only
    C8 = 8, ///< Purine only
    N9 = 9, ///< Purine only

    // === Nucleotide exocyclic atoms ===
    O6 = 10,  ///< Guanine carbonyl
    N6 = 11,  ///< Adenine amino
    O2 = 12,  ///< Uracil/cytosine carbonyl
    N2 = 13,  ///< Guanine amino
    O4 = 14,  ///< Uracil/thymine carbonyl
    N4 = 15,  ///< Cytosine amino
    C5M = 16, ///< Thymine methyl (C7 in some nomenclatures)
    C7 = 17,  ///< Alternative name for thymine methyl

    // === Nucleotide sugar atoms ===
    C1_PRIME = 20, ///< C1'
    C2_PRIME = 21, ///< C2'
    C3_PRIME = 22, ///< C3'
    C4_PRIME = 23, ///< C4'
    C5_PRIME = 24, ///< C5'
    O2_PRIME = 25, ///< O2' (RNA only)
    O3_PRIME = 26, ///< O3'
    O4_PRIME = 27, ///< O4'
    O5_PRIME = 28, ///< O5'

    // === Nucleotide backbone atoms ===
    P = 30,   ///< Phosphorus
    OP1 = 31, ///< Phosphate oxygen 1
    OP2 = 32, ///< Phosphate oxygen 2
    OP3 = 33, ///< Phosphate oxygen 3 (5' terminal)

    // === Amino acid backbone atoms ===
    N = 40,   ///< Backbone nitrogen
    CA = 41,  ///< Alpha carbon
    C = 42,   ///< Backbone carbonyl carbon
    O = 43,   ///< Backbone carbonyl oxygen
    OXT = 44, ///< C-terminal oxygen

    // === Common amino acid side chain atoms ===
    CB = 50, ///< Beta carbon
    CG = 51, ///< Gamma carbon
    CG1 = 52,
    CG2 = 53,
    CD = 54, ///< Delta carbon
    CD1 = 55,
    CD2 = 56,
    CE = 57, ///< Epsilon carbon
    CE1 = 58,
    CE2 = 59,
    CE3 = 60,
    CZ = 61, ///< Zeta carbon
    CZ2 = 62,
    CZ3 = 63,
    CH2 = 64,
    OG = 65, ///< Serine/threonine hydroxyl
    OG1 = 66,
    OD1 = 67, ///< Aspartate/asparagine
    OD2 = 68,
    OE1 = 69, ///< Glutamate/glutamine
    OE2 = 70,
    OH = 71,  ///< Tyrosine hydroxyl
    ND1 = 72, ///< Histidine
    ND2 = 73, ///< Asparagine
    NE = 74,  ///< Arginine
    NE1 = 75, ///< Tryptophan
    NE2 = 76, ///< Histidine/glutamine
    NH1 = 77, ///< Arginine
    NH2 = 78, ///< Arginine
    NZ = 79,  ///< Lysine
    SD = 80,  ///< Methionine sulfur
    SG = 81,  ///< Cysteine sulfur

    // === Water ===
    OW = 90, ///< Water oxygen

    // Marker for count
    COUNT = 100
};

/// Number of ring atoms (for array sizing)
constexpr size_t NUM_RING_ATOM_TYPES = 9;

/// Array of ring atom types in order (for indexed iteration)
constexpr AtomType RING_ATOM_TYPES[NUM_RING_ATOM_TYPES] = {AtomType::C4, AtomType::N3, AtomType::C2,
                                                           AtomType::N1, AtomType::C6, AtomType::C5,
                                                           AtomType::N7, AtomType::C8, AtomType::N9};

/// Check if atom type is a ring atom
[[nodiscard]] inline bool is_ring_atom(AtomType type) {
    auto val = static_cast<uint8_t>(type);
    return val >= 1 && val <= 9;
}

/// Check if atom type is a purine-only ring atom (N7, C8, N9)
[[nodiscard]] inline bool is_purine_ring_atom(AtomType type) {
    return type == AtomType::N7 || type == AtomType::C8 || type == AtomType::N9;
}

/// Check if atom type is a sugar atom (C1'-C5', O2'-O5', O4')
[[nodiscard]] inline bool is_sugar_atom(AtomType type) {
    auto val = static_cast<uint8_t>(type);
    return val >= 20 && val <= 28;
}

/// Check if atom type is a backbone atom (P, OP1, OP2, OP3)
[[nodiscard]] inline bool is_backbone_atom(AtomType type) {
    auto val = static_cast<uint8_t>(type);
    return val >= 30 && val <= 33;
}

/// Check if atom type is an exocyclic nucleobase atom (O6, N6, O2, N2, O4, N4, C5M, C7)
[[nodiscard]] inline bool is_exocyclic_atom(AtomType type) {
    auto val = static_cast<uint8_t>(type);
    return val >= 10 && val <= 17;
}

/// Convert AtomType to string name
[[nodiscard]] inline const char* to_string(AtomType type) {
    switch (type) {
        case AtomType::UNKNOWN:
            return "UNKNOWN";
        // Ring atoms
        case AtomType::C4:
            return "C4";
        case AtomType::N3:
            return "N3";
        case AtomType::C2:
            return "C2";
        case AtomType::N1:
            return "N1";
        case AtomType::C6:
            return "C6";
        case AtomType::C5:
            return "C5";
        case AtomType::N7:
            return "N7";
        case AtomType::C8:
            return "C8";
        case AtomType::N9:
            return "N9";
        // Exocyclic
        case AtomType::O6:
            return "O6";
        case AtomType::N6:
            return "N6";
        case AtomType::O2:
            return "O2";
        case AtomType::N2:
            return "N2";
        case AtomType::O4:
            return "O4";
        case AtomType::N4:
            return "N4";
        case AtomType::C5M:
            return "C5M";
        case AtomType::C7:
            return "C7";
        // Sugar atoms
        case AtomType::C1_PRIME:
            return "C1'";
        case AtomType::C2_PRIME:
            return "C2'";
        case AtomType::C3_PRIME:
            return "C3'";
        case AtomType::C4_PRIME:
            return "C4'";
        case AtomType::C5_PRIME:
            return "C5'";
        case AtomType::O2_PRIME:
            return "O2'";
        case AtomType::O3_PRIME:
            return "O3'";
        case AtomType::O4_PRIME:
            return "O4'";
        case AtomType::O5_PRIME:
            return "O5'";
        // Backbone atoms
        case AtomType::P:
            return "P";
        case AtomType::OP1:
            return "OP1";
        case AtomType::OP2:
            return "OP2";
        case AtomType::OP3:
            return "OP3";
        // Amino acid backbone
        case AtomType::N:
            return "N";
        case AtomType::CA:
            return "CA";
        case AtomType::C:
            return "C";
        case AtomType::O:
            return "O";
        case AtomType::OXT:
            return "OXT";
        // Common side chain atoms
        case AtomType::CB:
            return "CB";
        case AtomType::CG:
            return "CG";
        case AtomType::CG1:
            return "CG1";
        case AtomType::CG2:
            return "CG2";
        case AtomType::CD:
            return "CD";
        case AtomType::CD1:
            return "CD1";
        case AtomType::CD2:
            return "CD2";
        case AtomType::CE:
            return "CE";
        case AtomType::CE1:
            return "CE1";
        case AtomType::CE2:
            return "CE2";
        case AtomType::CE3:
            return "CE3";
        case AtomType::CZ:
            return "CZ";
        case AtomType::CZ2:
            return "CZ2";
        case AtomType::CZ3:
            return "CZ3";
        case AtomType::CH2:
            return "CH2";
        case AtomType::OG:
            return "OG";
        case AtomType::OG1:
            return "OG1";
        case AtomType::OD1:
            return "OD1";
        case AtomType::OD2:
            return "OD2";
        case AtomType::OE1:
            return "OE1";
        case AtomType::OE2:
            return "OE2";
        case AtomType::OH:
            return "OH";
        case AtomType::ND1:
            return "ND1";
        case AtomType::ND2:
            return "ND2";
        case AtomType::NE:
            return "NE";
        case AtomType::NE1:
            return "NE1";
        case AtomType::NE2:
            return "NE2";
        case AtomType::NH1:
            return "NH1";
        case AtomType::NH2:
            return "NH2";
        case AtomType::NZ:
            return "NZ";
        case AtomType::SD:
            return "SD";
        case AtomType::SG:
            return "SG";
        // Water
        case AtomType::OW:
            return "OW";
        case AtomType::COUNT:
            return "COUNT";
    }
    return "UNKNOWN";
}

/**
 * @enum HBondRole
 * @brief Classification of hydrogen bonding capability
 */
enum class HBondRole {
    UNKNOWN,
    DONOR,    ///< Can only donate H-bond (e.g., N-H)
    ACCEPTOR, ///< Can only accept H-bond (e.g., C=O)
    BOTH,     ///< Can donate and accept (e.g., O-H)
    NONE      ///< Cannot participate in H-bonding
};

/**
 * @brief Convert ElementType to string
 */
[[nodiscard]] inline const char* to_string(ElementType type) {
    switch (type) {
        case ElementType::UNKNOWN:
            return "UNKNOWN";
        case ElementType::CARBON:
            return "C";
        case ElementType::NITROGEN:
            return "N";
        case ElementType::OXYGEN:
            return "O";
        case ElementType::HYDROGEN:
            return "H";
        case ElementType::PHOSPHORUS:
            return "P";
        case ElementType::SULFUR:
            return "S";
        case ElementType::FLUORINE:
            return "F";
        case ElementType::CHLORINE:
            return "CL";
        case ElementType::BROMINE:
            return "BR";
        case ElementType::IODINE:
            return "I";
        case ElementType::SODIUM:
            return "NA";
        case ElementType::POTASSIUM:
            return "K";
        case ElementType::MAGNESIUM:
            return "MG";
        case ElementType::CALCIUM:
            return "CA";
        case ElementType::ZINC:
            return "ZN";
        case ElementType::IRON:
            return "FE";
        case ElementType::MANGANESE:
            return "MN";
        case ElementType::COPPER:
            return "CU";
        case ElementType::COBALT:
            return "CO";
        case ElementType::NICKEL:
            return "NI";
        case ElementType::SELENIUM:
            return "SE";
        case ElementType::CADMIUM:
            return "CD";
        case ElementType::BARIUM:
            return "BA";
        case ElementType::STRONTIUM:
            return "SR";
        case ElementType::RUBIDIUM:
            return "RB";
        case ElementType::CESIUM:
            return "CS";
        case ElementType::LITHIUM:
            return "LI";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert AtomLocation to string
 */
[[nodiscard]] inline const char* to_string(AtomLocation loc) {
    switch (loc) {
        case AtomLocation::UNKNOWN:
            return "UNKNOWN";
        case AtomLocation::BACKBONE:
            return "BACKBONE";
        case AtomLocation::SUGAR:
            return "SUGAR";
        case AtomLocation::NUCLEOBASE:
            return "NUCLEOBASE";
        case AtomLocation::MAINCHAIN:
            return "MAINCHAIN";
        case AtomLocation::SIDECHAIN:
            return "SIDECHAIN";
        case AtomLocation::SOLVENT:
            return "SOLVENT";
    }
    return "UNKNOWN";
}

/**
 * @brief Convert HBondRole to string
 */
[[nodiscard]] inline const char* to_string(HBondRole role) {
    switch (role) {
        case HBondRole::UNKNOWN:
            return "UNKNOWN";
        case HBondRole::DONOR:
            return "DONOR";
        case HBondRole::ACCEPTOR:
            return "ACCEPTOR";
        case HBondRole::BOTH:
            return "BOTH";
        case HBondRole::NONE:
            return "NONE";
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
        case ElementType::CARBON:
            return 1;
        case ElementType::OXYGEN:
            return 2;
        case ElementType::HYDROGEN:
            return 3;
        case ElementType::NITROGEN:
            return 4;
        case ElementType::SULFUR:
            return 5;
        case ElementType::PHOSPHORUS:
            return 6;
        default:
            return 0;
    }
}

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::AtomLocation;
using typing::AtomType;
using typing::ElementType;
using typing::HBondRole;
using typing::is_backbone_atom;
using typing::is_exocyclic_atom;
using typing::is_purine_ring_atom;
using typing::is_ring_atom;
using typing::is_sugar_atom;
using typing::NUM_RING_ATOM_TYPES;
using typing::RING_ATOM_TYPES;

// Backward compatibility alias
using StandardAtom = AtomType;

} // namespace core
} // namespace x3dna
