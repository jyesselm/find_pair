/**
 * @file hbond_types.hpp
 * @brief Core types for hydrogen bond representation and classification
 */

#pragma once

namespace x3dna {
namespace core {

/**
 * @brief Classification of H-bond validity based on donor/acceptor analysis
 *
 * Legacy mapping:
 *   STANDARD     -> '-'  (valid donor-acceptor relationship)
 *   NON_STANDARD -> '*'  (atoms can H-bond but role unclear)
 *   INVALID      -> ' '  (failed validation)
 */
enum class HBondClassification {
    UNKNOWN,      // Not yet classified
    STANDARD,     // '-' - Valid donor-acceptor confirmed
    NON_STANDARD, // '*' - Can form H-bond, role ambiguous
    INVALID       // ' ' - Failed validation or filtered
};

/**
 * @brief Context describing what structural elements the H-bond connects
 * Extended to support proteins and ligands
 */
enum class HBondContext {
    UNKNOWN,
    // Nucleic acid contexts
    BASE_BASE,         // Between nucleotide bases
    BASE_BACKBONE,     // Base to phosphate backbone
    BACKBONE_BACKBONE, // Between backbone atoms
    BASE_SUGAR,        // Base to ribose sugar
    SUGAR_SUGAR,       // Between sugar atoms
    // Protein contexts
    PROTEIN_MAINCHAIN, // Protein backbone N-H...O=C
    PROTEIN_SIDECHAIN, // Protein sidechain donors/acceptors
    // Cross-molecule contexts
    BASE_PROTEIN,     // Nucleic acid base to protein
    SUGAR_PROTEIN,    // Sugar to protein
    BACKBONE_PROTEIN, // NA backbone to protein
    // Ligand contexts
    BASE_LIGAND,    // Base to ligand
    PROTEIN_LIGAND, // Protein to ligand
    LIGAND_LIGAND   // Between ligands
};

/**
 * @brief State from conflict resolution algorithm
 *
 * When multiple H-bonds share the same atom, the shortest wins.
 * This tracks each bond's relationship to that process.
 */
enum class ConflictState {
    NO_CONFLICT,                 // 0 - Not involved in any conflict
    SHARES_DONOR_WITH_WINNER,    // 1 - Another bond using same donor won
    SHARES_ACCEPTOR_WITH_WINNER, // 2 - Another bond using same acceptor won
    SHARES_BOTH_WITH_WINNER,     // 3 - Shares both atoms (rare)
    IS_CONFLICT_WINNER           // 18 - This bond won the conflict
};

/**
 * @brief Role of an atom in H-bond
 */
enum class HBondAtomRole {
    DONOR,    // Has hydrogen to donate
    ACCEPTOR, // Has lone pair to accept
    EITHER,   // Can act as either
    UNKNOWN   // Not in lookup table
};

/**
 * @brief Type of molecular interaction for filtering
 */
enum class HBondInteractionType {
    BASE_BASE = 1 << 0,       // Nucleic acid base-base
    BASE_BACKBONE = 1 << 1,   // Base to NA backbone
    BASE_SUGAR = 1 << 2,      // Base to sugar
    BASE_PROTEIN = 1 << 3,    // Base to protein
    BASE_LIGAND = 1 << 4,     // Base to ligand
    PROTEIN_PROTEIN = 1 << 5, // Protein-protein
    PROTEIN_LIGAND = 1 << 6,  // Protein to ligand
    RNA_INTERNAL = 1 << 7,    // All within RNA (backbone, sugar, base)
    ANY = 0xFFFF              // All interactions
};

// Allow bitwise operations on HBondInteractionType
inline HBondInteractionType operator|(HBondInteractionType a, HBondInteractionType b) {
    return static_cast<HBondInteractionType>(static_cast<int>(a) | static_cast<int>(b));
}

inline bool operator&(HBondInteractionType a, HBondInteractionType b) {
    return (static_cast<int>(a) & static_cast<int>(b)) != 0;
}

// Conversion helpers
[[nodiscard]] inline const char* to_string(HBondClassification c) {
    switch (c) {
        case HBondClassification::UNKNOWN:
            return "UNKNOWN";
        case HBondClassification::STANDARD:
            return "STANDARD";
        case HBondClassification::NON_STANDARD:
            return "NON_STANDARD";
        case HBondClassification::INVALID:
            return "INVALID";
    }
    return "UNKNOWN";
}

[[nodiscard]] inline char to_legacy_char(HBondClassification c) {
    switch (c) {
        case HBondClassification::STANDARD:
            return '-';
        case HBondClassification::NON_STANDARD:
            return '*';
        default:
            return ' ';
    }
}

} // namespace core
} // namespace x3dna
