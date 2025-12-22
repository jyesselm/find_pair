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
 */
enum class HBondContext {
    UNKNOWN,
    BASE_BASE,         // Between nucleotide bases (Watson-Crick, etc.)
    BASE_BACKBONE,     // Base atom to phosphate/sugar backbone
    BACKBONE_BACKBONE, // Between backbone atoms
    BASE_SUGAR,        // Base to ribose sugar (O2', O3', O4')
    SUGAR_SUGAR        // Between sugar atoms
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
 * @brief Role of an atom in H-bond (from donor_acceptor lookup)
 */
enum class HBondAtomRole {
    DONOR,    // 'D' - Has hydrogen to donate
    ACCEPTOR, // 'A' - Has lone pair to accept
    EITHER,   // 'X' - Can act as either (ring N atoms)
    UNKNOWN   // '?' - Not in lookup table
};

} // namespace core
} // namespace x3dna
