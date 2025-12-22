/**
 * @file atom_classification.hpp
 * @brief AtomClassification struct and AtomClassifier for centralized atom type detection
 */

#pragma once

#include <string>
#include <x3dna/core/typing/atom_type.hpp>
#include <x3dna/core/typing/molecule_type.hpp>

namespace x3dna {
namespace core {
namespace typing {

/**
 * @struct AtomClassification
 * @brief Complete classification of an atom
 */
struct AtomClassification {
    ElementType element = ElementType::UNKNOWN;
    AtomLocation location = AtomLocation::UNKNOWN;
    HBondRole hbond_role = HBondRole::UNKNOWN;
    int legacy_element_index = 0;  ///< For backwards compatibility (0-6)
    bool is_ring_atom = false;     ///< Part of nucleobase ring system
};

/**
 * @class AtomClassifier
 * @brief Centralized atom classification logic
 *
 * Replaces scattered atom type checks throughout the codebase with
 * a single, unified classification system.
 */
class AtomClassifier {
public:
    // === Element detection ===

    /**
     * @brief Get ElementType from atom name
     * @param atom_name PDB-format atom name (4 characters, space-padded)
     * @return ElementType enum value
     */
    [[nodiscard]] static ElementType get_element(const std::string& atom_name);

    /**
     * @brief Get legacy element index for backwards compatibility
     * @param atom_name PDB-format atom name
     * @return Element index (0-6): 0=unknown, 1=C, 2=O, 3=H, 4=N, 5=S, 6=P
     */
    [[nodiscard]] static int get_legacy_element_index(const std::string& atom_name);

    // === Location detection (nucleotides) ===

    /**
     * @brief Check if atom is part of nucleotide backbone
     * @param atom_name PDB-format atom name
     * @return true for P, OP1, OP2, O1P, O2P, O5', O3'
     */
    [[nodiscard]] static bool is_backbone_atom(const std::string& atom_name);

    /**
     * @brief Check if atom is part of ribose sugar
     * @param atom_name PDB-format atom name
     * @return true for C1'-C5', O4', O2'
     */
    [[nodiscard]] static bool is_sugar_atom(const std::string& atom_name);

    /**
     * @brief Check if atom is part of the nucleobase
     * @param atom_name PDB-format atom name
     * @return true if not backbone or sugar
     */
    [[nodiscard]] static bool is_nucleobase_atom(const std::string& atom_name);

    /**
     * @brief Check if atom is a ring atom
     * @param atom_name Atom name (will be trimmed)
     * @return true for N1, C2, N3, C4, C5, C6, N7, C8, N9
     */
    [[nodiscard]] static bool is_ring_atom(const std::string& atom_name);

    // === Location detection (proteins) ===

    /**
     * @brief Check if atom is part of protein mainchain
     * @param atom_name PDB-format atom name
     * @return true for N, CA, C, O, OXT
     */
    [[nodiscard]] static bool is_mainchain_atom(const std::string& atom_name);

    /**
     * @brief Check if atom is part of protein sidechain
     * @param atom_name PDB-format atom name
     * @return true if not mainchain and not hydrogen
     */
    [[nodiscard]] static bool is_sidechain_atom(const std::string& atom_name);

    // === H-bond classification ===

    /**
     * @brief Check if atom can participate in hydrogen bonding
     * @param atom_name PDB-format atom name
     * @param allowed_elements Period-delimited elements (default ".O.N.")
     * @return true if atom element is in allowed set
     */
    [[nodiscard]] static bool can_form_hbond(
        const std::string& atom_name,
        const std::string& allowed_elements = ".O.N.");

    /**
     * @brief Check if two atoms can form an H-bond together
     * @param atom1 First atom name
     * @param atom2 Second atom name
     * @param allowed_elements Period-delimited elements (default ".O.N.")
     * @return true if both atoms can participate in H-bonding
     */
    [[nodiscard]] static bool can_form_hbond_pair(
        const std::string& atom1,
        const std::string& atom2,
        const std::string& allowed_elements = ".O.N.");

    /**
     * @brief Check if atom is a base atom for H-bond counting
     * @param atom_name PDB-format atom name
     * @return true if matches legacy pattern (C5M or XD where X not H/P)
     */
    [[nodiscard]] static bool is_base_atom_for_hbond(const std::string& atom_name);

    // === AtomType classification (for fast integer comparison) ===

    /**
     * @brief Get AtomType enum from atom name for fast comparison (without context)
     * @param atom_name Atom name (trimmed or PDB-format)
     * @return AtomType enum value (UNKNOWN if not a standard atom)
     *
     * WARNING: This version assigns AtomType without knowing molecule context.
     * Use get_atom_type(name, molecule_type) when context is available.
     */
    [[nodiscard]] static AtomType get_atom_type(const std::string& atom_name);

    /**
     * @brief Get AtomType enum from atom name with molecule context
     * @param atom_name Atom name (trimmed or PDB-format)
     * @param molecule_type The type of molecule this atom belongs to
     * @return AtomType enum value (UNKNOWN if atom doesn't belong to this molecule type)
     *
     * This is the preferred method - it only assigns nucleotide-specific atom types
     * (N7, C8, N9, sugar atoms, etc.) to nucleic acid atoms, and protein-specific
     * atom types (N, CA, C, O, side chains) to protein atoms.
     */
    [[nodiscard]] static AtomType get_atom_type(const std::string& atom_name, MoleculeType molecule_type);

    // Backward compatibility alias
    [[nodiscard]] static AtomType get_standard_atom(const std::string& atom_name) {
        return get_atom_type(atom_name);
    }

    // === Full classification ===

    /**
     * @brief Get complete classification for an atom with explicit molecule context
     * @param atom_name PDB-format atom name
     * @param molecule_type The type of molecule this atom belongs to
     * @return Full AtomClassification struct
     *
     * This is the preferred method - it requires knowing the context before classifying.
     */
    [[nodiscard]] static AtomClassification classify(
        const std::string& atom_name,
        MoleculeType molecule_type);

    /**
     * @brief Get complete classification for an atom in a nucleotide context
     * @param atom_name PDB-format atom name
     * @return Full AtomClassification struct
     */
    [[nodiscard]] static AtomClassification classify_nucleotide_atom(const std::string& atom_name);

    /**
     * @brief Get complete classification for an atom in a protein context
     * @param atom_name PDB-format atom name
     * @return Full AtomClassification struct
     */
    [[nodiscard]] static AtomClassification classify_protein_atom(const std::string& atom_name);
};

} // namespace typing

// Expose in core namespace for backwards compatibility
using typing::AtomClassification;
using typing::AtomClassifier;

} // namespace core
} // namespace x3dna
