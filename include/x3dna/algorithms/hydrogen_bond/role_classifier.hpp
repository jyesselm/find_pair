/**
 * @file role_classifier.hpp
 * @brief H-bond donor/acceptor role classification for all molecule types
 */

#pragma once

#include <string>
#include <vector>
#include <x3dna/core/hbond_types.hpp>
#include <x3dna/core/hbond.hpp>
#include <x3dna/core/typing/molecule_type.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @brief Classifies donor/acceptor roles for H-bond atoms
 *
 * Supports:
 * - Nucleotides (A, C, G, T, U, I) - legacy compatible
 * - Proteins (all 20 amino acids)
 * - Ligands (common patterns, element-based fallback)
 */
class HBondRoleClassifier {
public:
    // === Nucleotide classification (legacy compatible) ===

    /**
     * @brief Get atom role for a nucleotide base
     * @param base One-letter code (A, C, G, T, U, I)
     * @param atom_name Atom name (4-char or trimmed)
     * @return Role (DONOR, ACCEPTOR, EITHER, UNKNOWN)
     */
    [[nodiscard]] static core::HBondAtomRole get_nucleotide_atom_role(char base, const std::string& atom_name);

    /**
     * @brief Classify bond for nucleotide-nucleotide interaction
     * @param base1 First base one-letter code
     * @param base2 Second base one-letter code
     * @param atom1 First atom name
     * @param atom2 Second atom name
     * @return Classification (STANDARD, NON_STANDARD, INVALID)
     */
    [[nodiscard]] static core::HBondClassification classify_nucleotide_bond(char base1, char base2,
                                                                            const std::string& atom1,
                                                                            const std::string& atom2);

    // === Protein classification ===

    /**
     * @brief Get atom role for a protein residue
     * @param residue_name Three-letter code (ALA, ARG, etc.)
     * @param atom_name Atom name
     * @return Role (DONOR, ACCEPTOR, EITHER, UNKNOWN)
     */
    [[nodiscard]] static core::HBondAtomRole get_protein_atom_role(const std::string& residue_name,
                                                                   const std::string& atom_name);

    /**
     * @brief Check if protein mainchain atom
     */
    [[nodiscard]] static bool is_mainchain_atom(const std::string& atom_name);

    // === Element-based fallback ===

    /**
     * @brief Get atom role based on element (fallback for unknown atoms)
     * @param atom_name Atom name (element extracted from first 1-2 chars)
     * @return Role (EITHER for N/O/S, ACCEPTOR for F, UNKNOWN otherwise)
     */
    [[nodiscard]] static core::HBondAtomRole get_element_based_role(const std::string& atom_name);

    // === Ligand classification ===

    /**
     * @brief Get atom role for a ligand (element-based heuristic)
     * @param atom_name Atom name
     * @param element Element symbol (optional, derived from name if not provided)
     * @return Role (DONOR, ACCEPTOR, EITHER, UNKNOWN)
     */
    [[nodiscard]] static core::HBondAtomRole get_ligand_atom_role(const std::string& atom_name,
                                                                  const std::string& element = "");

    // === General classification ===

    /**
     * @brief Get atom role based on molecule type
     * @param molecule_type Type of molecule
     * @param residue_name Residue name (3-letter for protein, 1-letter for nucleotide)
     * @param atom_name Atom name
     * @return Role (DONOR, ACCEPTOR, EITHER, UNKNOWN)
     */
    [[nodiscard]] static core::HBondAtomRole get_atom_role(core::MoleculeType molecule_type,
                                                           const std::string& residue_name,
                                                           const std::string& atom_name);

    /**
     * @brief Classify H-bond between any two atoms
     * @param role1 Role of first atom
     * @param role2 Role of second atom
     * @return Classification based on donor-acceptor compatibility
     */
    [[nodiscard]] static core::HBondClassification classify_by_roles(core::HBondAtomRole role1,
                                                                     core::HBondAtomRole role2);

    // === Legacy compatibility (same as modern-backup) ===

    [[nodiscard]] static core::HBondAtomRole get_atom_role(char base, const std::string& atom_name) {
        return get_nucleotide_atom_role(base, atom_name);
    }

    [[nodiscard]] static core::HBondClassification classify_bond(char base1, char base2, const std::string& atom1,
                                                                 const std::string& atom2) {
        return classify_nucleotide_bond(base1, base2, atom1, atom2);
    }

    // === Utility methods ===

    [[nodiscard]] static bool is_good_hbond_distance(double distance, double min_dist = 2.5, double max_dist = 3.5);

    [[nodiscard]] static int count_good_hbonds(const std::vector<core::HBond>& bonds, double min_dist = 2.5,
                                               double max_dist = 3.5);
};

} // namespace algorithms
} // namespace x3dna
