/**
 * @file atom_capacity.hpp
 * @brief Donor and acceptor capacity tables for H-bond optimization
 *
 * Defines how many hydrogen atoms a donor can provide and how many
 * lone pairs an acceptor can use for H-bonding.
 */

#pragma once

#include <string>
#include <unordered_map>
#include <optional>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @class AtomCapacity
 * @brief Static utility class for looking up donor/acceptor capacities
 *
 * Provides capacity lookup for standard nucleotide atoms with fallback
 * to parent base type for modified residues.
 */
class AtomCapacity {
public:
    /**
     * @brief Get donor capacity (number of H atoms that can donate)
     * @param residue_code Residue code (e.g., "A", "G", "5MC")
     * @param atom_name Atom name (e.g., "N6", "N1")
     * @return Number of hydrogens that can donate (0 if not a donor)
     */
    [[nodiscard]] static int get_donor_capacity(const std::string& residue_code,
                                                 const std::string& atom_name);

    /**
     * @brief Get acceptor capacity (number of lone pairs that can accept)
     * @param residue_code Residue code (e.g., "A", "G", "5MC")
     * @param atom_name Atom name (e.g., "O6", "N1")
     * @return Number of lone pairs that can accept (0 if not an acceptor)
     */
    [[nodiscard]] static int get_acceptor_capacity(const std::string& residue_code,
                                                    const std::string& atom_name);

    /**
     * @brief Get parent base type for a residue code
     * @param residue_code Residue code (e.g., "5MC", "PSU")
     * @return Parent base type (A, G, C, U, T) or empty if unknown
     *
     * Used to look up capacities for modified bases by falling back
     * to their parent nucleotide.
     */
    [[nodiscard]] static std::optional<char> get_parent_base_type(const std::string& residue_code);

    /**
     * @brief Normalize atom name (handle OP1/O1P variants)
     * @param atom_name Raw atom name
     * @return Normalized atom name
     */
    [[nodiscard]] static std::string normalize_atom_name(const std::string& atom_name);

    /**
     * @brief Check if an atom is a backbone atom
     * @param atom_name Atom name
     * @return True if backbone (P, OP1, OP2, O3', O5')
     */
    [[nodiscard]] static bool is_backbone_atom(const std::string& atom_name);

private:
    // Implementation details are in the .cpp file
    AtomCapacity() = delete; // Static-only class
};

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
