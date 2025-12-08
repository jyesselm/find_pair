/**
 * @file residue_factory.hpp
 * @brief Factory for creating Residue objects with all properties initialized
 */

#pragma once

#include <string>
#include <vector>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/residue_type.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>

namespace x3dna {
namespace core {

/**
 * @brief Factory for creating Residue objects
 *
 * Centralizes residue creation logic and uses ModifiedNucleotideRegistry
 * to determine properties (one_letter_code, type, is_purine) at creation time.
 *
 * This ensures Residue objects are created with all properties set, rather
 * than computing them on-demand.
 */
class ResidueFactory {
public:
    /**
     * @brief Create a residue from PDB data
     *
     * @param name Three-letter residue name (e.g., "ALA", "ATP", "J48")
     * @param sequence_number PDB sequence number
     * @param chain_id Chain identifier
     * @param insertion_code Insertion code
     * @param atoms Vector of atoms in this residue
     * @return Residue with all properties initialized
     */
    static Residue create(const std::string& name, int sequence_number, char chain_id, char insertion_code,
                          const std::vector<Atom>& atoms);

private:
    /**
     * @brief Determine one-letter code from residue name
     * Uses ModifiedNucleotideRegistry for modified nucleotides
     */
    static char determine_one_letter_code(const std::string& name);

    /**
     * @brief Determine residue type from name and one-letter code
     */
    static ResidueType determine_type(const std::string& name, char one_letter_code);

    /**
     * @brief Determine if residue is a purine
     */
    static bool determine_is_purine(const std::string& name, ResidueType type);
};

} // namespace core
} // namespace x3dna
