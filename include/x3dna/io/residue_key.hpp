/**
 * @file residue_key.hpp
 * @brief Key for identifying unique residues during parsing
 */

#pragma once

#include <string>
#include <tuple>

namespace x3dna {
namespace io {

/**
 * @struct ResidueKey
 * @brief Unique identifier for a residue during parsing
 *
 * Used to group atoms into residues based on residue name, chain ID,
 * sequence number, and insertion code.
 */
struct ResidueKey {
    std::string residue_name;   ///< Residue name (e.g., "A", "G", "PSU")
    char chain_id;              ///< Chain identifier
    int residue_seq;            ///< Residue sequence number
    char insertion_code;        ///< Insertion code (usually ' ')

    /**
     * @brief Comparison operator for map ordering
     */
    bool operator<(const ResidueKey& other) const {
        return std::tie(chain_id, residue_seq, insertion_code, residue_name) <
               std::tie(other.chain_id, other.residue_seq, other.insertion_code, other.residue_name);
    }

    /**
     * @brief Equality comparison
     */
    bool operator==(const ResidueKey& other) const {
        return residue_name == other.residue_name &&
               chain_id == other.chain_id &&
               residue_seq == other.residue_seq &&
               insertion_code == other.insertion_code;
    }
};

} // namespace io
} // namespace x3dna
