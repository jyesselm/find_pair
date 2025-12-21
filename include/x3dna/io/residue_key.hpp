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
    std::string residue_name; ///< Residue name (e.g., "A", "G", "PSU")
    std::string chain_id;     ///< Chain identifier (string for CIF compatibility)
    int residue_seq;          ///< Residue sequence number
    std::string insertion_code; ///< Insertion code (usually empty)
    char record_type;         ///< PDB record type: 'A' for ATOM, 'H' for HETATM

    /**
     * @brief Comparison operator for map ordering
     */
    bool operator<(const ResidueKey& other) const {
        return std::tie(chain_id, residue_seq, insertion_code, residue_name, record_type) <
               std::tie(other.chain_id, other.residue_seq, other.insertion_code, other.residue_name, other.record_type);
    }

    /**
     * @brief Equality comparison
     */
    bool operator==(const ResidueKey& other) const {
        return residue_name == other.residue_name && chain_id == other.chain_id && residue_seq == other.residue_seq &&
               insertion_code == other.insertion_code && record_type == other.record_type;
    }
};

} // namespace io
} // namespace x3dna
