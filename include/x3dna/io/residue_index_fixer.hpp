/**
 * @file residue_index_fixer.hpp
 * @brief Helper to fix residue legacy indices by matching with legacy JSON
 */

#ifndef X3DNA_IO_RESIDUE_INDEX_FIXER_HPP
#define X3DNA_IO_RESIDUE_INDEX_FIXER_HPP

#include <string>
#include <map>
#include <tuple>
#include <x3dna/core/structure.hpp>
#include <nlohmann/json.hpp>

namespace x3dna {
namespace io {

/**
 * @brief Fix residue legacy indices by matching with legacy JSON
 * 
 * This function matches residues by PDB properties (residue_name, chain_id, 
 * residue_seq, insertion) and assigns legacy indices from the JSON file.
 * 
 * @param structure The structure to fix
 * @param legacy_json_file Path to legacy JSON file (e.g., base_frame_calc)
 * @return Number of residues that were matched and fixed
 */
int fix_residue_indices_from_json(core::Structure& structure, 
                                   const std::string& legacy_json_file);

} // namespace io
} // namespace x3dna

#endif // X3DNA_IO_RESIDUE_INDEX_FIXER_HPP

