/**
 * @file structure_legacy_order.hpp
 * @brief Utilities for ordering residues in legacy (PDB file) order
 *
 * Legacy counts residues by processing atoms in PDB file order and grouping
 * by (ResName, ChainID, ResSeq, insertion). This provides functions to get
 * residues in the same order as legacy.
 */

#pragma once

#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <x3dna/core/residue.hpp>

namespace x3dna {
namespace core {

// Forward declaration
class Structure;

/**
 * @brief Get all residues in legacy order (PDB file order)
 *
 * Legacy counts residues by:
 * 1. Processing atoms in PDB file order (by line_number)
 * 2. Grouping atoms with same (ResName, ChainID, ResSeq, insertion)
 * 3. Counting unique groups in order they first appear
 *
 * @param structure The structure to get residues from
 * @return Vector of residue pointers in legacy order (non-owning)
 */
std::vector<const Residue*> get_residues_in_legacy_order(const Structure& structure);

/**
 * @brief Get residue by legacy index (1-based)
 *
 * Finds the residue that would be at the given legacy index when counting
 * in legacy order (PDB file order).
 *
 * @param structure The structure to search
 * @param legacy_idx Legacy residue index (1-based)
 * @return Pointer to residue, or nullptr if not found
 */
const Residue* get_residue_by_legacy_idx(const Structure& structure, int legacy_idx);

/**
 * @brief Get legacy index for a residue
 *
 * Returns the legacy index (1-based) for the given residue when counting
 * in legacy order (PDB file order).
 *
 * @param structure The structure containing the residue
 * @param residue The residue to find the index for
 * @return Legacy residue index (1-based), or 0 if not found
 */
int get_legacy_idx_for_residue(const Structure& structure, const Residue* residue);

} // namespace core
} // namespace x3dna
