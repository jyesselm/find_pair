/**
 * @file residue_index_map.hpp
 * @brief Maps between legacy 1-based and modern 0-based residue indices
 */

#pragma once

#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <map>
#include <vector>
#include <optional>

namespace x3dna {
namespace algorithms {

/**
 * @class ResidueIndexMap
 * @brief Manages mapping between legacy 1-based and modern 0-based residue indices
 *
 * The legacy X3DNA code uses 1-based residue indices stored during PDB parsing.
 * This class provides a clean abstraction for:
 * - Looking up residues by legacy or modern index
 * - Converting between index systems
 * - Iterating over residues in legacy order (required for compatibility)
 *
 * Usage:
 * @code
 * ResidueIndexMap index_map;
 * index_map.build(structure);
 *
 * // Get residue by legacy index
 * const Residue* res = index_map.get_by_legacy_idx(42);
 *
 * // Iterate in legacy order
 * for (int legacy_idx : index_map.legacy_indices()) {
 *     const Residue* res = index_map.get_by_legacy_idx(legacy_idx);
 *     // ...
 * }
 * @endcode
 */
class ResidueIndexMap {
public:
    /**
     * @brief Build mapping from structure
     * @param structure Structure to index (residues must have legacy indices set via atoms)
     */
    void build(const core::Structure& structure);

    /**
     * @brief Clear all mappings
     */
    void clear();

    /**
     * @brief Check if mapping is empty
     */
    [[nodiscard]] bool empty() const { return by_legacy_.empty(); }

    /**
     * @brief Get number of mapped residues
     */
    [[nodiscard]] size_t size() const { return by_legacy_.size(); }

    // ==================== Lookups ====================

    /**
     * @brief Get residue by legacy 1-based index
     * @param legacy_idx Legacy index (from atom.legacy_residue_idx())
     * @return Pointer to residue, or nullptr if not found
     */
    [[nodiscard]] const core::Residue* get_by_legacy_idx(int legacy_idx) const;

    /**
     * @brief Get residue by modern 0-based index
     * @param modern_idx Modern index (position in structure)
     * @return Pointer to residue, or nullptr if not found
     */
    [[nodiscard]] const core::Residue* get_by_modern_idx(size_t modern_idx) const;

    /**
     * @brief Check if legacy index exists in map
     */
    [[nodiscard]] bool has_legacy_idx(int legacy_idx) const;

    /**
     * @brief Check if modern index exists in map
     */
    [[nodiscard]] bool has_modern_idx(size_t modern_idx) const;

    // ==================== Conversions ====================

    /**
     * @brief Convert legacy 1-based index to modern 0-based index
     * @param legacy_idx Legacy index
     * @return Modern index, or nullopt if not found
     */
    [[nodiscard]] std::optional<size_t> to_modern(int legacy_idx) const;

    /**
     * @brief Convert modern 0-based index to legacy 1-based index
     * @param modern_idx Modern index
     * @return Legacy index, or nullopt if not found
     */
    [[nodiscard]] std::optional<int> to_legacy(size_t modern_idx) const;

    // ==================== Range Info ====================

    /**
     * @brief Get maximum legacy index
     */
    [[nodiscard]] int max_legacy_idx() const { return max_legacy_idx_; }

    /**
     * @brief Get minimum legacy index (usually 1)
     */
    [[nodiscard]] int min_legacy_idx() const { return min_legacy_idx_; }

    // ==================== Iteration ====================

    /**
     * @brief Get all legacy indices in ascending order
     * @return Vector of legacy indices
     */
    [[nodiscard]] std::vector<int> legacy_indices() const;

    /**
     * @brief Get all legacy indices of nucleotide residues
     * @param nucleotide_checker Function to check if residue is nucleotide
     * @return Vector of legacy indices for nucleotides only
     */
    template <typename NucleotideChecker>
    [[nodiscard]] std::vector<int> nucleotide_legacy_indices(NucleotideChecker&& checker) const {
        std::vector<int> result;
        for (const auto& [legacy_idx, residue] : by_legacy_) {
            if (residue && checker(*residue)) {
                result.push_back(legacy_idx);
            }
        }
        return result;
    }

    /**
     * @brief Iterate over all (legacy_idx, residue) pairs
     */
    [[nodiscard]] const std::map<int, const core::Residue*>& all() const { return by_legacy_; }

private:
    std::map<int, const core::Residue*> by_legacy_;
    std::map<size_t, const core::Residue*> by_modern_;
    std::map<int, size_t> legacy_to_modern_;
    std::map<size_t, int> modern_to_legacy_;
    int max_legacy_idx_ = 0;
    int min_legacy_idx_ = 0;
};

}  // namespace algorithms
}  // namespace x3dna

