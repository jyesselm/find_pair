/**
 * @file ring_data_cache.hpp
 * @brief Cache for pre-computed ring atom data per residue
 *
 * Optimizes overlap calculation by caching ring atom indices and
 * exocyclic atom mapping, avoiding repeated O(n) lookups.
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {
namespace validation {

/**
 * @brief Pre-computed ring data for a single residue
 *
 * Stores atom indices (not pointers) for cache safety across residue copies.
 * Computed once per residue, reused for all pairs involving that residue.
 */
struct ResidueRingData {
    std::vector<size_t> ring_atom_indices;      // Indices of ring atoms in residue.atoms()
    std::vector<size_t> exocyclic_atom_indices; // For each ring atom, index of its exocyclic partner
                                                 // (same index as ring_atom if no exocyclic found)
    bool is_purine = false;                      // True if 9 ring atoms, false if 6
    bool is_valid = false;                       // True if at least 3 ring atoms found
};

/**
 * @brief Cache for ResidueRingData keyed by residue res_id
 *
 * Thread-safe caching of ring atom data. The cache key is res_id string
 * to handle residue copies safely.
 *
 * Usage:
 *   RingDataCache cache;
 *   const auto& data = cache.get_or_compute(residue);
 *   auto coords = cache.get_ring_coords(residue, oave);
 */
class RingDataCache {
public:
    RingDataCache() = default;

    /**
     * @brief Get or compute ring data for a residue
     * @param residue The residue to get ring data for
     * @return Cached or newly computed ResidueRingData
     */
    const ResidueRingData& get_or_compute(const core::Residue& residue);

    /**
     * @brief Get ring coordinates relative to oave
     * @param residue The residue
     * @param oave The average origin point
     * @return Vector of ring coordinates (using exocyclic atoms where available)
     *
     * This combines cached ring data with runtime oave to produce final coordinates.
     */
    std::vector<geometry::Vector3D> get_ring_coords(const core::Residue& residue,
                                                     const geometry::Vector3D& oave);

    /**
     * @brief Clear all cached data
     */
    void clear();

    /**
     * @brief Get number of cached entries
     */
    [[nodiscard]] size_t size() const { return cache_.size(); }

private:
    /**
     * @brief Compute ring data for a residue (internal)
     */
    static ResidueRingData compute_ring_data(const core::Residue& residue);

    std::unordered_map<std::string, ResidueRingData> cache_; // Keyed by res_id
};

} // namespace validation
} // namespace algorithms
} // namespace x3dna
