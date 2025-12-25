/**
 * @file slot_cache.hpp
 * @brief Per-residue caching of computed H and LP slots
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <x3dna/algorithms/hydrogen_bond/slot/slot.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/slot_predictor.hpp>
#include <x3dna/core/residue.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @class SlotCache
 * @brief Caches computed H and LP slots for a residue
 *
 * Slots are computed lazily on first access and cached for reuse.
 * The cache should be cleared between optimization runs to reset
 * slot usage tracking.
 */
class SlotCache {
public:
    /**
     * @brief Construct cache for a residue
     * @param residue Reference to the residue (must outlive cache)
     * @param base_type Single letter base type (A, G, C, U, T)
     */
    SlotCache(const core::Residue& residue, char base_type);

    /**
     * @brief Get H slots for a donor atom (computes if not cached)
     * @param atom_name Donor atom name
     * @return Reference to vector of HSlots
     */
    [[nodiscard]] std::vector<HSlot>& get_h_slots(const std::string& atom_name);

    /**
     * @brief Get LP slots for an acceptor atom (computes if not cached)
     * @param atom_name Acceptor atom name
     * @return Reference to vector of LPSlots
     */
    [[nodiscard]] std::vector<LPSlot>& get_lp_slots(const std::string& atom_name);

    /**
     * @brief Reset all slots to unused state (clears bond_directions)
     */
    void reset_slots();

    /**
     * @brief Clear entire cache (slots will be recomputed on next access)
     */
    void clear();

    /** @brief Get the base type */
    [[nodiscard]] char base_type() const { return base_type_; }

    /** @brief Get the residue reference */
    [[nodiscard]] const core::Residue& residue() const { return residue_; }

private:
    const core::Residue& residue_;
    char base_type_;
    geometry::Vector3D base_normal_;
    bool normal_computed_ = false;

    std::unordered_map<std::string, std::vector<HSlot>> h_slots_;
    std::unordered_map<std::string, std::vector<LPSlot>> lp_slots_;

    void ensure_base_normal();
};

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
