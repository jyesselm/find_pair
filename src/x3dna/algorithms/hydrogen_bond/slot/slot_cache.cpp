/**
 * @file slot_cache.cpp
 * @brief Implementation of SlotCache
 */

#include <x3dna/algorithms/hydrogen_bond/slot/slot_cache.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

SlotCache::SlotCache(const core::Residue& residue, char base_type)
    : residue_(residue)
    , base_type_(base_type) {
}

void SlotCache::ensure_base_normal() {
    if (!normal_computed_) {
        base_normal_ = SlotPredictor::compute_base_normal(residue_);
        normal_computed_ = true;
    }
}

std::vector<HSlot>& SlotCache::get_h_slots(const std::string& atom_name) {
    auto it = h_slots_.find(atom_name);
    if (it != h_slots_.end()) {
        return it->second;
    }

    // Compute and cache
    ensure_base_normal();
    auto slots = SlotPredictor::predict_h_slots(base_type_, atom_name, residue_, base_normal_);
    auto [inserted_it, _] = h_slots_.emplace(atom_name, std::move(slots));
    return inserted_it->second;
}

std::vector<LPSlot>& SlotCache::get_lp_slots(const std::string& atom_name) {
    auto it = lp_slots_.find(atom_name);
    if (it != lp_slots_.end()) {
        return it->second;
    }

    // Compute and cache
    ensure_base_normal();
    auto slots = SlotPredictor::predict_lp_slots(base_type_, atom_name, residue_, base_normal_);
    auto [inserted_it, _] = lp_slots_.emplace(atom_name, std::move(slots));
    return inserted_it->second;
}

void SlotCache::reset_slots() {
    for (auto& [name, slots] : h_slots_) {
        for (auto& slot : slots) {
            slot.reset();
        }
    }
    for (auto& [name, slots] : lp_slots_) {
        for (auto& slot : slots) {
            slot.reset();
        }
    }
}

void SlotCache::clear() {
    h_slots_.clear();
    lp_slots_.clear();
    normal_computed_ = false;
}

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
