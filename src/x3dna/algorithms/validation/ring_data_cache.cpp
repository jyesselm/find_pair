/**
 * @file ring_data_cache.cpp
 * @brief Implementation of RingDataCache
 */

#include <x3dna/algorithms/validation/ring_data_cache.hpp>
#include <x3dna/algorithms/validation_constants.hpp>
#include <x3dna/core/typing/atom_type.hpp>

namespace x3dna {
namespace algorithms {
namespace validation {

using core::StandardAtom;
using core::RING_ATOM_TYPES;
using core::NUM_RING_ATOM_TYPES;

const ResidueRingData& RingDataCache::get_or_compute(const core::Residue& residue) {
    const std::string& key = residue.res_id();

    auto it = cache_.find(key);
    if (it != cache_.end()) {
        return it->second;
    }

    // Compute and cache
    auto result = cache_.emplace(key, compute_ring_data(residue));
    return result.first->second;
}

std::vector<geometry::Vector3D> RingDataCache::get_ring_coords(const core::Residue& residue,
                                                                const geometry::Vector3D& oave) {
    const auto& data = get_or_compute(residue);

    std::vector<geometry::Vector3D> coords;
    coords.reserve(data.ring_atom_indices.size());

    const auto& atoms = residue.atoms();

    for (size_t i = 0; i < data.ring_atom_indices.size(); ++i) {
        size_t exo_idx = data.exocyclic_atom_indices[i];
        if (exo_idx < atoms.size()) {
            coords.push_back(atoms[exo_idx].position() - oave);
        }
    }

    return coords;
}

void RingDataCache::clear() {
    cache_.clear();
}

ResidueRingData RingDataCache::compute_ring_data(const core::Residue& residue) {
    ResidueRingData data;

    const auto& atoms = residue.atoms();

    // Find ring atom indices using enum comparison (O(1) per atom instead of string compare)
    for (size_t r = 0; r < NUM_RING_ATOM_TYPES; ++r) {
        StandardAtom target_type = RING_ATOM_TYPES[r];
        for (size_t i = 0; i < atoms.size(); ++i) {
            if (atoms[i].standard_atom() == target_type) {
                data.ring_atom_indices.push_back(i);
                break;
            }
        }
    }

    // Determine if purine (need at least N7, C8, N9 which are indices 6,7,8)
    data.is_purine = (data.ring_atom_indices.size() >= NUM_RING_ATOM_TYPES);
    data.is_valid = (data.ring_atom_indices.size() >= 3);

    if (!data.is_valid) {
        return data;
    }

    // For each ring atom, find exocyclic partner
    data.exocyclic_atom_indices.reserve(data.ring_atom_indices.size());

    for (size_t ring_idx : data.ring_atom_indices) {
        const auto& ring_atom = atoms[ring_idx];
        size_t best_exo_idx = ring_idx; // Default to ring atom itself
        double min_dist = validation_constants::BOND_DISTANCE;

        for (size_t i = 0; i < atoms.size(); ++i) {
            const auto& atom = atoms[i];

            // Skip ring atoms using enum check (O(1) instead of string comparison)
            if (core::is_ring_atom(atom.standard_atom())) {
                continue;
            }

            // Skip hydrogen atoms (check element type via first char - still fast)
            if (!atom.name().empty() && atom.name()[0] == 'H') {
                continue;
            }

            double dist = (atom.position() - ring_atom.position()).length();
            if (dist < min_dist && dist > validation_constants::MIN_ATOM_DISTANCE) {
                min_dist = dist;
                best_exo_idx = i;
            }
        }

        data.exocyclic_atom_indices.push_back(best_exo_idx);
    }

    return data;
}

} // namespace validation
} // namespace algorithms
} // namespace x3dna
