/**
 * @file ring_data_cache.cpp
 * @brief Implementation of RingDataCache
 */

#include <x3dna/algorithms/validation/ring_data_cache.hpp>
#include <x3dna/algorithms/validation_constants.hpp>

namespace x3dna {
namespace algorithms {
namespace validation {

// Ring atom names in order (purines use all 9, pyrimidines use first 6)
static const char* RING_ATOM_NAMES[] = {"C4", "N3", "C2", "N1", "C6", "C5", "N7", "C8", "N9"};
static constexpr size_t NUM_PURINE_ATOMS = 9;

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

    // Find ring atom indices
    for (size_t r = 0; r < NUM_PURINE_ATOMS; ++r) {
        for (size_t i = 0; i < atoms.size(); ++i) {
            if (atoms[i].name() == RING_ATOM_NAMES[r]) {
                data.ring_atom_indices.push_back(i);
                break;
            }
        }
    }

    // Determine if purine (need at least N7, C8, N9 which are indices 6,7,8)
    data.is_purine = (data.ring_atom_indices.size() >= NUM_PURINE_ATOMS);
    data.is_valid = (data.ring_atom_indices.size() >= 3);

    if (!data.is_valid) {
        return data;
    }

    // Build set of ring atom names for exclusion
    std::vector<std::string> ring_names;
    ring_names.reserve(data.ring_atom_indices.size());
    for (size_t idx : data.ring_atom_indices) {
        ring_names.push_back(atoms[idx].name());
    }

    // For each ring atom, find exocyclic partner
    data.exocyclic_atom_indices.reserve(data.ring_atom_indices.size());

    for (size_t ring_idx : data.ring_atom_indices) {
        const auto& ring_atom = atoms[ring_idx];
        size_t best_exo_idx = ring_idx; // Default to ring atom itself
        double min_dist = validation_constants::BOND_DISTANCE;

        for (size_t i = 0; i < atoms.size(); ++i) {
            const auto& atom = atoms[i];

            // Skip ring atoms
            bool is_ring = false;
            for (const auto& rn : ring_names) {
                if (atom.name() == rn) {
                    is_ring = true;
                    break;
                }
            }
            if (is_ring) continue;

            // Skip hydrogen atoms
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
