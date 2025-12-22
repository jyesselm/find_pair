/**
 * @file structure_legacy_order_impl.cpp
 * @brief Implementation of legacy order methods for Structure class
 *
 * Uses stored legacy_residue_idx values instead of recomputing from atoms.
 */

#include <x3dna/core/structure.hpp>
#include <x3dna/core/structure_legacy_order.hpp>
#include <algorithm>

namespace x3dna {
namespace core {

std::vector<const Residue*> Structure::residues_in_legacy_order() const {
    return ::x3dna::core::get_residues_in_legacy_order(*this);
}

const Residue* Structure::get_residue_by_legacy_idx(int legacy_idx) const {
    if (legacy_idx < 1) {
        return nullptr;
    }

    // Search for residue with matching legacy_residue_idx
    for (const auto& chain : chains_) {
        for (const auto& residue : chain.residues()) {
            if (residue.legacy_residue_idx() == legacy_idx) {
                return &residue;
            }
        }
    }

    return nullptr;
}

int Structure::get_legacy_idx_for_residue(const Residue* residue) const {
    if (!residue) {
        return 0;
    }

    // Return the stored legacy_residue_idx
    return residue->legacy_residue_idx();
}

} // namespace core
} // namespace x3dna
