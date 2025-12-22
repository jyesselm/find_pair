/**
 * @file structure_legacy_order_impl.cpp
 * @brief Implementation of legacy order methods for Structure class
 */

#include <x3dna/core/structure.hpp>
#include <x3dna/core/structure_legacy_order.hpp>

namespace x3dna {
namespace core {

std::vector<const Residue*> Structure::residues_in_legacy_order() const {
    // Use the standalone function from structure_legacy_order namespace
    // Need to fully qualify to avoid ambiguity with method name
    return ::x3dna::core::get_residues_in_legacy_order(*this);
}

const Residue* Structure::get_residue_by_legacy_idx(int legacy_idx) const {
    // Use the standalone function from structure_legacy_order namespace
    // Need to fully qualify to avoid ambiguity with method name
    return ::x3dna::core::get_residue_by_legacy_idx(*this, legacy_idx);
}

int Structure::get_legacy_idx_for_residue(const Residue* residue) const {
    // Use the standalone function from structure_legacy_order namespace
    // Need to fully qualify to avoid ambiguity with method name
    return ::x3dna::core::get_legacy_idx_for_residue(*this, residue);
}

} // namespace core
} // namespace x3dna
