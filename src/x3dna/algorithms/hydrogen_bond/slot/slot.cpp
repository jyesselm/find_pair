/**
 * @file slot.cpp
 * @brief Implementation of HSlot and LPSlot classes
 */

#include <x3dna/algorithms/hydrogen_bond/slot/slot.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

// ============================================================================
// HSlot Implementation
// ============================================================================

HSlot::HSlot(const geometry::Vector3D& direction, int max_bonds)
    : direction_(direction.normalized())
    , max_bonds_(max_bonds) {
}

bool HSlot::can_add_bond(const geometry::Vector3D& new_direction, double min_angle_deg) const {
    // If no bonds yet, always available
    if (bond_directions_.empty()) {
        return true;
    }

    // Check if we've reached max bonds
    if (static_cast<int>(bond_directions_.size()) >= max_bonds_) {
        return false;
    }

    // Check angular separation from all existing bonds (bifurcation check)
    for (const auto& existing : bond_directions_) {
        double angle = angle_between_degrees(new_direction.normalized(), existing);
        if (angle < min_angle_deg) {
            return false; // Too close to existing bond
        }
    }

    return true;
}

void HSlot::add_bond(const geometry::Vector3D& direction) {
    bond_directions_.push_back(direction.normalized());
}

void HSlot::reset() {
    bond_directions_.clear();
}

// ============================================================================
// LPSlot Implementation
// ============================================================================

LPSlot::LPSlot(const geometry::Vector3D& direction, int max_bonds)
    : direction_(direction.normalized())
    , max_bonds_(max_bonds) {
}

bool LPSlot::can_add_bond(const geometry::Vector3D& new_direction, double min_angle_deg) const {
    // If no bonds yet, always available
    if (bond_directions_.empty()) {
        return true;
    }

    // Check if we've reached max bonds
    if (static_cast<int>(bond_directions_.size()) >= max_bonds_) {
        return false;
    }

    // Check angular separation from all existing bonds (bifurcation check)
    for (const auto& existing : bond_directions_) {
        double angle = angle_between_degrees(new_direction.normalized(), existing);
        if (angle < min_angle_deg) {
            return false; // Too close to existing bond
        }
    }

    return true;
}

void LPSlot::add_bond(const geometry::Vector3D& direction) {
    bond_directions_.push_back(direction.normalized());
}

void LPSlot::reset() {
    bond_directions_.clear();
}

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
