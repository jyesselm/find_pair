/**
 * @file slot.hpp
 * @brief Hydrogen and lone pair slot definitions for H-bond optimization
 *
 * Slots represent the directional positions where hydrogen atoms (H slots)
 * or lone pairs (LP slots) can participate in hydrogen bonding. Each slot
 * tracks its direction and which bonds are currently using it.
 */

#pragma once

#include <vector>
#include <cmath>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @class HSlot
 * @brief Represents a hydrogen atom slot on a donor atom
 *
 * Tracks the direction of the H atom and supports bifurcation
 * (one H atom donating to two acceptors if angularly separated).
 */
class HSlot {
public:
    /**
     * @brief Construct an H slot with given direction
     * @param direction Unit vector pointing from donor to expected H position
     * @param max_bonds Maximum number of bonds this slot can participate in (default 2)
     */
    explicit HSlot(const geometry::Vector3D& direction, int max_bonds = 2);

    /** @brief Get the slot direction (unit vector) */
    [[nodiscard]] const geometry::Vector3D& direction() const { return direction_; }

    /** @brief Check if slot has not been used at all */
    [[nodiscard]] bool is_available() const { return bond_directions_.empty(); }

    /** @brief Get number of bonds currently using this slot */
    [[nodiscard]] size_t bond_count() const { return bond_directions_.size(); }

    /** @brief Get maximum bonds allowed for this slot */
    [[nodiscard]] int max_bonds() const { return max_bonds_; }

    /**
     * @brief Check if a new bond can be added to this slot
     * @param new_direction Direction of the new bond (donor to acceptor)
     * @param min_angle_deg Minimum angle between existing bonds (for bifurcation)
     * @return True if the bond can be added
     *
     * A bond can be added if:
     * 1. Slot has no bonds yet, OR
     * 2. Slot has bonds but new direction is sufficiently separated (bifurcation)
     *    AND we haven't exceeded max_bonds
     */
    [[nodiscard]] bool can_add_bond(const geometry::Vector3D& new_direction,
                                    double min_angle_deg = 60.0) const;

    /**
     * @brief Add a bond using this slot
     * @param direction Direction of the bond (donor to acceptor)
     */
    void add_bond(const geometry::Vector3D& direction);

    /** @brief Reset slot to unused state */
    void reset();

    /** @brief Get all bond directions currently using this slot */
    [[nodiscard]] const std::vector<geometry::Vector3D>& bond_directions() const {
        return bond_directions_;
    }

private:
    geometry::Vector3D direction_;                    ///< Expected H direction (unit vector)
    std::vector<geometry::Vector3D> bond_directions_; ///< Actual bond directions using this slot
    int max_bonds_;                                   ///< Maximum bonds allowed
};

/**
 * @class LPSlot
 * @brief Represents a lone pair slot on an acceptor atom
 *
 * Tracks the direction of the lone pair and supports bifurcation
 * (one LP accepting from two donors if angularly separated).
 */
class LPSlot {
public:
    /**
     * @brief Construct an LP slot with given direction
     * @param direction Unit vector pointing from acceptor in LP direction
     * @param max_bonds Maximum number of bonds this slot can participate in (default 2)
     */
    explicit LPSlot(const geometry::Vector3D& direction, int max_bonds = 2);

    /** @brief Get the slot direction (unit vector) */
    [[nodiscard]] const geometry::Vector3D& direction() const { return direction_; }

    /** @brief Check if slot has not been used at all */
    [[nodiscard]] bool is_available() const { return bond_directions_.empty(); }

    /** @brief Get number of bonds currently using this slot */
    [[nodiscard]] size_t bond_count() const { return bond_directions_.size(); }

    /** @brief Get maximum bonds allowed for this slot */
    [[nodiscard]] int max_bonds() const { return max_bonds_; }

    /**
     * @brief Check if a new bond can be added to this slot
     * @param new_direction Direction of the new bond (acceptor to donor, i.e., incoming)
     * @param min_angle_deg Minimum angle between existing bonds (for bifurcation)
     * @return True if the bond can be added
     */
    [[nodiscard]] bool can_add_bond(const geometry::Vector3D& new_direction,
                                    double min_angle_deg = 60.0) const;

    /**
     * @brief Add a bond using this slot
     * @param direction Direction of the bond (acceptor to donor)
     */
    void add_bond(const geometry::Vector3D& direction);

    /** @brief Reset slot to unused state */
    void reset();

    /** @brief Get all bond directions currently using this slot */
    [[nodiscard]] const std::vector<geometry::Vector3D>& bond_directions() const {
        return bond_directions_;
    }

private:
    geometry::Vector3D direction_;                    ///< Expected LP direction (unit vector)
    std::vector<geometry::Vector3D> bond_directions_; ///< Actual bond directions using this slot
    int max_bonds_;                                   ///< Maximum bonds allowed
};

/**
 * @brief Calculate angle between two vectors in degrees
 * @param v1 First vector
 * @param v2 Second vector
 * @return Angle in degrees [0, 180]
 */
[[nodiscard]] inline double angle_between_degrees(const geometry::Vector3D& v1,
                                                   const geometry::Vector3D& v2) {
    double dot = v1.dot(v2);
    // Clamp to [-1, 1] to handle numerical errors
    dot = std::max(-1.0, std::min(1.0, dot));
    return std::acos(dot) * 180.0 / M_PI;
}

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
