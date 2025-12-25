/**
 * @file slot_predictor.hpp
 * @brief Geometry-based prediction of H and LP slot positions
 *
 * Predicts the directions of hydrogen atoms and lone pairs based on
 * the local bonding geometry of donor/acceptor atoms.
 */

#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <x3dna/algorithms/hydrogen_bond/slot/slot.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/core/residue.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @class SlotPredictor
 * @brief Predicts H and LP slot positions from molecular geometry
 *
 * Uses connectivity information and base plane normal to predict
 * the directions of hydrogen atoms (for donors) and lone pairs
 * (for acceptors).
 */
class SlotPredictor {
public:
    /**
     * @brief Compute the normal vector to the base plane
     * @param residue Residue to compute normal for
     * @return Unit normal vector (or [0,0,1] if insufficient atoms)
     */
    [[nodiscard]] static geometry::Vector3D compute_base_normal(const core::Residue& residue);

    /**
     * @brief Predict hydrogen slot positions for a donor atom
     * @param base_type Single letter base type (A, G, C, U, T, P, I)
     * @param atom_name Donor atom name (e.g., "N6", "N1")
     * @param residue Residue containing the atom
     * @param base_normal Normal vector to the base plane
     * @return Vector of HSlot objects with predicted directions
     */
    [[nodiscard]] static std::vector<HSlot> predict_h_slots(
        char base_type,
        const std::string& atom_name,
        const core::Residue& residue,
        const geometry::Vector3D& base_normal);

    /**
     * @brief Predict lone pair slot positions for an acceptor atom
     * @param base_type Single letter base type (A, G, C, U, T, P, I)
     * @param atom_name Acceptor atom name (e.g., "O6", "N1")
     * @param residue Residue containing the atom
     * @param base_normal Normal vector to the base plane
     * @return Vector of LPSlot objects with predicted directions
     */
    [[nodiscard]] static std::vector<LPSlot> predict_lp_slots(
        char base_type,
        const std::string& atom_name,
        const core::Residue& residue,
        const geometry::Vector3D& base_normal);

private:
    /**
     * @brief Rotate a vector around an axis
     * @param v Vector to rotate
     * @param axis Rotation axis (must be unit vector)
     * @param angle_deg Rotation angle in degrees
     * @return Rotated vector
     */
    [[nodiscard]] static geometry::Vector3D rotate_around_axis(
        const geometry::Vector3D& v,
        const geometry::Vector3D& axis,
        double angle_deg);

    /**
     * @brief Get connected atom names for a base atom
     * @param base_type Base type (A, G, C, U, T)
     * @param atom_name Atom name
     * @return List of connected atom names
     */
    [[nodiscard]] static std::vector<std::string> get_connectivity(
        char base_type,
        const std::string& atom_name);

    /**
     * @brief Predict slots for sp2 NH2 amino group (2 hydrogens)
     */
    [[nodiscard]] static std::vector<HSlot> predict_sp2_amino_slots(
        const geometry::Vector3D& donor_pos,
        const geometry::Vector3D& bonded_pos,
        const geometry::Vector3D& base_normal);

    /**
     * @brief Predict slots for sp2 imino NH group (1 hydrogen)
     */
    [[nodiscard]] static std::vector<HSlot> predict_sp2_imino_slots(
        const geometry::Vector3D& donor_pos,
        const std::vector<geometry::Vector3D>& bonded_positions,
        const geometry::Vector3D& base_normal);

    /**
     * @brief Predict slots for sp2 carbonyl oxygen (2 lone pairs)
     */
    [[nodiscard]] static std::vector<LPSlot> predict_sp2_carbonyl_slots(
        const geometry::Vector3D& acceptor_pos,
        const geometry::Vector3D& bonded_pos,
        const geometry::Vector3D& base_normal);

    /**
     * @brief Predict slots for sp2 ring nitrogen (1 lone pair)
     */
    [[nodiscard]] static std::vector<LPSlot> predict_sp2_ring_nitrogen_slots(
        const geometry::Vector3D& acceptor_pos,
        const std::vector<geometry::Vector3D>& bonded_positions);

    /**
     * @brief Predict slots for sp3 hydroxyl (ribose O2', O3', etc.)
     */
    [[nodiscard]] static std::vector<HSlot> predict_sp3_hydroxyl_h_slots(
        const geometry::Vector3D& oxygen_pos,
        const geometry::Vector3D& bonded_carbon_pos);

    [[nodiscard]] static std::vector<LPSlot> predict_sp3_hydroxyl_lp_slots(
        const geometry::Vector3D& oxygen_pos,
        const geometry::Vector3D& bonded_carbon_pos);
};

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
