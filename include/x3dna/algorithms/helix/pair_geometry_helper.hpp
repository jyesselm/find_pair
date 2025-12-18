/**
 * @file pair_geometry_helper.hpp
 * @brief Geometry calculations for base pairs in helix organization
 *
 * Extracts pair geometry helper functions from helix_organizer for
 * calculating origins, z-axes, and strand residue mappings.
 */

#pragma once

#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/helix_organizer.hpp> // For StrandResidues
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna::algorithms::helix {

// Re-use types from helix_organizer
using algorithms::StrandResidues;

/**
 * @class PairGeometryHelper
 * @brief Calculates geometric properties of base pairs
 *
 * Provides utility functions for extracting geometric information
 * from base pairs during helix organization.
 */
class PairGeometryHelper {
public:
    /**
     * @brief Get the average origin of a base pair
     * @param pair Base pair with frames
     * @return Average of the two frame origins
     * @pre Both frames must have values
     */
    [[nodiscard]] static geometry::Vector3D get_pair_origin(const core::BasePair& pair);

    /**
     * @brief Get the average z-axis of a base pair
     *
     * Computes the bisector of the two z-axes, accounting for
     * whether they point in the same or opposite directions.
     *
     * @param pair Base pair with frames
     * @return Normalized average z-axis
     * @pre Both frames must have values
     */
    [[nodiscard]] static geometry::Vector3D get_pair_z_axis(const core::BasePair& pair);

    /**
     * @brief Get frame z-axis based on swap status
     * @param pair Base pair with frames
     * @param swapped If true, use frame2's z-axis; otherwise frame1's
     * @return Z-axis of selected frame
     * @pre Both frames must have values
     */
    [[nodiscard]] static geometry::Vector3D get_frame_z(const core::BasePair& pair, bool swapped);

    /**
     * @brief Get 1-based residue indices for strands
     *
     * Returns residue indices accounting for:
     * - The normalized (smaller, larger) storage order in BasePair
     * - The finding_order_swapped flag (original finding order)
     * - The five2three swap flag (strand assignment)
     *
     * @param pair Base pair
     * @param swapped Whether five2three algorithm swapped strands
     * @return StrandResidues with 1-based indices
     */
    [[nodiscard]] static StrandResidues get_strand_residues(const core::BasePair& pair, bool swapped);
};

} // namespace x3dna::algorithms::helix
