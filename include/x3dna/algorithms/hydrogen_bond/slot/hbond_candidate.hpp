/**
 * @file hbond_candidate.hpp
 * @brief Candidate hydrogen bond for slot-based optimization
 */

#pragma once

#include <string>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @struct HBondCandidate
 * @brief Represents a potential hydrogen bond with alignment scoring
 */
struct HBondCandidate {
    // Residue identifiers
    std::string donor_res_id;
    std::string acceptor_res_id;

    // Atom names
    std::string donor_atom;
    std::string acceptor_atom;

    // Positions
    geometry::Vector3D donor_pos;
    geometry::Vector3D acceptor_pos;

    // Distance
    double distance = 0.0;

    // Slot assignment (set during selection)
    int h_slot_idx = -1;
    int lp_slot_idx = -1;

    // Alignment score (0-2, higher is better)
    double alignment_score = 0.0;

    /**
     * @brief Calculate quality score for ranking
     * @return Quality score (higher is better)
     *
     * Score combines distance (shorter is better) with alignment
     * (better alignment can compensate for slightly longer distance).
     */
    [[nodiscard]] double quality_score() const {
        // Weight 0.4 allows good alignment to overcome ~0.25Å distance disadvantage
        return -distance + 0.4 * alignment_score;
    }

    /**
     * @brief Get direction from donor to acceptor
     * @return Unit vector
     */
    [[nodiscard]] geometry::Vector3D direction() const {
        return (acceptor_pos - donor_pos).normalized();
    }
};

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
