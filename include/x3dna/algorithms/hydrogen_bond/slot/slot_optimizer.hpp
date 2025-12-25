/**
 * @file slot_optimizer.hpp
 * @brief Slot-based hydrogen bond optimizer
 *
 * Implements greedy H-bond selection with slot saturation tracking,
 * alignment scoring, and bifurcation support.
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <x3dna/algorithms/hydrogen_bond/slot/slot.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/slot_cache.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/hbond_candidate.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/slot_optimizer_params.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/algorithms/hydrogen_bond/hbond.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

/**
 * @class SlotOptimizer
 * @brief Greedy H-bond optimizer with slot saturation tracking
 *
 * Selects optimal H-bonds by:
 * 1. Finding all candidate H-bonds within distance threshold
 * 2. Predicting H and LP slot positions from geometry
 * 3. Scoring alignment of each candidate with available slots
 * 4. Greedily selecting best candidates while respecting slot capacity
 *
 * In baseline mode, uses simpler distance-based selection matching
 * legacy C++ behavior.
 */
class SlotOptimizer {
public:
    /**
     * @brief Construct optimizer with given parameters
     * @param params Optimization parameters
     */
    explicit SlotOptimizer(const SlotOptimizerParams& params = SlotOptimizerParams::optimized());

    /**
     * @brief Optimize H-bonds between two residues
     * @param res1 First residue
     * @param res2 Second residue
     * @return Vector of selected H-bonds
     */
    [[nodiscard]] std::vector<core::HBond> optimize_pair(
        const core::Residue& res1,
        const core::Residue& res2);

    /**
     * @brief Get current parameters
     */
    [[nodiscard]] const SlotOptimizerParams& params() const { return params_; }

    /**
     * @brief Set new parameters
     */
    void set_params(const SlotOptimizerParams& params) { params_ = params; }

private:
    SlotOptimizerParams params_;

    /**
     * @brief Find all candidate H-bonds between two residues
     */
    [[nodiscard]] std::vector<HBondCandidate> find_candidates(
        const core::Residue& res1,
        const core::Residue& res2) const;

    /**
     * @brief Select optimal H-bonds using slot-based greedy algorithm
     */
    [[nodiscard]] std::vector<core::HBond> select_optimal(
        std::vector<HBondCandidate>& candidates,
        SlotCache& cache1,
        SlotCache& cache2);

    /**
     * @brief Select H-bonds using baseline (legacy-compatible) algorithm
     */
    [[nodiscard]] std::vector<core::HBond> select_baseline(
        std::vector<HBondCandidate>& candidates,
        const core::Residue& res1,
        const core::Residue& res2) const;

    /**
     * @brief Score alignment for a candidate and find best slots
     * @param candidate Candidate to score (h_slot_idx, lp_slot_idx, alignment_score updated)
     * @param h_slots Available H slots for the donor
     * @param lp_slots Available LP slots for the acceptor
     */
    void score_alignment(
        HBondCandidate& candidate,
        const std::vector<HSlot>& h_slots,
        const std::vector<LPSlot>& lp_slots) const;

    /**
     * @brief Try to find alternative available slots for a candidate
     * @return True if alternative slots found
     */
    [[nodiscard]] bool try_alternative_slots(
        HBondCandidate& candidate,
        std::vector<HSlot>& h_slots,
        std::vector<LPSlot>& lp_slots) const;

    /**
     * @brief Convert candidate to HBond
     */
    [[nodiscard]] core::HBond candidate_to_hbond(const HBondCandidate& candidate) const;

    /**
     * @brief Get base type character for a residue
     */
    [[nodiscard]] static char get_base_type(const core::Residue& residue);

    /**
     * @brief Check if both atoms are backbone atoms
     */
    [[nodiscard]] static bool is_backbone_backbone(
        const std::string& atom1,
        const std::string& atom2);
};

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
