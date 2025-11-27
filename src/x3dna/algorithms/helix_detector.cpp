/**
 * @file helix_detector.cpp
 * @brief HelixDetector implementation
 */

#include <x3dna/algorithms/helix_detector.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <algorithm>
#include <cmath>

namespace x3dna {
namespace algorithms {

HelixDetector::HelixDetector(double helix_break_distance)
    : helix_break_distance_(helix_break_distance) {
}

std::vector<Helix> HelixDetector::detect_helices(const std::vector<core::BasePair>& pairs) {
    std::vector<Helix> helices;

    if (pairs.empty()) {
        return helices;
    }

    if (pairs.size() == 1) {
        // Single pair forms a helix
        Helix helix(0, 0, false);
        helices.push_back(helix);
        return helices;
    }

    // Simple algorithm: group consecutive pairs that are neighbors
    // This is a simplified version - legacy code is more complex
    std::vector<bool> processed(pairs.size(), false);

    for (size_t i = 0; i < pairs.size(); ++i) {
        if (processed[i]) {
            continue;
        }

        // Start a new helix
        size_t start = i;
        size_t end = i;

        // Try to extend helix forward
        for (size_t j = i + 1; j < pairs.size(); ++j) {
            if (are_neighbors(pairs[end], pairs[j])) {
                end = j;
                processed[j] = true;
            } else {
                break; // Break in helix
            }
        }

        // Create helix
        Helix helix;
        helix.start_index = start;
        helix.end_index = end;
        for (size_t k = start; k <= end; ++k) {
            helix.base_pair_indices.push_back(k);
        }

        // Check if circular (first and last pairs are neighbors)
        helix.is_circular = is_circular(pairs, helix);

        helices.push_back(helix);
        processed[i] = true;
    }

    return helices;
}

std::vector<Helix> HelixDetector::detect_helices(const core::Structure& /* structure */) {
    // TODO: Extract base pairs from structure and call detect_helices(vector<BasePair>)
    // For now, return empty (this would require Structure to have base pairs)
    return std::vector<Helix>();
}

void HelixDetector::reorder_base_pairs(std::vector<core::BasePair>& pairs) {
    if (pairs.size() < 2) {
        return;
    }

    // Detect helices first
    auto helices = detect_helices(pairs);

    // Reorder each helix to 5'→3'
    for (const auto& helix : helices) {
        ensure_five_to_three_ordering(pairs, helix);
    }
}

void HelixDetector::ensure_five_to_three_ordering(std::vector<core::BasePair>& pairs, const Helix& helix) {
    // For now, this is a placeholder
    // Legacy code has complex logic for determining 5'→3' orientation
    // based on O3'→P linkage between consecutive residues
    // TODO: Implement full 5'→3' ordering logic

    // Basic check: ensure pairs are ordered by residue index
    // This is a simplified version
    (void)pairs;
    (void)helix;
}

double HelixDetector::calculate_pair_distance(const core::BasePair& pair1,
                                               const core::BasePair& pair2) const {
    // Calculate distance between base pair origins
    // Use frame origins if available, otherwise use residue centers

    auto frame1_opt = pair1.frame1();
    auto frame2_opt = pair2.frame1();

    if (frame1_opt.has_value() && frame2_opt.has_value()) {
        auto origin1 = frame1_opt->origin();
        auto origin2 = frame2_opt->origin();
        auto diff = origin2 - origin1;
        return diff.length();
    }

    // Fallback: return large distance if frames not available
    return helix_break_distance_ * 2.0;
}

bool HelixDetector::are_neighbors(const core::BasePair& pair1, const core::BasePair& pair2) const {
    double distance = calculate_pair_distance(pair1, pair2);
    return distance <= helix_break_distance_;
}

bool HelixDetector::is_circular(const std::vector<core::BasePair>& pairs, const Helix& helix) const {
    if (helix.base_pair_indices.size() < 2) {
        return false;
    }

    size_t first_idx = helix.base_pair_indices.front();
    size_t last_idx = helix.base_pair_indices.back();

    if (first_idx >= pairs.size() || last_idx >= pairs.size()) {
        return false;
    }

    return are_neighbors(pairs[first_idx], pairs[last_idx]);
}

std::vector<size_t> HelixDetector::find_neighbors(const std::vector<core::BasePair>& pairs,
                                                   size_t pair_index) const {
    std::vector<size_t> neighbors;

    if (pair_index >= pairs.size()) {
        return neighbors;
    }

    for (size_t i = 0; i < pairs.size(); ++i) {
        if (i == pair_index) {
            continue;
        }

        if (are_neighbors(pairs[pair_index], pairs[i])) {
            neighbors.push_back(i);
        }
    }

    return neighbors;
}

} // namespace algorithms
} // namespace x3dna

