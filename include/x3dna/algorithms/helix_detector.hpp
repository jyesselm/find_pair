/**
 * @file helix_detector.hpp
 * @brief HelixDetector class for detecting helices from base pairs
 */

#pragma once

#include <vector>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/reference_frame.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @struct Helix
 * @brief Represents a helix region containing base pairs
 */
struct Helix {
    /**
     * @brief Base pair indices in this helix (0-based, ordered 5' to 3')
     */
    std::vector<size_t> base_pair_indices;

    /**
     * @brief Start index of helix (0-based, inclusive)
     */
    size_t start_index = 0;

    /**
     * @brief End index of helix (0-based, inclusive)
     */
    size_t end_index = 0;

    /**
     * @brief Whether this helix is circular (closed structure)
     */
    bool is_circular = false;

    /**
     * @brief Default constructor
     */
    Helix() = default;

    /**
     * @brief Constructor with indices
     */
    Helix(size_t start, size_t end, bool circular = false) : start_index(start), end_index(end), is_circular(circular) {
        for (size_t i = start; i <= end; ++i) {
            base_pair_indices.push_back(i);
        }
    }
};

/**
 * @class HelixDetector
 * @brief Detects helices from base pairs and reorders them
 *
 * This class implements helix detection by:
 * 1. Analyzing base pair context (neighbors, distances)
 * 2. Grouping consecutive base pairs into helices
 * 3. Reordering pairs to 5'→3' orientation
 * 4. Handling circular structures
 */
class HelixDetector {
public:
    /**
     * @brief Constructor
     * @param helix_break_distance Distance threshold for helix breaks (default: 7.5 Angstroms)
     */
    explicit HelixDetector(double helix_break_distance = 7.5);

    /**
     * @brief Detect helices from base pairs
     * @param pairs Vector of base pairs
     * @return Vector of detected helices
     */
    [[nodiscard]] std::vector<Helix> detect_helices(const std::vector<core::BasePair>& pairs);

    /**
     * @brief Detect helices from structure (uses structure's base pairs)
     * @param structure Structure with base pairs
     * @return Vector of detected helices
     */
    [[nodiscard]] std::vector<Helix> detect_helices(const core::Structure& structure);

    /**
     * @brief Reorder base pairs to ensure 5'→3' orientation
     * @param pairs Vector of base pairs (modified in place)
     */
    void reorder_base_pairs(std::vector<core::BasePair>& pairs);

    /**
     * @brief Ensure 5'→3' ordering for pairs in a helix
     * @param pairs Vector of base pairs (modified in place)
     * @param helix Helix structure defining which pairs to reorder
     */
    void ensure_five_to_three_ordering(std::vector<core::BasePair>& pairs, const Helix& helix);

    /**
     * @brief Set helix break distance threshold
     */
    void set_helix_break_distance(double distance) {
        helix_break_distance_ = distance;
    }

    /**
     * @brief Get helix break distance threshold
     */
    [[nodiscard]] double helix_break_distance() const {
        return helix_break_distance_;
    }

private:
    /**
     * @brief Calculate distance between two base pairs
     * @param pair1 First base pair
     * @param pair2 Second base pair
     * @return Distance between pair origins (Angstroms)
     */
    [[nodiscard]] double calculate_pair_distance(const core::BasePair& pair1, const core::BasePair& pair2) const;

    /**
     * @brief Check if two base pairs are neighbors (within helix_break distance)
     * @param pair1 First base pair
     * @param pair2 Second base pair
     * @return True if pairs are neighbors
     */
    [[nodiscard]] bool are_neighbors(const core::BasePair& pair1, const core::BasePair& pair2) const;

    /**
     * @brief Check if helix is circular (first and last pairs are neighbors)
     * @param pairs Vector of base pairs
     * @param helix Helix to check
     * @return True if helix is circular
     */
    [[nodiscard]] bool is_circular(const std::vector<core::BasePair>& pairs, const Helix& helix) const;

    /**
     * @brief Analyze base pair context to find neighbors
     * @param pairs All base pairs
     * @param pair_index Index of pair to analyze
     * @return Vector of neighbor indices (within helix_break distance)
     */
    [[nodiscard]] std::vector<size_t> find_neighbors(const std::vector<core::BasePair>& pairs, size_t pair_index) const;

    double helix_break_distance_; // Distance threshold for helix breaks (Angstroms)
};

} // namespace algorithms
} // namespace x3dna
