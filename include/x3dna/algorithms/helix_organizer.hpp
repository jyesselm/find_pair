/**
 * @file helix_organizer.hpp
 * @brief Organizes base pairs by helical continuity
 * 
 * This replicates the legacy X3DNA bp_context/locate_helix/five2three algorithm
 * which orders base pairs so that consecutive pairs are spatially adjacent
 * within the same helix.
 */

#pragma once

#include <vector>
#include <optional>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna::algorithms {

/**
 * @brief Represents a helix segment (continuous run of base pairs)
 */
struct HelixSegment {
    size_t start_idx;       ///< Start index in the ordered pair list
    size_t end_idx;         ///< End index (inclusive) in the ordered pair list
    bool is_zdna = false;   ///< Z-DNA conformation detected
    bool has_break = false; ///< Broken O3'-P linkage within helix
};

/**
 * @brief Result of helix organization
 */
struct HelixOrdering {
    std::vector<size_t> pair_order;     ///< Indices into original pair list, in helix order
    std::vector<HelixSegment> helices;  ///< Helix segment boundaries
    std::vector<bool> strand_swapped;   ///< Whether strand assignment was swapped for each pair
};

/**
 * @brief Organizes base pairs by helical continuity
 * 
 * This class implements the legacy X3DNA algorithm for reordering base pairs
 * so that consecutive pairs in the output are spatially adjacent within
 * the same helix. This is essential for meaningful step parameter calculations.
 */
class HelixOrganizer {
public:
    /**
     * @brief Configuration parameters
     */
    struct Config {
        double helix_break;      ///< Max distance (Å) between adjacent pairs
        double neighbor_cutoff;  ///< Cutoff for neighbor detection
        
        Config() : helix_break(7.5), neighbor_cutoff(8.0) {}
    };
    
    explicit HelixOrganizer(const Config& config = Config());
    
    /**
     * @brief Organize base pairs by helical continuity
     * 
     * @param pairs Vector of base pairs (in selection order)
     * @return HelixOrdering with reordered indices and helix boundaries
     */
    HelixOrdering organize(const std::vector<core::BasePair>& pairs) const;
    
private:
    Config config_;
    
    /**
     * @brief Neighbor information for a base pair
     */
    struct PairContext {
        bool is_endpoint = true;        ///< True if at helix end (< 2 neighbors)
        std::optional<size_t> neighbor1; ///< Nearest neighbor
        std::optional<size_t> neighbor2; ///< 2nd neighbor (opposite z-side)
        double dist1 = 0.0;             ///< Distance to neighbor1
        double dist2 = 0.0;             ///< Distance to neighbor2
    };
    
    /**
     * @brief Calculate pair context (neighbors) for all pairs
     * Equivalent to legacy bp_context()
     */
    std::vector<PairContext> calculate_context(
        const std::vector<core::BasePair>& pairs) const;
    
    /**
     * @brief Find helix endpoints from context
     */
    std::vector<size_t> find_endpoints(
        const std::vector<PairContext>& context) const;
    
    /**
     * @brief Chain pairs into helices starting from endpoints
     * Equivalent to legacy locate_helix()
     */
    std::pair<std::vector<size_t>, std::vector<HelixSegment>> locate_helices(
        const std::vector<PairContext>& context,
        const std::vector<size_t>& endpoints,
        size_t num_pairs) const;
    
    /**
     * @brief Ensure 5'→3' direction within each helix
     * Equivalent to legacy five2three()
     */
    void ensure_five_to_three(
        const std::vector<core::BasePair>& pairs,
        std::vector<size_t>& pair_order,
        std::vector<HelixSegment>& helices,
        std::vector<bool>& strand_swapped) const;
    
    /**
     * @brief Get combined z-axis for a base pair (average of both base z-axes)
     */
    geometry::Vector3D get_pair_z_axis(const core::BasePair& pair) const;
    
    /**
     * @brief Get pair origin (midpoint of both base origins)
     */
    geometry::Vector3D get_pair_origin(const core::BasePair& pair) const;
};

} // namespace x3dna::algorithms

