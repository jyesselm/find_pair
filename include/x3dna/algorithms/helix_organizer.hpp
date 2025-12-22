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
#include <map>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna::algorithms {

/**
 * @brief Method for ordering base pairs within helices
 */
enum class OrderingMode {
    Legacy,    ///< Use legacy five_to_three algorithm (matches X3DNA)
    ChainBased ///< Use ChainDetector for backbone connectivity ordering
};

/**
 * @brief Direction of backbone linkage between residues
 */
enum class LinkDirection {
    None = 0,    ///< No O3'-P linkage detected
    Forward = 1, ///< i → j linkage (5'→3' direction)
    Reverse = -1 ///< j → i linkage (reverse direction)
};

/**
 * @brief Residue indices for a base pair's two strands
 *
 * Uses 1-based indexing for compatibility with legacy backbone data.
 */
struct StrandResidues {
    size_t strand1; ///< Residue index on strand 1 (1-based)
    size_t strand2; ///< Residue index on strand 2 (1-based)
};

/**
 * @brief Backbone atom coordinates for a residue
 */
struct BackboneAtoms {
    std::optional<geometry::Vector3D> O3_prime; ///< O3' atom coordinates
    std::optional<geometry::Vector3D> P;        ///< P atom coordinates
};

/**
 * @brief Map from residue index (1-based legacy) to backbone atoms
 */
using BackboneData = std::map<size_t, BackboneAtoms>;

/**
 * @brief Direction counts for backbone linkages in a helix
 */
struct DirectionCounts {
    int strand1_forward = 0; ///< i1 → i2 linkages (strand 1 forward)
    int strand1_reverse = 0; ///< i2 → i1 linkages (strand 1 reverse)
    int strand1_none = 0;    ///< No linkage on strand 1
    int strand2_forward = 0; ///< j1 → j2 linkages (strand 2 forward)
    int strand2_reverse = 0; ///< j2 → j1 linkages (strand 2 reverse)
    int strand2_none = 0;    ///< No linkage on strand 2
};

/**
 * @brief Represents a helix segment (continuous run of base pairs)
 */
struct HelixSegment {
    size_t start_idx;                 ///< Start index in the ordered pair list
    size_t end_idx;                   ///< End index (inclusive) in the ordered pair list
    bool is_zdna = false;             ///< Z-DNA conformation detected
    bool has_break = false;           ///< Broken O3'-P linkage within helix
    bool is_parallel = false;         ///< Parallel strand orientation (vs anti-parallel)
    bool has_mixed_direction = false; ///< Mixed strand directions detected
    DirectionCounts direction;        ///< Direction counts for backbone linkages (debug info)
};

/**
 * @brief Pair context for comparison with legacy bp_order
 */
struct PairContextInfo {
    bool is_endpoint = true;         ///< True if pair is at helix end
    std::optional<size_t> neighbor1; ///< First neighbor (closest)
    std::optional<size_t> neighbor2; ///< Second neighbor (opposite z-side)
};

/**
 * @brief Result of helix organization
 */
struct HelixOrdering {
    std::vector<size_t> pair_order;       ///< Indices into original pair list, in helix order
    std::vector<HelixSegment> helices;    ///< Helix segment boundaries
    std::vector<bool> strand_swapped;     ///< Whether strand assignment was swapped for each pair
    std::vector<bool> helix_breaks;       ///< True at positions that are helix boundaries (no backbone link)
    std::vector<PairContextInfo> context; ///< Neighbor context for each pair (for debugging)
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
        double helix_break;         ///< Max distance (Å) between adjacent pairs
        double neighbor_cutoff;     ///< Cutoff for neighbor detection
        double o3p_upper;           ///< Max O3'-P distance for backbone linkage (Å)
        double end_stack_xang;      ///< Max x-angle for stacked WC pairs (degrees)
        OrderingMode ordering_mode; ///< Method for ordering base pairs

        // Legacy uses helix_break=7.8 from $X3DNA/config/misc_3dna.par
        Config()
            : helix_break(7.8), neighbor_cutoff(8.5), o3p_upper(2.5), end_stack_xang(125.0),
              ordering_mode(OrderingMode::Legacy) {}
    };

    explicit HelixOrganizer(const Config& config = Config());

    /**
     * @brief Organize base pairs by helical continuity
     *
     * @param pairs Vector of base pairs (in selection order)
     * @param backbone Optional backbone data for 5'→3' direction checking
     * @param structure Optional structure for chain-based ordering (required if ordering_mode is ChainBased)
     * @return HelixOrdering with reordered indices and helix boundaries
     */
    [[nodiscard]] HelixOrdering organize(const std::vector<core::BasePair>& pairs, const BackboneData& backbone = {},
                                         const core::Structure* structure = nullptr) const;

private:
    Config config_;

    /**
     * @brief Neighbor information for a base pair
     */
    struct PairContext {
        bool is_endpoint = true;         ///< True if at helix end (< 2 neighbors)
        std::optional<size_t> neighbor1; ///< Nearest neighbor
        std::optional<size_t> neighbor2; ///< 2nd neighbor (opposite z-side)
        double dist1 = 0.0;              ///< Distance to neighbor1
        double dist2 = 0.0;              ///< Distance to neighbor2
        bool has_backbone_link1 = false; ///< Backbone connected to neighbor1
        bool has_backbone_link2 = false; ///< Backbone connected to neighbor2
    };

    // Context calculation
    [[nodiscard]] std::vector<PairContext> calculate_context(const std::vector<core::BasePair>& pairs,
                                                             const BackboneData& backbone) const;
    [[nodiscard]] std::vector<size_t> find_endpoints(const std::vector<PairContext>& context) const;
    [[nodiscard]] std::pair<std::vector<size_t>, std::vector<HelixSegment>> locate_helices(
        const std::vector<PairContext>& context, const std::vector<size_t>& endpoints, const BackboneData& backbone,
        size_t num_pairs) const;

    /** @brief Check if two pairs are backbone-connected (either strand) */
    [[nodiscard]] bool are_pairs_backbone_connected(const core::BasePair& pair1, const core::BasePair& pair2,
                                                    const BackboneData& backbone) const;

    // Five-to-three algorithm (legacy equivalent)
    void ensure_five_to_three(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                              std::vector<size_t>& pair_order, std::vector<HelixSegment>& helices,
                              std::vector<bool>& strand_swapped) const;

    // Chain-based ordering (new approach)
    void ensure_chain_order(const std::vector<core::BasePair>& pairs, const core::Structure& structure,
                            std::vector<size_t>& pair_order, std::vector<HelixSegment>& helices,
                            std::vector<bool>& strand_swapped) const;

    /** @brief Get residue indices based on swap status (returns 1-based indices) */
    [[nodiscard]] StrandResidues get_strand_residues(const core::BasePair& pair, bool swapped) const;

    // Five2three sub-functions (legacy equivalents)

    /** @brief Set initial strand assignment for first pair in helix */
    void first_step(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                    std::vector<size_t>& pair_order, const HelixSegment& helix, std::vector<bool>& swapped) const;

    /** @brief Check Watson-Crick base pair z-direction alignment */
    [[nodiscard]] bool wc_bporien(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m, bool swap_n,
                                  const BackboneData& backbone) const;

    /** @brief Check O3' distance patterns for swap indication */
    [[nodiscard]] bool check_o3dist(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                    bool swap_n, const BackboneData& backbone) const;

    /** @brief Check strand chain connectivity for swap indication */
    [[nodiscard]] bool check_schain(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                    bool swap_n, const BackboneData& backbone) const;

    /** @brief Check frame orientation alignment for swap indication */
    [[nodiscard]] bool check_others(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                    bool swap_n, const BackboneData& backbone) const;

    /** @brief Ensure strand 1 direction is consistent */
    [[nodiscard]] bool chain1dir(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m, bool swap_n,
                                 const BackboneData& backbone) const;

    /** @brief Count backbone linkage directions in a helix (may modify swapped and pair_order) */
    DirectionCounts check_direction(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                                    std::vector<size_t>& pair_order, HelixSegment& helix,
                                    std::vector<bool>& swapped) const;

    /** @brief Additional strand corrections based on direction counts */
    void check_strand2(const std::vector<core::BasePair>& pairs, const BackboneData& backbone,
                       const std::vector<size_t>& pair_order, HelixSegment& helix, std::vector<bool>& swapped,
                       const DirectionCounts& direction) const;

    /** @brief Check backbone linkage between residues */
    [[nodiscard]] LinkDirection check_linkage(size_t res_i, size_t res_j, const BackboneData& backbone) const;

    /** @brief Get O3'-O3' distance between residues */
    [[nodiscard]] double o3_distance(size_t i, size_t j, const BackboneData& backbone) const;

    // Geometry helpers
    [[nodiscard]] geometry::Vector3D get_pair_z_axis(const core::BasePair& pair) const;
    [[nodiscard]] geometry::Vector3D get_pair_origin(const core::BasePair& pair) const;

    /** @brief Get frame z-direction based on swap status */
    [[nodiscard]] geometry::Vector3D get_frame_z(const core::BasePair& pair, bool swapped) const;

    /** @brief Calculate angle between x-axes of two pair frames */
    [[nodiscard]] double wcbp_xang(const core::BasePair& pair_m, const core::BasePair& pair_n) const;

    /** @brief Calculate z-direction dot product for WC pairs */
    [[nodiscard]] double wcbp_zdir(const core::BasePair& pair_m, const core::BasePair& pair_n, bool swap_m,
                                   bool swap_n) const;
};

} // namespace x3dna::algorithms
