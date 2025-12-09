/**
 * @file helix_organizer.cpp
 * @brief Implementation of helix organization algorithm
 * 
 * This implements the legacy X3DNA five2three algorithm for ensuring
 * proper 5'→3' strand direction in base pair step calculations.
 */

#include <x3dna/algorithms/helix_organizer.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace x3dna::algorithms {

namespace {
    constexpr double PI = 3.14159265358979323846;
    
    // Convert dot product to angle in degrees
    double dot2ang(double d) {
        if (d > 1.0) d = 1.0;
        if (d < -1.0) d = -1.0;
        return std::acos(d) * 180.0 / PI;
    }
}

HelixOrganizer::HelixOrganizer(const Config& config) : config_(config) {}

// =============================================================================
// Geometry helpers
// =============================================================================

geometry::Vector3D HelixOrganizer::get_pair_origin(const core::BasePair& pair) const {
    if (!pair.frame1().has_value()) {
        return geometry::Vector3D(0, 0, 0);
    }
    return pair.frame1().value().origin();
}

geometry::Vector3D HelixOrganizer::get_pair_z_axis(const core::BasePair& pair) const {
    if (!pair.frame1().has_value()) {
        return geometry::Vector3D(0, 0, 1);
    }
    return pair.frame1().value().z_axis();
}

geometry::Vector3D HelixOrganizer::get_frame_z(const core::BasePair& pair, bool swapped) const {
    if (swapped) {
        if (pair.frame2().has_value()) {
            return pair.frame2().value().z_axis();
        }
    } else {
        if (pair.frame1().has_value()) {
            return pair.frame1().value().z_axis();
        }
    }
    return geometry::Vector3D(0, 0, 1);
}

void HelixOrganizer::get_ij(const core::BasePair& pair, bool swapped,
                            size_t& i, size_t& j) const {
    if (swapped) {
        i = pair.residue_idx2();
        j = pair.residue_idx1();
    } else {
        i = pair.residue_idx1();
        j = pair.residue_idx2();
    }
}

// =============================================================================
// Backbone connectivity
// =============================================================================

int HelixOrganizer::is_linked(size_t i, size_t j, const BackboneData& backbone) const {
    auto it_i = backbone.find(i);
    auto it_j = backbone.find(j);
    
    if (it_i == backbone.end() || it_j == backbone.end()) {
        return 0;
    }
    
    const auto& atoms_i = it_i->second;
    const auto& atoms_j = it_j->second;
    
    // Check O3'[i] → P[j] direction
    if (atoms_i.O3_prime.has_value() && atoms_j.P.has_value()) {
        double dist = (atoms_i.O3_prime.value() - atoms_j.P.value()).length();
        if (dist <= config_.o3p_upper) {
            return 1;  // i → j linkage (5'→3')
        }
    }
    
    // Check O3'[j] → P[i] direction (reverse)
    if (atoms_j.O3_prime.has_value() && atoms_i.P.has_value()) {
        double dist = (atoms_j.O3_prime.value() - atoms_i.P.value()).length();
        if (dist <= config_.o3p_upper) {
            return -1;  // j → i linkage (reverse)
        }
    }
    
    return 0;  // No linkage
}

double HelixOrganizer::o3_distance(size_t i, size_t j, const BackboneData& backbone) const {
    auto it_i = backbone.find(i);
    auto it_j = backbone.find(j);
    
    if (it_i == backbone.end() || it_j == backbone.end()) {
        return -1.0;
    }
    
    const auto& atoms_i = it_i->second;
    const auto& atoms_j = it_j->second;
    
    if (atoms_i.O3_prime.has_value() && atoms_j.O3_prime.has_value()) {
        return (atoms_i.O3_prime.value() - atoms_j.O3_prime.value()).length();
    }
    
    return -1.0;
}

// =============================================================================
// WC pair geometry checks
// =============================================================================

double HelixOrganizer::wcbp_xang(const core::BasePair& pair_m, const core::BasePair& pair_n) const {
    // Calculate angle between combined x-axes of the two pairs
    if (!pair_m.frame1().has_value() || !pair_m.frame2().has_value() ||
        !pair_n.frame1().has_value() || !pair_n.frame2().has_value()) {
        return 180.0;  // Return large angle if frames missing
    }
    
    // Sum of x-axes for each pair
    auto xm = pair_m.frame1().value().x_axis() + pair_m.frame2().value().x_axis();
    auto xn = pair_n.frame1().value().x_axis() + pair_n.frame2().value().x_axis();
    
    xm = xm.normalized();
    xn = xn.normalized();
    
    double dot = xm.dot(xn);
    return dot2ang(dot);
}

double HelixOrganizer::wcbp_zdir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                  bool swap_m, bool swap_n) const {
    // Get z-direction vectors based on swap status
    geometry::Vector3D zm, zn;
    
    if (swap_m) {
        if (pair_m.frame2().has_value() && pair_m.frame1().has_value()) {
            zm = pair_m.frame1().value().z_axis() - pair_m.frame2().value().z_axis();
        } else {
            return 0.0;
        }
    } else {
        if (pair_m.frame1().has_value() && pair_m.frame2().has_value()) {
            zm = pair_m.frame2().value().z_axis() - pair_m.frame1().value().z_axis();
        } else {
            return 0.0;
        }
    }
    
    if (swap_n) {
        if (pair_n.frame2().has_value() && pair_n.frame1().has_value()) {
            zn = pair_n.frame1().value().z_axis() - pair_n.frame2().value().z_axis();
        } else {
            return 0.0;
        }
    } else {
        if (pair_n.frame1().has_value() && pair_n.frame2().has_value()) {
            zn = pair_n.frame2().value().z_axis() - pair_n.frame1().value().z_axis();
        } else {
            return 0.0;
        }
    }
    
    zm = zm.normalized();
    zn = zn.normalized();
    
    return zm.dot(zn);
}

// =============================================================================
// Five2three sub-functions
// =============================================================================

void HelixOrganizer::first_step(const std::vector<core::BasePair>& pairs,
                                const BackboneData& backbone,
                                std::vector<size_t>& pair_order,
                                const HelixSegment& helix,
                                std::vector<bool>& swapped) const {
    // For single-pair helices, nothing to do
    if (helix.end_idx <= helix.start_idx) {
        return;
    }
    
    // Count backbone direction across the WHOLE helix
    int strand1_forward = 0;
    int strand1_reverse = 0;
    int strand2_forward = 0;
    int strand2_reverse = 0;
    
    for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
        size_t idx_m = pair_order[pos];
        size_t idx_n = pair_order[pos + 1];
        
        const auto& pair_m = pairs[idx_m];
        const auto& pair_n = pairs[idx_n];
        
        size_t i1, j1, i2, j2;
        get_ij(pair_m, swapped[idx_m], i1, j1);
        get_ij(pair_n, swapped[idx_n], i2, j2);
        
        int k1 = is_linked(i1, i2, backbone);
        int k2 = is_linked(j1, j2, backbone);
        
        if (k1 == 1) strand1_forward++;
        if (k1 == -1) strand1_reverse++;
        if (k2 == 1) strand2_forward++;
        if (k2 == -1) strand2_reverse++;
    }
    
    // Check if backbone direction is consistent
    bool backbone_consistent = (strand1_forward == 0 || strand1_reverse == 0) &&
                               (strand2_forward == 0 || strand2_reverse == 0);
    bool has_any_backbone = (strand1_forward + strand1_reverse + strand2_forward + strand2_reverse) > 0;
    
    if (has_any_backbone && backbone_consistent) {
        // Backbone is consistent - use it to determine direction
        // If strand 1 is reverse (-1) for most steps, that's the anti-parallel pattern
        // No reversal needed if pattern is consistent
        
        // Check first step for strand swap
        size_t idx_m = pair_order[helix.start_idx];
        size_t idx_n = pair_order[helix.start_idx + 1];
        
        size_t i1, j1, i2, j2;
        get_ij(pairs[idx_m], swapped[idx_m], i1, j1);
        get_ij(pairs[idx_n], swapped[idx_n], i2, j2);
        
        int k = is_linked(i1, i2, backbone);
        if (k == -1) {
            swapped[idx_m] = !swapped[idx_m];
        }
    } else {
        // Backbone is inconsistent or absent - use z-coordinate heuristic
        // First pair should have higher z than last pair for typical 5'->3' direction
        size_t first_idx = pair_order[helix.start_idx];
        size_t last_idx = pair_order[helix.end_idx];
        const auto& first_pair = pairs[first_idx];
        const auto& last_pair = pairs[last_idx];
        
        double first_z = 0, last_z = 0;
        if (first_pair.frame1().has_value()) {
            first_z = first_pair.frame1().value().origin().z();
        }
        if (last_pair.frame1().has_value()) {
            last_z = last_pair.frame1().value().origin().z();
        }
        
        // Reverse if first pair has lower z than last (going wrong direction)
        if (first_z < last_z) {
            std::reverse(pair_order.begin() + helix.start_idx,
                        pair_order.begin() + helix.end_idx + 1);
        }
    }
}

bool HelixOrganizer::wc_bporien(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                 bool swap_m, bool swap_n,
                                 const BackboneData& backbone) const {
    // Only check WC pairs (bp_type > 0 in legacy)
    // For now, check if both pairs have proper frames
    if (!pair_m.frame1().has_value() || !pair_m.frame2().has_value() ||
        !pair_n.frame1().has_value() || !pair_n.frame2().has_value()) {
        return false;
    }
    
    size_t i1, j1, i2, j2;
    get_ij(pair_m, swap_m, i1, j1);
    get_ij(pair_n, swap_n, i2, j2);
    
    // If x-angle is too large or backbone is linked, don't swap
    if (wcbp_xang(pair_m, pair_n) > config_.end_stack_xang ||
        is_linked(i1, i2, backbone) || is_linked(j1, j2, backbone)) {
        return false;
    }
    
    // Check z-direction alignment
    double zdir_normal = wcbp_zdir(pair_m, pair_n, swap_m, swap_n);
    double zdir_swapped = wcbp_zdir(pair_m, pair_n, swap_m, !swap_n);
    
    if (zdir_normal < 0.0 && zdir_swapped > 0.0) {
        return true;  // Need to swap pair_n
    }
    
    return false;
}

bool HelixOrganizer::check_o3dist(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const {
    size_t i1, j1, i2, j2;
    get_ij(pair_m, swap_m, i1, j1);
    get_ij(pair_n, swap_n, i2, j2);
    
    double di1_i2 = o3_distance(i1, i2, backbone);
    double di1_j2 = o3_distance(i1, j2, backbone);
    double dj1_i2 = o3_distance(j1, i2, backbone);
    double dj1_j2 = o3_distance(j1, j2, backbone);
    
    // If i1-i2 > i1-j2 AND j1-j2 > j1-i2, swap is indicated
    if ((di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2) &&
        (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2)) {
        return true;
    }
    
    return false;
}

bool HelixOrganizer::check_schain(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const {
    size_t i1, j1, i2, j2;
    get_ij(pair_m, swap_m, i1, j1);
    get_ij(pair_n, swap_n, i2, j2);
    
    // If no same-strand linkage but cross-strand linkage exists, swap is indicated
    if (!is_linked(i1, i2, backbone) && !is_linked(j1, j2, backbone) &&
        (is_linked(i1, j2, backbone) || is_linked(j1, i2, backbone))) {
        return true;
    }
    
    return false;
}

bool HelixOrganizer::check_others(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                   bool swap_m, bool swap_n,
                                   const BackboneData& backbone) const {
    size_t i1, j1, i2, j2;
    get_ij(pair_m, swap_m, i1, j1);
    get_ij(pair_n, swap_n, i2, j2);
    
    // If any backbone linkage exists, no swap needed
    if (is_linked(i1, i2, backbone) || is_linked(j1, j2, backbone) ||
        is_linked(i1, j2, backbone) || is_linked(j1, i2, backbone)) {
        return false;
    }
    
    // Check frame alignment patterns
    if (!pair_m.frame1().has_value() || !pair_m.frame2().has_value() ||
        !pair_n.frame1().has_value() || !pair_n.frame2().has_value()) {
        return false;
    }
    
    // Get frames based on swap status
    auto frame_m1 = swap_m ? pair_m.frame2().value() : pair_m.frame1().value();
    auto frame_m2 = swap_m ? pair_m.frame1().value() : pair_m.frame2().value();
    auto frame_n1 = swap_n ? pair_n.frame2().value() : pair_n.frame1().value();
    auto frame_n2 = swap_n ? pair_n.frame1().value() : pair_n.frame2().value();
    
    // Check axis alignment (similar to legacy dot product checks)
    double a1_x = frame_m1.x_axis().dot(frame_n1.x_axis());
    double a1_y = frame_m1.y_axis().dot(frame_n1.y_axis());
    double a1_z = frame_m1.z_axis().dot(frame_n1.z_axis());
    
    double a2_x = frame_m2.x_axis().dot(frame_n2.x_axis());
    double a2_y = frame_m2.y_axis().dot(frame_n2.y_axis());
    double a2_z = frame_m2.z_axis().dot(frame_n2.z_axis());
    
    bool aligned1 = (a1_x > 0.0 && a1_y > 0.0 && a1_z > 0.0);
    bool aligned2 = (a2_x > 0.0 && a2_y > 0.0 && a2_z > 0.0);
    
    if (aligned1 && aligned2) {
        return false;
    }
    
    // Check cross-alignment (m1 with n2, m2 with n1)
    double r1_x = frame_m1.x_axis().dot(frame_n2.x_axis());
    double r1_y = frame_m1.y_axis().dot(frame_n2.y_axis());
    double r1_z = frame_m1.z_axis().dot(frame_n2.z_axis());
    
    double r2_x = frame_m2.x_axis().dot(frame_n1.x_axis());
    double r2_y = frame_m2.y_axis().dot(frame_n1.y_axis());
    double r2_z = frame_m2.z_axis().dot(frame_n1.z_axis());
    
    bool cross1 = (r1_x > 0.0 && r1_y > 0.0 && r1_z > 0.0);
    bool cross2 = (r2_x > 0.0 && r2_y > 0.0 && r2_z > 0.0);
    
    if (!aligned1 && !aligned2) {
        if (cross1 || cross2) {
            return true;
        }
    }
    
    // Compare total angles for mixed cases
    if ((aligned1 || aligned2) && (cross1 || cross2)) {
        double sum_aligned = dot2ang(a1_x) + dot2ang(a1_y) + dot2ang(a1_z) +
                            dot2ang(a2_x) + dot2ang(a2_y) + dot2ang(a2_z);
        double sum_cross = dot2ang(r1_x) + dot2ang(r1_y) + dot2ang(r1_z) +
                          dot2ang(r2_x) + dot2ang(r2_y) + dot2ang(r2_z);
        
        if (sum_aligned > sum_cross) {
            return true;
        }
    }
    
    return false;
}

bool HelixOrganizer::chain1dir(const core::BasePair& pair_m, const core::BasePair& pair_n,
                                bool swap_m, bool swap_n,
                                const BackboneData& backbone) const {
    size_t i1, j1, i2, j2;
    get_ij(pair_m, swap_m, i1, j1);
    get_ij(pair_n, swap_n, i2, j2);
    
    int k = is_linked(i1, i2, backbone);
    return (k == -1);  // Reverse linkage indicates need to swap
}

DirectionCounts HelixOrganizer::check_direction(
    const std::vector<core::BasePair>& pairs,
    const BackboneData& backbone,
    const std::vector<size_t>& pair_order,
    const HelixSegment& helix,
    const std::vector<bool>& swapped) const {
    
    DirectionCounts dir;
    
    for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
        size_t idx_m = pair_order[pos];
        size_t idx_n = pair_order[pos + 1];
        
        const auto& pair_m = pairs[idx_m];
        const auto& pair_n = pairs[idx_n];
        
        size_t i1, j1, i2, j2;
        get_ij(pair_m, swapped[idx_m], i1, j1);
        get_ij(pair_n, swapped[idx_n], i2, j2);
        
        // Strand 1 direction
        int k = is_linked(i1, i2, backbone);
        if (k == 1) dir.strand1_forward++;
        else if (k == -1) dir.strand1_reverse++;
        else dir.strand1_none++;
        
        // Strand 2 direction
        k = is_linked(j1, j2, backbone);
        if (k == 1) dir.strand2_forward++;
        else if (k == -1) dir.strand2_reverse++;
        else dir.strand2_none++;
    }
    
    return dir;
}

void HelixOrganizer::check_strand2(const std::vector<core::BasePair>& pairs,
                                   const BackboneData& backbone,
                                   const std::vector<size_t>& pair_order,
                                   HelixSegment& helix,
                                   std::vector<bool>& swapped,
                                   const DirectionCounts& direction) const {
    
    bool mixed_direction = (direction.strand1_forward && direction.strand1_reverse) ||
                          (direction.strand2_forward && direction.strand2_reverse);
    
    if (!mixed_direction) {
        // Normal case - check for cross-strand swaps
        if (direction.strand1_forward + direction.strand1_reverse +
            direction.strand2_forward + direction.strand2_reverse == 0) {
            return;  // No linkages to check
        }
        
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];
            
            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];
            
            // Skip if WC orientation check would pass
            if (wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone)) {
                continue;
            }
            
            size_t i1, j1, i2, j2;
            get_ij(pair_m, swapped[idx_m], i1, j1);
            get_ij(pair_n, swapped[idx_n], i2, j2);
            
            // Check for cross-strand linkage indicating swap needed
            if (!is_linked(i1, i2, backbone) && !is_linked(j1, j2, backbone) &&
                ((is_linked(i1, j2, backbone) == 1) ||
                 (is_linked(i1, j2, backbone) && is_linked(j1, i2, backbone)))) {
                swapped[idx_n] = !swapped[idx_n];
            }
        }
    } else {
        // Mixed direction case - more complex handling
        bool anti_p = (direction.strand1_forward > direction.strand1_reverse) &&
                     (direction.strand2_forward < direction.strand2_reverse);
        bool parallel = (direction.strand1_forward > direction.strand1_reverse) &&
                       (direction.strand2_forward > direction.strand2_reverse);
        
        helix.is_parallel = parallel;
        
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];
            
            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];
            
            size_t i1, j1, i2, j2;
            get_ij(pair_m, swapped[idx_m], i1, j1);
            get_ij(pair_n, swapped[idx_n], i2, j2);
            
            int k = is_linked(j1, j2, backbone);
            if (!is_linked(i1, i2, backbone) && 
                ((anti_p && k == 1) || (parallel && k == -1))) {
                swapped[idx_n] = !swapped[idx_n];
            }
            
            // Re-get j2 after potential swap
            get_ij(pair_n, swapped[idx_n], i2, j2);
            
            if (!is_linked(i1, i2, backbone) && !is_linked(j1, j2, backbone)) {
                if ((anti_p && is_linked(j1, i2, backbone) == 1) ||
                    (parallel && is_linked(i1, j2, backbone) == -1)) {
                    swapped[idx_m] = !swapped[idx_m];
                } else if ((anti_p && is_linked(i1, j2, backbone) == 1) ||
                           (parallel && is_linked(j1, i2, backbone) == -1)) {
                    swapped[idx_n] = !swapped[idx_n];
                }
            }
        }
    }
}

// =============================================================================
// Context calculation (bp_context equivalent)
// =============================================================================

std::vector<HelixOrganizer::PairContext> HelixOrganizer::calculate_context(
    const std::vector<core::BasePair>& pairs) const {
    
    size_t n = pairs.size();
    std::vector<PairContext> context(n);
    
    if (n < 2) return context;
    
    for (size_t i = 0; i < n; ++i) {
        auto org_i = get_pair_origin(pairs[i]);
        auto z_i = get_pair_z_axis(pairs[i]);
        
        std::vector<std::pair<double, size_t>> neighbors;
        
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            
            auto org_j = get_pair_origin(pairs[j]);
            double dist = (org_j - org_i).length();
            
            if (dist <= config_.neighbor_cutoff) {
                neighbors.emplace_back(dist, j);
            }
        }
        
        std::sort(neighbors.begin(), neighbors.end());
        
        if (neighbors.empty()) {
            context[i].is_endpoint = true;
            continue;
        }
        
        context[i].neighbor1 = neighbors[0].second;
        context[i].dist1 = neighbors[0].first;
        
        if (context[i].dist1 > config_.helix_break) {
            context[i].is_endpoint = true;
            continue;
        }
        
        auto v1 = get_pair_origin(pairs[neighbors[0].second]) - org_i;
        double d1 = z_i.dot(v1);
        
        for (size_t k = 1; k < neighbors.size(); ++k) {
            if (neighbors[k].first > config_.helix_break) break;
            
            auto vk = get_pair_origin(pairs[neighbors[k].second]) - org_i;
            double dk = z_i.dot(vk);
            
            if (d1 * dk < 0) {
                context[i].neighbor2 = neighbors[k].second;
                context[i].dist2 = neighbors[k].first;
                context[i].is_endpoint = false;
                break;
            }
        }
        
        if (!context[i].neighbor2.has_value()) {
            context[i].is_endpoint = true;
        }
    }
    
    return context;
}

std::vector<size_t> HelixOrganizer::find_endpoints(
    const std::vector<PairContext>& context) const {
    
    std::vector<size_t> endpoints;
    
    for (size_t i = 0; i < context.size(); ++i) {
        if (context[i].is_endpoint) {
            endpoints.push_back(i);
        }
    }
    
    if (endpoints.empty() && !context.empty()) {
        endpoints.push_back(0);
    }
    
    return endpoints;
}

std::pair<std::vector<size_t>, std::vector<HelixSegment>> HelixOrganizer::locate_helices(
    const std::vector<PairContext>& context,
    const std::vector<size_t>& endpoints,
    size_t num_pairs) const {
    
    std::vector<size_t> pair_order;
    std::vector<HelixSegment> helices;
    std::vector<bool> visited(num_pairs, false);
    
    pair_order.reserve(num_pairs);
    
    for (size_t ep : endpoints) {
        if (visited[ep]) continue;
        
        HelixSegment helix;
        helix.start_idx = pair_order.size();
        
        size_t current = ep;
        std::optional<size_t> prev;
        
        while (!visited[current]) {
            visited[current] = true;
            pair_order.push_back(current);
            
            const auto& ctx = context[current];
            std::optional<size_t> next;
            
            if (ctx.neighbor1.has_value() && !visited[ctx.neighbor1.value()]) {
                if (!prev.has_value() || ctx.neighbor1.value() != prev.value()) {
                    next = ctx.neighbor1;
                }
            }
            
            if (!next.has_value() && ctx.neighbor2.has_value() && 
                !visited[ctx.neighbor2.value()]) {
                if (!prev.has_value() || ctx.neighbor2.value() != prev.value()) {
                    next = ctx.neighbor2;
                }
            }
            
            if (!next.has_value()) break;
            
            prev = current;
            current = next.value();
        }
        
        helix.end_idx = pair_order.size() - 1;
        if (helix.end_idx >= helix.start_idx) {
            helices.push_back(helix);
        }
    }
    
    for (size_t i = 0; i < num_pairs; ++i) {
        if (!visited[i]) {
            HelixSegment helix;
            helix.start_idx = pair_order.size();
            pair_order.push_back(i);
            helix.end_idx = pair_order.size() - 1;
            helices.push_back(helix);
        }
    }
    
    return {pair_order, helices};
}

// =============================================================================
// Main five2three algorithm
// =============================================================================

void HelixOrganizer::ensure_five_to_three(
    const std::vector<core::BasePair>& pairs,
    const BackboneData& backbone,
    std::vector<size_t>& pair_order,
    std::vector<HelixSegment>& helices,
    std::vector<bool>& swapped) const {
    
    swapped.resize(pairs.size(), false);
    
    if (backbone.empty()) {
        return;
    }
    
    // Process each helix
    for (auto& helix : helices) {
        if (helix.start_idx > helix.end_idx) continue;
        
        // STEP 1: first_step - set initial strand assignment
        // This also handles helix reversal through backbone connectivity
        first_step(pairs, backbone, pair_order, helix, swapped);
        
        // STEP 3: First pass through steps - check each consecutive pair
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];
            
            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];
            
            bool rev_wc = wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_o3d = check_o3dist(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_csc = check_schain(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            bool rev_oth = check_others(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            
            // Apply swap based on checks
            if (rev_wc) {
                swapped[idx_n] = !swapped[idx_n];
            } else if (rev_o3d || rev_csc || rev_oth) {
                swapped[idx_n] = !swapped[idx_n];
            }
            
            // Check strand 1 direction
            bool rev_s1 = chain1dir(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (rev_s1) {
                swapped[idx_n] = !swapped[idx_n];
            }
        }
        
        // STEP 4: Second pass - re-check WC orientation
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];
            
            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];
            
            bool rev_wc = wc_bporien(pair_m, pair_n, swapped[idx_m], swapped[idx_n], backbone);
            if (rev_wc) {
                swapped[idx_m] = !swapped[idx_m];
            }
        }
        
        // STEP 5: check_direction - count backbone linkage directions
        DirectionCounts direction = check_direction(pairs, backbone, pair_order, helix, swapped);
        
        // STEP 6: check_strand2 - additional corrections based on direction
        check_strand2(pairs, backbone, pair_order, helix, swapped, direction);
    }
}

// =============================================================================
// Main organize function
// =============================================================================

HelixOrdering HelixOrganizer::organize(const std::vector<core::BasePair>& pairs,
                                        const BackboneData& backbone) const {
    HelixOrdering result;
    
    if (pairs.empty()) {
        return result;
    }
    
    if (pairs.size() == 1) {
        result.pair_order = {0};
        result.helices = {{0, 0, false, false, false}};
        result.strand_swapped = {false};
        return result;
    }
    
    // Step 1: Calculate neighbor context
    auto context = calculate_context(pairs);
    
    // Step 2: Find helix endpoints
    auto endpoints = find_endpoints(context);
    
    // Step 3: Chain pairs into helices
    auto [pair_order, helices] = locate_helices(context, endpoints, pairs.size());
    
    // Step 4: Ensure 5'→3' direction (full five2three algorithm)
    std::vector<bool> strand_swapped;
    ensure_five_to_three(pairs, backbone, pair_order, helices, strand_swapped);
    
    result.pair_order = std::move(pair_order);
    result.helices = std::move(helices);
    result.strand_swapped = std::move(strand_swapped);
    
    return result;
}

} // namespace x3dna::algorithms
