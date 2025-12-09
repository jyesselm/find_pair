/**
 * @file helix_organizer.cpp
 * @brief Implementation of helix organization algorithm
 */

#include <x3dna/algorithms/helix_organizer.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace x3dna::algorithms {

HelixOrganizer::HelixOrganizer(const Config& config) : config_(config) {}

geometry::Vector3D HelixOrganizer::get_pair_origin(const core::BasePair& pair) const {
    if (!pair.frame1().has_value()) {
        return geometry::Vector3D(0, 0, 0);
    }
    // Use frame1 origin (strand 1 base) - matches legacy behavior
    return pair.frame1().value().origin();
}

geometry::Vector3D HelixOrganizer::get_pair_z_axis(const core::BasePair& pair) const {
    if (!pair.frame1().has_value()) {
        return geometry::Vector3D(0, 0, 1);
    }
    // Use frame1 z-axis
    return pair.frame1().value().z_axis();
}

int HelixOrganizer::is_linked(size_t i, size_t j, const BackboneData& backbone) const {
    // Check if O3' of residue i is linked to P of residue j
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

std::vector<HelixOrganizer::PairContext> HelixOrganizer::calculate_context(
    const std::vector<core::BasePair>& pairs) const {
    
    size_t n = pairs.size();
    std::vector<PairContext> context(n);
    
    if (n < 2) return context;
    
    // For each pair, find its neighbors
    for (size_t i = 0; i < n; ++i) {
        auto org_i = get_pair_origin(pairs[i]);
        auto z_i = get_pair_z_axis(pairs[i]);
        
        // Find all neighbors within cutoff, sorted by distance
        std::vector<std::pair<double, size_t>> neighbors;
        
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            
            auto org_j = get_pair_origin(pairs[j]);
            double dist = (org_j - org_i).length();
            
            if (dist <= config_.neighbor_cutoff) {
                neighbors.emplace_back(dist, j);
            }
        }
        
        // Sort by distance
        std::sort(neighbors.begin(), neighbors.end());
        
        if (neighbors.empty()) {
            // Isolated pair
            context[i].is_endpoint = true;
            continue;
        }
        
        // Nearest neighbor
        context[i].neighbor1 = neighbors[0].second;
        context[i].dist1 = neighbors[0].first;
        
        // Check if nearest is too far (helix break)
        if (context[i].dist1 > config_.helix_break) {
            context[i].is_endpoint = true;
            continue;
        }
        
        // Look for 2nd neighbor on opposite z-side
        auto v1 = get_pair_origin(pairs[neighbors[0].second]) - org_i;
        double d1 = z_i.dot(v1);  // Projection onto z-axis
        
        for (size_t k = 1; k < neighbors.size(); ++k) {
            if (neighbors[k].first > config_.helix_break) break;
            
            auto vk = get_pair_origin(pairs[neighbors[k].second]) - org_i;
            double dk = z_i.dot(vk);
            
            // Check if on opposite z-side
            if (d1 * dk < 0) {
                context[i].neighbor2 = neighbors[k].second;
                context[i].dist2 = neighbors[k].first;
                context[i].is_endpoint = false;
                break;
            }
        }
        
        // If no opposite-side neighbor found, this is an endpoint
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
    
    // If no endpoints found (circular structure), start from pair 0
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
        
        // Start chain from this endpoint
        size_t current = ep;
        std::optional<size_t> prev;
        
        while (!visited[current]) {
            visited[current] = true;
            pair_order.push_back(current);
            
            // Find next unvisited neighbor
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
    
    // Add any remaining unvisited pairs (isolated or in separate components)
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

void HelixOrganizer::ensure_five_to_three(
    const std::vector<core::BasePair>& pairs,
    const BackboneData& backbone,
    std::vector<size_t>& pair_order,
    std::vector<HelixSegment>& helices,
    std::vector<bool>& strand_swapped) const {
    
    strand_swapped.resize(pairs.size(), false);
    
    // If no backbone data provided, skip 5'→3' checking
    if (backbone.empty()) {
        return;
    }
    
    // For each helix, check backbone connectivity and determine:
    // 1. Whether the helix segment needs to be reversed (going wrong direction)
    // 2. Which pairs need strand swapping
    for (auto& helix : helices) {
        if (helix.start_idx >= helix.end_idx) continue;
        
        // Count forward and reverse backbone linkages to determine direction
        int forward_links = 0;   // i1 → i2 or j1 → j2 (continuing same strands)
        int cross_links = 0;     // i1 → j2 or j1 → i2 (cross-strand, needs swap)
        
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];
            
            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];
            
            size_t i1 = pair_m.residue_idx1();
            size_t j1_m = pair_m.residue_idx2();
            size_t i2 = pair_n.residue_idx1();
            size_t j2 = pair_n.residue_idx2();
            
            // Check all linkage patterns
            if (is_linked(i1, i2, backbone) == 1 || is_linked(j1_m, j2, backbone) == 1) {
                forward_links++;
            }
            if (is_linked(i1, j2, backbone) == 1 || is_linked(j1_m, i2, backbone) == 1) {
                cross_links++;
            }
        }
        
        // If most links are cross-strand, all pairs in this helix need swapping
        bool swap_all = (cross_links > forward_links);
        
        // Apply strand swapping
        for (size_t pos = helix.start_idx; pos <= helix.end_idx; ++pos) {
            size_t idx = pair_order[pos];
            strand_swapped[idx] = swap_all;
        }
        
        // After swapping, do fine-grained check for remaining issues
        for (size_t pos = helix.start_idx; pos < helix.end_idx; ++pos) {
            size_t idx_m = pair_order[pos];
            size_t idx_n = pair_order[pos + 1];
            
            const auto& pair_m = pairs[idx_m];
            const auto& pair_n = pairs[idx_n];
            
            // Get residue indices based on current swap status
            size_t i1 = strand_swapped[idx_m] ? pair_m.residue_idx2() : pair_m.residue_idx1();
            size_t j1_m = strand_swapped[idx_m] ? pair_m.residue_idx1() : pair_m.residue_idx2();
            size_t i2_default = pair_n.residue_idx1();
            size_t j2_default = pair_n.residue_idx2();
            
            // If current swap status doesn't give proper linkage, toggle it
            int link_i1_i2 = is_linked(i1, i2_default, backbone);
            int link_j1_j2 = is_linked(j1_m, j2_default, backbone);
            int link_i1_j2 = is_linked(i1, j2_default, backbone);
            int link_j1_i2 = is_linked(j1_m, i2_default, backbone);
            
            // Check if swap is needed for pair n relative to current assignment
            bool needs_swap_toggle = false;
            
            if (!strand_swapped[idx_n]) {
                // Currently not swapped - check if we should swap
                if ((link_i1_j2 == 1 || link_j1_i2 == 1) && 
                    !(link_i1_i2 == 1 || link_j1_j2 == 1)) {
                    needs_swap_toggle = true;
                }
            } else {
                // Currently swapped - check if we should un-swap
                // When swapped, i2 becomes j2 and vice versa
                if ((link_i1_i2 == 1 || link_j1_j2 == 1) && 
                    !(link_i1_j2 == 1 || link_j1_i2 == 1)) {
                    needs_swap_toggle = true;
                }
            }
            
            if (needs_swap_toggle) {
                strand_swapped[idx_n] = !strand_swapped[idx_n];
            }
        }
    }
}

HelixOrdering HelixOrganizer::organize(const std::vector<core::BasePair>& pairs,
                                        const BackboneData& backbone) const {
    HelixOrdering result;
    
    if (pairs.empty()) {
        return result;
    }
    
    if (pairs.size() == 1) {
        result.pair_order = {0};
        result.helices = {{0, 0, false, false}};
        result.strand_swapped = {false};
        return result;
    }
    
    // Step 1: Calculate neighbor context for each pair
    auto context = calculate_context(pairs);
    
    // Step 2: Find helix endpoints
    auto endpoints = find_endpoints(context);
    
    // Step 3: Chain pairs into helices
    auto [pair_order, helices] = locate_helices(context, endpoints, pairs.size());
    
    // Step 4: Ensure 5'→3' direction using backbone connectivity
    std::vector<bool> strand_swapped;
    ensure_five_to_three(pairs, backbone, pair_order, helices, strand_swapped);
    
    result.pair_order = std::move(pair_order);
    result.helices = std::move(helices);
    result.strand_swapped = std::move(strand_swapped);
    
    return result;
}

} // namespace x3dna::algorithms
