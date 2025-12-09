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
    std::vector<size_t>& /*pair_order*/,
    std::vector<HelixSegment>& /*helices*/,
    std::vector<bool>& strand_swapped) const {
    
    strand_swapped.resize(pairs.size(), false);
    
    // For now, skip the 5'→3' direction check
    // This requires O3'-P atom connectivity information which we don't have yet
    // The basic helix ordering should still work for step parameter calculation
    
    // TODO: Implement proper 5'→3' direction check using backbone atoms
}

HelixOrdering HelixOrganizer::organize(const std::vector<core::BasePair>& pairs) const {
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
    
    // Step 4: Ensure 5'→3' direction
    std::vector<bool> strand_swapped;
    ensure_five_to_three(pairs, pair_order, helices, strand_swapped);
    
    result.pair_order = std::move(pair_order);
    result.helices = std::move(helices);
    result.strand_swapped = std::move(strand_swapped);
    
    return result;
}

} // namespace x3dna::algorithms

