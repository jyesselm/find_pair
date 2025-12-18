/**
 * @file helix_context_calculator.cpp
 * @brief Implementation of helix context calculation
 */

#include <x3dna/algorithms/helix/helix_context_calculator.hpp>
#include <x3dna/config/config_manager.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace x3dna::algorithms::helix {

namespace {
/**
 * @brief Check if five2three debugging is enabled
 */
bool is_debug_enabled() {
    static bool initialized = false;
    static bool debug_enabled = false;
    if (!initialized) {
        auto& cfg = config::ConfigManager::instance();
        cfg.init_debug_from_environment();
        debug_enabled = cfg.debug_config().debug_five2three;
        initialized = true;
    }
    return debug_enabled;
}
} // namespace

std::vector<PairContext> HelixContextCalculator::calculate_context(const std::vector<core::BasePair>& pairs,
                                                                   const BackboneData& backbone) const {

    size_t n = pairs.size();
    std::vector<PairContext> context(n);

    if (n < 2)
        return context;

    for (size_t i = 0; i < n; ++i) {
        auto org_i = PairGeometryHelper::get_pair_origin(pairs[i]);
        auto z_i = PairGeometryHelper::get_pair_z_axis(pairs[i]);

        std::vector<std::pair<double, size_t>> neighbors;

        for (size_t j = 0; j < n; ++j) {
            if (j == i)
                continue;

            auto org_j = PairGeometryHelper::get_pair_origin(pairs[j]);
            double dist = (org_j - org_i).length();

            if (dist <= config_.neighbor_cutoff) {
                neighbors.emplace_back(dist, j);
            }
        }

        std::sort(neighbors.begin(), neighbors.end());

        // Legacy behavior: if no neighbors within helix_break, pair is an isolated endpoint
        // with NO stored neighbors (end_list only stores the endpoint itself)
        if (neighbors.empty() || neighbors[0].first > config_.helix_break) {
            context[i].is_endpoint = true;
            continue; // Don't store neighbor1 - matches legacy line 963-965
        }

        context[i].neighbor1 = neighbors[0].second;
        context[i].dist1 = neighbors[0].first;

        // Check backbone connectivity to neighbor1
        context[i].has_backbone_link1 = linkage_checker_.are_pairs_connected(pairs[i], pairs[neighbors[0].second],
                                                                             backbone);

        auto v1 = PairGeometryHelper::get_pair_origin(pairs[neighbors[0].second]) - org_i;
        double d1 = z_i.dot(v1);

        // Legacy lines 931-941: If 2nd and 3rd closest are both on opposite z-side,
        // swap them if 2nd has larger |z-distance| (prefer smaller |z-distance|)
        if (neighbors.size() >= 3 && neighbors[1].first <= config_.helix_break &&
            neighbors[2].first <= config_.helix_break) {
            auto v2 = PairGeometryHelper::get_pair_origin(pairs[neighbors[1].second]) - org_i;
            auto v3 = PairGeometryHelper::get_pair_origin(pairs[neighbors[2].second]) - org_i;
            double d2 = z_i.dot(v2);
            double d3 = z_i.dot(v3);

            // Both on opposite z-side from n1, and 2nd has larger |z-dist|
            const bool both_opposite_from_n1 = are_on_opposite_z_sides(d1, d2) && are_on_opposite_z_sides(d1, d3);
            const bool second_has_larger_z_dist = std::abs(d2) > std::abs(d3);
            if (both_opposite_from_n1 && second_has_larger_z_dist) {
                std::swap(neighbors[1], neighbors[2]);
            }
        }

        for (size_t k = 1; k < neighbors.size(); ++k) {
            if (neighbors[k].first > config_.helix_break)
                break;

            auto vk = PairGeometryHelper::get_pair_origin(pairs[neighbors[k].second]) - org_i;
            double dk = z_i.dot(vk);

            if (are_on_opposite_z_sides(d1, dk)) {
                context[i].neighbor2 = neighbors[k].second;
                context[i].dist2 = neighbors[k].first;
                // Check backbone connectivity to neighbor2
                context[i].has_backbone_link2 = linkage_checker_.are_pairs_connected(pairs[i],
                                                                                     pairs[neighbors[k].second],
                                                                                     backbone);
                context[i].is_endpoint = false;
                break;
            }
        }

        if (!context[i].neighbor2.has_value()) {
            context[i].is_endpoint = true;

            // Legacy special case (find_pair.c lines 971-976):
            // Even for endpoints, try to find n2 through an indirect check:
            // If vector from n1 to 2nd closest is on opposite z-side AND within helix_break
            if (neighbors.size() >= 2) {
                size_t n2_idx = neighbors[1].second;
                auto org_n1 = PairGeometryHelper::get_pair_origin(pairs[neighbors[0].second]);
                auto org_n2 = PairGeometryHelper::get_pair_origin(pairs[n2_idx]);
                // Legacy ddxyz(n1, n2) = n1 - n2, so we compute org_n1 - org_n2
                auto v_n2_n1 = org_n1 - org_n2;
                double dist_n1_n2 = v_n2_n1.length();
                double d2 = z_i.dot(v_n2_n1);

                // If n2->n1 is on opposite z-side from i->n1 AND within helix_break
                if (are_on_opposite_z_sides(d1, d2) && dist_n1_n2 <= config_.helix_break) {
                    context[i].neighbor2 = n2_idx;
                    context[i].dist2 = neighbors[1].first;
                    context[i].has_backbone_link2 = linkage_checker_.are_pairs_connected(pairs[i], pairs[n2_idx],
                                                                                         backbone);
                    // Still an endpoint, but now has a neighbor2
                }
            }
        }
    }

    return context;
}

std::vector<size_t> HelixContextCalculator::find_endpoints(const std::vector<PairContext>& context) const {

    bool debug = is_debug_enabled();

    std::vector<size_t> endpoints;

    for (size_t i = 0; i < context.size(); ++i) {
        if (context[i].is_endpoint) {
            endpoints.push_back(i);
        }
    }

    if (endpoints.empty() && !context.empty()) {
        endpoints.push_back(0);
    }

    if (debug) {
        std::cerr << "[find_endpoints] Endpoints found: ";
        for (size_t ep : endpoints) {
            std::cerr << ep << " ";
        }
        std::cerr << std::endl;
        for (size_t i = 0; i < context.size(); ++i) {
            std::cerr << "[context] pair " << i << ": ep=" << context[i].is_endpoint << " n1="
                      << (context[i].neighbor1.has_value() ? std::to_string(context[i].neighbor1.value()) : "-")
                      << " n2="
                      << (context[i].neighbor2.has_value() ? std::to_string(context[i].neighbor2.value()) : "-")
                      << std::endl;
        }
    }

    return endpoints;
}

std::pair<std::vector<size_t>, std::vector<HelixSegment>> HelixContextCalculator::locate_helices(
    const std::vector<PairContext>& context, const std::vector<size_t>& endpoints, const BackboneData& backbone,
    size_t num_pairs) const {

    bool debug = is_debug_enabled();

    std::vector<size_t> pair_order;
    std::vector<HelixSegment> helices;
    std::vector<bool> visited(num_pairs, false);

    (void)backbone; // Not used in traverse - backbone connectivity checked in calculate_context
    pair_order.reserve(num_pairs);

    // Legacy locate_helix algorithm (find_pair.c lines 1029-1082):
    // For each endpoint, add endpoint + neighbors from end_list, then traverse
    for (size_t ep : endpoints) {
        if (debug) {
            std::cerr << "[locate_helices] Processing endpoint " << ep << std::endl;
        }
        // Skip if all items in this endpoint's "end_list" are already matched
        const auto& ep_ctx = context[ep];
        int matched_count = 0;
        int item_count = 1; // The endpoint itself
        if (visited[ep])
            matched_count++;
        if (ep_ctx.neighbor1.has_value()) {
            item_count++;
            if (visited[ep_ctx.neighbor1.value()])
                matched_count++;
        }
        if (ep_ctx.neighbor2.has_value()) {
            item_count++;
            if (visited[ep_ctx.neighbor2.value()])
                matched_count++;
        }
        if (matched_count == item_count)
            continue;

        HelixSegment helix;
        helix.start_idx = pair_order.size();

        // Legacy lines 1039-1045: Add endpoint and its neighbors from end_list
        // Add endpoint
        if (!visited[ep]) {
            pair_order.push_back(ep);
            visited[ep] = true;
        }
        // Add neighbor1
        if (ep_ctx.neighbor1.has_value() && !visited[ep_ctx.neighbor1.value()]) {
            pair_order.push_back(ep_ctx.neighbor1.value());
            visited[ep_ctx.neighbor1.value()] = true;
        }
        // Add neighbor2
        if (ep_ctx.neighbor2.has_value() && !visited[ep_ctx.neighbor2.value()]) {
            pair_order.push_back(ep_ctx.neighbor2.value());
            visited[ep_ctx.neighbor2.value()] = true;
        }

        // Legacy lines 1046-1068: Traverse from the last added pair
        // Continue adding neighbors until we can't anymore
        while (pair_order.size() >= helix.start_idx + 1) {
            size_t ip = pair_order.size() - 1;
            size_t current = pair_order[ip];
            const auto& ctx = context[current];

            // If this is an endpoint (bp_order[k][1] == 0 in legacy means endpoint)
            if (ctx.is_endpoint) {
                // Legacy lines 1050-1055: Only add neighbor1 if it exists, not matched, and no neighbor2
                if (ctx.neighbor1.has_value() && !visited[ctx.neighbor1.value()] && !ctx.neighbor2.has_value()) {
                    pair_order.push_back(ctx.neighbor1.value());
                    visited[ctx.neighbor1.value()] = true;
                }
                break; // Stop traversal for this helix
            }

            // Not an endpoint - continue traversal
            bool n1_matched = !ctx.neighbor1.has_value() || visited[ctx.neighbor1.value()];
            bool n2_matched = !ctx.neighbor2.has_value() || visited[ctx.neighbor2.value()];

            // Legacy line 1057-1059: If both or neither neighbor matched, stop
            if ((n1_matched && n2_matched) ||
                (!n1_matched && !n2_matched && ctx.neighbor1.has_value() && ctx.neighbor2.has_value())) {
                // Both matched or neither matched - stop
                // For "neither matched" case, only apply if both neighbors actually exist
                break;
            }

            // Get previous pair
            std::optional<size_t> prev;
            if (ip > helix.start_idx) {
                prev = pair_order[ip - 1];
            }

            // Legacy lines 1060-1067: Go to the other neighbor
            std::optional<size_t> next;
            if (ctx.neighbor1.has_value() && prev.has_value() && ctx.neighbor1.value() == prev.value()) {
                // Previous was neighbor1, add neighbor2
                if (ctx.neighbor2.has_value() && !visited[ctx.neighbor2.value()]) {
                    next = ctx.neighbor2;
                }
            } else if (ctx.neighbor2.has_value() && prev.has_value() && ctx.neighbor2.value() == prev.value()) {
                // Previous was neighbor2, add neighbor1
                if (ctx.neighbor1.has_value() && !visited[ctx.neighbor1.value()]) {
                    next = ctx.neighbor1;
                }
            }

            if (!next.has_value()) {
                break;
            }

            pair_order.push_back(next.value());
            visited[next.value()] = true;
        }

        helix.end_idx = pair_order.size() - 1;
        if (helix.end_idx >= helix.start_idx) {
            helices.push_back(helix);
            if (debug) {
                std::cerr << "[locate_helices] Created helix " << helices.size() << " (pos " << helix.start_idx << "-"
                          << helix.end_idx << "): ";
                for (size_t p = helix.start_idx; p <= helix.end_idx; ++p) {
                    std::cerr << pair_order[p] << " ";
                }
                std::cerr << std::endl;
            }
        }
    }

    // Handle any leftover pairs not reached from endpoints
    // Legacy behavior (lines 1120-1128): Put ALL leftover pairs into ONE helix region
    // This handles "complicated structures" with isolated pairs
    std::vector<size_t> leftover;
    for (size_t i = 0; i < num_pairs; ++i) {
        if (!visited[i]) {
            leftover.push_back(i);
        }
    }

    if (!leftover.empty()) {
        // All leftover pairs go into a single helix (matching legacy behavior)
        HelixSegment helix;
        helix.start_idx = pair_order.size();
        for (size_t idx : leftover) {
            pair_order.push_back(idx);
        }
        helix.end_idx = pair_order.size() - 1;
        helices.push_back(helix);
    }

    return {pair_order, helices};
}

std::vector<PairContextInfo> HelixContextCalculator::to_public_context(const std::vector<PairContext>& context) {

    std::vector<PairContextInfo> result;
    result.reserve(context.size());

    for (const auto& ctx : context) {
        PairContextInfo info;
        info.is_endpoint = ctx.is_endpoint;
        info.neighbor1 = ctx.neighbor1;
        info.neighbor2 = ctx.neighbor2;
        result.push_back(info);
    }

    return result;
}

} // namespace x3dna::algorithms::helix
