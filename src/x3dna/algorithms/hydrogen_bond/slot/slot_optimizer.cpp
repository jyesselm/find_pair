/**
 * @file slot_optimizer.cpp
 * @brief Implementation of slot-based H-bond optimizer
 */

#include <x3dna/algorithms/hydrogen_bond/slot/slot_optimizer.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/atom_capacity.hpp>
#include <algorithm>
#include <unordered_map>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

SlotOptimizer::SlotOptimizer(const SlotOptimizerParams& params)
    : params_(params) {
}

char SlotOptimizer::get_base_type(const core::Residue& residue) {
    // Try to get single-letter base type
    std::string name = residue.name();
    if (name.length() == 1) {
        return std::toupper(name[0]);
    }

    // Try parent base type for modified residues
    auto parent = AtomCapacity::get_parent_base_type(name);
    if (parent) {
        return *parent;
    }

    // Common 2-3 letter codes
    if (name == "DA" || name == "ADE") return 'A';
    if (name == "DG" || name == "GUA") return 'G';
    if (name == "DC" || name == "CYT") return 'C';
    if (name == "DT" || name == "THY") return 'T';
    if (name == "URA") return 'U';

    return '?';
}

bool SlotOptimizer::is_backbone_backbone(const std::string& atom1, const std::string& atom2) {
    return AtomCapacity::is_backbone_atom(atom1) && AtomCapacity::is_backbone_atom(atom2);
}

std::vector<core::HBond> SlotOptimizer::optimize_pair(
    const core::Residue& res1,
    const core::Residue& res2) {

    // Find all candidates
    auto candidates = find_candidates(res1, res2);

    if (candidates.empty()) {
        return {};
    }

    if (params_.baseline_mode) {
        return select_baseline(candidates, res1, res2);
    }

    // Create slot caches
    char base1 = get_base_type(res1);
    char base2 = get_base_type(res2);
    SlotCache cache1(res1, base1);
    SlotCache cache2(res2, base2);

    return select_optimal(candidates, cache1, cache2);
}

std::vector<HBondCandidate> SlotOptimizer::find_candidates(
    const core::Residue& res1,
    const core::Residue& res2) const {

    std::vector<HBondCandidate> candidates;
    std::string code1 = res1.name();
    std::string code2 = res2.name();

    // Check res1 donors -> res2 acceptors
    for (const auto& donor_atom : res1.atoms()) {
        std::string donor_name = AtomCapacity::normalize_atom_name(donor_atom.name());
        if (AtomCapacity::get_donor_capacity(code1, donor_name) == 0) {
            continue;
        }

        for (const auto& acceptor_atom : res2.atoms()) {
            std::string acceptor_name = AtomCapacity::normalize_atom_name(acceptor_atom.name());
            if (AtomCapacity::get_acceptor_capacity(code2, acceptor_name) == 0) {
                continue;
            }

            double dist = (acceptor_atom.position() - donor_atom.position()).length();
            if (dist <= params_.max_distance) {
                HBondCandidate c;
                c.donor_res_id = res1.res_id();
                c.acceptor_res_id = res2.res_id();
                c.donor_atom = donor_name;
                c.acceptor_atom = acceptor_name;
                c.donor_pos = donor_atom.position();
                c.acceptor_pos = acceptor_atom.position();
                c.distance = dist;
                candidates.push_back(c);
            }
        }
    }

    // Check res2 donors -> res1 acceptors
    for (const auto& donor_atom : res2.atoms()) {
        std::string donor_name = AtomCapacity::normalize_atom_name(donor_atom.name());
        if (AtomCapacity::get_donor_capacity(code2, donor_name) == 0) {
            continue;
        }

        for (const auto& acceptor_atom : res1.atoms()) {
            std::string acceptor_name = AtomCapacity::normalize_atom_name(acceptor_atom.name());
            if (AtomCapacity::get_acceptor_capacity(code1, acceptor_name) == 0) {
                continue;
            }

            double dist = (acceptor_atom.position() - donor_atom.position()).length();
            if (dist <= params_.max_distance) {
                HBondCandidate c;
                c.donor_res_id = res2.res_id();
                c.acceptor_res_id = res1.res_id();
                c.donor_atom = donor_name;
                c.acceptor_atom = acceptor_name;
                c.donor_pos = donor_atom.position();
                c.acceptor_pos = acceptor_atom.position();
                c.distance = dist;
                candidates.push_back(c);
            }
        }
    }

    return candidates;
}

void SlotOptimizer::score_alignment(
    HBondCandidate& candidate,
    const std::vector<HSlot>& h_slots,
    const std::vector<LPSlot>& lp_slots) const {

    if (h_slots.empty() || lp_slots.empty()) {
        candidate.alignment_score = 0.0;
        return;
    }

    geometry::Vector3D d2a = candidate.direction();
    geometry::Vector3D a2d = -d2a;

    double best_score = -999.0;
    int best_h = 0;
    int best_lp = 0;

    for (size_t hi = 0; hi < h_slots.size(); ++hi) {
        double h_align = h_slots[hi].direction().dot(d2a);

        for (size_t li = 0; li < lp_slots.size(); ++li) {
            double lp_align = lp_slots[li].direction().dot(a2d);
            double total = h_align + lp_align;

            if (total > best_score) {
                best_score = total;
                best_h = static_cast<int>(hi);
                best_lp = static_cast<int>(li);
            }
        }
    }

    candidate.h_slot_idx = best_h;
    candidate.lp_slot_idx = best_lp;
    candidate.alignment_score = best_score;
}

bool SlotOptimizer::try_alternative_slots(
    HBondCandidate& candidate,
    std::vector<HSlot>& h_slots,
    std::vector<LPSlot>& lp_slots) const {

    geometry::Vector3D d2a = candidate.direction();
    geometry::Vector3D a2d = -d2a;

    for (size_t hi = 0; hi < h_slots.size(); ++hi) {
        if (!h_slots[hi].can_add_bond(d2a, params_.min_bifurcation_angle)) {
            continue;
        }

        for (size_t li = 0; li < lp_slots.size(); ++li) {
            if (!lp_slots[li].can_add_bond(a2d, params_.min_bifurcation_angle)) {
                continue;
            }

            // Compute alignment for this pair
            double h_align = h_slots[hi].direction().dot(d2a);
            double lp_align = lp_slots[li].direction().dot(a2d);
            double score = h_align + lp_align;

            // Check alignment threshold
            bool is_bifurcated = !h_slots[hi].is_available() || !lp_slots[li].is_available();
            double min_align = is_bifurcated ? params_.min_bifurcation_alignment : params_.min_alignment;

            if (score >= min_align) {
                candidate.h_slot_idx = static_cast<int>(hi);
                candidate.lp_slot_idx = static_cast<int>(li);
                candidate.alignment_score = score;
                return true;
            }
        }
    }

    return false;
}

std::vector<core::HBond> SlotOptimizer::select_optimal(
    std::vector<HBondCandidate>& candidates,
    SlotCache& cache1,
    SlotCache& cache2) {

    // Reset slots
    cache1.reset_slots();
    cache2.reset_slots();

    // Score alignments
    for (auto& c : candidates) {
        SlotCache& donor_cache = (c.donor_res_id == cache1.residue().res_id()) ? cache1 : cache2;
        SlotCache& acceptor_cache = (c.acceptor_res_id == cache1.residue().res_id()) ? cache1 : cache2;

        auto& h_slots = donor_cache.get_h_slots(c.donor_atom);
        auto& lp_slots = acceptor_cache.get_lp_slots(c.acceptor_atom);

        score_alignment(c, h_slots, lp_slots);
    }

    // Sort by quality (best first)
    std::sort(candidates.begin(), candidates.end(),
              [](const HBondCandidate& a, const HBondCandidate& b) {
                  return a.quality_score() > b.quality_score();
              });

    std::vector<core::HBond> selected;

    for (auto& c : candidates) {
        SlotCache& donor_cache = (c.donor_res_id == cache1.residue().res_id()) ? cache1 : cache2;
        SlotCache& acceptor_cache = (c.acceptor_res_id == cache1.residue().res_id()) ? cache1 : cache2;

        auto& h_slots = donor_cache.get_h_slots(c.donor_atom);
        auto& lp_slots = acceptor_cache.get_lp_slots(c.acceptor_atom);

        if (h_slots.empty() || lp_slots.empty()) {
            continue;
        }

        // Check if assigned slots are still available
        int hi = c.h_slot_idx;
        int li = c.lp_slot_idx;

        if (hi < 0 || li < 0 || hi >= static_cast<int>(h_slots.size()) ||
            li >= static_cast<int>(lp_slots.size())) {
            continue;
        }

        geometry::Vector3D d2a = c.direction();
        geometry::Vector3D a2d = -d2a;

        bool h_ok = h_slots[hi].can_add_bond(d2a, params_.min_bifurcation_angle);
        bool lp_ok = lp_slots[li].can_add_bond(a2d, params_.min_bifurcation_angle);

        if (!h_ok || !lp_ok) {
            // Try alternative slots
            if (!try_alternative_slots(c, h_slots, lp_slots)) {
                continue;
            }
            hi = c.h_slot_idx;
            li = c.lp_slot_idx;
        }

        // Check alignment threshold (skip for short distances)
        if (c.distance >= params_.short_distance_threshold) {
            bool is_bifurcated = !h_slots[hi].is_available() || !lp_slots[li].is_available();
            double min_align = is_bifurcated ? params_.min_bifurcation_alignment : params_.min_alignment;

            if (c.alignment_score < min_align) {
                continue;
            }
        }

        // Accept this H-bond
        h_slots[hi].add_bond(d2a);
        lp_slots[li].add_bond(a2d);

        selected.push_back(candidate_to_hbond(c));
    }

    return selected;
}

std::vector<core::HBond> SlotOptimizer::select_baseline(
    std::vector<HBondCandidate>& candidates,
    const core::Residue& res1,
    const core::Residue& res2) const {

    // Filter by distance and backbone
    std::vector<HBondCandidate> valid;
    for (const auto& c : candidates) {
        if (is_backbone_backbone(c.donor_atom, c.acceptor_atom)) {
            continue;
        }
        if (c.distance < params_.baseline_min_distance ||
            c.distance > params_.baseline_max_distance) {
            continue;
        }
        valid.push_back(c);
    }

    // Sort by distance
    std::sort(valid.begin(), valid.end(),
              [](const HBondCandidate& a, const HBondCandidate& b) {
                  return a.distance < b.distance;
              });

    // Track usage by (res_id, atom_name)
    std::unordered_map<std::string, int> donor_usage;
    std::unordered_map<std::string, int> acceptor_usage;

    auto make_key = [](const std::string& res_id, const std::string& atom) {
        return res_id + ":" + atom;
    };

    std::vector<core::HBond> selected;

    for (const auto& c : valid) {
        std::string donor_key = make_key(c.donor_res_id, c.donor_atom);
        std::string acceptor_key = make_key(c.acceptor_res_id, c.acceptor_atom);

        // Get capacities
        std::string donor_code = (c.donor_res_id == res1.res_id()) ? res1.name() : res2.name();
        std::string acceptor_code = (c.acceptor_res_id == res1.res_id()) ? res1.name() : res2.name();

        int donor_cap = AtomCapacity::get_donor_capacity(donor_code, c.donor_atom);
        int acceptor_cap = AtomCapacity::get_acceptor_capacity(acceptor_code, c.acceptor_atom);

        // Check capacity
        if (donor_usage[donor_key] >= donor_cap) {
            continue;
        }
        if (acceptor_usage[acceptor_key] >= acceptor_cap) {
            continue;
        }

        // Accept
        donor_usage[donor_key]++;
        acceptor_usage[acceptor_key]++;

        selected.push_back(candidate_to_hbond(c));
    }

    return selected;
}

core::HBond SlotOptimizer::candidate_to_hbond(const HBondCandidate& candidate) const {
    core::HBond hb;
    hb.donor_atom_name = candidate.donor_atom;
    hb.acceptor_atom_name = candidate.acceptor_atom;
    hb.distance = candidate.distance;
    // Note: Additional fields like donor_residue_idx would need to be set by caller
    return hb;
}

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
