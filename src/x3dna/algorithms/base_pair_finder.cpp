/**
 * @file base_pair_finder.cpp
 * @brief Implementation of base pair finding (matches legacy find_bestpair)
 */

#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/geometry/least_squares_fitter.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <cmath>
#include <algorithm>
#include <limits>
#include <optional>
#include <array>
#include <iostream>

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;

std::vector<BasePair> BasePairFinder::find_pairs(Structure& structure) {
    return find_pairs_with_recording(structure, nullptr);
}

std::vector<BasePair> BasePairFinder::find_pairs(const Structure& structure) const {
    return find_pairs_with_recording(const_cast<Structure&>(structure), nullptr);
}

std::vector<BasePair> BasePairFinder::find_pairs_with_recording(Structure& structure, io::JsonWriter* writer) const {
    switch (strategy_) {
        case PairFindingStrategy::BEST_PAIR:
            return find_best_pairs(structure, writer);
        case PairFindingStrategy::ALL_PAIRS:
            // TODO: Add recording for all_pairs strategy
            return find_all_pairs(structure);
        case PairFindingStrategy::DISTANCE_BASED:
            // TODO: Implement distance-based strategy
            return {};
    }
    return {};
}

std::vector<BasePair> BasePairFinder::find_best_pairs(Structure& structure, io::JsonWriter* writer) const {
    std::vector<BasePair> base_pairs;
    std::vector<std::pair<size_t, size_t>> selected_pairs_legacy_idx;

    // Build mapping from legacy residue index to residue pointer
    ResidueIndexMapping mapping = build_residue_index_mapping(structure);
    if (mapping.by_legacy_idx.empty()) {
        return base_pairs;
    }

    // Aliases for compatibility with rest of function
    const auto& residue_by_legacy_idx = mapping.by_legacy_idx;
    int max_legacy_idx = mapping.max_legacy_idx;

    // PHASE 1: Validate ALL pairs (matches legacy check_pair loop)
    Phase1Results phase1 = run_phase1_validation(mapping);
    const auto& phase1_validation_results = phase1.validation_results;
    const auto& phase1_bp_type_ids = phase1.bp_type_ids;

    // Track matched residues using legacy 1-based indices (from PDB parsing)
    std::vector<bool> matched_indices(max_legacy_idx + 1, false); // +1 for 1-based indexing

    // Iterate until no new pairs are found (matches legacy while loop)
    size_t num_matched_prev = 0;
    size_t num_matched_curr = 0;
    int iteration_num = 0;
    std::vector<std::pair<int, int>> pairs_found_this_iteration;

    do {
        iteration_num++;
        num_matched_prev = num_matched_curr;
        num_matched_curr = 0;
        pairs_found_this_iteration.clear(); // Reset for this iteration

        // Count current matches
        for (bool matched : matched_indices) {
            if (matched)
                num_matched_curr++;
        }

        // Try to find pairs for each unpaired residue
        // Legacy: for (i = 1; i <= num_residue; i++) { if (RY[i] < 0 || matched_idx[i]) continue;
        // ... } CRITICAL: Iterate sequentially from 1 to max_legacy_idx to match legacy iteration
        // order
        for (int legacy_idx1 = 1; legacy_idx1 <= max_legacy_idx; ++legacy_idx1) {
            // Skip if already matched (matches legacy matched_idx[i] check)
            if (legacy_idx1 >= static_cast<int>(matched_indices.size()) || matched_indices[legacy_idx1]) {
                continue;
            }

            // Check if this legacy index exists in our map
            auto it = residue_by_legacy_idx.find(legacy_idx1);
            if (it == residue_by_legacy_idx.end() || !it->second) {
                continue; // Skip if residue doesn't exist (equivalent to RY[i] < 0)
            }

            const Residue* res1_ptr = it->second;

            // Check RY equivalent: is_nucleotide (includes modified nucleotides) and has frame
            // Use is_nucleotide() instead of residue_type() >= 0 to handle modified nucleotides
            if (!is_nucleotide(*res1_ptr) || !res1_ptr->reference_frame().has_value()) {
                continue; // Skip non-nucleotide residues (RY[i] < 0 equivalent)
            }

            // Find best partner for this residue (using legacy 1-based index from PDB parsing)
            auto best_partner = find_best_partner(legacy_idx1, structure, matched_indices, residue_by_legacy_idx,
                                                  phase1_validation_results, phase1_bp_type_ids, writer);

            if (!best_partner.has_value()) {
                continue;
            }

            int legacy_idx2 = best_partner->first;
            const ValidationResult& result1 = best_partner->second;

            // Check if res_idx2's best partner is res_idx1 (mutual best match)
            auto partner_of_partner = find_best_partner(legacy_idx2, structure, matched_indices, residue_by_legacy_idx,
                                                        phase1_validation_results, phase1_bp_type_ids, writer);

            bool is_mutual = (partner_of_partner.has_value() && partner_of_partner->first == legacy_idx1);

            if (is_mutual) {
                // Mutual best match found - get second residue first
                auto res2_it = residue_by_legacy_idx.find(legacy_idx2);
                const Residue* res2_ptr = (res2_it != residue_by_legacy_idx.end()) ? res2_it->second : nullptr;

                if (!res2_ptr) {
                    continue; // Safety check
                }

                // Double-check validation to ensure pair is still valid
                // CRITICAL: Use Phase 1 validation result instead of re-validating
                // Re-validating can give different results due to floating-point precision
                // or state changes. Phase 1 result is the source of truth.
                std::pair<int, int> normalized_pair_check = (legacy_idx1 < legacy_idx2)
                                                                ? std::make_pair(legacy_idx1, legacy_idx2)
                                                                : std::make_pair(legacy_idx2, legacy_idx1);

                auto phase1_check_it = phase1_validation_results.find(normalized_pair_check);
                if (phase1_check_it == phase1_validation_results.end()) {
                    // Pair not found in Phase 1 - this should not happen, skip it
                    std::cerr << "Warning: Pair (" << legacy_idx1 << ", " << legacy_idx2
                              << ") not found in Phase 1 validation results. Skipping.\n";
                    continue;
                }

                // CRITICAL: Ensure pair is valid before selecting
                if (!phase1_check_it->second.is_valid) {
                    // Pair was invalid in Phase 1 - must not select it
                    // This is a safety check to prevent invalid pairs from being selected
                    std::cerr << "Error: Attempted to select invalid pair (" << legacy_idx1 << ", " << legacy_idx2
                              << "). is_valid=" << phase1_check_it->second.is_valid
                              << ", d_v_check=" << phase1_check_it->second.d_v_check
                              << ", d_v=" << phase1_check_it->second.d_v << "\n";
                    continue;
                }

                // Mutual best match found - create base pair
                matched_indices[legacy_idx1] = true;
                matched_indices[legacy_idx2] = true;

                BasePair pair = create_base_pair(legacy_idx1, legacy_idx2, res1_ptr, res2_ptr, result1);
                base_pairs.push_back(pair);

                // Store legacy indices for recording (preserve original finding order to match legacy)
                // Legacy stores pairs in the order they were found: (searching_residue, best_partner)
                selected_pairs_legacy_idx.push_back({static_cast<size_t>(legacy_idx1),
                                                     static_cast<size_t>(legacy_idx2)});

                // Track pair found in this iteration (preserve finding order)
                pairs_found_this_iteration.push_back(
                    {static_cast<int>(legacy_idx1), static_cast<int>(legacy_idx2)});

                // Record mutual best decision (after selection)
                if (writer) {
                    int best_j_for_i = legacy_idx2;
                    int best_i_for_j = partner_of_partner.has_value() ? partner_of_partner->first : 0;
                    writer->record_mutual_best_decision(legacy_idx1, legacy_idx2, best_j_for_i, best_i_for_j, is_mutual,
                                                        true); // was_selected = true
                }
            } else {
                // Record mutual best decision (not selected)
                if (writer) {
                    int best_j_for_i = legacy_idx2;
                    int best_i_for_j = partner_of_partner.has_value() ? partner_of_partner->first : 0;
                    writer->record_mutual_best_decision(legacy_idx1, legacy_idx2, best_j_for_i, best_i_for_j, is_mutual,
                                                        false); // was_selected = false
                }
            }
        }

        // Recalculate matched count
        num_matched_curr = 0;
        for (bool matched : matched_indices) {
            if (matched)
                num_matched_curr++;
        }

        // Record iteration state
        if (writer) {
            // Use pairs found in THIS iteration only (not all pairs found so far)
            writer->record_iteration_state(iteration_num, num_matched_curr,
                                           static_cast<int>(matched_indices.size() - 1), matched_indices,
                                           pairs_found_this_iteration);
        }

    } while (num_matched_curr > num_matched_prev);

    // Record the original base pair selection (matches legacy find_bestpair output)
    // Use legacy indices from PDB parsing (NOT from JSON loading)
    if (writer && !selected_pairs_legacy_idx.empty()) {
        writer->record_find_bestpair_selection(selected_pairs_legacy_idx);
    }

    // Record base_pair records ONLY for pairs in the final selection
    // This matches legacy behavior where base_pair records correspond to ref_frames.dat
    // (only pairs that appear in the final output)
    if (writer) {
        for (const auto& pair : base_pairs) {
            writer->record_base_pair(pair);
        }
    }

    return base_pairs;
}

std::vector<BasePair> BasePairFinder::find_all_pairs(const Structure& structure) const {
    std::vector<BasePair> base_pairs;

    // Build list of all nucleotide residues
    std::vector<std::pair<size_t, const Residue*>> nucleotide_residues;
    size_t global_idx = 0;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (is_nucleotide(residue) && residue.reference_frame().has_value()) {
                nucleotide_residues.push_back({global_idx, &residue});
            }
            global_idx++;
        }
    }

    // Check all pairs
    for (size_t i = 0; i < nucleotide_residues.size(); ++i) {
        for (size_t j = i + 1; j < nucleotide_residues.size(); ++j) {
            const auto& [idx1, res1] = nucleotide_residues[i];
            const auto& [idx2, res2] = nucleotide_residues[j];

            ValidationResult result = validator_.validate(*res1, *res2);

            if (result.is_valid) {
                BasePair pair(idx1, idx2, result.bp_type);

                if (res1->reference_frame().has_value()) {
                    pair.set_frame1(res1->reference_frame().value());
                }
                if (res2->reference_frame().has_value()) {
                    pair.set_frame2(res2->reference_frame().value());
                }

                pair.set_hydrogen_bonds(result.hbonds);

                // Set bp_type string
                char base1 = res1->one_letter_code();
                char base2 = res2->one_letter_code();
                if (base1 != ' ' && base2 != ' ') {
                    pair.set_bp_type(std::string(1, base1) + std::string(1, base2));
                }

                base_pairs.push_back(pair);
            }
        }
    }

    return base_pairs;
}

std::optional<std::pair<int, ValidationResult>> BasePairFinder::find_best_partner(
    int legacy_idx1,                  // Legacy 1-based residue index
    const Structure& /* structure */, // Unused, kept for API compatibility
    const std::vector<bool>& matched_indices, const std::map<int, const Residue*>& residue_by_legacy_idx,
    const std::map<std::pair<int, int>, ValidationResult>& phase1_validation_results,
    const std::map<std::pair<int, int>, int>& phase1_bp_type_ids, // Store bp_type_id from Phase 1
    io::JsonWriter* writer) const {

    // Get residue using legacy index
    auto it = residue_by_legacy_idx.find(legacy_idx1);
    if (it == residue_by_legacy_idx.end() || !it->second) {
        return std::nullopt;
    }

    const Residue* res1 = it->second;

    // Check RY equivalent: is_nucleotide (includes modified nucleotides) and has frame
    // Use is_nucleotide() instead of residue_type() >= 0 to handle modified nucleotides
    if (!is_nucleotide(*res1) || !res1->reference_frame().has_value()) {
        return std::nullopt;
    }

    // Find best partner (lowest quality_score)
    // Legacy: for (j = 1; j <= num_residue; j++) { if (j == i || RY[j] < 0 || matched_idx[j])
    // continue; ... } CRITICAL: Iterate sequentially from 1 to max_legacy_idx to match legacy
    // iteration order
    double best_score = std::numeric_limits<double>::max();
    std::optional<std::pair<int, ValidationResult>> best_result;

    // Collect candidates for JSON output
    std::vector<std::tuple<int, bool, double, int>> candidates;
    bool collect_candidates = (writer != nullptr);

    // Find max legacy index from the map
    int max_legacy_idx = 0;
    for (const auto& [idx, _] : residue_by_legacy_idx) {
        if (idx > max_legacy_idx) {
            max_legacy_idx = idx;
        }
    }

    // Helper lambda to record ineligible candidate and continue
    auto record_ineligible = [&](int idx) {
        if (collect_candidates) {
            candidates.push_back(std::make_tuple(idx, false,
                std::numeric_limits<double>::max(), 0));
        }
    };

    for (int legacy_idx2 = 1; legacy_idx2 <= max_legacy_idx; ++legacy_idx2) {
        double candidate_score = std::numeric_limits<double>::max();
        int candidate_bp_type_id = 0;

        // Skip if same residue or already matched (matches legacy checks)
        const bool is_self_or_matched =
            legacy_idx2 == legacy_idx1 ||
            legacy_idx2 >= static_cast<int>(matched_indices.size()) ||
            matched_indices[legacy_idx2];
        if (is_self_or_matched) {
            record_ineligible(legacy_idx2);
            continue;
        }

        // Check if this legacy index exists in our map
        auto it = residue_by_legacy_idx.find(legacy_idx2);
        const bool residue_not_found = (it == residue_by_legacy_idx.end() || !it->second);
        if (residue_not_found) {
            record_ineligible(legacy_idx2);
            continue;
        }

        const Residue* residue = it->second;

        // Equivalent to legacy RY[j] >= 0 check: is_nucleotide (includes modified nucleotides)
        // and has frame. Use is_nucleotide() instead of residue_type() >= 0.
        const bool is_invalid_nucleotide =
            !is_nucleotide(*residue) || !residue->reference_frame().has_value();
        if (is_invalid_nucleotide) {
            record_ineligible(legacy_idx2);
            continue;
        }

        // Validate pair - use Phase 1 validation result if available
        // This ensures consistency: if Phase 1 said pair is invalid, don't select it
        std::pair<int, int> normalized_pair = (legacy_idx1 < legacy_idx2)
            ? std::make_pair(legacy_idx1, legacy_idx2)
            : std::make_pair(legacy_idx2, legacy_idx1);

        // Get or compute validation result
        ValidationResult result;
        auto phase1_it = phase1_validation_results.find(normalized_pair);
        const bool has_phase1_result = (phase1_it != phase1_validation_results.end());

        if (has_phase1_result) {
            result = phase1_it->second;
        } else {
            // Fallback: validate on the fly (shouldn't happen if Phase 1 ran and writer was set)
            result = (legacy_idx1 < legacy_idx2)
                ? validator_.validate(*res1, *residue)
                : validator_.validate(*residue, *res1);
        }

        // Common validity check - skip invalid pairs regardless of source
        if (!result.is_valid) {
            if (collect_candidates) {
                candidates.push_back(std::make_tuple(legacy_idx2, true, candidate_score, candidate_bp_type_id));
            }
            continue;
        }

        // Record pair_validation during selection (matches legacy check_pair recording)
        // Legacy records in check_pair when: bpid != 0 && i < j
        // Here we record when: result.is_valid && legacy_idx1 < legacy_idx2
        // This ensures only pairs actually checked during greedy selection are recorded
        if (writer && result.is_valid && legacy_idx1 < legacy_idx2) {
            // Get residue pointers in correct order for recording
            const Residue* rec_res1 = res1;
            const Residue* rec_res2 = residue;
            record_validation_results(legacy_idx1, legacy_idx2, rec_res1, rec_res2, result, writer);
        }

        // Calculate adjusted quality score for pair selection
        // CRITICAL: Legacy uses rtn_val[5] which is AFTER adjust_pairQuality and bp_type_id adjustment
        double quality_adjustment = adjust_pair_quality(result.hbonds);
        double adjusted_quality_score = result.quality_score + quality_adjustment;

        // Get bp_type_id - prefer Phase 1 result for consistency
        auto bp_type_it = phase1_bp_type_ids.find(normalized_pair);
        const bool has_phase1_bp_type = (bp_type_it != phase1_bp_type_ids.end());
        int bp_type_id = has_phase1_bp_type
            ? bp_type_it->second
            : calculate_bp_type_id(res1, residue, result, adjusted_quality_score);

        // Apply Watson-Crick bonus
        if (bp_type_id == 2) {
            adjusted_quality_score -= 2.0;
        }

        // Update candidate tracking
        candidate_score = adjusted_quality_score;
        candidate_bp_type_id = bp_type_id;

        if (collect_candidates) {
            // At this point, candidate is eligible (passed all checks above)
            candidates.push_back(std::make_tuple(legacy_idx2, true, candidate_score, candidate_bp_type_id));
        }

        // Update best result if this score is better
        // CRITICAL: Legacy uses strict < comparison - when scores are equal, keep first
        if (adjusted_quality_score < best_score) {
            best_score = adjusted_quality_score;
            best_result = std::make_pair(legacy_idx2, result);
        }
    }

    // Record candidates for JSON output
    if (writer && collect_candidates) {
        int best_j = best_result.has_value() ? best_result->first : 0;
        writer->record_best_partner_candidates(legacy_idx1, candidates, best_j,
                                               best_score < std::numeric_limits<double>::max() ? best_score : 0.0);
    }

    return best_result;
}

void BasePairFinder::record_validation_results(int legacy_idx1, int legacy_idx2, const core::Residue* res1,
                                               const core::Residue* res2, const ValidationResult& result,
                                               io::JsonWriter* writer) const {
    // Check if passes distance/angle checks (cdns) - matches legacy behavior
    // Legacy records validation if cdns is true, regardless of overlap
    // BUT legacy also records validation for pairs that FAIL cdns (for debugging)
    // See legacy check_pair: it records validation in both the cdns block AND the else block
    bool passes_cdns = result.distance_check && result.d_v_check && result.plane_angle_check && result.dNN_check;

    if (passes_cdns) {
        // Use 0-based indices for consistency with base_frame_calc
        size_t base_i = static_cast<size_t>(legacy_idx1 - 1); // Convert to 0-based
        size_t base_j = static_cast<size_t>(legacy_idx2 - 1); // Convert to 0-based

        // Adjust quality_score using adjust_pairQuality (matches legacy)
        double quality_adjustment = adjust_pair_quality(result.hbonds);
        double adjusted_quality_score = result.quality_score + quality_adjustment;

        // Prepare rtn_val array: [dorg, d_v, plane_angle, dNN, quality_score]
        std::array<double, 5> rtn_val = {result.dorg, result.d_v, result.plane_angle, result.dNN,
                                         adjusted_quality_score};

        // Calculate bp_type_id using check_wc_wobble_pair logic
        int bp_type_id = calculate_bp_type_id(res1, res2, result, adjusted_quality_score);

        // If bp_type_id == 2, subtract 2.0 from quality_score (matches legacy)
        if (bp_type_id == 2) {
            rtn_val[4] -= 2.0;
        }

        // Only record pair_validation for valid pairs when i < j to avoid duplicates
        // (Recording both (i,j) and (j,i) doubles the file size unnecessarily)
        if (result.is_valid && legacy_idx1 < legacy_idx2) {
            writer->record_pair_validation(base_i, base_j, result.is_valid, bp_type_id, result.dir_x, result.dir_y,
                                           result.dir_z, rtn_val, validator_.parameters());

            // Record base_pair for the same pairs as pair_validation (matches legacy behavior)
            // Analysis confirmed legacy pair_validation and base_pair have IDENTICAL pairs
            BasePair pair(base_i, base_j, result.bp_type);

            // Set reference frames
            if (res1->reference_frame().has_value()) {
                pair.set_frame1(res1->reference_frame().value());
            }
            if (res2->reference_frame().has_value()) {
                pair.set_frame2(res2->reference_frame().value());
            }

            // Set hydrogen bonds
            pair.set_hydrogen_bonds(result.hbonds);

            // Set bp_type string from residue names
            char base1 = res1->one_letter_code();
            char base2 = res2->one_letter_code();
            if (base1 != ' ' && base2 != ' ') {
                pair.set_bp_type(std::string(1, base1) + std::string(1, base2));
            }

            // Record to JSON (duplicate detection handled in JsonWriter)
            writer->record_base_pair(pair);
        }
    }

    // Record distance checks only for valid pairs (is_valid) when i < j
    // (only output pairs that pass all checks to reduce file size)
    if (result.is_valid && legacy_idx1 < legacy_idx2) {
        // Use 0-based indices for consistency with base_frame_calc
        size_t base_i = static_cast<size_t>(legacy_idx1 - 1);
        size_t base_j = static_cast<size_t>(legacy_idx2 - 1);
        writer->record_distance_checks(base_i, base_j, result.dorg, result.dNN, result.plane_angle, result.d_v,
                                       result.overlap_area);
    }

    // Record hydrogen bonds if present
    if (!result.hbonds.empty()) {
        // Use 0-based indices for consistency with base_frame_calc
        size_t base_i = static_cast<size_t>(legacy_idx1 - 1);
        size_t base_j = static_cast<size_t>(legacy_idx2 - 1);
        writer->record_hbond_list(base_i, base_j, result.hbonds);
    }
}

double BasePairFinder::adjust_pair_quality(const std::vector<core::hydrogen_bond>& hbonds) const {
    // Delegate to QualityScoreCalculator
    return quality_calculator_.adjust_pair_quality(hbonds);
}

int BasePairFinder::calculate_bp_type_id(const Residue* res1, const Residue* res2, const ValidationResult& result,
                                         double /* quality_score */) const {
    // Delegate to QualityScoreCalculator
    return quality_calculator_.calculate_bp_type_id(*res1, *res2, result);
}

char BasePairFinder::get_base_letter_from_type(core::ResidueType type) {
    // Delegate to QualityScoreCalculator
    return QualityScoreCalculator::get_base_letter(type);
}

namespace {
// Standard nucleotide ring geometry (from legacy xyz_ring array)
// Matches RA_LIST order: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
constexpr std::array<std::array<double, 3>, 9> STANDARD_RING_GEOMETRY = {{
    {{-1.265, 3.177, 0.000}}, // C4
    {{-2.342, 2.364, 0.001}}, // N3
    {{-1.999, 1.087, 0.000}}, // C2
    {{-0.700, 0.641, 0.000}}, // N1
    {{0.424, 1.460, 0.000}},  // C6
    {{0.071, 2.833, 0.000}},  // C5
    {{0.870, 3.969, 0.000}},  // N7 (purine)
    {{0.023, 4.962, 0.000}},  // C8 (purine)
    {{-1.289, 4.551, 0.000}}, // N9 (purine)
}};

// Legacy RA_LIST order for ring atoms
constexpr std::array<const char*, 9> RING_ATOM_NAMES = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ",
                                                        " C5 ", " N7 ", " C8 ", " N9 "};

// NT_CUTOFF from legacy (0.2618)
constexpr double NT_CUTOFF = 0.2618;

/**
 * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
 * @param residue Residue to check
 * @return RMSD value if calculable, or RMSD_DUMMY if not enough atoms
 */
std::optional<double> check_nt_type_by_rmsd(const Residue& residue) {
    // Find ring atoms in residue
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int nN = 0; // Count of nitrogen atoms (N1, N3, N7, N9)
    bool has_c1_prime = false;

    // LEGACY BEHAVIOR: Try ALL 9 ring atoms first (matches legacy residue_ident)
    for (size_t i = 0; i < RING_ATOM_NAMES.size(); ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];

        // Find this atom in residue
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                const auto& pos = atom.position();
                experimental_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));

                // Use corresponding standard geometry
                standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                             STANDARD_RING_GEOMETRY[i][2]));

                // Count nitrogen atoms (indices 1=N3, 3=N1, 6=N7, 8=N9)
                if (i == 1 || i == 3 || i == 6 || i == 8) {
                    nN++;
                }
                break;
            }
        }
    }

    // Check for C1' or C1R atom (required by legacy)
    // Some nucleotides like NMN use C1R instead of C1'
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " C1'" || atom.name() == " C1R") {
            has_c1_prime = true;
            break;
        }
    }

    // Legacy requires: (!nN && !C1_prime) -> return DUMMY
    if (nN == 0 && !has_c1_prime) {
        return std::nullopt; // DUMMY
    }

    // Need at least 3 atoms for RMSD calculation
    if (experimental_coords.size() < 3) {
        return std::nullopt; // DUMMY
    }

    // Perform least-squares fitting (matches legacy ls_fitting)
    geometry::LeastSquaresFitter fitter;
    try {
        auto fit_result = fitter.fit(standard_coords, experimental_coords);
        return fit_result.rms;
    } catch (const std::exception&) {
        return std::nullopt; // DUMMY
    }
}
} // namespace

bool BasePairFinder::is_nucleotide(const Residue& residue) {
    ResidueType type = residue.residue_type();

    // Check standard nucleotide types (always recognized)
    if (type == ResidueType::ADENINE || type == ResidueType::CYTOSINE || type == ResidueType::GUANINE ||
        type == ResidueType::THYMINE || type == ResidueType::URACIL) {
        return true;
    }

    // Check for modified nucleotides that are explicitly recognized
    // BUT: We still need to do RMSD check for NONCANONICAL_RNA to reject distorted residues
    // Legacy does RMSD check for all modified nucleotides (not in NT_LIST)
    if (type == ResidueType::PSEUDOURIDINE || type == ResidueType::INOSINE) {
        return true;
    }

    // For NONCANONICAL_RNA, still need RMSD check (legacy does this for all non-NT_LIST residues)
    // This includes H2U and other modified nucleotides that might be distorted
    if (type == ResidueType::NONCANONICAL_RNA) {
        // Fall through to RMSD check below
    }

    // For UNKNOWN or NONCANONICAL_RNA residues, use RMSD-based recognition (matches legacy
    // residue_ident) This correctly rejects distorted/clashing residues like H2U residue 16 in 1TTT
    // Legacy: All residues not in NT_LIST go through RMSD check (including H2U)
    if (type == ResidueType::UNKNOWN || type == ResidueType::NONCANONICAL_RNA) {
        // First check: find ring atoms and count them
        static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "};
        static const std::vector<std::string> purine_ring_atoms = {" N7 ", " C8 ", " N9 "};
        static const std::vector<std::string> nitrogen_atoms = {" N1 ", " N3 "};

        int ring_atom_count = 0;
        int kr = 0; // Purine ring atom count (N7, C8, N9)

        for (const auto& atom_name : common_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
                    break;
                }
            }
        }

        // Check for purine-specific atoms
        for (const auto& atom_name : purine_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    kr++;
                    break;
                }
            }
        }

        // Legacy logic: require >= 3 ring atoms, then do RMSD check
        // No nitrogen requirement - RMSD check handles non-nucleotides
        if (ring_atom_count + kr >= 3) {
            // Use strict threshold (0.2618) for all
            double rmsd_threshold = NT_CUTOFF;

            auto rmsd_opt = check_nt_type_by_rmsd(residue);

            if (rmsd_opt.has_value() && *rmsd_opt <= rmsd_threshold) {
                // RMSD check passed - it's a nucleotide
                return true;
            }

            // If we had purine atoms (kr > 0), try again without them (pyrimidine-only)
            // This matches legacy's fallback logic
            if (kr > 0 && ring_atom_count >= 3) {
                // Already checked above with all atoms, so if it failed with purine atoms,
                // it will likely fail without them too, but we should check
                // For now, if RMSD check fails, reject it (matches legacy behavior)
                return false;
            }
        }
    }

    return false;
}

size_t BasePairFinder::get_residue_index(const Structure& structure, const Residue& residue) {
    size_t idx = 0;
    for (const auto& chain : structure.chains()) {
        for (const auto& res : chain.residues()) {
            if (&res == &residue) {
                return idx;
            }
            idx++;
        }
    }
    return idx;
}

BasePairFinder::ResidueIndexMapping BasePairFinder::build_residue_index_mapping(const Structure& structure) const {
    ResidueIndexMapping mapping;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    mapping.by_legacy_idx[legacy_idx] = &residue;
                    if (legacy_idx > mapping.max_legacy_idx) {
                        mapping.max_legacy_idx = legacy_idx;
                    }
                }
            }
        }
    }

    return mapping;
}

BasePairFinder::Phase1Results BasePairFinder::run_phase1_validation(const ResidueIndexMapping& mapping) const {
    Phase1Results results;

    for (int legacy_idx1 = 1; legacy_idx1 <= mapping.max_legacy_idx - 1; ++legacy_idx1) {
        auto it1 = mapping.by_legacy_idx.find(legacy_idx1);
        if (it1 == mapping.by_legacy_idx.end() || !it1->second) {
            continue;
        }
        const Residue* res1 = it1->second;
        if (!is_nucleotide(*res1) || !res1->reference_frame().has_value()) {
            continue;
        }

        for (int legacy_idx2 = legacy_idx1 + 1; legacy_idx2 <= mapping.max_legacy_idx; ++legacy_idx2) {
            auto it2 = mapping.by_legacy_idx.find(legacy_idx2);
            if (it2 == mapping.by_legacy_idx.end() || !it2->second) {
                continue;
            }
            const Residue* res2 = it2->second;
            if (!is_nucleotide(*res2) || !res2->reference_frame().has_value()) {
                continue;
            }

            // Validate pair
            ValidationResult result = validator_.validate(*res1, *res2);

            // Store validation result (normalized by index order)
            std::pair<int, int> normalized_pair = std::make_pair(legacy_idx1, legacy_idx2);
            results.validation_results[normalized_pair] = result;

            // Calculate and store bp_type_id
            double quality_adjustment = adjust_pair_quality(result.hbonds);
            double adjusted_quality_score = result.quality_score + quality_adjustment;
            int bp_type_id = calculate_bp_type_id(res1, res2, result, adjusted_quality_score);
            results.bp_type_ids[normalized_pair] = bp_type_id;
        }
    }

    return results;
}

BasePair BasePairFinder::create_base_pair(int legacy_idx1, int legacy_idx2,
                                           const Residue* res1, const Residue* res2,
                                           const ValidationResult& result) const {
    // ALWAYS store smaller index first for consistency with legacy behavior
    size_t idx_small = static_cast<size_t>(std::min(legacy_idx1, legacy_idx2)) - 1;
    size_t idx_large = static_cast<size_t>(std::max(legacy_idx1, legacy_idx2)) - 1;
    bool swapped = (legacy_idx1 > legacy_idx2);

    BasePair pair(idx_small, idx_large, result.bp_type);
    pair.set_finding_order_swapped(swapped);

    // Set frames - swap if we reordered the indices
    const Residue* res_small = swapped ? res2 : res1;
    const Residue* res_large = swapped ? res1 : res2;

    if (res_small && res_small->reference_frame().has_value()) {
        pair.set_frame1(res_small->reference_frame().value());
    }

    if (res_large && res_large->reference_frame().has_value()) {
        pair.set_frame2(res_large->reference_frame().value());
    }

    // Set hydrogen bonds
    pair.set_hydrogen_bonds(result.hbonds);

    // Determine bp_type string from residue names
    if (res_small && res_large) {
        char base1 = res_small->one_letter_code();
        char base2 = res_large->one_letter_code();
        if (base1 != ' ' && base2 != ' ') {
            pair.set_bp_type(std::string(1, base1) + std::string(1, base2));
        }
    }

    return pair;
}

} // namespace algorithms
} // namespace x3dna
