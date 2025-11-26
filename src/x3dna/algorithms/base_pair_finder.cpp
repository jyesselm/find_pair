/**
 * @file base_pair_finder.cpp
 * @brief Implementation of base pair finding (matches legacy find_bestpair)
 */

#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/io/json_writer.hpp>
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

std::vector<BasePair> BasePairFinder::find_pairs_with_recording(Structure& structure,
                                                                io::JsonWriter* writer) const {
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

std::vector<BasePair> BasePairFinder::find_best_pairs(Structure& structure,
                                                       io::JsonWriter* writer) const {
    std::vector<BasePair> base_pairs;
    std::vector<std::pair<size_t, size_t>> selected_pairs_legacy_idx; // Store legacy indices for recording

    // Build mapping from legacy residue index to residue pointer
    // Use the legacy_residue_idx already stored on atoms (set during PDB parsing, NOT from JSON)
    std::map<int, const Residue*> residue_by_legacy_idx;
    int max_legacy_idx = 0;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            // Get legacy residue index from first atom in residue (set during PDB parsing)
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residue_by_legacy_idx[legacy_idx] = &residue;
                    if (legacy_idx > max_legacy_idx) {
                        max_legacy_idx = legacy_idx;
                    }
                }
            }
        }
    }

    if (residue_by_legacy_idx.empty()) {
        return base_pairs;
    }

    // PHASE 1: Validate ALL pairs and record base_pair for valid pairs (matches legacy check_pair loop)
    // Legacy: for (i = 1; i < num_residue; i++) { for (j = i + 1; j <= num_residue; j++) { check_pair(i, j, ...); } }
    // This validates all pairs BEFORE best partner selection
    // NOTE: Legacy uses i < num_residue, which means i goes from 1 to num_residue-1 (inclusive)
    // So we need legacy_idx1 to go from 1 to max_legacy_idx-1 (inclusive), which is < max_legacy_idx
    if (writer) {
        for (int legacy_idx1 = 1; legacy_idx1 <= max_legacy_idx - 1; ++legacy_idx1) {
            auto it1 = residue_by_legacy_idx.find(legacy_idx1);
            if (it1 == residue_by_legacy_idx.end() || !it1->second) {
                continue;
            }
            const Residue* res1 = it1->second;
            if (!is_nucleotide(*res1) || !res1->reference_frame().has_value()) {
                continue;
            }

            for (int legacy_idx2 = legacy_idx1 + 1; legacy_idx2 <= max_legacy_idx; ++legacy_idx2) {
                auto it2 = residue_by_legacy_idx.find(legacy_idx2);
                if (it2 == residue_by_legacy_idx.end() || !it2->second) {
                    continue;
                }
                const Residue* res2 = it2->second;
                if (!is_nucleotide(*res2) || !res2->reference_frame().has_value()) {
                    continue;
                }

                // Validate pair (matches legacy check_pair)
                ValidationResult result = validator_.validate(*res1, *res2);

                // Record validation and base_pair for all pairs that pass validation
                // (matches legacy check_pair -> calculate_more_bppars -> json_writer_record_base_pair)
                record_validation_results(legacy_idx1, legacy_idx2, res1, res2, result, writer);
            }
        }
    }

    // Track matched residues using legacy 1-based indices (from PDB parsing)
    std::vector<bool> matched_indices(max_legacy_idx + 1, false); // +1 for 1-based indexing

    // Iterate until no new pairs are found (matches legacy while loop)
    size_t num_matched_prev = 0;
    size_t num_matched_curr = 0;

    do {
        num_matched_prev = num_matched_curr;
        num_matched_curr = 0;

        // Count current matches
        for (bool matched : matched_indices) {
            if (matched)
                num_matched_curr++;
        }

        // Try to find pairs for each unpaired residue
        // Legacy: for (i = 1; i <= num_residue; i++) { if (RY[i] < 0 || matched_idx[i]) continue; ... }
        // CRITICAL: Iterate sequentially from 1 to max_legacy_idx to match legacy iteration order
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
            auto best_partner = find_best_partner(legacy_idx1, structure, matched_indices,
                                                  residue_by_legacy_idx, writer);

            if (!best_partner.has_value()) {
                continue;
            }

            int legacy_idx2 = best_partner->first;
            const ValidationResult& result1 = best_partner->second;

            // Check if res_idx2's best partner is res_idx1 (mutual best match)
            auto partner_of_partner = find_best_partner(legacy_idx2, structure, matched_indices,
                                                        residue_by_legacy_idx, writer);

            if (partner_of_partner.has_value() && partner_of_partner->first == legacy_idx1) {
                // Mutual best match found - create base pair
                matched_indices[legacy_idx1] = true;
                matched_indices[legacy_idx2] = true;

                // Get second residue using legacy index
                auto res2_it = residue_by_legacy_idx.find(legacy_idx2);
                const Residue* res2_ptr = (res2_it != residue_by_legacy_idx.end()) ? res2_it->second : nullptr;

                // Create BasePair object (using 0-based indices for BasePair)
                // Convert from legacy 1-based indices to 0-based by subtracting 1
                // This matches the method used in validation_pair creation
                BasePair pair(static_cast<size_t>(legacy_idx1) - 1, 
                             static_cast<size_t>(legacy_idx2) - 1, result1.bp_type);

                // Set frames if available
                if (res1_ptr->reference_frame().has_value()) {
                    pair.set_frame1(res1_ptr->reference_frame().value());
                }

                if (res2_ptr && res2_ptr->reference_frame().has_value()) {
                    pair.set_frame2(res2_ptr->reference_frame().value());
                }

                // Set hydrogen bonds
                pair.set_hydrogen_bonds(result1.hbonds);

                // Determine bp_type string from residue names
                std::string bp_type_str;
                if (res1_ptr && res2_ptr) {
                    char base1 = res1_ptr->one_letter_code();
                    char base2 = res2_ptr->one_letter_code();
                    if (base1 != ' ' && base2 != ' ') {
                        bp_type_str = std::string(1, base1) + std::string(1, base2);
                        pair.set_bp_type(bp_type_str);
                    }
                }

                base_pairs.push_back(pair);
                
                // Store legacy indices for recording (from PDB parsing, NOT from JSON)
                selected_pairs_legacy_idx.push_back({static_cast<size_t>(legacy_idx1), static_cast<size_t>(legacy_idx2)});
            }
        }

        // Recalculate matched count
        num_matched_curr = 0;
        for (bool matched : matched_indices) {
            if (matched)
                num_matched_curr++;
        }

    } while (num_matched_curr > num_matched_prev);

    // Record the original base pair selection (matches legacy find_bestpair output)
    // Use legacy indices from PDB parsing (NOT from JSON loading)
    if (writer && !selected_pairs_legacy_idx.empty()) {
        writer->record_find_bestpair_selection(selected_pairs_legacy_idx);
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
    const std::vector<bool>& matched_indices,
    const std::map<int, const Residue*>& residue_by_legacy_idx, io::JsonWriter* /* writer */) const {

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
    // Legacy: for (j = 1; j <= num_residue; j++) { if (j == i || RY[j] < 0 || matched_idx[j]) continue; ... }
    // CRITICAL: Iterate sequentially from 1 to max_legacy_idx to match legacy iteration order
    double best_score = std::numeric_limits<double>::max();
    std::optional<std::pair<int, ValidationResult>> best_result;

    // Find max legacy index from the map
    int max_legacy_idx = 0;
    for (const auto& [idx, _] : residue_by_legacy_idx) {
        if (idx > max_legacy_idx) {
            max_legacy_idx = idx;
        }
    }

    for (int legacy_idx2 = 1; legacy_idx2 <= max_legacy_idx; ++legacy_idx2) {
        // Skip if same residue or already matched (matches legacy checks)
        if (legacy_idx2 == legacy_idx1 || 
            legacy_idx2 >= static_cast<int>(matched_indices.size()) ||
            matched_indices[legacy_idx2]) {
            continue;
        }

        // Check if this legacy index exists in our map
        auto it = residue_by_legacy_idx.find(legacy_idx2);
        if (it == residue_by_legacy_idx.end() || !it->second) {
            continue; // Skip if residue doesn't exist (equivalent to RY[j] < 0)
        }

        const Residue* residue = it->second;

        // Equivalent to legacy RY[j] >= 0 check: is_nucleotide (includes modified nucleotides) and has frame
        // Use is_nucleotide() instead of residue_type() >= 0 to handle modified nucleotides
        if (!is_nucleotide(*residue) || !residue->reference_frame().has_value()) {
            continue; // Skip non-nucleotide residues (RY[j] < 0 equivalent)
        }

        // Validate pair
        ValidationResult result = validator_.validate(*res1, *residue);

        // NOTE: Validation results are now recorded in Phase 1 (validate_all_pairs)
        // We only use the result here for best partner selection

        // Check if valid and has better score
        // CRITICAL: Legacy uses rtn_val[5] which is AFTER adjust_pairQuality and bp_type_id adjustment
        // We need to calculate the adjusted quality_score here for pair selection
        if (result.is_valid) {
            // Apply adjust_pairQuality adjustment
            double quality_adjustment = adjust_pair_quality(result.hbonds);
            double adjusted_quality_score = result.quality_score + quality_adjustment;
            
            // Apply bp_type_id == 2 adjustment (if applicable)
            // We need to calculate bp_type_id first to know if we should subtract 2.0
            int bp_type_id = calculate_bp_type_id(res1, residue, result, adjusted_quality_score);
            if (bp_type_id == 2) {
                adjusted_quality_score -= 2.0;
            }
            
            // Use adjusted quality_score for pair selection (matches legacy rtn_val[5])
            if (adjusted_quality_score < best_score) {
                best_score = adjusted_quality_score;
                best_result = std::make_pair(legacy_idx2, result);
            }
        }
    }

    return best_result;
}

void BasePairFinder::record_validation_results(int legacy_idx1, int legacy_idx2,
                                                const core::Residue* res1, const core::Residue* res2,
                                                const ValidationResult& result,
                                                io::JsonWriter* writer) const {
    // Check if passes distance/angle checks (cdns) - matches legacy behavior
    // Legacy records validation if cdns is true, regardless of overlap
    // BUT legacy also records validation for pairs that FAIL cdns (for debugging)
    // See legacy check_pair: it records validation in both the cdns block AND the else block
    bool passes_cdns = result.distance_check && result.d_v_check &&
                       result.plane_angle_check && result.dNN_check;

    // Record pair validation for ALL pairs (matches legacy - records even for failed pairs)
    // Legacy: records validation in both cdns block (with is_valid) and else block (with is_valid=0)
    size_t base_i = static_cast<size_t>(legacy_idx1);
    size_t base_j = static_cast<size_t>(legacy_idx2);

    if (passes_cdns) {
        // Use legacy residue indices directly (already 1-based from atoms)
        size_t base_i = static_cast<size_t>(legacy_idx1);
        size_t base_j = static_cast<size_t>(legacy_idx2);

        // Adjust quality_score using adjust_pairQuality (matches legacy)
        double quality_adjustment = adjust_pair_quality(result.hbonds);
        double adjusted_quality_score = result.quality_score + quality_adjustment;

        // Prepare rtn_val array: [dorg, d_v, plane_angle, dNN, quality_score]
        std::array<double, 5> rtn_val = {result.dorg, result.d_v, result.plane_angle,
                                         result.dNN, adjusted_quality_score};

        // Calculate bp_type_id using check_wc_wobble_pair logic
        int bp_type_id = calculate_bp_type_id(res1, res2, result, adjusted_quality_score);

        // If bp_type_id == 2, subtract 2.0 from quality_score (matches legacy)
        if (bp_type_id == 2) {
            rtn_val[4] -= 2.0;
        }

        // Record pair validation (legacy records this for all pairs passing cdns+overlap)
        writer->record_pair_validation(base_i, base_j, result.is_valid, bp_type_id,
                                       result.dir_x, result.dir_y, result.dir_z, rtn_val,
                                       validator_.parameters());
    } else {
        // Legacy also records validation for pairs that FAIL cdns (for debugging)
        // See legacy check_pair else block: json_writer_record_pair_validation(i, j, 0, 0, ...)
        // We need to calculate rtn_val even for failed pairs
        std::array<double, 5> rtn_val = {result.dorg, result.d_v, result.plane_angle,
                                         result.dNN, result.quality_score};
        writer->record_pair_validation(base_i, base_j, 0, 0,
                                       result.dir_x, result.dir_y, result.dir_z, rtn_val,
                                       validator_.parameters());
    }

    // Record base_pair for ALL pairs that pass validation (matches legacy calculate_more_bppars)
    // Legacy records base_pair inside calculate_more_bppars for all pairs that pass validation
    if (result.is_valid) {
        // Create BasePair object with validation results
        BasePair validation_pair(static_cast<size_t>(legacy_idx1) - 1, 
                                static_cast<size_t>(legacy_idx2) - 1, result.bp_type);
        
        // Set frames if available
        if (res1->reference_frame().has_value()) {
            validation_pair.set_frame1(res1->reference_frame().value());
        }
        if (res2->reference_frame().has_value()) {
            validation_pair.set_frame2(res2->reference_frame().value());
        }
        
        // Set hydrogen bonds
        validation_pair.set_hydrogen_bonds(result.hbonds);
        
        // Set bp_type string
        char base1 = res1->one_letter_code();
        char base2 = res2->one_letter_code();
        if (base1 != ' ' && base2 != ' ') {
            validation_pair.set_bp_type(std::string(1, base1) + std::string(1, base2));
        }
        
        // Record base_pair (will assign indices automatically)
        writer->record_base_pair(validation_pair);
    }
    
    // Record distance checks only if also passes hydrogen bond check
    // (legacy only records distance_checks for pairs that pass hbond check)
    if (result.hbond_check) {
        size_t base_i = static_cast<size_t>(legacy_idx1);
        size_t base_j = static_cast<size_t>(legacy_idx2);
        writer->record_distance_checks(base_i, base_j, result.dorg, result.dNN,
                                       result.plane_angle, result.d_v,
                                       result.overlap_area);
    }

    // Record hydrogen bonds if present
    if (!result.hbonds.empty()) {
        size_t base_i = static_cast<size_t>(legacy_idx1);
        size_t base_j = static_cast<size_t>(legacy_idx2);
        writer->record_hbond_list(base_i, base_j, result.hbonds);
    }
}

double BasePairFinder::adjust_pair_quality(const std::vector<core::hydrogen_bond>& hbonds) const {
    // Match legacy adjust_pairQuality logic
    // Count good hydrogen bonds (distance in [2.5, 3.5] Angstroms)
    // Legacy: dval = num_list[k][3] / MFACTOR, then checks if in [2.5, 3.5]
    // Modern: hbonds already have distance calculated
    int num_good_hb = 0;
    for (const auto& hbond : hbonds) {
        // Skip non-standard hydrogen bonds (type != '-')
        if (hbond.type != '-') {
            continue;
        }
        // Check if distance is in good range [2.5, 3.5]
        if (hbond.distance >= 2.5 && hbond.distance <= 3.5) {
            num_good_hb++;
        }
    }
    
    // Legacy: if (num_good_hb >= 2) return -3.0; else return -num_good_hb;
    if (num_good_hb >= 2) {
        return -3.0;
    } else {
        return -static_cast<double>(num_good_hb);
    }
}

int BasePairFinder::calculate_bp_type_id(const Residue* res1, const Residue* res2,
                                         const ValidationResult& result,
                                         double /* quality_score */) const {
    // Match legacy check_wc_wobble_pair logic
    // Legacy: *bpid = -1 initially in calculate_more_bppars
    // Then: if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
    //           check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
    //           if (*bpid == 2) rtn_val[5] -= 2.0;
    //       }
    // check_wc_wobble_pair sets:
    //   *bpid = 1 if fabs(shear) in [1.8, 2.8]
    //   *bpid = 2 if fabs(shear) <= 1.8 && base pair matches WC_LIST
    
    // For now, use simplified version (full implementation requires bpstep_par from Stage 6)
    // Start with -1 (legacy initial value)
    int bp_type_id = -1;
    
    // Check direction vector conditions (legacy requirement)
    if (result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) {
        // Get base pair string
        char base1 = res1->one_letter_code();
        char base2 = res2->one_letter_code();
        std::string bp_str = std::string(1, std::toupper(base1)) + std::string(1, std::toupper(base2));
        
        // WC_LIST: "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
        // Simplified: if base pair matches WC_LIST and is valid, set to 2
        // Full implementation would require shear/stretch/opening from bpstep_par
        // For now, if base pair is in WC_LIST and direction vectors are correct, assume bp_type_id = 2
        // This is a heuristic - full accuracy requires bpstep_par (Stage 6)
        if (result.is_valid) {
            if (bp_str == "AT" || bp_str == "TA" || bp_str == "AU" || bp_str == "UA" ||
                bp_str == "GC" || bp_str == "CG" || bp_str == "IC" || bp_str == "CI") {
                // Base pair is in WC_LIST - set to 2 (Watson-Crick)
                // Note: Full implementation would check fabs(shear) <= 1.8
                bp_type_id = 2;
            } else if (result.bp_type == core::BasePairType::WOBBLE) {
                // Wobble pair (not in WC_LIST) - set to 1
                // Note: Full implementation would check fabs(shear) in [1.8, 2.8]
                bp_type_id = 1;
            }
        }
    }
    
    // Legacy behavior: bp_type_id = -1 is kept for valid pairs when:
    // 1. Direction vectors don't match (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0)
    // 2. OR check_wc_wobble_pair doesn't set it to 1 or 2
    // Only set to 0 for invalid pairs
    if (!result.is_valid) {
        bp_type_id = 0;
    }
    // DO NOT convert -1 to 0 for valid pairs - legacy keeps -1 in some cases
    // If bp_type_id is still -1 and is_valid is true, keep it as -1
    
    return bp_type_id;
}

bool BasePairFinder::is_nucleotide(const Residue& residue) {
    ResidueType type = residue.residue_type();
    
    // Check standard nucleotide types
    if (type == ResidueType::ADENINE || type == ResidueType::CYTOSINE ||
        type == ResidueType::GUANINE || type == ResidueType::THYMINE ||
        type == ResidueType::URACIL) {
        return true;
    }
    
    // Check for modified nucleotides (like XGR, XCR, XTR, XAR in 3KNC)
    // These have ResidueType::UNKNOWN but have ring atoms
    // Legacy code treats residues with ring atoms as nucleotides (RY >= 0)
    if (type == ResidueType::UNKNOWN) {
        // Check for ring atoms (similar to legacy residue_ident and generate_modern_json)
        static const std::vector<std::string> common_ring_atoms = {
            " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "};
        int ring_atom_count = 0;
        for (const auto& atom_name : common_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
                    break;
                }
            }
        }
        // If has >= 3 ring atoms, treat as nucleotide (matches legacy and generate_modern_json)
        if (ring_atom_count >= 3) {
            return true;
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

} // namespace algorithms
} // namespace x3dna
