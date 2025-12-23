/**
 * @file base_pair_finder.cpp
 * @brief Implementation of base pair finding (matches legacy find_bestpair)
 */

#include <x3dna/algorithms/pair_identification/base_pair_finder.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/typing/atom_type.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/geometry/least_squares_fitter.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <cmath>
#include <algorithm>
#include <limits>
#include <optional>
#include <array>
#include <iostream>
#include <chrono>
#include <iomanip>

// Timing support for profiling (controlled by environment variable or compile flag)
namespace {
    bool g_profile_pair_finding = false;

    class ScopedTimer {
    public:
        ScopedTimer(const char* name, bool& enabled) : name_(name), enabled_(enabled) {
            if (enabled_) start_ = std::chrono::high_resolution_clock::now();
        }
        ~ScopedTimer() {
            if (enabled_) {
                auto end = std::chrono::high_resolution_clock::now();
                auto ms = std::chrono::duration<double, std::milli>(end - start_).count();
                std::cout << "[PAIR_TIMING] " << std::setw(30) << std::left << name_
                          << std::fixed << std::setprecision(1) << ms << " ms\n";
            }
        }
    private:
        const char* name_;
        bool& enabled_;
        std::chrono::high_resolution_clock::time_point start_;
    };
}

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;
using core::AtomType;
using core::RING_ATOM_TYPES;
using core::NUM_RING_ATOM_TYPES;

// ============================================================================
// Helper methods - small, focused functions
// ============================================================================

bool BasePairFinder::can_participate_in_pairing(const Residue* res) {
    if (!res)
        return false;
    return is_nucleotide(*res) && res->reference_frame().has_value();
}

bool BasePairFinder::is_matched(int legacy_idx, const std::vector<bool>& matched) {
    if (legacy_idx < 0 || legacy_idx >= static_cast<int>(matched.size())) {
        return false;
    }
    return matched[legacy_idx];
}

double BasePairFinder::calculate_adjusted_score(const ValidationResult& result, int bp_type_id) const {
    double quality_adjustment = adjust_pair_quality(result.hbonds);
    double score = result.quality_score + quality_adjustment;

    // Watson-Crick pairs get a bonus (lower is better)
    if (bp_type_id == 2) {
        score -= validation_constants::WC_QUALITY_BONUS;
    }
    return score;
}

// ============================================================================
// Public interface
// ============================================================================

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
            return find_all_pairs(structure);
        case PairFindingStrategy::DISTANCE_BASED:
            return {};
    }
    return {};
}

bool BasePairFinder::try_select_mutual_pair(int legacy_idx1, int legacy_idx2, const Residue* res1, const Residue* res2,
                                            const ValidationResult& result, const PartnerSearchContext& ctx,
                                            PairSelectionState& state) const {
    // Verify pair is valid in Phase 1 results
    const auto* phase1_result = ctx.phase1.get_result(legacy_idx1, legacy_idx2);
    if (!phase1_result) {
        std::cerr << "Warning: Pair (" << legacy_idx1 << ", " << legacy_idx2
                  << ") not found in Phase 1 validation results. Skipping.\n";
        return false;
    }

    if (!phase1_result->is_valid) {
        std::cerr << "Error: Attempted to select invalid pair (" << legacy_idx1 << ", " << legacy_idx2
                  << "). is_valid=" << phase1_result->is_valid << ", d_v_check=" << phase1_result->d_v_check
                  << ", d_v=" << phase1_result->d_v << "\n";
        return false;
    }

    // Create and store the pair
    state.mark_matched(legacy_idx1, legacy_idx2);
    state.base_pairs.push_back(create_base_pair(legacy_idx1, legacy_idx2, res1, res2, result));
    state.selected_pairs_legacy_idx.push_back({static_cast<size_t>(legacy_idx1), static_cast<size_t>(legacy_idx2)});
    state.pairs_found_this_iteration.push_back({legacy_idx1, legacy_idx2});

    return true;
}

std::vector<BasePair> BasePairFinder::find_best_pairs(Structure& structure, io::JsonWriter* writer) const {
    // Check for profiling environment variable
    static bool profile_checked = false;
    if (!profile_checked) {
        if (const char* env = std::getenv("X3DNA_PROFILE_PAIRS")) {
            g_profile_pair_finding = (std::string(env) == "1");
        }
        profile_checked = true;
    }

    ResidueIndexMapping mapping = [&]() {
        ScopedTimer t("Build residue mapping", g_profile_pair_finding);
        return build_residue_index_mapping(structure);
    }();

    if (mapping.empty()) {
        return {};
    }

    if (g_profile_pair_finding) {
        std::cout << "[PAIR_TIMING] Nucleotide count: " << mapping.by_legacy_idx.size()
                  << ", max_legacy_idx: " << mapping.max_legacy_idx << "\n";
    }

    Phase1Results phase1 = [&]() {
        ScopedTimer t("Phase 1 validation", g_profile_pair_finding);
        return run_phase1_validation(mapping);
    }();

    if (g_profile_pair_finding) {
        std::cout << "[PAIR_TIMING] Phase 1 pairs validated: " << phase1.validation_results.size() << "\n";
    }

    PairSelectionState state(mapping.max_legacy_idx);
    PartnerSearchContext ctx{state.matched_indices, mapping, phase1, writer};

    int iteration_num = 0;
    size_t prev_matched = 0;

    auto iteration_start = std::chrono::high_resolution_clock::now();

    // Iterate until no new pairs found
    do {
        iteration_num++;
        prev_matched = state.count_matched();
        state.pairs_found_this_iteration.clear();

        for (int idx1 = 1; idx1 <= mapping.max_legacy_idx; ++idx1) {
            if (is_matched(idx1, state.matched_indices))
                continue;

            const Residue* res1 = mapping.get(idx1);
            if (!can_participate_in_pairing(res1))
                continue;

            auto best = find_best_partner(idx1, ctx);
            if (!best.has_value())
                continue;

            int idx2 = best->first;
            const ValidationResult& result = best->second;

            // Check for mutual best match
            auto reverse = find_best_partner(idx2, ctx);
            const bool is_mutual = reverse.has_value() && reverse->first == idx1;

            if (is_mutual) {
                const Residue* res2 = mapping.get(idx2);
                if (res2) {
                    (void)try_select_mutual_pair(idx1, idx2, res1, res2, result, ctx, state);
                }
            }

            // Record decision for JSON output
            if (writer) {
                int best_j_for_i = idx2;
                int best_i_for_j = reverse.has_value() ? reverse->first : 0;
                writer->record_mutual_best_decision(idx1, idx2, best_j_for_i, best_i_for_j, is_mutual, is_mutual);
            }
        }

        if (writer) {
            writer->record_iteration_state(iteration_num, static_cast<int>(state.count_matched()),
                                           mapping.max_legacy_idx, state.matched_indices,
                                           state.pairs_found_this_iteration);
        }
    } while (state.count_matched() > prev_matched);

    if (g_profile_pair_finding) {
        auto iteration_end = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration<double, std::milli>(iteration_end - iteration_start).count();
        std::cout << "[PAIR_TIMING] Mutual best matching      " << std::fixed << std::setprecision(1) << ms << " ms"
                  << " (" << iteration_num << " iterations, " << state.base_pairs.size() << " pairs found)\n";
    }

    // Record final results
    if (writer) {
        if (!state.selected_pairs_legacy_idx.empty()) {
            writer->record_find_bestpair_selection(state.selected_pairs_legacy_idx);
        }
        for (const auto& pair : state.base_pairs) {
            writer->record_base_pair(pair);
        }
    }

    return state.base_pairs;
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
                // Validation already ensures both residues have frames, so we can access them directly
                // Store by value since reference_frame().value() returns a temporary
                auto frame1 = res1->reference_frame().value();
                auto frame2 = res2->reference_frame().value();
                BasePair pair(idx1, idx2, frame1, frame2, result.bp_type);

                // Set unique residue identifiers
                pair.set_res_id1(res1->res_id());
                pair.set_res_id2(res2->res_id());

                pair.set_hydrogen_bonds(result.hbonds);

                // Set bp_type string
                char base1 = core::one_letter_code(*res1);
                char base2 = core::one_letter_code(*res2);
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
    int legacy_idx1, const PartnerSearchContext& ctx) const {

    const Residue* res1 = ctx.mapping.get(legacy_idx1);
    if (!can_participate_in_pairing(res1)) {
        return std::nullopt;
    }

    // Early rejection threshold (squared to avoid sqrt)
    const double max_origin_distance_sq = validator_.parameters().max_dorg * validator_.parameters().max_dorg;
    const auto& origin1 = res1->reference_frame()->origin();

    double best_score = std::numeric_limits<double>::max();
    std::optional<std::pair<int, ValidationResult>> best_result;
    std::vector<std::tuple<int, bool, double, int>> candidates;
    const bool collect = (ctx.writer != nullptr);

    for (int idx2 = 1; idx2 <= ctx.mapping.max_legacy_idx; ++idx2) {
        // Skip self or already matched
        if (idx2 == legacy_idx1 || is_matched(idx2, ctx.matched_indices)) {
            if (collect)
                candidates.emplace_back(idx2, false, std::numeric_limits<double>::max(), 0);
            continue;
        }

        const Residue* res2 = ctx.mapping.get(idx2);
        if (!can_participate_in_pairing(res2)) {
            if (collect)
                candidates.emplace_back(idx2, false, std::numeric_limits<double>::max(), 0);
            continue;
        }

        // Early distance rejection - skip pairs that are too far apart
        const auto& origin2 = res2->reference_frame()->origin();
        double dx = origin2.x() - origin1.x();
        double dy = origin2.y() - origin1.y();
        double dz = origin2.z() - origin1.z();
        double dist_sq = dx * dx + dy * dy + dz * dz;

        if (dist_sq > max_origin_distance_sq) {
            if (collect)
                candidates.emplace_back(idx2, false, std::numeric_limits<double>::max(), 0);
            continue;  // Skip - too far apart
        }

        // Get validation result from Phase 1 or compute on the fly
        const auto* phase1_result = ctx.phase1.get_result(legacy_idx1, idx2);
        ValidationResult fallback_result; // Stays in scope for the loop iteration

        if (!phase1_result) {
            fallback_result = (legacy_idx1 < idx2) ? validator_.validate(*res1, *res2)
                                                   : validator_.validate(*res2, *res1);
            phase1_result = &fallback_result;
        }

        const ValidationResult& result = *phase1_result;
        if (!result.is_valid) {
            if (collect)
                candidates.emplace_back(idx2, true, std::numeric_limits<double>::max(), 0);
            continue;
        }

        // Record validation for JSON output
        if (ctx.writer && legacy_idx1 < idx2) {
            record_validation_results(legacy_idx1, idx2, res1, res2, result, ctx.writer);
        }

        // Calculate score
        int bp_type_id = ctx.phase1.get_bp_type_id(legacy_idx1, idx2);
        if (bp_type_id == 0) {
            double adj_score = result.quality_score + adjust_pair_quality(result.hbonds);
            bp_type_id = calculate_bp_type_id(res1, res2, result, adj_score);
        }
        double score = calculate_adjusted_score(result, bp_type_id);

        if (collect) {
            candidates.emplace_back(idx2, true, score, bp_type_id);
        }

        if (score < best_score) {
            best_score = score;
            best_result = std::make_pair(idx2, result);
        }
    }

    if (ctx.writer && collect) {
        int best_j = best_result.has_value() ? best_result->first : 0;
        double final_score = (best_score < std::numeric_limits<double>::max()) ? best_score : 0.0;
        ctx.writer->record_best_partner_candidates(legacy_idx1, candidates, best_j, final_score);
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

        // If bp_type_id == 2, subtract WC bonus from quality_score (matches legacy)
        if (bp_type_id == 2) {
            rtn_val[4] -= validation_constants::WC_QUALITY_BONUS;
        }

        // Only record pair_validation for valid pairs when i < j to avoid duplicates
        // (Recording both (i,j) and (j,i) doubles the file size unnecessarily)
        if (result.is_valid && legacy_idx1 < legacy_idx2) {
            writer->record_pair_validation(base_i, base_j, result.is_valid, bp_type_id, result.dir_x, result.dir_y,
                                           result.dir_z, rtn_val, validator_.parameters(),
                                           res1->res_id(), res2->res_id());

            // Record base_pair for the same pairs as pair_validation (matches legacy behavior)
            // Analysis confirmed legacy pair_validation and base_pair have IDENTICAL pairs
            // Validation already ensures both residues have frames
            // Store by value since reference_frame().value() returns a temporary
            auto frame1 = res1->reference_frame().value();
            auto frame2 = res2->reference_frame().value();
            BasePair pair(base_i, base_j, frame1, frame2, result.bp_type);

            // Set unique residue identifiers
            pair.set_res_id1(res1->res_id());
            pair.set_res_id2(res2->res_id());

            // Set hydrogen bonds
            pair.set_hydrogen_bonds(result.hbonds);

            // Set bp_type string from residue names
            char base1 = core::one_letter_code(*res1);
            char base2 = core::one_letter_code(*res2);
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
                                       result.overlap_area, res1->res_id(), res2->res_id());
    }

    // Record hydrogen bonds if present
    if (!result.hbonds.empty()) {
        // Use 0-based indices for consistency with base_frame_calc
        size_t base_i = static_cast<size_t>(legacy_idx1 - 1);
        size_t base_j = static_cast<size_t>(legacy_idx2 - 1);
        writer->record_hbond_list(base_i, base_j, result.hbonds, res1->res_id(), res2->res_id());
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

// Legacy RA_LIST order for ring atoms (using trimmed names)
// Ring atom names array removed - now using RING_ATOM_TYPES from atom_type.hpp

// Use constant from validation_constants.hpp
using validation_constants::NT_RMSD_CUTOFF;

/**
 * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
 * @param residue Residue to check
 * @return RMSD value if calculable, or RMSD_DUMMY if not enough atoms
 */
std::optional<double> check_nt_type_by_rmsd(const Residue& residue) {
    // Find ring atoms in residue using AtomType for O(1) comparison
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int nN = 0; // Count of nitrogen atoms (N1, N3, N7, N9)
    bool has_c1_prime = false;

    // LEGACY BEHAVIOR: Try ALL 9 ring atoms first (matches legacy residue_ident)
    for (size_t i = 0; i < NUM_RING_ATOM_TYPES; ++i) {
        AtomType target_type = RING_ATOM_TYPES[i];
        const Atom* atom = residue.find_atom_by_type(target_type);
        if (atom != nullptr) {
            const auto& pos = atom->position();
            experimental_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));

            // Use corresponding standard geometry
            standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                         STANDARD_RING_GEOMETRY[i][2]));

            // Count nitrogen atoms (indices 1=N3, 3=N1, 6=N7, 8=N9)
            if (i == 1 || i == 3 || i == 6 || i == 8) {
                nN++;
            }
        }
    }

    // Check for C1' using AtomType, or C1R via string (alternative name)
    if (residue.has_atom_type(AtomType::C1_PRIME)) {
        has_c1_prime = true;
    } else {
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == "C1R") {
                has_c1_prime = true;
                break;
            }
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

// Additional helpers for is_nucleotide using AtomType

// Common ring atom types (pyrimidine ring)
constexpr std::array<AtomType, 6> COMMON_RING_ATOM_TYPES = {
    AtomType::C4, AtomType::N3, AtomType::C2, AtomType::N1, AtomType::C6, AtomType::C5
};
// Purine-only ring atom types
constexpr std::array<AtomType, 3> PURINE_RING_ATOM_TYPES = {AtomType::N7, AtomType::C8, AtomType::N9};

bool is_standard_nucleotide(core::typing::BaseType type) {
    return type == core::typing::BaseType::ADENINE || type == core::typing::BaseType::CYTOSINE ||
           type == core::typing::BaseType::GUANINE || type == core::typing::BaseType::THYMINE ||
           type == core::typing::BaseType::URACIL;
}

bool is_recognized_modified_nucleotide(core::typing::BaseType type) {
    return type == core::typing::BaseType::PSEUDOURIDINE || type == core::typing::BaseType::INOSINE;
}

bool needs_rmsd_validation(const Residue& residue) {
    auto mol_type = residue.molecule_type();
    auto base_type = residue.base_type();
    return mol_type == core::typing::MoleculeType::UNKNOWN ||
           (mol_type == core::typing::MoleculeType::NUCLEIC_ACID && base_type == core::typing::BaseType::UNKNOWN);
}

// Count matching atoms using AtomType array
template <size_t N> int count_matching_atom_types(const Residue& residue, const std::array<AtomType, N>& atom_types) {
    int count = 0;
    for (auto type : atom_types) {
        if (residue.has_atom_type(type))
            count++;
    }
    return count;
}

bool passes_rmsd_nucleotide_check(const Residue& residue) {
    const int common_count = count_matching_atom_types(residue, COMMON_RING_ATOM_TYPES);
    const int purine_count = count_matching_atom_types(residue, PURINE_RING_ATOM_TYPES);
    const int total_ring_atoms = common_count + purine_count;

    if (total_ring_atoms < 3)
        return false;

    auto rmsd = check_nt_type_by_rmsd(residue);
    return rmsd.has_value() && *rmsd <= NT_RMSD_CUTOFF;
}

} // namespace

bool BasePairFinder::is_nucleotide(const Residue& residue) {
    const auto base_type = residue.base_type();

    // Standard nucleotides (A, C, G, T, U)
    if (is_standard_nucleotide(base_type))
        return true;

    // Explicitly recognized modified nucleotides
    if (is_recognized_modified_nucleotide(base_type))
        return true;

    // Unknown or noncanonical residues need RMSD validation
    if (needs_rmsd_validation(residue)) {
        return passes_rmsd_nucleotide_check(residue);
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
            int legacy_idx = residue.legacy_residue_idx();
            if (legacy_idx > 0) {
                mapping.by_legacy_idx[legacy_idx] = &residue;
                if (legacy_idx > mapping.max_legacy_idx) {
                    mapping.max_legacy_idx = legacy_idx;
                }
            }
        }
    }

    return mapping;
}

BasePairFinder::Phase1Results BasePairFinder::run_phase1_validation(const ResidueIndexMapping& mapping) const {
    Phase1Results results;

    // Early rejection threshold - pairs with origin distance > this are skipped
    // This matches max_dorg in ValidationParameters (default 15.0)
    const double max_origin_distance_sq = validator_.parameters().max_dorg * validator_.parameters().max_dorg;

    for (int legacy_idx1 = 1; legacy_idx1 <= mapping.max_legacy_idx - 1; ++legacy_idx1) {
        auto it1 = mapping.by_legacy_idx.find(legacy_idx1);
        if (it1 == mapping.by_legacy_idx.end() || !it1->second) {
            continue;
        }
        const Residue* res1 = it1->second;
        if (!is_nucleotide(*res1) || !res1->reference_frame().has_value()) {
            continue;
        }

        // Cache origin for res1 to avoid repeated access
        const auto& origin1 = res1->reference_frame()->origin();

        for (int legacy_idx2 = legacy_idx1 + 1; legacy_idx2 <= mapping.max_legacy_idx; ++legacy_idx2) {
            auto it2 = mapping.by_legacy_idx.find(legacy_idx2);
            if (it2 == mapping.by_legacy_idx.end() || !it2->second) {
                continue;
            }
            const Residue* res2 = it2->second;
            if (!is_nucleotide(*res2) || !res2->reference_frame().has_value()) {
                continue;
            }

            // Early distance rejection - skip pairs that are too far apart
            // Uses squared distance to avoid sqrt overhead
            const auto& origin2 = res2->reference_frame()->origin();
            double dx = origin2.x() - origin1.x();
            double dy = origin2.y() - origin1.y();
            double dz = origin2.z() - origin1.z();
            double dist_sq = dx * dx + dy * dy + dz * dz;

            if (dist_sq > max_origin_distance_sq) {
                continue;  // Skip - too far apart to form a base pair
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

BasePair BasePairFinder::create_base_pair(int legacy_idx1, int legacy_idx2, const Residue* res1, const Residue* res2,
                                          const ValidationResult& result) const {
    // ALWAYS store smaller index first for consistency with legacy behavior
    size_t idx_small = static_cast<size_t>(std::min(legacy_idx1, legacy_idx2)) - 1;
    size_t idx_large = static_cast<size_t>(std::max(legacy_idx1, legacy_idx2)) - 1;
    bool swapped = (legacy_idx1 > legacy_idx2);

    // Set frames - swap if we reordered the indices
    const Residue* res_small = swapped ? res2 : res1;
    const Residue* res_large = swapped ? res1 : res2;

    // Get frames (use identity if not available, though validation should ensure they exist)
    core::ReferenceFrame frame1, frame2;
    if (res_small && res_small->reference_frame().has_value()) {
        frame1 = res_small->reference_frame().value();
    }
    if (res_large && res_large->reference_frame().has_value()) {
        frame2 = res_large->reference_frame().value();
    }

    BasePair pair(idx_small, idx_large, frame1, frame2, result.bp_type);
    pair.set_finding_order_swapped(swapped);

    // Set unique residue identifiers
    if (res_small) {
        pair.set_res_id1(res_small->res_id());
    }
    if (res_large) {
        pair.set_res_id2(res_large->res_id());
    }

    // Set hydrogen bonds
    pair.set_hydrogen_bonds(result.hbonds);

    // Determine bp_type string from residue names
    if (res_small && res_large) {
        char base1 = core::one_letter_code(*res_small);
        char base2 = core::one_letter_code(*res_large);
        if (base1 != ' ' && base2 != ' ') {
            pair.set_bp_type(std::string(1, base1) + std::string(1, base2));
        }
    }

    return pair;
}

} // namespace algorithms
} // namespace x3dna
