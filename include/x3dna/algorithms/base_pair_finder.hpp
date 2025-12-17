/**
 * @file base_pair_finder.hpp
 * @brief Base pair finding algorithm (matches legacy find_bestpair)
 */

#pragma once

#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/quality_score_calculator.hpp>
#include <x3dna/algorithms/pair_candidate_cache.hpp>
#include <x3dna/algorithms/pair_finding_observer.hpp>
#include <x3dna/algorithms/pair_selection_strategy.hpp>
#include <x3dna/io/json_writer.hpp>
#include <vector>
#include <memory>
#include <optional>

namespace x3dna {
namespace algorithms {

/**
 * @enum PairFindingStrategy
 * @brief Strategy for finding base pairs
 */
enum class PairFindingStrategy {
    BEST_PAIR,     // Greedy mutual best match (matches legacy find_bestpair)
    ALL_PAIRS,     // Exhaustive search (matches legacy all_pairs)
    DISTANCE_BASED // Simple distance-based
};

/**
 * @class BasePairFinder
 * @brief Finds base pairs in a structure using various strategies
 *
 * This class serves as a facade over several specialized components:
 * - BasePairValidator: Validates individual base pairs
 * - QualityScoreCalculator: Calculates adjusted quality scores
 * - PairCandidateCache: Caches Phase 1 validation results
 * - IPairSelectionStrategy: Implements pair selection algorithm
 * - IPairFindingObserver: Records events during pair finding
 *
 * Implements the legacy find_bestpair algorithm which uses a greedy mutual
 * best match strategy: for each unpaired residue, find its best partner,
 * then check if that partner's best partner is the original residue.
 */
class BasePairFinder {
public:
    /**
     * @brief Constructor
     * @param params Validation parameters
     */
    explicit BasePairFinder(const ValidationParameters& params = ValidationParameters::defaults())
        : validator_(params), strategy_(PairFindingStrategy::BEST_PAIR) {}

    /**
     * @brief Find base pairs in a structure
     * @param structure Structure to search (residues must have frames calculated)
     * @return Vector of found base pairs
     */
    [[nodiscard]] std::vector<core::BasePair> find_pairs(core::Structure& structure);

    /**
     * @brief Find base pairs (const version)
     */
    [[nodiscard]] std::vector<core::BasePair> find_pairs(const core::Structure& structure) const;

    /**
     * @brief Find base pairs and record validation results to JSON
     * @param structure Structure to search
     * @param writer JsonWriter to record validation results (can be nullptr to skip recording)
     * @return Vector of found base pairs
     */
    [[nodiscard]] std::vector<core::BasePair> find_pairs_with_recording(core::Structure& structure,
                                                          io::JsonWriter* writer = nullptr) const;

    /**
     * @brief Set finding strategy
     */
    void set_strategy(PairFindingStrategy strategy) {
        strategy_ = strategy;
    }

    /**
     * @brief Get finding strategy
     */
    [[nodiscard]] PairFindingStrategy strategy() const {
        return strategy_;
    }

    /**
     * @brief Set validation parameters
     */
    void set_parameters(const ValidationParameters& params) {
        validator_.set_parameters(params);
    }

    /**
     * @brief Get validation parameters
     */
    [[nodiscard]] const ValidationParameters& parameters() const {
        return validator_.parameters();
    }

    /**
     * @brief Check if residue is a nucleotide
     * @param residue Residue to check
     * @return True if residue is a nucleotide (standard or modified)
     */
    [[nodiscard]] static bool is_nucleotide(const core::Residue& residue);

private:
    // Core components
    BasePairValidator validator_;
    QualityScoreCalculator quality_calculator_;
    PairFindingStrategy strategy_;

    // Mutable cache for const methods
    mutable PairCandidateCache cache_;

    /**
     * @brief Find best pairs using greedy mutual best match (legacy find_bestpair)
     * @param writer Optional JsonWriter to record validation results
     */
    [[nodiscard]] std::vector<core::BasePair> find_best_pairs(core::Structure& structure, io::JsonWriter* writer = nullptr) const;

    /**
     * @brief Find all valid pairs (exhaustive search)
     */
    [[nodiscard]] std::vector<core::BasePair> find_all_pairs(const core::Structure& structure) const;

    /**
     * @brief Find best partner for a residue
     * @param residue_idx Index of residue (0-based)
     * @param structure Structure to search
     * @param matched_indices Set of already matched residue indices
     * @param writer Optional JsonWriter to record validation results
     * @return Best partner index and validation result, or nullopt if none found
     */
    [[nodiscard]] std::optional<std::pair<int, ValidationResult>> find_best_partner(
        int legacy_idx1, // Legacy 1-based residue index (from atom.legacy_residue_idx() set during
                         // PDB parsing)
        const core::Structure& structure, const std::vector<bool>& matched_indices,
        const std::map<int, const core::Residue*>& residue_by_legacy_idx,
        const std::map<std::pair<int, int>, ValidationResult>& phase1_validation_results,
        const std::map<std::pair<int, int>, int>& phase1_bp_type_ids, // Store bp_type_id from Phase 1
        io::JsonWriter* writer = nullptr) const;

    /**
     * @brief Adjust quality score based on hydrogen bonds (matches legacy adjust_pairQuality)
     * @param hbonds Vector of hydrogen bonds
     * @return Adjustment value to add to quality_score (negative value)
     */
    [[nodiscard]] double adjust_pair_quality(const std::vector<core::hydrogen_bond>& hbonds) const;

    /**
     * @brief Record validation results for a pair (matches legacy check_pair ->
     * calculate_more_bppars)
     * @param legacy_idx1 First residue legacy index (1-based)
     * @param legacy_idx2 Second residue legacy index (1-based)
     * @param res1 First residue pointer
     * @param res2 Second residue pointer
     * @param result Validation result
     * @param writer JSON writer (can be nullptr)
     */
    void record_validation_results(int legacy_idx1, int legacy_idx2, const core::Residue* res1,
                                   const core::Residue* res2, const ValidationResult& result,
                                   io::JsonWriter* writer) const;

    /**
     * @brief Calculate bp_type_id using check_wc_wobble_pair logic
     * @param res1 First residue
     * @param res2 Second residue
     * @param result Validation result
     * @param quality_score Adjusted quality score
     * @return bp_type_id (-1, 0, 1, or 2)
     *
     * Note: Full implementation requires bpstep_par (shear, stretch, opening) from Stage 6.
     * This is a simplified version that uses base pair type and direction vectors.
     */
    [[nodiscard]] int calculate_bp_type_id(const core::Residue* res1, const core::Residue* res2, const ValidationResult& result,
                             double quality_score) const;

    /**
     * @brief Get residue index in structure (0-based, counting all residues)
     */
    [[nodiscard]] static size_t get_residue_index(const core::Structure& structure, const core::Residue& residue);

    /**
     * @brief Get base letter from ResidueType (matches legacy bseq character)
     * @param type ResidueType enum value
     * @return One-letter code (A, C, G, T, U) or '?' if unknown
     */
    [[nodiscard]] static char get_base_letter_from_type(core::ResidueType type);
};

} // namespace algorithms
} // namespace x3dna
