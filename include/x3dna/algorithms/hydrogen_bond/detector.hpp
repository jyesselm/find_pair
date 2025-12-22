/**
 * @file detector.hpp
 * @brief General-purpose H-bond detector supporting multiple molecule types
 */

#pragma once

#include <vector>
#include <x3dna/core/hbond.hpp>
#include <x3dna/core/hbond_types.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/algorithms/hydrogen_bond/detection_params.hpp>

namespace x3dna {

// Forward declaration
namespace core {
class Structure;
}

namespace algorithms {
namespace hydrogen_bond {

/**
 * @brief Result from detailed H-bond detection pipeline
 */
struct HBondPipelineResult {
    std::vector<core::HBond> initial_candidates;        // All potential H-bonds found
    std::vector<core::HBond> after_conflict_resolution; // After atom-sharing conflicts resolved
    std::vector<core::HBond> all_classified_bonds;      // All bonds with classification (including INVALID)
    std::vector<core::HBond> final_bonds;               // Only valid bonds (STANDARD + NON_STANDARD)

    int standard_bond_count = 0; // Count of STANDARD classification
    int good_bond_count = 0;     // Count in good distance range [2.5, 3.5]
};

/**
 * @brief General-purpose H-bond detector with configurable parameters
 *
 * Key design principles:
 * - Stateless algorithms (all methods const)
 * - Always calculate angles (no conditional logic)
 * - Context-aware distance thresholds
 * - Pipeline architecture: detect -> resolve -> classify -> filter
 * - Supports nucleic acids, proteins, and ligands
 *
 * For legacy base-pair H-bond detection, use legacy_compatible() params.
 */
class HBondDetector {
public:
    explicit HBondDetector(const HBondDetectionParams& params = HBondDetectionParams::legacy_compatible());

    // === Base-Pair H-Bond Detection (legacy behavior) ===

    /**
     * @brief Detect H-bonds between two residues (base atoms only)
     * @param residue1 First residue
     * @param residue2 Second residue
     * @return Vector of validated H-bonds
     *
     * Full pipeline: detect -> resolve conflicts -> classify -> filter
     */
    [[nodiscard]] std::vector<core::HBond> detect_base_hbonds(const core::Residue& residue1,
                                                              const core::Residue& residue2) const;

    /**
     * @brief Detect H-bonds with detailed pipeline results (for debugging)
     * @param residue1 First residue
     * @param residue2 Second residue
     * @return Detailed results including all pipeline stages
     */
    [[nodiscard]] HBondPipelineResult detect_base_hbonds_detailed(const core::Residue& residue1,
                                                                  const core::Residue& residue2) const;

    // === General H-Bond Detection ===

    /**
     * @brief Detect all H-bonds between two residues (any atom types)
     * @param residue1 First residue
     * @param residue2 Second residue
     * @param mol1_type Molecule type for residue1
     * @param mol2_type Molecule type for residue2
     * @return Vector of all H-bonds found
     */
    [[nodiscard]] std::vector<core::HBond> detect_all_hbonds_between(
        const core::Residue& residue1, const core::Residue& residue2,
        core::typing::MoleculeType mol1_type = core::typing::MoleculeType::NUCLEIC_ACID,
        core::typing::MoleculeType mol2_type = core::typing::MoleculeType::NUCLEIC_ACID) const;

    /**
     * @brief Detect all H-bonds with detailed pipeline results
     * @param residue1 First residue
     * @param residue2 Second residue
     * @param mol1_type Molecule type for residue1
     * @param mol2_type Molecule type for residue2
     * @return Detailed results including all pipeline stages
     */
    [[nodiscard]] HBondPipelineResult detect_all_hbonds_detailed(
        const core::Residue& residue1, const core::Residue& residue2,
        core::typing::MoleculeType mol1_type = core::typing::MoleculeType::NUCLEIC_ACID,
        core::typing::MoleculeType mol2_type = core::typing::MoleculeType::NUCLEIC_ACID) const;

    // === Counting (for validation checks) ===

    /**
     * @brief Count potential H-bonds without full detection
     * @param residue1 First residue
     * @param residue2 Second residue
     * @param base_hbond_count Output: count of base-base H-bonds
     * @param o2_prime_hbond_count Output: count of O2' H-bonds
     *
     * Fast counting for validation checks (matches legacy check_pair).
     */
    void count_potential_hbonds(const core::Residue& residue1, const core::Residue& residue2, int& base_hbond_count,
                                int& o2_prime_hbond_count) const;

    [[nodiscard]] const HBondDetectionParams& params() const {
        return params_;
    }

private:
    HBondDetectionParams params_;

    // === Pipeline Implementation ===

    /**
     * @brief Internal detection with explicit options
     */
    [[nodiscard]] HBondPipelineResult detect_internal(const core::Residue& residue1, const core::Residue& residue2,
                                                       bool base_atoms_only, core::typing::MoleculeType mol1_type,
                                                       core::typing::MoleculeType mol2_type) const;

    /**
     * @brief Find candidate H-bonds based on distance and element criteria
     * @param residue1 First residue
     * @param residue2 Second residue
     * @param base_atoms_only If true, only check base atoms
     * @param mol1_type Molecule type for residue1
     * @param mol2_type Molecule type for residue2
     * @return Vector of candidate H-bonds
     */
    [[nodiscard]] std::vector<core::HBond> find_candidate_bonds(const core::Residue& residue1,
                                                                const core::Residue& residue2, bool base_atoms_only,
                                                                core::typing::MoleculeType mol1_type,
                                                                core::typing::MoleculeType mol2_type) const;

    /**
     * @brief Resolve conflicts when same atom participates in multiple H-bonds
     * @param bonds H-bonds to resolve (modified in place)
     *
     * 3-phase algorithm:
     * 1. Mark conflicts by finding shortest bond for each shared atom
     * 2. Calculate linkage types
     * 3. Apply filtering based on linkage types
     */
    void resolve_atom_sharing_conflicts(std::vector<core::HBond>& bonds) const;

    /**
     * @brief Classify H-bonds based on donor/acceptor roles
     * @param bonds H-bonds to classify (modified in place)
     * @param base1_type One-letter code for residue1
     * @param base2_type One-letter code for residue2
     */
    void classify_bonds(std::vector<core::HBond>& bonds, char base1_type, char base2_type) const;

    /**
     * @brief Calculate angles for all H-bonds
     * @param bonds H-bonds to process (modified in place)
     * @param residue1 First residue
     * @param residue2 Second residue
     */
    void calculate_angles(std::vector<core::HBond>& bonds, const core::Residue& residue1,
                          const core::Residue& residue2) const;

    /**
     * @brief Apply post-validation filtering (marks bonds as INVALID but doesn't remove them)
     * @param bonds H-bonds to filter (modified in place)
     */
    void apply_post_validation_filtering(std::vector<core::HBond>& bonds) const;

    /**
     * @brief Get base type for H-bond detection (handles modified nucleotides)
     * @param residue Residue to get base type for
     * @return One-letter base type (A, C, G, T, U, or '?')
     */
    [[nodiscard]] static char get_base_type_for_hbond(const core::Residue& residue);

    // === Phase 1-3: Conflict Resolution ===

    void resolve_conflicts_phase1(std::vector<core::HBond>& bonds) const;
    void resolve_conflicts_phase2(std::vector<core::HBond>& bonds) const;
    void resolve_conflicts_phase3(std::vector<core::HBond>& bonds) const;
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
