/**
 * @file hydrogen_bond_finder.hpp
 * @brief Hydrogen bond finder - matches legacy get_hbond_ij and hb_numlist
 */

#pragma once

#include <vector>
#include <string>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/structure/residue.hpp>  // Polymorphic types
#include <x3dna/core/atom.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @struct HydrogenBondResult
 * @brief Result of hydrogen bond finding (for comparison with legacy)
 */
struct HydrogenBondResult {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;        // '-' for standard, '*' for non-standard, ' ' for invalid
    int linkage_type; // Not yet implemented, defaults to 0

    // For comparison
    bool operator==(const HydrogenBondResult& other) const {
        return donor_atom == other.donor_atom && acceptor_atom == other.acceptor_atom &&
               std::abs(distance - other.distance) < 0.001 && type == other.type;
    }
};

/**
 * @struct DetailedHBondResult
 * @brief Detailed results from H-bond finding including all steps
 */
struct DetailedHBondResult {
    std::vector<HydrogenBondResult> initial_hbonds; // Before conflict resolution
    std::vector<HydrogenBondResult> after_conflict_resolution;
    std::vector<HydrogenBondResult> after_validation; // ALL H-bonds after validation (including
                                                      // type=' ') - matches legacy JSON recording
    std::vector<HydrogenBondResult> final_hbonds;     // Only H-bonds with type != ' ' (for quality adjustment counting)
    int num_good_hb;                                  // Count of H-bonds with type='-' and distance in [2.5, 3.5]
};

/**
 * @class HydrogenBondFinder
 * @brief Finds hydrogen bonds between two residues - matches legacy get_hbond_ij
 *
 * This class is designed to match legacy's get_hbond_ij function exactly,
 * making it easy to compare and debug differences.
 */
class HydrogenBondFinder {
public:
    /**
     * @brief Count H-bonds simply (no validation) - matches legacy check_pair
     * @param res1 First residue
     * @param res2 Second residue
     * @param hb_lower Lower distance limit
     * @param hb_dist1 Upper distance limit
     * @param hb_atoms H-bond atom list (default ".O.N")
     * @param num_base_hb Output: count of base-base H-bonds
     * @param num_o2_hb Output: count of O2' H-bonds
     *
     * Matches legacy check_pair H-bond counting (lines 4605-4614 in cmn_fncs.c)
     * Counts H-bonds BEFORE validation - this is the key difference from validated counting
     */
    static void count_simple(const core::Residue& res1, const core::Residue& res2, double hb_lower, double hb_dist1,
                             const std::string& hb_atoms, int& num_base_hb, int& num_o2_hb);

    /**
     * @brief Find hydrogen bonds between two residues
     * @param res1 First residue
     * @param res2 Second residue
     * @param hb_lower Lower distance limit
     * @param hb_dist1 Upper distance limit
     * @return Vector of hydrogen bond results
     *
     * Matches legacy get_hbond_ij flow:
     * 1. Find all potential H-bonds (good_hbatoms + within_limits)
     * 2. Resolve conflicts (hb_atompair)
     * 3. Validate H-bonds (validate_hbonds)
     * 4. Return only H-bonds with type != ' '
     */
    [[nodiscard]] static std::vector<HydrogenBondResult> find_hydrogen_bonds(const core::Residue& res1,
                                                                             const core::Residue& res2, double hb_lower,
                                                                             double hb_dist1);

    /**
     * @brief Find hydrogen bonds and return detailed comparison info
     * @param res1 First residue
     * @param res2 Second residue
     * @param hb_lower Lower distance limit
     * @param hb_dist1 Upper distance limit
     * @return Detailed results including all steps
     */
    [[nodiscard]] static DetailedHBondResult find_hydrogen_bonds_detailed(const core::Residue& res1,
                                                                          const core::Residue& res2, double hb_lower,
                                                                          double hb_dist1, double hb_dist2 = 4.5);

    // === Polymorphic overloads ===

    /**
     * @brief Count H-bonds simply for polymorphic residues
     */
    static void count_simple(const core::structure::IResidue& res1, const core::structure::IResidue& res2, double hb_lower,
                             double hb_dist1, const std::string& hb_atoms, int& num_base_hb, int& num_o2_hb);

    /**
     * @brief Find hydrogen bonds between two polymorphic residues
     */
    [[nodiscard]] static std::vector<HydrogenBondResult> find_hydrogen_bonds(const core::structure::IResidue& res1,
                                                                             const core::structure::IResidue& res2,
                                                                             double hb_lower, double hb_dist1);

    /**
     * @brief Find hydrogen bonds and return detailed comparison info for polymorphic residues
     */
    [[nodiscard]] static DetailedHBondResult find_hydrogen_bonds_detailed(const core::structure::IResidue& res1,
                                                                          const core::structure::IResidue& res2,
                                                                          double hb_lower, double hb_dist1,
                                                                          double hb_dist2 = 4.5);

private:
    /**
     * @brief Resolve conflicts when same atom has multiple H-bonds
     * Matches legacy hb_atompair logic with full iterative algorithm and linkage type calculation
     * @param hbonds H-bonds to resolve (distances may be negated to mark conflicts)
     * @param hb_lower Lower distance limit
     * @param hb_dist2 Upper distance limit for linkage type checking
     */
    static void resolve_conflicts(std::vector<HydrogenBondResult>& hbonds, double hb_lower, double hb_dist2);

    /**
     * @brief Validate H-bonds based on donor-acceptor relationship
     * Matches legacy validate_hbonds logic
     * Only processes H-bonds with positive distance (conflicts marked by negative distance)
     */
    static void validate_hbonds(std::vector<HydrogenBondResult>& hbonds, char base1, char base2);

    /**
     * @brief Get base type for H-bond detection (handles modified nucleotides)
     * Matches legacy behavior: uses one_letter_code() if available, otherwise uses residue_type()
     * For modified nucleotides, returns appropriate base type (A, C, G, T, U)
     */
    [[nodiscard]] static char get_base_type_for_hbond(const core::Residue& residue);
};

} // namespace algorithms
} // namespace x3dna
