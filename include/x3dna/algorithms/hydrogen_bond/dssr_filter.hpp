/**
 * @file dssr_filter.hpp
 * @brief DSSR-style H-bond filtering to match DSSR's selection criteria
 *
 * DSSR uses tiered distance thresholds based on element pairs:
 * - N-O, O-N, N-N bonds: standard donor-acceptor, up to 3.5Å
 * - O-O bonds involving O2': sugar hydroxyl, up to 3.3Å
 * - Other O-O bonds: acceptor-acceptor, up to 3.0Å (stricter)
 */

#pragma once

#include <string>
#include <vector>
#include <x3dna/algorithms/hydrogen_bond/hbond.hpp>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

// Forward declaration
struct StructureHBondResult;

/**
 * @brief Parameters for DSSR-style H-bond filtering
 */
struct DSSRFilterParams {
    // Distance thresholds by element pair type
    // Note: DSSR uses 4.0Å detection threshold for all bonds
    double n_containing_max_distance = 4.0;    // N-O, O-N, N-N bonds
    double o2prime_oo_max_distance = 3.7;      // O-O bonds involving O2' (sugar hydroxyl)
    double other_oo_max_distance = 3.5;        // Other O-O bonds

    // Minimum distance for all bonds
    double min_distance = 2.0;

    // Whether to include specific bond types that DSSR reports
    bool include_intra_residue_o2prime = true;  // O2'-N3, O2'-O2 within residue
    bool include_sequential_backbone = true;    // O5'-OP2 between adjacent residues

    // Whether to filter chemically unlikely pairs (amino-amino, carbonyl-carbonyl)
    bool filter_unlikely_pairs = true;

    /**
     * @brief Default DSSR-compatible parameters
     */
    static DSSRFilterParams defaults() {
        return DSSRFilterParams{};
    }

    /**
     * @brief Tighter thresholds - balance between count reduction and match rate
     */
    static DSSRFilterParams tight() {
        DSSRFilterParams params;
        params.n_containing_max_distance = 3.6;   // 3.6Å for N-O, N-N
        params.o2prime_oo_max_distance = 3.4;     // 3.4Å for O2'-O
        params.other_oo_max_distance = 3.2;       // 3.2Å for other O-O
        return params;
    }

    /**
     * @brief Stricter parameters for higher precision (may lose some matches)
     */
    static DSSRFilterParams strict() {
        DSSRFilterParams params;
        params.n_containing_max_distance = 3.4;
        params.o2prime_oo_max_distance = 3.2;
        params.other_oo_max_distance = 2.9;
        return params;
    }
};

/**
 * @brief DSSR-style H-bond filter
 *
 * Applies tiered distance thresholds based on element pairs to match
 * DSSR's H-bond selection behavior. Key insight: DSSR is stricter on
 * O-O bonds (acceptor-acceptor) than on N-containing bonds (donor-acceptor).
 */
class DSSRStyleFilter {
public:
    /**
     * @brief Check if an atom pair is chemically unlikely to form an H-bond
     *
     * Filters out:
     * - Amino-amino pairs (N6-N4, N6-N2, N4-N2) - both are donors
     * - Carbonyl-carbonyl pairs (O6-O4, O6-O2, O4-O2) - both are acceptors
     *
     * @param atom1_name First atom name
     * @param atom2_name Second atom name
     * @return true if pair is chemically unlikely
     */
    [[nodiscard]] static bool is_chemically_unlikely_pair(
        const std::string& atom1_name,
        const std::string& atom2_name);

    /**
     * @brief Check if an H-bond should be kept based on DSSR criteria
     * @param hb The H-bond to evaluate
     * @param params Filter parameters
     * @return true if bond passes DSSR-style filtering
     */
    [[nodiscard]] static bool should_keep(const core::HBond& hb,
                                          const DSSRFilterParams& params = DSSRFilterParams::defaults());

    /**
     * @brief Filter a vector of H-bonds, returning only those that pass
     * @param hbonds Input H-bonds
     * @param params Filter parameters
     * @return Filtered H-bonds
     */
    [[nodiscard]] static std::vector<core::HBond> filter(
        const std::vector<core::HBond>& hbonds,
        const DSSRFilterParams& params = DSSRFilterParams::defaults());

    /**
     * @brief Filter StructureHBondResult in place
     * @param result Structure H-bond result to filter (modified in place)
     * @param params Filter parameters
     */
    static void filter_in_place(StructureHBondResult& result,
                                const DSSRFilterParams& params = DSSRFilterParams::defaults());

    /**
     * @brief Apply scoring-based occupancy filter
     * @param result Structure H-bond result to filter (modified in place)
     * @param max_bonds_per_atom Maximum bonds per atom (keeps highest scoring)
     *
     * Unlike simple distance-based occupancy, this scores each bond using
     * geometric criteria (distance + angles) and keeps the highest-scoring
     * bonds for each atom. This better matches DSSR's selection behavior.
     */
    static void apply_scored_occupancy_filter(StructureHBondResult& result,
                                              int max_bonds_per_atom = 2);

    /**
     * @brief Get element from atom name (N, O, C, etc.)
     * @param atom_name Atom name (e.g., "N7", "O2'", "OP1")
     * @return Single character element symbol
     */
    [[nodiscard]] static char get_element(const std::string& atom_name);

    /**
     * @brief Check if atom is O2' (ribose 2' hydroxyl)
     * @param atom_name Atom name
     * @return true if atom is O2'
     */
    [[nodiscard]] static bool is_o2_prime(const std::string& atom_name);

    /**
     * @brief Get the appropriate distance threshold for an atom pair
     * @param atom1_name First atom name
     * @param atom2_name Second atom name
     * @param params Filter parameters
     * @return Maximum allowed distance for this atom pair
     */
    [[nodiscard]] static double get_distance_threshold(
        const std::string& atom1_name,
        const std::string& atom2_name,
        const DSSRFilterParams& params = DSSRFilterParams::defaults());
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
