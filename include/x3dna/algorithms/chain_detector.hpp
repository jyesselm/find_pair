/**
 * @file chain_detector.hpp
 * @brief Chain detection based on physical backbone connectivity
 */

#pragma once

#include <string>
#include <vector>
#include <optional>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {

/**
 * @struct BackboneConnectivity
 * @brief Stores backbone atom positions for connectivity checking
 */
struct BackboneConnectivity {
    std::optional<geometry::Vector3D> O3_prime; // RNA: O3' position
    std::optional<geometry::Vector3D> P;        // RNA: P position
    std::optional<geometry::Vector3D> PA;       // RNA: PA position (triphosphate)
    std::optional<geometry::Vector3D> C;        // Protein: carbonyl C
    std::optional<geometry::Vector3D> N;        // Protein: amide N
};

/**
 * @struct ConnectedChain
 * @brief Represents a chain of residues connected by backbone bonds
 */
struct ConnectedChain {
    std::vector<const core::Residue*> residues; // Residues in physical order (5'->3' or N->C)
    std::string chain_id;                       // PDB chain identifier
    bool is_rna = false;                        // True if RNA chain
    bool is_protein = false;                    // True if protein chain
};

/**
 * @class ChainDetector
 * @brief Detects physical chains in RNA and proteins based on backbone connectivity
 *
 * This class implements chain detection that orders residues based on physical
 * backbone connections rather than PDB sequence numbers:
 * - RNA/DNA: O3'(res1) → P(res2) distance < 2.75Å means connected in 5'→3' direction
 * - Protein: C(res1) → N(res2) distance < 2.0Å means connected in N-term→C-term direction
 *
 * The algorithm:
 * 1. Filter residues by type (RNA or protein)
 * 2. Sort by (chain_id, seq_num) to establish initial ordering
 * 3. Build chains by extending from unassigned residues using connectivity
 * 4. Optionally merge adjacent chains if they are close in space
 *
 * Example usage:
 * @code
 * // Parse a structure
 * Structure structure = PdbParser::parse("structure.pdb");
 *
 * // Detect RNA chains with default settings
 * ChainDetector detector;
 * auto rna_chains = detector.detect_rna_chains(structure);
 *
 * // Process each chain
 * for (const auto& chain : rna_chains) {
 *     std::cout << "Chain " << chain.chain_id << " has "
 *               << chain.residues.size() << " residues\n";
 * }
 *
 * // Custom configuration for stricter connectivity
 * ChainDetector::Config config;
 * config.rna_connectivity_cutoff = 2.5;
 * config.merge_adjacent_chains = false;
 * ChainDetector strict_detector(config);
 * auto strict_chains = strict_detector.detect_rna_chains(structure);
 * @endcode
 */
class ChainDetector {
public:
    /**
     * @struct Config
     * @brief Configuration parameters for chain detection
     */
    struct Config {
        double rna_connectivity_cutoff;     // O3'-P distance cutoff (Angstroms)
        double protein_connectivity_cutoff; // C-N distance cutoff (Angstroms)
        double chain_merge_distance;        // Sugar-sugar distance for merging chains (Angstroms)
        bool merge_adjacent_chains;         // Enable chain merging

        Config()
            : rna_connectivity_cutoff(2.75), protein_connectivity_cutoff(2.0), chain_merge_distance(8.0),
              merge_adjacent_chains(true) {}
    };

    /**
     * @brief Constructor with configuration
     * @param config Configuration parameters
     */
    explicit ChainDetector(const Config& config = Config());

    // Main detection methods

    /**
     * @brief Detect RNA/DNA chains based on backbone connectivity
     * @param structure Structure containing residues
     * @return Vector of connected RNA chains
     */
    [[nodiscard]] std::vector<ConnectedChain> detect_rna_chains(const core::Structure& structure) const;

    /**
     * @brief Detect protein chains based on peptide bond connectivity
     * @param structure Structure containing residues
     * @return Vector of connected protein chains
     */
    [[nodiscard]] std::vector<ConnectedChain> detect_protein_chains(const core::Structure& structure) const;

    /**
     * @brief Detect all chains (both RNA and protein)
     * @param structure Structure containing residues
     * @return Vector of all connected chains
     */
    [[nodiscard]] std::vector<ConnectedChain> detect_all_chains(const core::Structure& structure) const;

    // Connectivity checking

    /**
     * @brief Check if two RNA residues are connected
     * @param res1 First residue
     * @param res2 Second residue
     * @return +1 if res1->res2 (5'->3'), -1 if res2->res1 (3'->5'), 0 if not connected
     */
    [[nodiscard]] int are_rna_residues_connected(const core::Residue& res1, const core::Residue& res2) const;

    /**
     * @brief Check if two protein residues are connected
     * @param res1 First residue
     * @param res2 Second residue
     * @return +1 if res1->res2 (N->C), -1 if res2->res1 (C->N), 0 if not connected
     */
    [[nodiscard]] int are_protein_residues_connected(const core::Residue& res1, const core::Residue& res2) const;

private:
    Config config_;

    // Helper methods

    /**
     * @brief Extract backbone atom positions from a residue
     * @param residue Residue to extract from
     * @return BackboneConnectivity with available atom positions
     */
    [[nodiscard]] BackboneConnectivity extract_backbone(const core::Residue& residue) const;

    /**
     * @brief Filter residues to RNA/DNA only
     * @param structure Source structure
     * @return Vector of RNA residue pointers
     */
    [[nodiscard]] std::vector<const core::Residue*> filter_rna_residues(const core::Structure& structure) const;

    /**
     * @brief Filter residues to proteins only
     * @param structure Source structure
     * @return Vector of protein residue pointers
     */
    [[nodiscard]] std::vector<const core::Residue*> filter_protein_residues(const core::Structure& structure) const;

    /**
     * @brief Sort residues by chain_id and seq_num
     * @param residues Residues to sort (modifies in place)
     */
    void sort_by_chain_and_num(std::vector<const core::Residue*>& residues) const;

    /**
     * @brief Build chains from sorted residues using connectivity function
     * @param residues Sorted residues to process
     * @param connectivity_func Function that returns +1, -1, or 0 for connectivity
     * @return Vector of connected chains
     */
    using ConnectivityFunc = std::function<int(const core::Residue&, const core::Residue&)>;
    [[nodiscard]] std::vector<ConnectedChain> build_chains(std::vector<const core::Residue*>& residues,
                                                           const ConnectivityFunc& connectivity_func,
                                                           bool is_rna) const;

    /**
     * @brief Merge adjacent chains if they are close in space
     * @param chains Chains to merge
     * @return Merged chains
     */
    [[nodiscard]] std::vector<ConnectedChain> merge_adjacent_chains(std::vector<ConnectedChain>& chains) const;

    /**
     * @brief Calculate mean sugar atom position for an RNA residue
     * @param residue RNA residue
     * @return Mean position of sugar atoms, or nullopt if no sugar atoms
     */
    [[nodiscard]] std::optional<geometry::Vector3D> calculate_sugar_center(const core::Residue& residue) const;

    /**
     * @brief Get distance between residue positions
     * @param res1 First residue
     * @param res2 Second residue
     * @return Distance between sugar centers (RNA) or CA atoms (protein), or infinity if unavailable
     */
    [[nodiscard]] double get_residue_distance(const core::Residue& res1, const core::Residue& res2) const;
};

} // namespace algorithms
} // namespace x3dna
