/**
 * @file detection_params.hpp
 * @brief Parameters for hydrogen bond detection algorithm
 */

#pragma once

#include <string>
#include <x3dna/core/hbond_types.hpp>

// Forward declaration
namespace x3dna {
namespace config {
struct HBondParameters;
}
}

namespace x3dna {
namespace algorithms {

/**
 * @brief Context-specific distance thresholds for H-bond detection
 */
struct HBondDistanceThresholds {
    // Nucleic acid contexts
    double base_base_max = 4.0;
    double base_backbone_max = 3.5;
    double backbone_backbone_max = 3.5;
    double base_sugar_max = 3.5;
    double sugar_sugar_max = 3.5;

    // Protein contexts
    double protein_mainchain_max = 3.5;
    double protein_sidechain_max = 3.5;

    // Cross-molecule contexts
    double base_protein_max = 3.5;
    double protein_ligand_max = 3.5;
    double base_ligand_max = 3.5;

    // General limits
    double min_distance = 1.8;
    double conflict_filter_distance = 4.5;

    [[nodiscard]] double max_for_context(core::HBondContext ctx) const;
};

/**
 * @brief Parameters for H-bond detection algorithm
 */
struct HBondDetectionParams {
    HBondDistanceThresholds distances;

    // Element filter
    std::string allowed_elements = ".O.N.";

    // Quality thresholds
    double good_bond_min_distance = 2.5;
    double good_bond_max_distance = 3.5;

    // Post-validation filtering
    double post_validation_max_distance = 3.6;
    double nonstandard_min_distance = 2.6;
    double nonstandard_max_distance = 3.2;

    // Minimum H-bonds for valid base pair
    int min_base_hbonds_required = 1;

    // Interaction type filter (default: all)
    core::HBondInteractionType interaction_filter = core::HBondInteractionType::ANY;

    // === Optional angle-based filtering (off by default for legacy compatibility) ===
    bool enable_angle_filtering = false;    // Filter H-bonds by angle quality
    double min_donor_angle = 90.0;          // Below = impossible geometry
    double min_acceptor_angle = 70.0;       // Below = impossible geometry

    // === Quality scoring (optional - applied after detection) ===
    bool enable_quality_scoring = false;    // Calculate quality scores for H-bonds
    bool filter_invalid_scores = false;     // Remove H-bonds with INVALID tier

    // === Unlikely chemistry filtering ===
    // AA (acceptor-acceptor) and DD (donor-donor) pairs are chemically unlikely
    // without unusual protonation or tautomerization.
    bool include_unlikely_chemistry = false;  // Include AA/DD bonds for exploration

    // === Backbone H-bond filtering ===
    // Legacy mode filters out bonds where BOTH atoms are backbone (PO list)
    // Set to true to include backbone-backbone H-bonds (like DSSR)
    bool include_backbone_backbone = false;

    // === Intra-residue H-bond detection ===
    // Set to true to detect H-bonds within a single residue
    // (e.g., O2'-N3 within the same nucleotide)
    bool include_intra_residue = false;

    // Presets (load from config file)
    [[nodiscard]] static HBondDetectionParams legacy_compatible();
    [[nodiscard]] static HBondDetectionParams modern();
    [[nodiscard]] static HBondDetectionParams general();
    [[nodiscard]] static HBondDetectionParams dssr_like();  // DSSR-compatible thresholds

    /**
     * @brief Create detection params from unified config
     * @param config Loaded HBondParameters configuration
     * @return Detection params populated from config
     */
    [[nodiscard]] static HBondDetectionParams from_config(const config::HBondParameters& config);
};

} // namespace algorithms
} // namespace x3dna
