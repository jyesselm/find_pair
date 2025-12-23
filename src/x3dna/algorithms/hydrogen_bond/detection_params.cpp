/**
 * @file detection_params.cpp
 * @brief Implementation of H-bond detection parameters
 */

#include <x3dna/algorithms/hydrogen_bond/detection_params.hpp>

namespace x3dna {
namespace algorithms {

double HBondDistanceThresholds::max_for_context(core::HBondContext ctx) const {
    switch (ctx) {
        case core::HBondContext::BASE_BASE:
            return base_base_max;
        case core::HBondContext::BASE_BACKBONE:
            return base_backbone_max;
        case core::HBondContext::BACKBONE_BACKBONE:
            return backbone_backbone_max;
        case core::HBondContext::BASE_SUGAR:
            return base_sugar_max;
        case core::HBondContext::SUGAR_SUGAR:
            return sugar_sugar_max;
        case core::HBondContext::PROTEIN_MAINCHAIN:
            return protein_mainchain_max;
        case core::HBondContext::PROTEIN_SIDECHAIN:
            return protein_sidechain_max;
        case core::HBondContext::BASE_PROTEIN:
            return base_protein_max;
        case core::HBondContext::SUGAR_PROTEIN:
        case core::HBondContext::BACKBONE_PROTEIN:
            return base_protein_max;
        case core::HBondContext::BASE_LIGAND:
            return base_ligand_max;
        case core::HBondContext::PROTEIN_LIGAND:
            return protein_ligand_max;
        case core::HBondContext::LIGAND_LIGAND:
            return protein_ligand_max;
        default:
            return base_base_max;
    }
}

HBondDetectionParams HBondDetectionParams::legacy_compatible() {
    HBondDetectionParams params;
    // Legacy uses hb_dist1=4.0 for initial detection (from base_pair_validator.hpp)
    params.distances.base_base_max = 4.0;
    params.distances.min_distance = 2.0;  // Legacy hb_lower=2.0
    // Production uses hb_dist2=0.0, which means Phase 3 doesn't promote non-winners
    params.distances.conflict_filter_distance = 0.0;
    params.allowed_elements = ".O.N.";
    params.good_bond_min_distance = 2.5;
    params.good_bond_max_distance = 3.5;
    params.post_validation_max_distance = 3.6;
    params.nonstandard_min_distance = 2.6;
    params.nonstandard_max_distance = 3.2;
    params.interaction_filter = core::HBondInteractionType::BASE_BASE;
    return params;
}

HBondDetectionParams HBondDetectionParams::modern() {
    HBondDetectionParams params;
    // Tighter thresholds for better accuracy
    params.distances.base_base_max = 3.5;
    params.distances.base_backbone_max = 3.3;
    params.distances.backbone_backbone_max = 3.3;
    params.distances.min_distance = 2.0;
    params.allowed_elements = ".O.N.";
    params.interaction_filter = core::HBondInteractionType::RNA_INTERNAL;
    return params;
}

HBondDetectionParams HBondDetectionParams::general() {
    HBondDetectionParams params;
    // Support all molecule types
    params.distances.base_base_max = 3.5;
    params.distances.protein_mainchain_max = 3.5;
    params.distances.protein_sidechain_max = 3.5;
    params.distances.base_protein_max = 3.5;
    params.distances.protein_ligand_max = 3.5;
    params.distances.min_distance = 2.0;
    params.allowed_elements = ".O.N.S."; // Include sulfur for Cys
    params.interaction_filter = core::HBondInteractionType::ANY;
    return params;
}

HBondDetectionParams HBondDetectionParams::dssr_like() {
    HBondDetectionParams params;
    // DSSR-compatible thresholds
    // DSSR reports all H-bonds: RNA internal, protein, and RNA-protein
    // Use 3.5Å cutoff to match DSSR's coverage
    params.distances.base_base_max = 3.5;
    params.distances.base_backbone_max = 3.5;
    params.distances.backbone_backbone_max = 3.5;
    params.distances.base_sugar_max = 3.5;
    params.distances.sugar_sugar_max = 3.5;
    params.distances.protein_mainchain_max = 3.5;
    params.distances.protein_sidechain_max = 3.5;
    params.distances.base_protein_max = 3.5;
    params.distances.protein_ligand_max = 3.5;
    params.distances.min_distance = 2.0;
    params.distances.conflict_filter_distance = 4.5;  // Allow phase 3 promotion
    params.allowed_elements = ".O.N.";
    params.good_bond_min_distance = 2.5;
    params.good_bond_max_distance = 3.5;
    params.post_validation_max_distance = 3.6;
    params.nonstandard_min_distance = 2.6;
    params.nonstandard_max_distance = 3.2;
    // Include ALL interactions (RNA, protein, RNA-protein)
    params.interaction_filter = core::HBondInteractionType::ANY;
    return params;
}

} // namespace algorithms
} // namespace x3dna
