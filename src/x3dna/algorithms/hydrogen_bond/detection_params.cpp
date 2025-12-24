/**
 * @file detection_params.cpp
 * @brief Implementation of H-bond detection parameters
 */

#include <x3dna/algorithms/hydrogen_bond/detection_params.hpp>
#include <x3dna/config/hbond_parameters.hpp>
#include <x3dna/config/hbond_parameters_loader.hpp>

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
    // Try to load from config preset if available
    try {
        if (config::HBondParametersLoader::has_preset("legacy_compatible")) {
            auto params = from_config(config::HBondParametersLoader::load_preset("legacy_compatible"));
            // Use ANY filter - base_atoms_only flag already filters to base atoms
            params.interaction_filter = core::HBondInteractionType::ANY;
            return params;
        }
    } catch (...) {
        // Fall through to hardcoded values
    }
    // Fallback to hardcoded values if config not available
    HBondDetectionParams params;
    params.distances.base_base_max = 4.0;
    params.distances.min_distance = 2.0;
    params.distances.conflict_filter_distance = 0.0;
    params.allowed_elements = ".O.N.";
    params.good_bond_min_distance = 2.5;
    params.good_bond_max_distance = 3.5;
    params.post_validation_max_distance = 3.6;
    params.nonstandard_min_distance = 2.6;
    params.nonstandard_max_distance = 3.2;
    params.interaction_filter = core::HBondInteractionType::ANY;
    return params;
}

HBondDetectionParams HBondDetectionParams::modern() {
    // Try to load from config preset if available
    try {
        if (config::HBondParametersLoader::has_preset("modern")) {
            auto params = from_config(config::HBondParametersLoader::load_preset("modern"));
            params.interaction_filter = core::HBondInteractionType::RNA_INTERNAL;
            return params;
        }
    } catch (...) {
        // Fall through to hardcoded values
    }
    // Fallback to hardcoded values if config not available
    HBondDetectionParams params;
    params.distances.base_base_max = 3.5;
    params.distances.base_backbone_max = 3.3;
    params.distances.backbone_backbone_max = 3.3;
    params.distances.min_distance = 2.0;
    params.allowed_elements = ".O.N.";
    params.interaction_filter = core::HBondInteractionType::RNA_INTERNAL;
    return params;
}

HBondDetectionParams HBondDetectionParams::general() {
    // Try to load from config preset if available
    try {
        if (config::HBondParametersLoader::has_preset("general")) {
            auto params = from_config(config::HBondParametersLoader::load_preset("general"));
            params.interaction_filter = core::HBondInteractionType::ANY;
            return params;
        }
    } catch (...) {
        // Fall through to hardcoded values
    }
    // Fallback to hardcoded values if config not available
    HBondDetectionParams params;
    params.distances.base_base_max = 3.5;
    params.distances.protein_mainchain_max = 3.5;
    params.distances.protein_sidechain_max = 3.5;
    params.distances.base_protein_max = 3.5;
    params.distances.protein_ligand_max = 3.5;
    params.distances.min_distance = 2.0;
    params.allowed_elements = ".O.N.S.";
    params.interaction_filter = core::HBondInteractionType::ANY;
    return params;
}

HBondDetectionParams HBondDetectionParams::dssr_like() {
    // Try to load from config preset if available
    try {
        if (config::HBondParametersLoader::has_preset("dssr_like")) {
            return from_config(config::HBondParametersLoader::load_preset("dssr_like"));
        }
    } catch (...) {
        // Fall through to hardcoded values
    }
    // Fallback to hardcoded values if config not available
    HBondDetectionParams params;
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
    params.distances.conflict_filter_distance = 4.5;
    params.allowed_elements = ".O.N.";
    params.good_bond_min_distance = 2.5;
    params.good_bond_max_distance = 3.5;
    params.post_validation_max_distance = 3.6;
    params.nonstandard_min_distance = 2.6;
    params.nonstandard_max_distance = 3.2;
    params.interaction_filter = core::HBondInteractionType::ANY;
    params.include_backbone_backbone = true;
    return params;
}

HBondDetectionParams HBondDetectionParams::from_config(const config::HBondParameters& config) {
    HBondDetectionParams params;

    // Distance thresholds
    params.distances.min_distance = config.detection.distance.min;
    params.distances.base_base_max = config.detection.distance.base_base_max;
    params.distances.base_backbone_max = config.detection.distance.base_backbone_max;
    params.distances.backbone_backbone_max = config.detection.distance.backbone_backbone_max;
    params.distances.base_sugar_max = config.detection.distance.base_sugar_max;
    params.distances.sugar_sugar_max = config.detection.distance.sugar_sugar_max;
    params.distances.protein_mainchain_max = config.detection.distance.protein_mainchain_max;
    params.distances.protein_sidechain_max = config.detection.distance.protein_sidechain_max;
    params.distances.base_protein_max = config.detection.distance.base_protein_max;
    params.distances.protein_ligand_max = config.detection.distance.protein_ligand_max;
    params.distances.base_ligand_max = config.detection.distance.base_ligand_max;
    params.distances.conflict_filter_distance = config.detection.distance.conflict_filter;

    // Element filter
    params.allowed_elements = config.detection.elements.allowed;

    // Quality thresholds
    params.good_bond_min_distance = config.detection.thresholds.good_bond.min;
    params.good_bond_max_distance = config.detection.thresholds.good_bond.max;
    params.post_validation_max_distance = config.detection.thresholds.post_validation_max;
    params.nonstandard_min_distance = config.detection.thresholds.nonstandard.min;
    params.nonstandard_max_distance = config.detection.thresholds.nonstandard.max;

    // Validation
    params.min_base_hbonds_required = config.detection.validation.min_base_hbonds;

    // Options
    params.enable_angle_filtering = config.detection.options.enable_angle_filtering;
    params.min_donor_angle = config.geometry.donor_angle.min;
    params.min_acceptor_angle = config.geometry.acceptor_angle.min;
    params.enable_quality_scoring = config.detection.options.enable_quality_scoring;
    params.filter_invalid_scores = config.detection.options.filter_invalid_scores;
    params.include_unlikely_chemistry = config.detection.options.include_unlikely_chemistry;
    params.include_backbone_backbone = config.detection.options.include_backbone_backbone;

    return params;
}

} // namespace algorithms
} // namespace x3dna
