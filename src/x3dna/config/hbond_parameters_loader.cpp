/**
 * @file hbond_parameters_loader.cpp
 * @brief Implementation of H-bond parameters loading from JSON
 */

#include <x3dna/config/hbond_parameters_loader.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace x3dna {
namespace config {

// Static member initialization
std::optional<HBondParameters> HBondParametersLoader::cached_params_;
nlohmann::json HBondParametersLoader::cached_json_;

HBondParameters HBondParameters::defaults() {
    return HBondParameters{};
}

std::filesystem::path HBondParametersLoader::default_config_path() {
    if (ResourceLocator::is_initialized()) {
        return ResourceLocator::config_file("hbond_parameters.json");
    }
    // Fallback to relative path
    return std::filesystem::path("resources/config/hbond_parameters.json");
}

HBondParameters HBondParametersLoader::load() {
    auto path = default_config_path();
    if (std::filesystem::exists(path)) {
        return load_from_file(path);
    }
    // Return defaults if file not found
    return HBondParameters::defaults();
}

HBondParameters HBondParametersLoader::load_from_file(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("H-bond config file not found: " + path.string());
    }

    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Failed to open H-bond config file: " + path.string());
    }

    nlohmann::json json;
    try {
        file >> json;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("Failed to parse H-bond config file: " + std::string(e.what()));
    }

    // Cache the JSON for preset loading
    cached_json_ = json;

    return load_from_json(json);
}

HBondParameters HBondParametersLoader::load_from_json(const nlohmann::json& json) {
    HBondParameters params;

    if (json.contains("detection")) {
        load_detection(params.detection, json["detection"]);
    }
    if (json.contains("geometry")) {
        load_geometry(params.geometry, json["geometry"]);
    }
    if (json.contains("scoring")) {
        load_scoring(params.scoring, json["scoring"]);
    }
    if (json.contains("quality_tiers")) {
        load_quality_tiers(params.quality_tiers, json["quality_tiers"]);
    }

    return params;
}

void HBondParametersLoader::load_detection(HBondDetectionConfig& config, const nlohmann::json& json) {
    // Distance config
    if (json.contains("distance")) {
        const auto& dist = json["distance"];
        if (dist.contains("min")) config.distance.min = dist["min"];
        if (dist.contains("base_base_max")) config.distance.base_base_max = dist["base_base_max"];
        if (dist.contains("base_backbone_max")) config.distance.base_backbone_max = dist["base_backbone_max"];
        if (dist.contains("backbone_backbone_max")) config.distance.backbone_backbone_max = dist["backbone_backbone_max"];
        if (dist.contains("base_sugar_max")) config.distance.base_sugar_max = dist["base_sugar_max"];
        if (dist.contains("sugar_sugar_max")) config.distance.sugar_sugar_max = dist["sugar_sugar_max"];
        if (dist.contains("protein_mainchain_max")) config.distance.protein_mainchain_max = dist["protein_mainchain_max"];
        if (dist.contains("protein_sidechain_max")) config.distance.protein_sidechain_max = dist["protein_sidechain_max"];
        if (dist.contains("base_protein_max")) config.distance.base_protein_max = dist["base_protein_max"];
        if (dist.contains("protein_ligand_max")) config.distance.protein_ligand_max = dist["protein_ligand_max"];
        if (dist.contains("base_ligand_max")) config.distance.base_ligand_max = dist["base_ligand_max"];
        if (dist.contains("conflict_filter")) config.distance.conflict_filter = dist["conflict_filter"];
    }

    // Elements config
    if (json.contains("elements")) {
        const auto& elem = json["elements"];
        if (elem.contains("allowed")) config.elements.allowed = elem["allowed"];
    }

    // Thresholds config
    if (json.contains("thresholds")) {
        const auto& thresh = json["thresholds"];
        if (thresh.contains("good_bond")) {
            const auto& good = thresh["good_bond"];
            if (good.contains("min")) config.thresholds.good_bond.min = good["min"];
            if (good.contains("max")) config.thresholds.good_bond.max = good["max"];
        }
        if (thresh.contains("post_validation_max")) {
            config.thresholds.post_validation_max = thresh["post_validation_max"];
        }
        if (thresh.contains("nonstandard")) {
            const auto& ns = thresh["nonstandard"];
            if (ns.contains("min")) config.thresholds.nonstandard.min = ns["min"];
            if (ns.contains("max")) config.thresholds.nonstandard.max = ns["max"];
        }
    }

    // Validation config
    if (json.contains("validation")) {
        const auto& val = json["validation"];
        if (val.contains("min_base_hbonds")) config.validation.min_base_hbonds = val["min_base_hbonds"];
    }

    // Options config
    if (json.contains("options")) {
        const auto& opt = json["options"];
        if (opt.contains("enable_angle_filtering")) config.options.enable_angle_filtering = opt["enable_angle_filtering"];
        if (opt.contains("enable_quality_scoring")) config.options.enable_quality_scoring = opt["enable_quality_scoring"];
        if (opt.contains("filter_invalid_scores")) config.options.filter_invalid_scores = opt["filter_invalid_scores"];
        if (opt.contains("include_unlikely_chemistry")) config.options.include_unlikely_chemistry = opt["include_unlikely_chemistry"];
        if (opt.contains("include_backbone_backbone")) config.options.include_backbone_backbone = opt["include_backbone_backbone"];
        if (opt.contains("include_intra_residue")) config.options.include_intra_residue = opt["include_intra_residue"];
    }
}

void HBondParametersLoader::load_geometry(HBondGeometryConfig& config, const nlohmann::json& json) {
    if (json.contains("donor_angle")) {
        const auto& donor = json["donor_angle"];
        if (donor.contains("min")) config.donor_angle.min = donor["min"];
        if (donor.contains("ideal")) config.donor_angle.ideal = donor["ideal"];
    }

    if (json.contains("acceptor_angle")) {
        const auto& acc = json["acceptor_angle"];
        if (acc.contains("min")) config.acceptor_angle.min = acc["min"];
        if (acc.contains("ideal_sp2")) config.acceptor_angle.ideal_sp2 = acc["ideal_sp2"];
        if (acc.contains("ideal_sp3")) config.acceptor_angle.ideal_sp3 = acc["ideal_sp3"];
    }
}

void HBondParametersLoader::load_scoring(HBondScoringConfig& config, const nlohmann::json& json) {
    if (json.contains("distance")) {
        const auto& dist = json["distance"];
        if (dist.contains("ideal")) config.distance.ideal = dist["ideal"];
        if (dist.contains("sigma")) config.distance.sigma = dist["sigma"];
        if (dist.contains("min")) config.distance.min = dist["min"];
        if (dist.contains("max")) config.distance.max = dist["max"];
    }

    if (json.contains("weights")) {
        const auto& w = json["weights"];
        if (w.contains("distance")) config.weights.distance = w["distance"];
        if (w.contains("donor_angle")) config.weights.donor_angle = w["donor_angle"];
        if (w.contains("acceptor_angle")) config.weights.acceptor_angle = w["acceptor_angle"];
    }

    if (json.contains("resolution")) {
        const auto& res = json["resolution"];
        if (res.contains("apply_penalty")) config.resolution.apply_penalty = res["apply_penalty"];
        if (res.contains("high_res_threshold")) config.resolution.high_res_threshold = res["high_res_threshold"];
        if (res.contains("low_res_threshold")) config.resolution.low_res_threshold = res["low_res_threshold"];
    }
}

void HBondParametersLoader::load_quality_tiers(QualityTiersConfig& config, const nlohmann::json& json) {
    if (json.contains("excellent")) {
        const auto& tier = json["excellent"];
        if (tier.contains("min_score")) config.excellent_min = tier["min_score"];
    }
    if (json.contains("standard")) {
        const auto& tier = json["standard"];
        if (tier.contains("min_score")) config.standard_min = tier["min_score"];
    }
    if (json.contains("acceptable")) {
        const auto& tier = json["acceptable"];
        if (tier.contains("min_score")) config.acceptable_min = tier["min_score"];
    }
    if (json.contains("questionable")) {
        const auto& tier = json["questionable"];
        if (tier.contains("min_score")) config.questionable_min = tier["min_score"];
    }
}

HBondParameters HBondParametersLoader::load_preset(const std::string& preset_name) {
    // Ensure we have the JSON loaded
    if (cached_json_.empty()) {
        load();  // This populates cached_json_
    }

    if (!cached_json_.contains("presets")) {
        throw std::runtime_error("No presets defined in H-bond config");
    }

    const auto& presets = cached_json_["presets"];
    if (!presets.contains(preset_name)) {
        throw std::runtime_error("Unknown H-bond preset: " + preset_name);
    }

    // Start with defaults
    HBondParameters params = HBondParameters::defaults();

    // Apply preset overrides
    apply_preset(params, presets[preset_name]);

    return params;
}

void HBondParametersLoader::apply_preset(HBondParameters& base, const nlohmann::json& preset_json) {
    if (preset_json.contains("detection")) {
        load_detection(base.detection, preset_json["detection"]);
    }
    if (preset_json.contains("geometry")) {
        load_geometry(base.geometry, preset_json["geometry"]);
    }
    if (preset_json.contains("scoring")) {
        load_scoring(base.scoring, preset_json["scoring"]);
    }
    if (preset_json.contains("quality_tiers")) {
        load_quality_tiers(base.quality_tiers, preset_json["quality_tiers"]);
    }
}

const HBondParameters& HBondParametersLoader::instance() {
    if (!cached_params_.has_value()) {
        cached_params_ = load();
    }
    return *cached_params_;
}

void HBondParametersLoader::reload() {
    cached_params_.reset();
    cached_json_.clear();
    cached_params_ = load();
}

std::vector<std::string> HBondParametersLoader::available_presets() {
    // Ensure we have the JSON loaded
    if (cached_json_.empty()) {
        load();
    }

    std::vector<std::string> names;
    if (cached_json_.contains("presets")) {
        for (auto& [key, value] : cached_json_["presets"].items()) {
            // Skip description fields
            if (key[0] != '_') {
                names.push_back(key);
            }
        }
    }
    return names;
}

bool HBondParametersLoader::has_preset(const std::string& name) {
    // Ensure we have the JSON loaded
    if (cached_json_.empty()) {
        load();
    }

    return cached_json_.contains("presets") && cached_json_["presets"].contains(name);
}

} // namespace config
} // namespace x3dna
