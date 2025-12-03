/**
 * @file config_manager.cpp
 * @brief Configuration management implementation
 */

#include <x3dna/config/config_manager.hpp>
#include <fstream>
#include <iostream>
#include <cstdlib>

namespace x3dna {
namespace config {

ConfigManager& ConfigManager::instance() {
    static ConfigManager instance;
    return instance;
}

void ConfigManager::load_from_file(const std::filesystem::path& config_path) {
    if (!std::filesystem::exists(config_path)) {
        std::cerr << "Warning: Config file not found: " << config_path << "\n";
        set_defaults();
        return;
    }

    try {
        std::ifstream file(config_path);
        nlohmann::json json;
        file >> json;
        load_from_json(json);
    } catch (const std::exception& e) {
        std::cerr << "Error loading config file: " << e.what() << "\n";
        set_defaults();
    }
}

void ConfigManager::load_from_json(const nlohmann::json& json) {
    // Load thresholds if present
    if (json.contains("thresholds")) {
        const auto& thresh = json["thresholds"];
        if (thresh.contains("min_dorg"))
            thresholds_.min_dorg = thresh["min_dorg"];
        if (thresh.contains("max_dorg"))
            thresholds_.max_dorg = thresh["max_dorg"];
        if (thresh.contains("min_dv"))
            thresholds_.min_dv = thresh["min_dv"];
        if (thresh.contains("max_dv"))
            thresholds_.max_dv = thresh["max_dv"];
        if (thresh.contains("min_dNN"))
            thresholds_.min_dNN = thresh["min_dNN"];
        if (thresh.contains("max_dNN"))
            thresholds_.max_dNN = thresh["max_dNN"];
        if (thresh.contains("min_plane_angle"))
            thresholds_.min_plane_angle = thresh["min_plane_angle"];
        if (thresh.contains("max_plane_angle"))
            thresholds_.max_plane_angle = thresh["max_plane_angle"];
        if (thresh.contains("min_base_hb"))
            thresholds_.min_base_hb = thresh["min_base_hb"];
        if (thresh.contains("hb_lower"))
            thresholds_.hb_lower = thresh["hb_lower"];
        if (thresh.contains("hb_dist1"))
            thresholds_.hb_dist1 = thresh["hb_dist1"];
        if (thresh.contains("hb_dist2"))
            thresholds_.hb_dist2 = thresh["hb_dist2"];
        if (thresh.contains("hb_atoms"))
            thresholds_.hb_atoms = thresh["hb_atoms"];
        if (thresh.contains("overlap_threshold"))
            thresholds_.overlap_threshold = thresh["overlap_threshold"];
        if (thresh.contains("helix_break"))
            thresholds_.helix_break = thresh["helix_break"];
    }

    // Load paths
    if (json.contains("x3dna_home")) {
        x3dna_home_ = json["x3dna_home"].get<std::string>();
    }

    // Load options
    if (json.contains("include_hetatm")) {
        include_hetatm_ = json["include_hetatm"];
    }
    if (json.contains("include_waters")) {
        include_waters_ = json["include_waters"];
    }
    if (json.contains("legacy_mode")) {
        legacy_mode_ = json["legacy_mode"];
    }
}

void ConfigManager::set_defaults() {
    // Reset thresholds to defaults
    thresholds_ = ParameterThresholds();

    // Reset options
    include_hetatm_ = false;
    include_waters_ = false;
    legacy_mode_ = false;

    // Set x3dna_home if not already set
    if (x3dna_home_.empty()) {
        // Try to get from environment variable
        const char* env_home = std::getenv("X3DNA_HOMEDIR");
        if (env_home) {
            x3dna_home_ = env_home;
        }
    }
}

void ConfigManager::set_x3dna_home(const std::filesystem::path& path) {
    x3dna_home_ = path;
}

std::filesystem::path ConfigManager::standard_base_path() const {
    if (x3dna_home_.empty()) {
        return std::filesystem::path("data/templates");
    }
    return x3dna_home_ / "templates";
}

} // namespace config
} // namespace x3dna
