/**
 * @file config_manager.hpp
 * @brief Configuration management for X3DNA
 */

#pragma once

#include <filesystem>
#include <nlohmann/json.hpp>

namespace x3dna {
namespace config {

/**
 * @struct ParameterThresholds
 * @brief Validation and algorithm parameters (matches legacy miscPars)
 */
struct ParameterThresholds {
    // Distance constraints
    double min_dorg = 0.0;
    double max_dorg = 15.0;
    double min_dv = 0.0;
    double max_dv = 2.5;
    double min_dNN = 4.5;
    double max_dNN = 1e18; // XBIG

    // Angle constraints
    double min_plane_angle = 0.0;
    double max_plane_angle = 65.0;

    // Hydrogen bond constraints
    int min_base_hb = 1;
    double hb_lower = 1.8;
    double hb_dist1 = 4.0;
    double hb_dist2 = 0.0; // CRITICAL: Must be 0.0 for exact legacy match

    // H-bond atom list (default ".O.N" - matches legacy default)
    std::string hb_atoms = ".O.N";

    // Overlap threshold (matches legacy OVERLAP = 0.01)
    double overlap_threshold = 0.01;

    // Helix parameters
    double helix_break = 7.5;

    // Other parameters
    std::string alt_list = "A1";
    double std_curved = 0.6;
    double water_dist = 3.2;
    double water_dlow = 0.0;
    std::string water_atoms = ".O.N";
    double o3p_dist = 4.5;
};

/**
 * @class ConfigManager
 * @brief Singleton configuration manager
 */
class ConfigManager {
public:
    /**
     * @brief Get singleton instance
     */
    static ConfigManager& instance();

    // Configuration loading
    void load_from_file(const std::filesystem::path& config_path);
    void load_from_json(const nlohmann::json& json);
    void set_defaults();

    // Parameter thresholds
    const ParameterThresholds& thresholds() const {
        return thresholds_;
    }
    ParameterThresholds& thresholds() {
        return thresholds_;
    }

    // Paths
    void set_x3dna_home(const std::filesystem::path& path);
    std::filesystem::path x3dna_home() const {
        return x3dna_home_;
    }
    std::filesystem::path standard_base_path() const;

    // Options
    bool include_hetatm() const {
        return include_hetatm_;
    }
    void set_include_hetatm(bool value) {
        include_hetatm_ = value;
    }
    bool include_waters() const {
        return include_waters_;
    }
    void set_include_waters(bool value) {
        include_waters_ = value;
    }

    // Legacy mode (for exact compatibility with legacy code)
    // When enabled, breaks some OOP principles for exact matching
    bool legacy_mode() const {
        return legacy_mode_;
    }
    void set_legacy_mode(bool value) {
        legacy_mode_ = value;
    }

    // Delete copy and move
    ConfigManager(const ConfigManager&) = delete;
    ConfigManager& operator=(const ConfigManager&) = delete;
    ConfigManager(ConfigManager&&) = delete;
    ConfigManager& operator=(ConfigManager&&) = delete;

private:
    ConfigManager() = default;
    ~ConfigManager() = default;

    ParameterThresholds thresholds_;
    std::filesystem::path x3dna_home_;
    bool include_hetatm_ = false;
    bool include_waters_ = false;
    bool legacy_mode_ = false; // Enable legacy compatibility mode
};

} // namespace config
} // namespace x3dna
