/**
 * @file hbond_parameters_loader.hpp
 * @brief Loader for H-bond parameters from JSON configuration
 */

#pragma once

#include <x3dna/config/hbond_parameters.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <string>
#include <optional>

namespace x3dna {
namespace config {

/**
 * @class HBondParametersLoader
 * @brief Loads and manages H-bond parameters from JSON configuration
 *
 * Provides static methods to load parameters from:
 * - Default config file (resources/config/hbond_parameters.json)
 * - Custom file path
 * - JSON object
 * - Named presets
 *
 * Example usage:
 * @code
 *   // Load defaults
 *   auto params = HBondParametersLoader::load();
 *
 *   // Load a preset
 *   auto legacy = HBondParametersLoader::load_preset("legacy_compatible");
 *
 *   // Load from custom file
 *   auto custom = HBondParametersLoader::load_from_file("my_config.json");
 * @endcode
 */
class HBondParametersLoader {
public:
    /**
     * @brief Load parameters from default config file
     * @return Loaded parameters, or defaults if file not found
     */
    static HBondParameters load();

    /**
     * @brief Load parameters from a specific file
     * @param path Path to JSON configuration file
     * @return Loaded parameters
     * @throws std::runtime_error if file cannot be read or parsed
     */
    static HBondParameters load_from_file(const std::filesystem::path& path);

    /**
     * @brief Load parameters from JSON object
     * @param json JSON object containing parameters
     * @return Loaded parameters
     */
    static HBondParameters load_from_json(const nlohmann::json& json);

    /**
     * @brief Load a named preset with overrides applied to defaults
     * @param preset_name Name of the preset (e.g., "legacy_compatible", "modern")
     * @return Parameters with preset overrides applied
     * @throws std::runtime_error if preset not found
     */
    static HBondParameters load_preset(const std::string& preset_name);

    /**
     * @brief Get singleton instance of loaded parameters
     * @return Reference to global parameters (loaded once on first call)
     */
    static const HBondParameters& instance();

    /**
     * @brief Reload singleton instance from file
     *
     * Forces a reload of the singleton parameters. Useful for testing
     * or when config file has changed.
     */
    static void reload();

    /**
     * @brief Get list of available preset names
     * @return Vector of preset names
     */
    static std::vector<std::string> available_presets();

    /**
     * @brief Check if a preset exists
     * @param name Preset name to check
     * @return true if preset exists
     */
    static bool has_preset(const std::string& name);

private:
    /**
     * @brief Apply preset overrides to base parameters
     * @param base Base parameters to modify
     * @param preset_json JSON object containing preset overrides
     */
    static void apply_preset(HBondParameters& base, const nlohmann::json& preset_json);

    /**
     * @brief Load detection config from JSON
     */
    static void load_detection(HBondDetectionConfig& config, const nlohmann::json& json);

    /**
     * @brief Load geometry config from JSON
     */
    static void load_geometry(HBondGeometryConfig& config, const nlohmann::json& json);

    /**
     * @brief Load scoring config from JSON
     */
    static void load_scoring(HBondScoringConfig& config, const nlohmann::json& json);

    /**
     * @brief Load quality tiers from JSON
     */
    static void load_quality_tiers(QualityTiersConfig& config, const nlohmann::json& json);

    /**
     * @brief Get path to default config file
     */
    static std::filesystem::path default_config_path();

    // Singleton storage
    static std::optional<HBondParameters> cached_params_;
    static nlohmann::json cached_json_;
};

} // namespace config
} // namespace x3dna
