/**
 * @file resource_locator.hpp
 * @brief Centralized resource path management for portable library usage
 *
 * This class provides a single point of configuration for locating library
 * resources (templates, config files). It enables the library to be embedded
 * in external projects by allowing explicit resource path configuration.
 *
 * Usage:
 *   // Option 1: Explicit initialization (recommended for embedded use)
 *   x3dna::config::ResourceLocator::initialize("path/to/x3dna/resources");
 *
 *   // Option 2: Auto-detect from environment (fallback)
 *   // Uses X3DNA_HOMEDIR or X3DNA environment variables
 */

#pragma once

#include <filesystem>
#include <optional>
#include <string>

namespace x3dna {
namespace config {

/**
 * @class ResourceLocator
 * @brief Singleton for centralized resource path management
 *
 * This replaces the scattered path-finding logic throughout the codebase
 * with a single, configurable source of truth.
 */
class ResourceLocator {
public:
    /**
     * @brief Initialize the resource locator with explicit path
     *
     * This should be called once at application startup before using
     * any library functionality. The path should point to the directory
     * containing 'templates/' and 'config/' subdirectories.
     *
     * @param resources_path Path to resources directory
     * @throws std::runtime_error if path doesn't exist or is invalid
     */
    static void initialize(const std::filesystem::path& resources_path);

    /**
     * @brief Initialize with auto-detection (fallback behavior)
     *
     * Searches in this order:
     * 1. Common relative paths from current working directory
     * 2. X3DNA_HOMEDIR environment variable
     * 3. X3DNA environment variable
     *
     * @return true if resources were found, false otherwise
     */
    static bool initialize_from_environment();

    /**
     * @brief Reset to uninitialized state (primarily for testing)
     */
    static void reset();

    /**
     * @brief Check if resources have been located
     */
    [[nodiscard]] static bool is_initialized();

    /**
     * @brief Get the base resources directory
     * @throws std::runtime_error if not initialized
     */
    [[nodiscard]] static const std::filesystem::path& resources_path();

    /**
     * @brief Get the templates directory (resources/templates/)
     * @throws std::runtime_error if not initialized
     */
    [[nodiscard]] static std::filesystem::path templates_dir();

    /**
     * @brief Get the config directory (resources/config/)
     * @throws std::runtime_error if not initialized
     */
    [[nodiscard]] static std::filesystem::path config_dir();

    /**
     * @brief Get path to a specific template file
     * @param filename Template filename (e.g., "Atomic_A.pdb")
     * @throws std::runtime_error if not initialized
     */
    [[nodiscard]] static std::filesystem::path template_file(const std::string& filename);

    /**
     * @brief Get path to a specific config file
     * @param filename Config filename (e.g., "atomlist.dat")
     * @throws std::runtime_error if not initialized
     */
    [[nodiscard]] static std::filesystem::path config_file(const std::string& filename);

    /**
     * @brief Check if a template file exists
     * @param filename Template filename
     */
    [[nodiscard]] static bool template_exists(const std::string& filename);

    /**
     * @brief Check if a config file exists
     * @param filename Config filename
     */
    [[nodiscard]] static bool config_exists(const std::string& filename);

    // Delete copy/move - singleton pattern
    ResourceLocator(const ResourceLocator&) = delete;
    ResourceLocator& operator=(const ResourceLocator&) = delete;
    ResourceLocator(ResourceLocator&&) = delete;
    ResourceLocator& operator=(ResourceLocator&&) = delete;

private:
    ResourceLocator() = default;
    ~ResourceLocator() = default;

    static ResourceLocator& instance();

    /**
     * @brief Validate that a path contains expected resource structure
     * @param path Path to validate
     * @return true if path contains templates/ and config/ subdirs
     */
    [[nodiscard]] static bool validate_resources_path(const std::filesystem::path& path);

    /**
     * @brief Search for resources in common locations
     * @return Path if found, empty optional otherwise
     */
    [[nodiscard]] static std::optional<std::filesystem::path> find_resources_auto();

    /**
     * @brief Throw if not initialized
     */
    static void require_initialized();

    std::filesystem::path resources_path_;
    bool initialized_ = false;
};

} // namespace config
} // namespace x3dna
