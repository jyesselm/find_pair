/**
 * @file x3dna.hpp
 * @brief Main entry point for the x3dna library
 * @version 1.0.0
 *
 * This header provides the library initialization API and version information.
 * Include this header and call x3dna::init() before using any x3dna functionality.
 *
 * Example usage:
 * @code
 *   #include <x3dna/x3dna.hpp>
 *
 *   int main() {
 *       // Initialize with auto-detection
 *       if (!x3dna::init()) {
 *           std::cerr << "Failed to initialize x3dna\n";
 *           return 1;
 *       }
 *
 *       // Or initialize with explicit path
 *       x3dna::init("/path/to/resources");
 *
 *       // Use x3dna...
 *
 *       return 0;
 *   }
 * @endcode
 */

#pragma once

#include <x3dna/version.hpp>
#include <x3dna/forward_declarations.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <filesystem>

/**
 * @namespace x3dna
 * @brief Main namespace for the X3DNA library
 */
namespace x3dna {

/**
 * @brief Get library version string
 * @return Version string (e.g., "1.0.0")
 */
[[nodiscard]] inline const char* version() {
    return X3DNA_VERSION_STRING;
}

/**
 * @brief Initialize the x3dna library with explicit resources path
 *
 * Must be called before using any x3dna functionality.
 *
 * @param resources_path Path to the resources directory containing
 *                       templates/ and config/ subdirectories
 * @return true if initialization succeeded
 * @return false if resources path is invalid
 */
inline bool init(const std::filesystem::path& resources_path) {
    try {
        config::ResourceLocator::initialize(resources_path);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Initialize the x3dna library with auto-detection
 *
 * Attempts to find resources directory by searching:
 * 1. Common relative paths ("resources", "../resources", etc.)
 * 2. X3DNA_HOMEDIR environment variable
 * 3. X3DNA environment variable
 *
 * @return true if resources were found and initialization succeeded
 * @return false if resources could not be found
 */
inline bool init() {
    return config::ResourceLocator::initialize_from_environment();
}

/**
 * @brief Check if the library has been initialized
 *
 * @return true if init() was called successfully
 */
[[nodiscard]] inline bool is_initialized() {
    return config::ResourceLocator::is_initialized();
}

/**
 * @brief Reset the library to uninitialized state
 *
 * After calling this, init() must be called again before using the library.
 */
inline void shutdown() {
    config::ResourceLocator::reset();
}

/**
 * @brief Get the path to the resources directory
 *
 * @return Path to resources directory
 * @throws std::runtime_error if library not initialized
 */
[[nodiscard]] inline const std::filesystem::path& resources_path() {
    return config::ResourceLocator::resources_path();
}

} // namespace x3dna
