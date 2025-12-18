/**
 * @file resource_locator.cpp
 * @brief Implementation of centralized resource path management
 */

#include <x3dna/config/resource_locator.hpp>
#include <cstdlib>
#include <stdexcept>

namespace x3dna {
namespace config {

ResourceLocator& ResourceLocator::instance() {
    static ResourceLocator inst;
    return inst;
}

void ResourceLocator::initialize(const std::filesystem::path& resources_path) {
    if (!std::filesystem::exists(resources_path)) {
        throw std::runtime_error("ResourceLocator: Path does not exist: " + resources_path.string());
    }

    if (!validate_resources_path(resources_path)) {
        throw std::runtime_error("ResourceLocator: Invalid resources directory (missing templates/ or config/): " +
                                 resources_path.string());
    }

    auto& inst = instance();
    inst.resources_path_ = std::filesystem::canonical(resources_path);
    inst.initialized_ = true;
}

bool ResourceLocator::initialize_from_environment() {
    auto found_path = find_resources_auto();
    if (found_path.has_value()) {
        auto& inst = instance();
        inst.resources_path_ = std::filesystem::canonical(found_path.value());
        inst.initialized_ = true;
        return true;
    }
    return false;
}

void ResourceLocator::reset() {
    auto& inst = instance();
    inst.resources_path_.clear();
    inst.initialized_ = false;
}

bool ResourceLocator::is_initialized() {
    return instance().initialized_;
}

const std::filesystem::path& ResourceLocator::resources_path() {
    require_initialized();
    return instance().resources_path_;
}

std::filesystem::path ResourceLocator::templates_dir() {
    require_initialized();
    return instance().resources_path_ / "templates";
}

std::filesystem::path ResourceLocator::config_dir() {
    require_initialized();
    return instance().resources_path_ / "config";
}

std::filesystem::path ResourceLocator::template_file(const std::string& filename) {
    return templates_dir() / filename;
}

std::filesystem::path ResourceLocator::config_file(const std::string& filename) {
    return config_dir() / filename;
}

bool ResourceLocator::template_exists(const std::string& filename) {
    if (!is_initialized()) {
        return false;
    }
    return std::filesystem::exists(template_file(filename));
}

bool ResourceLocator::config_exists(const std::string& filename) {
    if (!is_initialized()) {
        return false;
    }
    return std::filesystem::exists(config_file(filename));
}

bool ResourceLocator::validate_resources_path(const std::filesystem::path& path) {
    // Check for expected subdirectories
    const auto templates = path / "templates";
    const auto config = path / "config";

    return std::filesystem::exists(templates) && std::filesystem::is_directory(templates) &&
           std::filesystem::exists(config) && std::filesystem::is_directory(config);
}

std::optional<std::filesystem::path> ResourceLocator::find_resources_auto() {
    // Priority 1: Common relative paths from current working directory
    const std::vector<std::filesystem::path> search_paths = {
        "resources",
        "../resources",
        "../../resources",
        "../../../resources",
    };

    for (const auto& path : search_paths) {
        if (std::filesystem::exists(path) && validate_resources_path(path)) {
            return path;
        }
    }

    // Priority 2: X3DNA_HOMEDIR environment variable
    if (const char* home_dir = std::getenv("X3DNA_HOMEDIR")) {
        std::filesystem::path home_path(home_dir);
        // Check if it's the resources directory itself
        if (validate_resources_path(home_path)) {
            return home_path;
        }
        // Check if resources is a subdirectory
        auto resources_subdir = home_path / "resources";
        if (std::filesystem::exists(resources_subdir) && validate_resources_path(resources_subdir)) {
            return resources_subdir;
        }
    }

    // Priority 3: X3DNA environment variable (legacy)
    if (const char* x3dna = std::getenv("X3DNA")) {
        std::filesystem::path x3dna_path(x3dna);
        // Legacy X3DNA has config/ directly (templates are in config/)
        // Create a mapping: config/ -> templates/, config/ -> config/
        if (std::filesystem::exists(x3dna_path / "config")) {
            // For legacy X3DNA, the config directory serves as both templates and config
            // We'll return the x3dna path and handle the mapping in the calling code
            // Actually, for simplicity, if they have the new structure, use it
            auto resources_subdir = x3dna_path / "resources";
            if (std::filesystem::exists(resources_subdir) && validate_resources_path(resources_subdir)) {
                return resources_subdir;
            }
        }
    }

    return std::nullopt;
}

void ResourceLocator::require_initialized() {
    if (!instance().initialized_) {
        // Try auto-initialization as a convenience
        if (!initialize_from_environment()) {
            throw std::runtime_error("ResourceLocator: Not initialized. Call ResourceLocator::initialize() "
                                     "or set X3DNA_HOMEDIR environment variable before using x3dna library.");
        }
    }
}

} // namespace config
} // namespace x3dna
