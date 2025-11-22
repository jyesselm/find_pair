/**
 * @file standard_base_templates.cpp
 * @brief Implementation of standard base template loader
 */

#include <x3dna/algorithms/standard_base_templates.hpp>
#include <filesystem>
#include <cstdlib>
#include <stdexcept>

namespace x3dna {
namespace algorithms {

StandardBaseTemplates::StandardBaseTemplates() {
    template_path_ = find_template_path();
    if (template_path_.empty()) {
        // Default to local data/templates directory
        template_path_ = std::filesystem::path("data/templates");
    }
}

StandardBaseTemplates::StandardBaseTemplates(const std::filesystem::path& template_path)
    : template_path_(template_path) {
    if (!std::filesystem::exists(template_path_)) {
        throw std::runtime_error("Template path does not exist: " + template_path_.string());
    }
}

std::filesystem::path StandardBaseTemplates::find_template_path() {
    // Check X3DNA_HOMEDIR environment variable
    const char* x3dna_home = std::getenv("X3DNA_HOMEDIR");
    if (x3dna_home) {
        std::filesystem::path config_path = std::filesystem::path(x3dna_home) / "config";
        if (std::filesystem::exists(config_path)) {
            return config_path;
        }
    }
    
    // Check local data/templates directory
    std::filesystem::path local_path = "data/templates";
    if (std::filesystem::exists(local_path)) {
        return local_path;
    }
    
    // Return empty path if not found
    return std::filesystem::path();
}

std::string StandardBaseTemplates::type_to_filename(core::ResidueType type) {
    switch (type) {
        case core::ResidueType::ADENINE:
            return "Atomic_A.pdb";
        case core::ResidueType::CYTOSINE:
            return "Atomic_C.pdb";
        case core::ResidueType::GUANINE:
            return "Atomic_G.pdb";
        case core::ResidueType::THYMINE:
            return "Atomic_T.pdb";
        case core::ResidueType::URACIL:
            return "Atomic_U.pdb";
        default:
            throw std::invalid_argument("Invalid residue type for template loading");
    }
}

std::filesystem::path StandardBaseTemplates::get_template_path(core::ResidueType type) const {
    std::string filename = type_to_filename(type);
    return template_path_ / filename;
}

bool StandardBaseTemplates::template_exists(core::ResidueType type) const {
    std::filesystem::path template_file = get_template_path(type);
    return std::filesystem::exists(template_file) && std::filesystem::is_regular_file(template_file);
}

core::Structure StandardBaseTemplates::load_template(core::ResidueType type) {
    // Check cache first
    auto it = cache_.find(type);
    if (it != cache_.end() && it->second) {
        return *(it->second);
    }
    
    // Get template file path
    std::filesystem::path template_file = get_template_path(type);
    
    if (!std::filesystem::exists(template_file)) {
        throw std::runtime_error("Template file not found: " + template_file.string());
    }
    
    // Load template using PDB parser
    io::PdbParser parser;
    core::Structure template_structure = parser.parse_file(template_file);
    
    // Cache the template
    cache_[type] = std::make_shared<core::Structure>(template_structure);
    
    return template_structure;
}

void StandardBaseTemplates::set_template_path(const std::filesystem::path& template_path) {
    if (!std::filesystem::exists(template_path)) {
        throw std::runtime_error("Template path does not exist: " + template_path.string());
    }
    template_path_ = template_path;
    // Clear cache when path changes
    clear_cache();
}

void StandardBaseTemplates::clear_cache() {
    cache_.clear();
}

bool StandardBaseTemplates::file_exists(const std::filesystem::path& path) {
    return std::filesystem::exists(path) && std::filesystem::is_regular_file(path);
}

} // namespace algorithms
} // namespace x3dna

