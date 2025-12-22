/**
 * @file standard_base_templates.cpp
 * @brief Implementation of standard base template loader
 */

#include <x3dna/algorithms/standard_base_templates.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <filesystem>
#include <stdexcept>

namespace x3dna {
namespace algorithms {

StandardBaseTemplates::StandardBaseTemplates() {
    // Use ResourceLocator for centralized path management
    // ResourceLocator will auto-initialize from environment if needed
    template_path_ = config::ResourceLocator::templates_dir();
}

StandardBaseTemplates::StandardBaseTemplates(const std::filesystem::path& template_path)
    : template_path_(template_path) {
    if (!std::filesystem::exists(template_path_)) {
        throw std::runtime_error("Template path does not exist: " + template_path_.string());
    }
}

std::string StandardBaseTemplates::type_to_filename(core::ResidueType type, bool is_modified) {
    // Legacy: uppercase one_letter_code -> Atomic_X.pdb, lowercase -> Atomic.x.pdb
    // is_modified=true uses lowercase template (for modified nucleotides)
    char base_char;
    switch (type) {
        case core::ResidueType::ADENINE:
            base_char = 'a';
            break;
        case core::ResidueType::CYTOSINE:
            base_char = 'c';
            break;
        case core::ResidueType::GUANINE:
            base_char = 'g';
            break;
        case core::ResidueType::THYMINE:
            base_char = 't';
            break;
        case core::ResidueType::URACIL:
            base_char = 'u';
            break;
        case core::ResidueType::PSEUDOURIDINE:
            base_char = 'p';
            break;
        case core::ResidueType::INOSINE:
            base_char = 'i';
            break;
        default:
            throw std::invalid_argument("Invalid residue type for template loading");
    }

    if (is_modified) {
        // Modified nucleotide: Atomic.x.pdb (lowercase)
        return std::string("Atomic.") + base_char + ".pdb";
    } else {
        // Standard nucleotide: Atomic_X.pdb (uppercase)
        return std::string("Atomic_") + static_cast<char>(std::toupper(base_char)) + ".pdb";
    }
}

// Backwards compatible version for existing code
std::string StandardBaseTemplates::type_to_filename(core::ResidueType type) {
    return type_to_filename(type, false); // Default to standard (uppercase) template
}

std::filesystem::path StandardBaseTemplates::get_template_path(core::ResidueType type, bool is_modified) const {
    std::string filename = type_to_filename(type, is_modified);
    return template_path_ / filename;
}

std::filesystem::path StandardBaseTemplates::get_template_path(core::ResidueType type) const {
    return get_template_path(type, false); // Default to standard (uppercase) template
}

bool StandardBaseTemplates::template_exists(core::ResidueType type) const {
    std::filesystem::path template_file = get_template_path(type);
    return std::filesystem::exists(template_file) && std::filesystem::is_regular_file(template_file);
}

core::Structure StandardBaseTemplates::load_template(core::ResidueType type, bool is_modified) {
    // Create cache key that includes is_modified
    // Use a simple encoding: type * 2 + is_modified
    auto cache_key = static_cast<core::ResidueType>(static_cast<int>(type) * 2 + (is_modified ? 1 : 0));

    // Check cache first
    auto it = cache_.find(cache_key);
    if (it != cache_.end() && it->second) {
        return *(it->second);
    }

    // Get template file path
    std::filesystem::path template_file = get_template_path(type, is_modified);

    if (!std::filesystem::exists(template_file)) {
        throw std::runtime_error("Template file not found: " + template_file.string());
    }

    // Load template using PDB parser
    io::PdbParser parser;
    core::Structure template_structure = parser.parse_file(template_file);

    // Cache the template
    cache_[cache_key] = std::make_shared<core::Structure>(template_structure);

    return template_structure;
}

// Backwards compatible version
core::Structure StandardBaseTemplates::load_template(core::ResidueType type) {
    return load_template(type, false); // Default to standard (uppercase) template
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
