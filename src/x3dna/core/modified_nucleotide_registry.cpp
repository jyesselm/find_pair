#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <mutex>

using json = nlohmann::json;

namespace x3dna {
namespace core {

// Helper to convert string to ResidueType
static ResidueType string_to_residue_type(const std::string& type_str) {
    if (type_str == "ADENINE")
        return ResidueType::ADENINE;
    if (type_str == "CYTOSINE")
        return ResidueType::CYTOSINE;
    if (type_str == "GUANINE")
        return ResidueType::GUANINE;
    if (type_str == "THYMINE")
        return ResidueType::THYMINE;
    if (type_str == "URACIL")
        return ResidueType::URACIL;
    if (type_str == "INOSINE")
        return ResidueType::INOSINE;
    if (type_str == "PSEUDOURIDINE")
        return ResidueType::PSEUDOURIDINE;
    throw std::runtime_error("Unknown residue type: " + type_str);
}

// Standard nucleotides - always available, no file loading needed
static const std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo> STANDARD_NUCLEOTIDES = {
    // RNA nucleotides
    {"A",   {'A', ResidueType::ADENINE, true, "Adenosine"}},
    {"ADE", {'A', ResidueType::ADENINE, true, "Adenosine"}},
    {"C",   {'C', ResidueType::CYTOSINE, false, "Cytidine"}},
    {"CYT", {'C', ResidueType::CYTOSINE, false, "Cytidine"}},
    {"G",   {'G', ResidueType::GUANINE, true, "Guanosine"}},
    {"GUA", {'G', ResidueType::GUANINE, true, "Guanosine"}},
    {"U",   {'U', ResidueType::URACIL, false, "Uridine"}},
    {"URA", {'U', ResidueType::URACIL, false, "Uridine"}},
    // DNA nucleotides
    {"DA",  {'A', ResidueType::ADENINE, true, "Deoxyadenosine"}},
    {"DC",  {'C', ResidueType::CYTOSINE, false, "Deoxycytidine"}},
    {"DG",  {'G', ResidueType::GUANINE, true, "Deoxyguanosine"}},
    {"DT",  {'T', ResidueType::THYMINE, false, "Deoxythymidine"}},
    {"T",   {'T', ResidueType::THYMINE, false, "Thymidine"}},
    {"THY", {'T', ResidueType::THYMINE, false, "Thymidine"}},
    {"DU",  {'U', ResidueType::URACIL, false, "Deoxyuridine"}},
    // Special nucleotides
    {"I",   {'I', ResidueType::INOSINE, true, "Inosine"}},
    {"INO", {'I', ResidueType::INOSINE, true, "Inosine"}},
    {"P",   {'P', ResidueType::PSEUDOURIDINE, false, "Pseudouridine"}},
    {"PSU", {'P', ResidueType::PSEUDOURIDINE, false, "Pseudouridine"}},
};

// Thread-safe lazy-loaded registry singleton
static std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo>& get_registry() {
    static std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo> registry;
    static std::once_flag init_flag;

    std::call_once(init_flag, []() {
        // Start with standard nucleotides
        registry = STANDARD_NUCLEOTIDES;

        // ResourceLocator MUST be initialized before using the registry
        if (!config::ResourceLocator::is_initialized()) {
            throw std::runtime_error(
                "ModifiedNucleotideRegistry: ResourceLocator not initialized. "
                "Call x3dna::config::ResourceLocator::set_resources_path() before using the library.");
        }

        // Load modified nucleotides from JSON - file MUST exist
        std::filesystem::path config_file = config::ResourceLocator::config_file("modified_nucleotides.json");

        std::ifstream file(config_file);
        if (!file.is_open()) {
            throw std::runtime_error(
                "ModifiedNucleotideRegistry: Cannot open config file: " + config_file.string() +
                ". Ensure the resources directory contains modified_nucleotides.json");
        }

        try {
            json j = json::parse(file);

            // Load all categories
            for (const auto& [category, nucleotides] : j["modified_nucleotides"].items()) {
                for (const auto& [name, info] : nucleotides.items()) {
                    ModifiedNucleotideRegistry::NucleotideInfo nucinfo;

                    std::string code_str = info["code"].get<std::string>();
                    nucinfo.one_letter_code = code_str.empty() ? '?' : code_str[0];
                    nucinfo.base_type = string_to_residue_type(info["type"].get<std::string>());
                    nucinfo.is_purine = info["is_purine"].get<bool>();
                    nucinfo.description = info["description"].get<std::string>();

                    registry[name] = nucinfo;
                }
            }
        } catch (const json::exception& e) {
            throw std::runtime_error(
                "ModifiedNucleotideRegistry: Error parsing modified_nucleotides.json: " + std::string(e.what()));
        }
    });

    return registry;
}

// For backward compatibility - REGISTRY is now a reference to lazy-loaded data
// Note: This maintains the same interface but delays loading until first access
const std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo>&
    ModifiedNucleotideRegistry::REGISTRY = get_registry();

std::optional<ModifiedNucleotideRegistry::NucleotideInfo> ModifiedNucleotideRegistry::get_info(
    const std::string& residue_name) {
    // Trim whitespace
    std::string trimmed = residue_name;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

    auto it = REGISTRY.find(trimmed);
    if (it != REGISTRY.end()) {
        return it->second;
    }
    return std::nullopt;
}

char ModifiedNucleotideRegistry::get_one_letter_code(const std::string& residue_name) {
    auto info = get_info(residue_name);
    return info.has_value() ? info->one_letter_code : '?';
}

std::optional<ResidueType> ModifiedNucleotideRegistry::get_base_type(const std::string& residue_name) {
    auto info = get_info(residue_name);
    return info.has_value() ? std::optional<ResidueType>(info->base_type) : std::nullopt;
}

std::optional<bool> ModifiedNucleotideRegistry::is_purine(const std::string& residue_name) {
    auto info = get_info(residue_name);
    return info.has_value() ? std::optional<bool>(info->is_purine) : std::nullopt;
}

bool ModifiedNucleotideRegistry::contains(const std::string& residue_name) {
    return get_info(residue_name).has_value();
}

} // namespace core
} // namespace x3dna
