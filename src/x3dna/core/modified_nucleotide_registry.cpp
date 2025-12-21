#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <mutex>
#include <set>
#include <iostream>

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

// Thread-safe lazy-loaded registry singleton
// Auto-initializes ResourceLocator if not already initialized
static std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo>& get_registry() {
    static std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo> registry;
    static std::once_flag init_flag;

    std::call_once(init_flag, []() {
        // config_file() will auto-initialize ResourceLocator if possible
        // (searches relative paths and environment variables)
        std::filesystem::path config_file = config::ResourceLocator::config_file("modified_nucleotides.json");

        std::ifstream file(config_file);
        if (!file.is_open()) {
            throw std::runtime_error("ModifiedNucleotideRegistry: Cannot open config file: " + config_file.string() +
                                     ". Ensure the resources directory contains modified_nucleotides.json");
        }

        try {
            json j = json::parse(file);

            // Load all categories (including standard_nucleotides)
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
            throw std::runtime_error("ModifiedNucleotideRegistry: Error parsing modified_nucleotides.json: " +
                                     std::string(e.what()));
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
    // Residue names are stored trimmed - direct lookup
    auto it = REGISTRY.find(residue_name);
    if (it != REGISTRY.end()) {
        return it->second;
    }
    return std::nullopt;
}

// Track which residues we've already warned about to avoid spam
static std::set<std::string>& get_warned_residues() {
    static std::set<std::string> warned;
    return warned;
}

// Check if residue name is a known non-nucleotide (water, ion, ligand, amino acid)
// These should not trigger warnings since they're not expected to be in the registry
static bool is_known_non_nucleotide(const std::string& name) {
    // Water molecules
    static const std::set<std::string> WATER = {"HOH", "WAT", "DOD", "H2O", "OH2"};
    if (WATER.count(name)) return true;

    // Common ions
    static const std::set<std::string> IONS = {
        "MG", "CA", "NA", "K", "CL", "ZN", "FE", "MN", "CO", "NI", "CU",
        "CD", "HG", "PB", "SR", "BA", "LI", "RB", "CS", "BR", "I", "F",
        "MG2", "CA2", "ZN2", "FE2", "FE3", "MN2", "CO2", "NI2", "CU2", "CU1"
    };
    if (IONS.count(name)) return true;

    // Standard amino acids (3-letter codes)
    static const std::set<std::string> AMINO_ACIDS = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        // Common modified amino acids
        "MSE", "SEP", "TPO", "PTR", "CSO", "OCS", "CME", "CSD", "CSX", "MLY"
    };
    if (AMINO_ACIDS.count(name)) return true;

    // Common small molecule ligands/cofactors (not nucleotide-derived)
    static const std::set<std::string> LIGANDS = {
        "SO4", "PO4", "GOL", "EDO", "ACT", "ACE", "NH4", "NO3", "CO3",
        "CIT", "TRS", "BME", "DMS", "IMD", "EPE", "PEG", "MPD", "HEM"
    };
    if (LIGANDS.count(name)) return true;

    return false;
}

char ModifiedNucleotideRegistry::get_one_letter_code(const std::string& residue_name) {
    auto info = get_info(residue_name);
    if (info.has_value()) {
        return info->one_letter_code;
    }

    // Only warn for residues that might be nucleotides (not water, ions, amino acids, etc.)
    if (!is_known_non_nucleotide(residue_name)) {
        auto& warned = get_warned_residues();
        if (warned.find(residue_name) == warned.end()) {
            warned.insert(residue_name);
            std::cerr << "Warning: Unknown residue '" << residue_name
                      << "' not found in modified_nucleotides.json registry\n";
        }
    }
    return '?';
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
