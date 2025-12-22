#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <mutex>
#include <set>

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

// registry_ is a reference to lazy-loaded data
// Note: This maintains the same interface but delays loading until first access
const std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo>&
    ModifiedNucleotideRegistry::registry_ = get_registry();

std::optional<ModifiedNucleotideRegistry::NucleotideInfo> ModifiedNucleotideRegistry::get_info(
    const std::string& residue_name) {
    // Trim whitespace
    std::string trimmed = residue_name;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

    auto it = registry_.find(trimmed);
    if (it != registry_.end()) {
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

ResidueClassification ModifiedNucleotideRegistry::classify(const std::string& residue_name) {
    ResidueClassification result;
    result.residue_name = residue_name;

    // Check for water
    if (residue_name == "HOH" || residue_name == "WAT" || residue_name == "H2O") {
        result.molecule_type = MoleculeType::WATER;
        return result;
    }

    // Look up in nucleotide registry first (takes priority over ions)
    // This handles I (inosine) correctly instead of treating it as iodine
    auto info = get_info(residue_name);
    if (info) {
        result.molecule_type = MoleculeType::NUCLEIC_ACID;

        // Determine RNA vs DNA from residue name
        // DNA: starts with 'D' (DA, DC, DG, DT) or is T/THY
        bool is_dna = (residue_name.size() >= 2 && residue_name[0] == 'D') || residue_name == "T" ||
                      residue_name == "THY";
        result.nucleic_acid_type = is_dna ? NucleicAcidType::DNA : NucleicAcidType::RNA;

        // Map base type
        switch (info->base_type) {
            case ResidueType::ADENINE:
                result.base_type = BaseType::ADENINE;
                result.canonical_code = 'A';
                break;
            case ResidueType::GUANINE:
                result.base_type = BaseType::GUANINE;
                result.canonical_code = 'G';
                break;
            case ResidueType::CYTOSINE:
                result.base_type = BaseType::CYTOSINE;
                result.canonical_code = 'C';
                break;
            case ResidueType::THYMINE:
                result.base_type = BaseType::THYMINE;
                result.canonical_code = 'T';
                break;
            case ResidueType::URACIL:
                result.base_type = BaseType::URACIL;
                result.canonical_code = 'U';
                break;
            case ResidueType::INOSINE:
                result.base_type = BaseType::INOSINE;
                result.canonical_code = 'I';
                break;
            case ResidueType::PSEUDOURIDINE:
                result.base_type = BaseType::PSEUDOURIDINE;
                result.canonical_code = 'U';
                break;
            default:
                break;
        }

        // Set category
        result.base_category = get_base_category(result.base_type);

        // Determine if modified: lowercase code = modified
        result.is_modified = std::islower(static_cast<unsigned char>(info->one_letter_code));

        return result;
    }

    // Check for amino acids (common 3-letter codes)
    static const std::set<std::string> amino_acids = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                                                      "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                                                      "TYR", "VAL", "SEC", "PYL", "ASX", "GLX", "XLE", "UNK"};
    if (amino_acids.count(residue_name) > 0) {
        result.molecule_type = MoleculeType::PROTEIN;
        return result;
    }

    // Check for common ions (after nucleotide check to avoid conflicts with I=inosine)
    static const std::set<std::string> ions = {"MG", "NA", "K",  "CA", "ZN", "FE", "MN", "CL", "BR", "I",
                                               "CO", "NI", "CU", "CD", "BA", "SR", "RB", "CS", "LI"};
    if (ions.count(residue_name) > 0) {
        result.molecule_type = MoleculeType::ION;
        return result;
    }

    // Unknown - treat as ligand
    result.molecule_type = MoleculeType::LIGAND;
    return result;
}

} // namespace core
} // namespace x3dna
