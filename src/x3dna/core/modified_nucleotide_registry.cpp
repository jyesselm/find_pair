#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
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

        // Load modified nucleotides from JSON
        std::filesystem::path config_file;
        try {
            config_file = config::ResourceLocator::config_file("modified_nucleotides.json");
        } catch (const std::runtime_error&) {
            // ResourceLocator not initialized - use only standard nucleotides
            return;
        }

        std::ifstream file(config_file);
        if (!file.is_open()) {
            return;
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
        } catch (const std::exception& e) {
            std::cerr << "ERROR loading modified_nucleotides.json: " << e.what() << "\n";
        }
    });

    return registry;
}

// For backward compatibility - REGISTRY is now a reference to lazy-loaded data
// Note: This maintains the same interface but delays loading until first access
const std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo>&
    ModifiedNucleotideRegistry::REGISTRY = get_registry();

// Fallback hardcoded registry (in case JSON file not found)
// This was the old approach - now just for reference/backup
static const std::map<std::string, ModifiedNucleotideRegistry::NucleotideInfo> FALLBACK_REGISTRY = {
    // ========== MODIFIED ADENINES ==========
    {"A2M", {'a', ResidueType::ADENINE, true, "2-methyladenosine"}},
    {"1MA", {'a', ResidueType::ADENINE, true, "1-methyladenosine"}},
    {"2MA", {'a', ResidueType::ADENINE, true, "2-methyladenosine"}},
    {"OMA", {'a', ResidueType::ADENINE, true, "O-methyladenosine"}},
    {"6MA", {'a', ResidueType::ADENINE, true, "N6-methyladenosine"}},
    {"MIA", {'a', ResidueType::ADENINE, true, "2-methylthio-N6-isopentenyladenosine"}},
    {"I6A", {'a', ResidueType::ADENINE, true, "N6-isopentenyladenosine"}},
    {"T6A", {'a', ResidueType::ADENINE, true, "N6-threonylcarbamoyladenosine"}},
    {"M6A", {'a', ResidueType::ADENINE, true, "N6-methyladenosine"}},
    {"LCA", {'a', ResidueType::ADENINE, true, "Locked adenosine (LNA)"}},
    {"A23", {'a', ResidueType::ADENINE, true, "2'-deoxy-2'-fluoroadenosine"}},
    {"ATP", {'a', ResidueType::ADENINE, true, "Adenosine triphosphate"}},
    {"ADP", {'a', ResidueType::ADENINE, true, "Adenosine diphosphate"}},
    {"AMP", {'a', ResidueType::ADENINE, true, "Adenosine monophosphate"}},
    {"SAM", {'a', ResidueType::ADENINE, true, "S-adenosylmethionine"}},
    {"SAH", {'a', ResidueType::ADENINE, true, "S-adenosylhomocysteine"}},
    {"APC", {'a', ResidueType::ADENINE, true, "Adenosine phosphate carrier"}},
    {"AF2", {'a', ResidueType::ADENINE, true, "Modified adenosine"}},
    {"0A", {'a', ResidueType::ADENINE, true, "Numbered adenine variant"}},

    // ========== MODIFIED GUANINES ==========
    {"OMG", {'g', ResidueType::GUANINE, true, "O-methylguanosine"}},
    {"1MG", {'g', ResidueType::GUANINE, true, "1-methylguanosine"}},
    {"2MG", {'g', ResidueType::GUANINE, true, "2-methylguanosine"}},
    {"7MG", {'g', ResidueType::GUANINE, true, "7-methylguanosine"}},
    {"M2G", {'g', ResidueType::GUANINE, true, "N2-methylguanosine"}},
    {"YYG", {'g', ResidueType::GUANINE, true, "Modified guanosine"}},
    {"YG", {'g', ResidueType::GUANINE, true, "Modified guanosine"}},
    {"QUO", {'g', ResidueType::GUANINE, true, "Queuosine"}},
    {"LCG", {'g', ResidueType::GUANINE, true, "Locked guanosine (LNA)"}},
    {"GTP", {'g', ResidueType::GUANINE, true, "Guanosine triphosphate"}},
    {"GDP", {'g', ResidueType::GUANINE, true, "Guanosine diphosphate"}},
    {"5GP", {'g', ResidueType::GUANINE, true, "5'-guanosine monophosphate"}},
    {"PRF", {'g', ResidueType::GUANINE, true, "Prefluorinated guanine"}},
    {"BGM", {'g', ResidueType::GUANINE, true, "Beta-D-glucosyl-guanine"}},
    {"G48", {'g', ResidueType::GUANINE, true, "Modified guanine"}},
    {"0G", {'g', ResidueType::GUANINE, true, "Numbered guanine variant"}},
    {"IGU", {'g', ResidueType::GUANINE, true, "Isoguanosine"}},
    {"DGP", {'g', ResidueType::GUANINE, true, "Deoxyguanosine phosphate"}},
    {"G7M", {'g', ResidueType::GUANINE, true, "N7-methyl-guanosine-5'-monophosphate"}},

    // ========== INOSINES ==========
    {"I", {'I', ResidueType::INOSINE, true, "Inosine"}},
    {"DI", {'I', ResidueType::INOSINE, true, "2'-deoxyinosine"}},
    {"IMP", {'I', ResidueType::INOSINE, true, "Inosine monophosphate"}},

    // ========== MODIFIED CYTOSINES ==========
    {"5MC", {'c', ResidueType::CYTOSINE, false, "5-methylcytidine"}},
    {"OMC", {'c', ResidueType::CYTOSINE, false, "O-methylcytidine"}},
    {"S4C", {'c', ResidueType::CYTOSINE, false, "4-thiocytidine"}},
    {"5IC", {'c', ResidueType::CYTOSINE, false, "5-iodocytidine"}},
    {"5FC", {'c', ResidueType::CYTOSINE, false, "5-fluorocytidine"}},
    {"CBR", {'c', ResidueType::CYTOSINE, false, "5-bromocytidine"}},
    {"LCC", {'c', ResidueType::CYTOSINE, false, "Locked cytosine (LNA)"}},
    {"EPE", {'c', ResidueType::CYTOSINE, false, "Modified cytosine analog"}},
    {"2YR", {'c', ResidueType::CYTOSINE, false, "2'-O-ribosylcytidine"}},
    {"CTP", {'c', ResidueType::CYTOSINE, false, "Cytidine triphosphate"}},
    {"CSL", {'c', ResidueType::CYTOSINE, false, "Modified cytosine"}},
    {"CBV", {'c', ResidueType::CYTOSINE, false, "Carbovir"}},
    {"CCC", {'c', ResidueType::CYTOSINE, false, "Modified cytosine"}},
    {"CFZ", {'c', ResidueType::CYTOSINE, false, "Modified cytosine"}},
    {"C43", {'c', ResidueType::CYTOSINE, false, "Modified cytosine"}},
    {"A6C", {'c', ResidueType::CYTOSINE, false, "N6-acetyl-cytosine"}},
    {"CFL", {'c', ResidueType::CYTOSINE, false, "Modified cytosine"}},
    {"RSQ", {'c', ResidueType::CYTOSINE, false, "R-squared cytosine"}},
    {"0C", {'c', ResidueType::CYTOSINE, false, "Numbered cytosine variant"}},

    // ========== MODIFIED THYMINES ==========
    {"5MU", {'t', ResidueType::THYMINE, false, "5-methyluridine (ribothymidine)"}},
    {"RT", {'t', ResidueType::THYMINE, false, "Ribothymidine"}},

    // ========== MODIFIED URACILS ==========
    {"70U", {'u', ResidueType::URACIL, false, "7-methyluridine"}},
    {"H2U", {'u', ResidueType::URACIL, false, "Dihydrouridine"}},
    {"DHU", {'u', ResidueType::URACIL, false, "Dihydrouridine"}},
    {"OMU", {'u', ResidueType::URACIL, false, "O-methyluridine"}},
    {"4SU", {'u', ResidueType::URACIL, false, "4-thiouridine"}},
    {"S4U", {'u', ResidueType::URACIL, false, "4-thiouridine"}},
    {"5BU", {'u', ResidueType::URACIL, false, "5-bromouridine"}},
    {"2MU", {'u', ResidueType::URACIL, false, "2-methyluridine"}},
    {"UR3", {'u', ResidueType::URACIL, false, "Modified uridine"}},
    {"TLN", {'u', ResidueType::URACIL, false, "Thymidine LNA"}},
    {"J48", {'u', ResidueType::URACIL, false, "Hypermodified nucleotide (wybutosine precursor)"}},
    {"NMN", {'u', ResidueType::URACIL, false, "Nicotinamide mononucleotide"}},
    {"NNR", {'u', ResidueType::URACIL, false, "Nicotinamide riboside"}},
    {"WVQ", {'u', ResidueType::URACIL, false, "Modified uracil"}},
    {"US5", {'u', ResidueType::URACIL, false, "5-hydroxymethyluridine"}},
    {"UTP", {'u', ResidueType::URACIL, false, "Uridine triphosphate"}},
    {"UFT", {'u', ResidueType::URACIL, false, "Modified uracil"}},
    {"0U", {'u', ResidueType::URACIL, false, "Numbered uracil variant"}},
    {"KIR", {'u', ResidueType::URACIL, false, "Kinetin riboside"}},

    // ========== PSEUDOURIDINES ==========
    {"B8H", {'P', ResidueType::PSEUDOURIDINE, false, "Pseudouridine derivative"}},
}; // End of FALLBACK_REGISTRY (not used if JSON loads successfully)

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
