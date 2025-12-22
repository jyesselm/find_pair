/**
 * @file type_registry.cpp
 * @brief Implementation of TypeRegistry singleton
 */

#include "x3dna/core/typing/type_registry.hpp"
#include "x3dna/config/resource_locator.hpp"

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
namespace typing {

namespace {

// Track which residues we've already warned about to avoid spam
std::set<std::string>& get_warned_residues() {
    static std::set<std::string> warned;
    return warned;
}

// Helper to convert string to BaseType
BaseType string_to_base_type(const std::string& type_str) {
    if (type_str == "ADENINE") return BaseType::ADENINE;
    if (type_str == "CYTOSINE") return BaseType::CYTOSINE;
    if (type_str == "GUANINE") return BaseType::GUANINE;
    if (type_str == "THYMINE") return BaseType::THYMINE;
    if (type_str == "URACIL") return BaseType::URACIL;
    if (type_str == "INOSINE") return BaseType::INOSINE;
    if (type_str == "PSEUDOURIDINE") return BaseType::PSEUDOURIDINE;
    return BaseType::UNKNOWN;
}

} // anonymous namespace

const TypeRegistry& TypeRegistry::instance() {
    static TypeRegistry instance;
    return instance;
}

TypeRegistry::TypeRegistry() {
    load_nucleotides();
    load_amino_acids();
    load_waters();
    load_ions();
}

void TypeRegistry::load_nucleotides() {
    try {
        std::filesystem::path config_file = config::ResourceLocator::config_file("modified_nucleotides.json");
        std::ifstream file(config_file);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open modified_nucleotides.json");
        }

        json j = json::parse(file);

        for (const auto& [category, nucleotides] : j["modified_nucleotides"].items()) {
            for (const auto& [name, info] : nucleotides.items()) {
                NucleotideInfo nucinfo;

                std::string code_str = info["code"].get<std::string>();
                nucinfo.one_letter_code = code_str.empty() ? '?' : code_str[0];
                nucinfo.base_type = string_to_base_type(info["type"].get<std::string>());
                nucinfo.is_purine = info["is_purine"].get<bool>();
                nucinfo.description = info["description"].get<std::string>();
                nucinfo.is_modified = std::islower(static_cast<unsigned char>(nucinfo.one_letter_code));

                nucleotides_[name] = nucinfo;
            }
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("TypeRegistry: Error loading nucleotides: " + std::string(e.what()));
    }
}

void TypeRegistry::load_amino_acids() {
    // Standard 20 amino acids
    struct AAData {
        const char* name;
        char code;
        AminoAcidType type;
        AminoAcidCategory category;
    };

    static const AAData standard_aas[] = {
        {"ALA", 'A', AminoAcidType::ALA, AminoAcidCategory::HYDROPHOBIC},
        {"ARG", 'R', AminoAcidType::ARG, AminoAcidCategory::POSITIVE},
        {"ASN", 'N', AminoAcidType::ASN, AminoAcidCategory::POLAR},
        {"ASP", 'D', AminoAcidType::ASP, AminoAcidCategory::NEGATIVE},
        {"CYS", 'C', AminoAcidType::CYS, AminoAcidCategory::POLAR},
        {"GLN", 'Q', AminoAcidType::GLN, AminoAcidCategory::POLAR},
        {"GLU", 'E', AminoAcidType::GLU, AminoAcidCategory::NEGATIVE},
        {"GLY", 'G', AminoAcidType::GLY, AminoAcidCategory::HYDROPHOBIC},
        {"HIS", 'H', AminoAcidType::HIS, AminoAcidCategory::POSITIVE},
        {"ILE", 'I', AminoAcidType::ILE, AminoAcidCategory::HYDROPHOBIC},
        {"LEU", 'L', AminoAcidType::LEU, AminoAcidCategory::HYDROPHOBIC},
        {"LYS", 'K', AminoAcidType::LYS, AminoAcidCategory::POSITIVE},
        {"MET", 'M', AminoAcidType::MET, AminoAcidCategory::HYDROPHOBIC},
        {"PHE", 'F', AminoAcidType::PHE, AminoAcidCategory::HYDROPHOBIC},
        {"PRO", 'P', AminoAcidType::PRO, AminoAcidCategory::HYDROPHOBIC},
        {"SER", 'S', AminoAcidType::SER, AminoAcidCategory::POLAR},
        {"THR", 'T', AminoAcidType::THR, AminoAcidCategory::POLAR},
        {"TRP", 'W', AminoAcidType::TRP, AminoAcidCategory::HYDROPHOBIC},
        {"TYR", 'Y', AminoAcidType::TYR, AminoAcidCategory::POLAR},
        {"VAL", 'V', AminoAcidType::VAL, AminoAcidCategory::HYDROPHOBIC},
        // Non-standard
        {"SEC", 'U', AminoAcidType::SEC, AminoAcidCategory::POLAR},
        {"PYL", 'O', AminoAcidType::PYL, AminoAcidCategory::POSITIVE},
        // Ambiguous
        {"ASX", 'B', AminoAcidType::ASX, AminoAcidCategory::UNKNOWN},
        {"GLX", 'Z', AminoAcidType::GLX, AminoAcidCategory::UNKNOWN},
        {"XLE", 'J', AminoAcidType::XLE, AminoAcidCategory::HYDROPHOBIC},
        {"UNK", 'X', AminoAcidType::UNK, AminoAcidCategory::UNKNOWN},
    };

    for (const auto& aa : standard_aas) {
        AminoAcidInfo info;
        info.one_letter_code = aa.code;
        info.type = aa.type;
        info.category = aa.category;
        info.is_modified = false;
        amino_acids_[aa.name] = info;
    }
}

void TypeRegistry::load_waters() {
    water_names_ = {"HOH", "WAT", "H2O", "OH2", "SOL", "DOD"};
}

void TypeRegistry::load_ions() {
    // Map residue names to ion types
    ion_types_ = {
        // Alkali metals
        {"LI", IonType::LITHIUM},
        {"NA", IonType::SODIUM},
        {"K", IonType::POTASSIUM},
        {"RB", IonType::RUBIDIUM},
        {"CS", IonType::CESIUM},
        // Alkaline earth
        {"MG", IonType::MAGNESIUM},
        {"CA", IonType::CALCIUM},
        {"SR", IonType::STRONTIUM},
        {"BA", IonType::BARIUM},
        // Transition metals
        {"MN", IonType::MANGANESE},
        {"FE", IonType::IRON},
        {"CO", IonType::COBALT},
        {"NI", IonType::NICKEL},
        {"CU", IonType::COPPER},
        {"ZN", IonType::ZINC},
        {"CD", IonType::CADMIUM},
        // Halogens
        {"F", IonType::FLUORIDE},
        {"CL", IonType::CHLORIDE},
        {"BR", IonType::BROMIDE},
        // Note: I (iodide) conflicts with inosine, handled specially
    };
}

ResidueClassification TypeRegistry::classify_residue(const std::string& residue_name) const {
    // Residue names are stored trimmed - direct use
    ResidueClassification result;
    result.residue_name = residue_name;

    // Check for water first
    if (water_names_.count(residue_name) > 0) {
        result.molecule_type = MoleculeType::WATER;
        result.solvent_type = SolventType::WATER;
        return result;
    }

    // Check nucleotides (takes priority over ions to handle I=inosine correctly)
    auto nuc_it = nucleotides_.find(residue_name);
    if (nuc_it != nucleotides_.end()) {
        const auto& info = nuc_it->second;
        result.molecule_type = MoleculeType::NUCLEIC_ACID;
        result.one_letter_code = info.one_letter_code;
        result.base_type = info.base_type;
        result.is_modified_nucleotide = info.is_modified;

        // Determine RNA vs DNA
        bool is_dna = (residue_name.size() >= 2 && residue_name[0] == 'D') ||
                      residue_name == "T" || residue_name == "THY";
        result.nucleic_acid_type = is_dna ? NucleicAcidType::DNA : NucleicAcidType::RNA;

        // Set category
        result.base_category = info.is_purine ? BaseCategory::PURINE : BaseCategory::PYRIMIDINE;

        // Set canonical code
        switch (info.base_type) {
            case BaseType::ADENINE: result.canonical_code = 'A'; break;
            case BaseType::GUANINE: result.canonical_code = 'G'; break;
            case BaseType::CYTOSINE: result.canonical_code = 'C'; break;
            case BaseType::THYMINE: result.canonical_code = 'T'; break;
            case BaseType::URACIL: result.canonical_code = 'U'; break;
            case BaseType::INOSINE: result.canonical_code = 'I'; break;
            case BaseType::PSEUDOURIDINE: result.canonical_code = 'U'; break;
            default: result.canonical_code = '?'; break;
        }

        return result;
    }

    // Check amino acids
    auto aa_it = amino_acids_.find(residue_name);
    if (aa_it != amino_acids_.end()) {
        const auto& info = aa_it->second;
        result.molecule_type = MoleculeType::PROTEIN;
        result.amino_acid_type = info.type;
        result.amino_acid_category = info.category;
        result.one_letter_code = info.one_letter_code;
        result.canonical_code = info.one_letter_code;
        result.is_modified_amino_acid = info.is_modified;
        return result;
    }

    // Check ions
    auto ion_it = ion_types_.find(residue_name);
    if (ion_it != ion_types_.end()) {
        result.molecule_type = MoleculeType::ION;
        result.ion_type = ion_it->second;
        return result;
    }

    // Unknown - treat as ligand
    result.molecule_type = MoleculeType::LIGAND;
    return result;
}

bool TypeRegistry::is_water(const std::string& residue_name) const {
    return water_names_.count(residue_name) > 0;
}

bool TypeRegistry::is_ion(const std::string& residue_name) const {
    return ion_types_.count(residue_name) > 0;
}

bool TypeRegistry::is_amino_acid(const std::string& residue_name) const {
    return amino_acids_.count(residue_name) > 0;
}

bool TypeRegistry::is_nucleotide(const std::string& residue_name) const {
    return nucleotides_.count(residue_name) > 0;
}

std::optional<NucleotideInfo> TypeRegistry::get_nucleotide_info(const std::string& residue_name) const {
    auto it = nucleotides_.find(residue_name);
    if (it != nucleotides_.end()) {
        return it->second;
    }
    return std::nullopt;
}

char TypeRegistry::get_one_letter_code(const std::string& residue_name) const {
    // Residue names are stored trimmed - direct lookup
    auto nuc_it = nucleotides_.find(residue_name);
    if (nuc_it != nucleotides_.end()) {
        return nuc_it->second.one_letter_code;
    }

    auto aa_it = amino_acids_.find(residue_name);
    if (aa_it != amino_acids_.end()) {
        return aa_it->second.one_letter_code;
    }

    // Only warn for residues that might be nucleotides (not water, ions, amino acids, ligands)
    // Check if it's a known non-nucleotide type
    if (!is_water(residue_name) && !is_ion(residue_name) && !is_amino_acid(residue_name)) {
        auto& warned = get_warned_residues();
        if (warned.find(residue_name) == warned.end()) {
            warned.insert(residue_name);
            std::cerr << "Warning: Unknown residue '" << residue_name
                      << "' not found in modified_nucleotides.json registry\n";
        }
    }

    return '?';
}

std::optional<bool> TypeRegistry::is_purine(const std::string& residue_name) const {
    auto it = nucleotides_.find(residue_name);
    if (it != nucleotides_.end()) {
        return it->second.is_purine;
    }
    return std::nullopt;
}

std::optional<AminoAcidInfo> TypeRegistry::get_amino_acid_info(const std::string& residue_name) const {
    auto it = amino_acids_.find(residue_name);
    if (it != amino_acids_.end()) {
        return it->second;
    }
    return std::nullopt;
}

IonType TypeRegistry::get_ion_type(const std::string& residue_name) const {
    auto it = ion_types_.find(residue_name);
    if (it != ion_types_.end()) {
        return it->second;
    }
    return IonType::UNKNOWN;
}

} // namespace typing
} // namespace core
} // namespace x3dna
