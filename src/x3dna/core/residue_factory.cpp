#include <x3dna/core/residue_factory.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace core {

// Helper to trim whitespace
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t");
    if (start == std::string::npos)
        return "";
    size_t end = str.find_last_not_of(" \t");
    return str.substr(start, end - start + 1);
}

Residue ResidueFactory::create(const std::string& name, int sequence_number, char chain_id, char insertion_code,
                               const std::vector<Atom>& atoms) {
    // Determine properties using registry
    char one_letter = determine_one_letter_code(name);
    ResidueType type = determine_type(name, one_letter);
    bool purine = determine_is_purine(name, type);

    // Create residue using Builder pattern
    return Residue::create(name, sequence_number, chain_id)
        .insertion(insertion_code)
        .one_letter_code(one_letter)
        .type(type)
        .is_purine(purine)
        .atoms(atoms)
        .build();
}

char ResidueFactory::determine_one_letter_code(const std::string& name) {
    std::string trimmed = trim(name);

    // Standard nucleotides - uppercase
    if (trimmed == "A" || trimmed == "ADE" || trimmed == "DA")
        return 'A';
    if (trimmed == "C" || trimmed == "CYT" || trimmed == "DC")
        return 'C';
    if (trimmed == "G" || trimmed == "GUA" || trimmed == "DG")
        return 'G';
    if (trimmed == "T" || trimmed == "THY" || trimmed == "DT")
        return 'T';
    if (trimmed == "U" || trimmed == "URA" || trimmed == "DU")
        return 'U';
    if (trimmed == "I" || trimmed == "INO")
        return 'I';
    if (trimmed == "P" || trimmed == "PSU")
        return 'P';

    // Check modified nucleotide registry
    char code = ModifiedNucleotideRegistry::get_one_letter_code(trimmed);
    if (code != '?') {
        return code;
    }

    return '?';
}

ResidueType ResidueFactory::determine_type(const std::string& name, char one_letter_code) {
    // Standard nucleotides
    if (one_letter_code == 'A')
        return ResidueType::ADENINE;
    if (one_letter_code == 'C')
        return ResidueType::CYTOSINE;
    if (one_letter_code == 'G')
        return ResidueType::GUANINE;
    if (one_letter_code == 'T')
        return ResidueType::THYMINE;
    if (one_letter_code == 'U')
        return ResidueType::URACIL;
    if (one_letter_code == 'I')
        return ResidueType::INOSINE;
    if (one_letter_code == 'P')
        return ResidueType::PSEUDOURIDINE;

    // Modified nucleotides (lowercase) - map to parent type
    if (one_letter_code == 'a')
        return ResidueType::ADENINE;
    if (one_letter_code == 'c')
        return ResidueType::CYTOSINE;
    if (one_letter_code == 'g')
        return ResidueType::GUANINE;
    if (one_letter_code == 't')
        return ResidueType::THYMINE;
    if (one_letter_code == 'u')
        return ResidueType::URACIL;

    // Check registry for additional info
    auto type_opt = ModifiedNucleotideRegistry::get_base_type(trim(name));
    if (type_opt.has_value()) {
        return type_opt.value();
    }

    // Check for water
    std::string trimmed = trim(name);
    if (trimmed == "HOH" || trimmed == "WAT") {
        return ResidueType::WATER;
    }

    // Check for common ions
    std::string upper = trimmed;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    if (upper == "MG" || upper == "NA" || upper == "CL" || upper == "K" || upper == "CA" || upper == "ZN" ||
        upper == "FE" || upper == "MN") {
        return ResidueType::ION;
    }

    return ResidueType::UNKNOWN;
}

bool ResidueFactory::determine_is_purine(const std::string& name, ResidueType type) {
    // Check registry first
    auto is_purine_opt = ModifiedNucleotideRegistry::is_purine(trim(name));
    if (is_purine_opt.has_value()) {
        return is_purine_opt.value();
    }

    // Fallback to type-based determination
    return (type == ResidueType::ADENINE || type == ResidueType::GUANINE || type == ResidueType::INOSINE);
}

} // namespace core
} // namespace x3dna
