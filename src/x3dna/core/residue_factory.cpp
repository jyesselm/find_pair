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
    // Delegate to the centralized registry (handles both standard and modified nucleotides)
    return ModifiedNucleotideRegistry::get_one_letter_code(trim(name));
}

ResidueType ResidueFactory::determine_type(const std::string& name, char /* one_letter_code */) {
    std::string trimmed = trim(name);

    // Check registry first (handles both standard and modified nucleotides)
    auto type_opt = ModifiedNucleotideRegistry::get_base_type(trimmed);
    if (type_opt.has_value()) {
        return type_opt.value();
    }

    // Check for water
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
    // Check registry first (handles both standard and modified nucleotides)
    auto is_purine_opt = ModifiedNucleotideRegistry::is_purine(trim(name));
    if (is_purine_opt.has_value()) {
        return is_purine_opt.value();
    }

    // Fallback to type-based determination for unknown residues
    return (type == ResidueType::ADENINE || type == ResidueType::GUANINE || type == ResidueType::INOSINE);
}

} // namespace core
} // namespace x3dna
