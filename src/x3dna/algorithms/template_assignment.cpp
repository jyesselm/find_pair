#include <x3dna/algorithms/template_assignment.hpp>

namespace x3dna {
namespace algorithms {

// Modified purines - hardcoded to match legacy template selection
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PURINES = {
    {"A23", core::ResidueType::ADENINE}, // 2'-deoxy-2'-fluoroadenosine
    // Add more as discovered during validation
};

// Modified pyrimidines - hardcoded to match legacy template selection
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PYRIMIDINES = {
    {"5MU", core::ResidueType::THYMINE},  // 5-methyluridine - has C5M
    {"TLN", core::ResidueType::THYMINE},  // 3-methyluridine - has C5M
    {"70U", core::ResidueType::URACIL},   // 7-methyluridine
    {"H2U", core::ResidueType::URACIL},   // dihydrouridine
    {"DHU", core::ResidueType::URACIL},   // dihydrouridine
    {"OMU", core::ResidueType::URACIL},   // O-methyluridine
    {"4SU", core::ResidueType::URACIL},   // 4-thiouridine
    {"S4U", core::ResidueType::URACIL},   // 4-thiouridine
    {"2MU", core::ResidueType::URACIL},   // 2-methyluridine
    {"KIR", core::ResidueType::URACIL},   // kinetin riboside
    {"EPE", core::ResidueType::CYTOSINE}, // Modified cytosine analog - has N4
    // Add more as discovered during validation
};

std::optional<core::ResidueType> TemplateAssignment::get_type_for_modified(
    const std::string& residue_name, bool /* is_purine */ // Unused - we check both tables
) {
    // Check BOTH lookup tables regardless of is_purine flag
    // Reason: Legacy fallback logic can misclassify purines as pyrimidines
    // (e.g., A23 has C8+N9 but pyrimidine-only RMSD check passes, so has_purine_atoms becomes
    // false)

    // Try purines first
    {
        auto it = MODIFIED_PURINES.find(residue_name);
        if (it != MODIFIED_PURINES.end()) {
            return it->second;
        }
    }

    // Then try pyrimidines
    {
        auto it = MODIFIED_PYRIMIDINES.find(residue_name);
        if (it != MODIFIED_PYRIMIDINES.end()) {
            return it->second;
        }
    }

    return std::nullopt;
}

std::optional<std::vector<std::string>> TemplateAssignment::get_matching_atoms(
    const std::string& residue_name) {
    // Hardcoded atom lists for modified residues that need exact legacy matching
    // Format: residue_name -> list of atom names to match (in order)
    static const std::map<std::string, std::vector<std::string>> MATCHING_ATOMS = {
        // A23: 2'-deoxy-2'-fluoroadenosine - use 9 purine ring atoms
        {"A23", {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "}},
        // KIR: kinetin riboside - use 5 pyrimidine atoms (no N3, has C3 instead)
        {"KIR", {" C4 ", " C2 ", " N1 ", " C6 ", " C5 "}},
        // Add more as needed during validation
    };

    auto it = MATCHING_ATOMS.find(residue_name);
    if (it != MATCHING_ATOMS.end()) {
        return it->second;
    }

    return std::nullopt;
}

} // namespace algorithms
} // namespace x3dna
