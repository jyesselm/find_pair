#include <x3dna/algorithms/template_assignment.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>

namespace x3dna {
namespace algorithms {

// Legacy approach - now using ModifiedNucleotideRegistry instead
// Keep these for backwards compatibility but they're no longer populated
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PURINES = {
    // Previously added
    {"A23", core::ResidueType::ADENINE}, // 2'-deoxy-2'-fluoroadenosine
    {"GTP", core::ResidueType::GUANINE}, // Guanosine triphosphate
    {"IGU", core::ResidueType::GUANINE}, // Isoguanosine
    {"DGP", core::ResidueType::GUANINE}, // Deoxyguanosine phosphate

    // Common modified adenines (>50 occurrences in failures)
    {"ATP", core::ResidueType::ADENINE}, // Adenosine triphosphate
    {"ADP", core::ResidueType::ADENINE}, // Adenosine diphosphate
    {"AMP", core::ResidueType::ADENINE}, // Adenosine monophosphate
    {"SAM", core::ResidueType::ADENINE}, // S-adenosylmethionine
    {"SAH", core::ResidueType::ADENINE}, // S-adenosylhomocysteine
    {"APC", core::ResidueType::ADENINE}, // Adenosine phosphate carrier
    {"AF2", core::ResidueType::ADENINE}, // Modified adenosine
    {"0A", core::ResidueType::ADENINE},  // Numbered adenine variant

    // Common modified guanines (>50 occurrences)
    {"GDP", core::ResidueType::GUANINE}, // Guanosine diphosphate
    {"5GP", core::ResidueType::GUANINE}, // 5'-guanosine monophosphate
    {"PRF", core::ResidueType::GUANINE}, // Prefluorinated guanine
    {"BGM", core::ResidueType::GUANINE}, // Beta-D-glucosyl-guanine
    {"G48", core::ResidueType::GUANINE}, // Modified guanine
    {"0G", core::ResidueType::GUANINE},  // Numbered guanine variant

    // Inosines
    {"IMP", core::ResidueType::INOSINE}, // Inosine monophosphate
};

// Modified pyrimidines - hardcoded to match legacy template selection
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PYRIMIDINES = {
    // Modified thymines
    {"5MU", core::ResidueType::THYMINE}, // 5-methyluridine - has C5M
    {"TLN", core::ResidueType::THYMINE}, // 3-methyluridine - has C5M

    // Modified uracils
    {"70U", core::ResidueType::URACIL}, // 7-methyluridine
    {"H2U", core::ResidueType::URACIL}, // dihydrouridine
    {"DHU", core::ResidueType::URACIL}, // dihydrouridine
    {"OMU", core::ResidueType::URACIL}, // O-methyluridine
    {"4SU", core::ResidueType::URACIL}, // 4-thiouridine
    {"S4U", core::ResidueType::URACIL}, // 4-thiouridine
    {"2MU", core::ResidueType::URACIL}, // 2-methyluridine
    {"KIR", core::ResidueType::URACIL}, // kinetin riboside
    {"US5", core::ResidueType::URACIL}, // 5-hydroxymethyluridine
    {"UTP", core::ResidueType::URACIL}, // Uridine triphosphate
    {"J48", core::ResidueType::URACIL}, // Hypermodified nucleotide (wybutosine precursor)
    {"NMN", core::ResidueType::URACIL}, // Nicotinamide mononucleotide (treat as uracil)
    {"NNR", core::ResidueType::URACIL}, // Nicotinamide riboside (treat as uracil)
    {"WVQ", core::ResidueType::URACIL}, // Unknown modified uracil
    {"UFT", core::ResidueType::URACIL}, // Modified uracil
    {"0U", core::ResidueType::URACIL},  // Numbered uracil variant

    // Modified cytosines
    {"EPE", core::ResidueType::CYTOSINE}, // Modified cytosine analog - has N4
    {"2YR", core::ResidueType::CYTOSINE}, // 2'-O-ribosylcytidine
    {"CTP", core::ResidueType::CYTOSINE}, // Cytidine triphosphate
    {"CSL", core::ResidueType::CYTOSINE}, // Unknown modified cytosine
    {"CBV", core::ResidueType::CYTOSINE}, // Carbovir (modified cytosine)
    {"CCC", core::ResidueType::CYTOSINE}, // Modified cytosine
    {"CFZ", core::ResidueType::CYTOSINE}, // Modified cytosine
    {"C43", core::ResidueType::CYTOSINE}, // Modified cytosine
    {"A6C", core::ResidueType::CYTOSINE}, // N6-acetyl-cytosine
    {"CFL", core::ResidueType::CYTOSINE}, // Modified cytosine
    {"RSQ", core::ResidueType::CYTOSINE}, // R-squared cytosine
    {"0C", core::ResidueType::CYTOSINE},  // Numbered cytosine variant

    // Pseudouridines
    {"B8H", core::ResidueType::PSEUDOURIDINE}, // Pseudouridine derivative
};

std::optional<core::ResidueType> TemplateAssignment::get_type_for_modified(
    const std::string& residue_name, bool /* is_purine */ // Unused - we check registry
) {
    // Use the centralized registry instead of hardcoded maps
    return core::ModifiedNucleotideRegistry::get_base_type(residue_name);
}

std::optional<std::vector<std::string>>
TemplateAssignment::get_matching_atoms(const std::string& residue_name) {
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
