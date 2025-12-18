#include <x3dna/algorithms/template_assignment.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>

namespace x3dna {
namespace algorithms {

// These maps are now obsolete - all lookups go through ModifiedNucleotideRegistry
// Kept as empty maps to satisfy linker requirements from header declaration
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PURINES = {};
const std::map<std::string, core::ResidueType> TemplateAssignment::MODIFIED_PYRIMIDINES = {};

std::optional<core::ResidueType> TemplateAssignment::get_type_for_modified(const std::string& residue_name,
                                                                           bool /* is_purine */
) {
    // Use the centralized registry
    return core::ModifiedNucleotideRegistry::get_base_type(residue_name);
}

std::optional<std::vector<std::string>> TemplateAssignment::get_matching_atoms(const std::string& residue_name) {
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
