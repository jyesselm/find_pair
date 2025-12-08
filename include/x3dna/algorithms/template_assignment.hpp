#pragma once

#include <x3dna/core/residue.hpp> // Includes ResidueType enum
#include <string>
#include <map>
#include <optional>

namespace x3dna {
namespace algorithms {

/**
 * @brief Handles template assignment for modified nucleotides
 *
 * This class provides lookup tables to override automatic detection
 * for modified nucleotides where atom-based classification would fail
 * or differ from legacy behavior.
 */
class TemplateAssignment {
public:
    /**
     * @brief Get the residue type for a modified nucleotide
     *
     * @param residue_name The 3-letter residue name
     * @param is_purine Whether automatic detection identified it as purine
     * @return ResidueType if found in lookup, nullopt otherwise
     */
    static std::optional<core::ResidueType> get_type_for_modified(const std::string& residue_name, bool is_purine);

    /**
     * @brief Get specific atoms to use for matching a modified nucleotide
     *
     * Some modified nucleotides need specific atom lists to match legacy exactly.
     * This avoids differences in auto-detection logic.
     *
     * @param residue_name The 3-letter residue name
     * @return Vector of atom names if specific list exists, nullopt otherwise
     */
    static std::optional<std::vector<std::string>> get_matching_atoms(const std::string& residue_name);

private:
    // Modified purines that need explicit template assignment
    static const std::map<std::string, core::ResidueType> MODIFIED_PURINES;

    // Modified pyrimidines that need explicit template assignment
    static const std::map<std::string, core::ResidueType> MODIFIED_PYRIMIDINES;
};

} // namespace algorithms
} // namespace x3dna
