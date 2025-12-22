#pragma once

#include <string>
#include <map>
#include <optional>
#include <x3dna/core/residue_type.hpp>
#include <x3dna/core/residue_classification.hpp>

namespace x3dna {
namespace core {

/**
 * @brief Registry for modified nucleotide properties
 *
 * Centralized lookup table for all modified nucleotides, providing:
 * - One-letter code mapping
 * - Base type (A, C, G, U, T, I, P)
 * - Purine/Pyrimidine classification
 *
 * This replaces scattered if-statements with a clean, data-driven approach.
 */
class ModifiedNucleotideRegistry {
public:
    struct NucleotideInfo {
        char one_letter_code;    // 'a', 'c', 'g', 'u', 't', 'I', 'P'
        ResidueType base_type;   // ADENINE, CYTOSINE, etc.
        bool is_purine;          // true for A/G/I, false for C/U/T/P
        std::string description; // Human-readable description
    };

    /**
     * @brief Get information for a modified nucleotide
     * @param residue_name Three-letter residue name (e.g., "ATP", "SAM")
     * @return NucleotideInfo if found, nullopt otherwise
     */
    [[nodiscard]] static std::optional<NucleotideInfo> get_info(const std::string& residue_name);

    /**
     * @brief Get one-letter code for a residue name
     * @param residue_name Three-letter residue name
     * @return One-letter code if found, '?' otherwise
     */
    [[nodiscard]] static char get_one_letter_code(const std::string& residue_name);

    /**
     * @brief Get base type for a modified nucleotide
     * @param residue_name Three-letter residue name
     * @return ResidueType if found, nullopt otherwise
     */
    [[nodiscard]] static std::optional<ResidueType> get_base_type(const std::string& residue_name);

    /**
     * @brief Check if a residue is a purine derivative
     * @param residue_name Three-letter residue name
     * @return true if purine, false if pyrimidine, nullopt if not found
     */
    [[nodiscard]] static std::optional<bool> is_purine(const std::string& residue_name);

    /**
     * @brief Check if a residue is in the registry
     * @param residue_name Three-letter residue name
     * @return true if residue is known/registered, false otherwise
     */
    [[nodiscard]] static bool contains(const std::string& residue_name);

    /**
     * @brief Classify a residue by name
     * @param residue_name Three-letter residue name (e.g., "5MC", "DA", "HOH")
     * @return Full ResidueClassification with all type information
     */
    [[nodiscard]] static ResidueClassification classify(const std::string& residue_name);

private:
    // The central registry - lazy-loaded on first access via get_registry()
    // This allows ResourceLocator to be initialized before the registry is loaded
    static const std::map<std::string, NucleotideInfo>& registry_;
};

} // namespace core
} // namespace x3dna
