/**
 * @file type_registry.hpp
 * @brief Unified registry for all molecular entity type lookups
 */

#pragma once

#include <string>
#include <map>
#include <set>
#include <optional>
#include <x3dna/core/typing/molecule_type.hpp>
#include <x3dna/core/typing/nucleotide_type.hpp>
#include <x3dna/core/typing/protein_type.hpp>
#include <x3dna/core/typing/solvent_type.hpp>
#include <x3dna/core/typing/atom_classification.hpp>
#include <x3dna/core/typing/residue_classification.hpp>

namespace x3dna {
namespace core {
namespace typing {

/**
 * @struct NucleotideInfo
 * @brief Information about a nucleotide residue
 */
struct NucleotideInfo {
    char one_letter_code = '?';
    BaseType base_type = BaseType::UNKNOWN;
    bool is_purine = false;
    bool is_modified = false;
    std::string description;
};

/**
 * @struct AminoAcidInfo
 * @brief Information about an amino acid residue
 */
struct AminoAcidInfo {
    char one_letter_code = '?';
    AminoAcidType type = AminoAcidType::UNKNOWN;
    AminoAcidCategory category = AminoAcidCategory::UNKNOWN;
    bool is_modified = false;
    std::string description;
};

/**
 * @class TypeRegistry
 * @brief Singleton registry for all molecular entity type lookups
 *
 * Provides a single source of truth for:
 * - Residue classification (nucleotide, protein, water, ion, ligand)
 * - Atom classification (backbone, sugar, nucleobase, mainchain, sidechain)
 * - H-bond role determination
 *
 * Replaces scattered type checks throughout the codebase with
 * a centralized, data-driven approach.
 */
class TypeRegistry {
public:
    /**
     * @brief Get the singleton instance
     * @return Reference to the TypeRegistry singleton
     */
    [[nodiscard]] static const TypeRegistry& instance();

    // === Residue classification ===

    /**
     * @brief Classify a residue by name
     * @param residue_name Three-letter residue code (e.g., "ALA", "5MC", "HOH")
     * @return Full ResidueClassification with all type information
     */
    [[nodiscard]] ResidueClassification classify_residue(const std::string& residue_name) const;

    /**
     * @brief Check if a residue name is a water molecule
     * @param residue_name Residue name to check
     * @return true for HOH, WAT, H2O, OH2, SOL
     */
    [[nodiscard]] bool is_water(const std::string& residue_name) const;

    /**
     * @brief Check if a residue name is an ion
     * @param residue_name Residue name to check
     * @return true for MG, NA, K, CA, ZN, etc.
     */
    [[nodiscard]] bool is_ion(const std::string& residue_name) const;

    /**
     * @brief Check if a residue name is a standard amino acid
     * @param residue_name Residue name to check
     * @return true for the standard 20 amino acids
     */
    [[nodiscard]] bool is_amino_acid(const std::string& residue_name) const;

    /**
     * @brief Check if a residue name is a nucleotide
     * @param residue_name Residue name to check
     * @return true for known nucleotides (canonical and modified)
     */
    [[nodiscard]] bool is_nucleotide(const std::string& residue_name) const;

    // === Nucleotide-specific lookups ===

    /**
     * @brief Get nucleotide info if available
     * @param residue_name Residue name
     * @return NucleotideInfo if found, nullopt otherwise
     */
    [[nodiscard]] std::optional<NucleotideInfo> get_nucleotide_info(const std::string& residue_name) const;

    /**
     * @brief Get one-letter code for a residue
     * @param residue_name Residue name
     * @return One-letter code or '?' if not found
     */
    [[nodiscard]] char get_one_letter_code(const std::string& residue_name) const;

    /**
     * @brief Check if a nucleotide is a purine
     * @param residue_name Residue name
     * @return true if purine, false otherwise
     */
    [[nodiscard]] std::optional<bool> is_purine(const std::string& residue_name) const;

    // === Amino acid-specific lookups ===

    /**
     * @brief Get amino acid info if available
     * @param residue_name Residue name
     * @return AminoAcidInfo if found, nullopt otherwise
     */
    [[nodiscard]] std::optional<AminoAcidInfo> get_amino_acid_info(const std::string& residue_name) const;

    // === Ion-specific lookups ===

    /**
     * @brief Get ion type for a residue name
     * @param residue_name Residue name
     * @return IonType if ion, IonType::UNKNOWN otherwise
     */
    [[nodiscard]] IonType get_ion_type(const std::string& residue_name) const;

private:
    TypeRegistry();
    ~TypeRegistry() = default;

    // Non-copyable
    TypeRegistry(const TypeRegistry&) = delete;
    TypeRegistry& operator=(const TypeRegistry&) = delete;

    void load_nucleotides();
    void load_amino_acids();
    void load_waters();
    void load_ions();

    // Data storage
    std::map<std::string, NucleotideInfo> nucleotides_;
    std::map<std::string, AminoAcidInfo> amino_acids_;
    std::set<std::string> water_names_;
    std::map<std::string, IonType> ion_types_;
};

} // namespace typing

// Expose in core namespace for convenience
using typing::TypeRegistry;
using typing::NucleotideInfo;
using typing::AminoAcidInfo;

} // namespace core
} // namespace x3dna
