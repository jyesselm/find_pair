/**
 * @file residue_classification.hpp
 * @brief Enhanced residue classification with full protein support
 */

#pragma once

#include <string>
#include <x3dna/core/typing/molecule_type.hpp>
#include <x3dna/core/typing/nucleotide_type.hpp>
#include <x3dna/core/typing/protein_type.hpp>
#include <x3dna/core/typing/solvent_type.hpp>

namespace x3dna {
namespace core {

// Forward declaration for legacy type
enum class ResidueType;

namespace typing {

/**
 * @struct ResidueClassification
 * @brief Complete classification of a residue with hierarchical type information
 *
 * Provides a unified way to query residue properties:
 * - Molecule type (nucleic acid, protein, water, ion, ligand)
 * - For nucleic acids: RNA vs DNA, canonical vs modified, base type
 * - For proteins: amino acid type, category, standard vs modified
 * - For solvents: water vs ion type
 * - Backwards compatible with legacy ResidueType enum
 */
struct ResidueClassification {
    // === Level 1: Molecule type ===
    MoleculeType molecule_type = MoleculeType::UNKNOWN;

    // === Level 2a: For nucleic acids ===
    NucleicAcidType nucleic_acid_type = NucleicAcidType::UNKNOWN;
    BaseType base_type = BaseType::UNKNOWN;
    BaseCategory base_category = BaseCategory::UNKNOWN;
    bool is_modified_nucleotide = false;

    // === Level 2b: For proteins ===
    AminoAcidType amino_acid_type = AminoAcidType::UNKNOWN;
    AminoAcidCategory amino_acid_category = AminoAcidCategory::UNKNOWN;
    bool is_modified_amino_acid = false;

    // === Level 2c: For solvents ===
    SolventType solvent_type = SolventType::UNKNOWN;
    IonType ion_type = IonType::UNKNOWN;

    // === Common fields ===
    std::string residue_name;       ///< Original 3-letter code from PDB
    char one_letter_code = '?';     ///< Single letter representation
    char canonical_code = '?';      ///< Canonical base/AA code

    // === Query methods: Molecule type ===

    [[nodiscard]] bool is_nucleotide() const {
        return molecule_type == MoleculeType::NUCLEIC_ACID;
    }

    [[nodiscard]] bool is_protein() const {
        return molecule_type == MoleculeType::PROTEIN;
    }

    [[nodiscard]] bool is_water() const {
        return molecule_type == MoleculeType::WATER;
    }

    [[nodiscard]] bool is_ion() const {
        return molecule_type == MoleculeType::ION;
    }

    [[nodiscard]] bool is_ligand() const {
        return molecule_type == MoleculeType::LIGAND;
    }

    // === Query methods: Nucleic acids ===

    [[nodiscard]] bool is_rna() const {
        return molecule_type == MoleculeType::NUCLEIC_ACID &&
               nucleic_acid_type == NucleicAcidType::RNA;
    }

    [[nodiscard]] bool is_dna() const {
        return molecule_type == MoleculeType::NUCLEIC_ACID &&
               nucleic_acid_type == NucleicAcidType::DNA;
    }

    [[nodiscard]] bool is_purine() const {
        return base_category == BaseCategory::PURINE;
    }

    [[nodiscard]] bool is_pyrimidine() const {
        return base_category == BaseCategory::PYRIMIDINE;
    }

    [[nodiscard]] bool is_canonical_nucleotide() const {
        return is_nucleotide() && !is_modified_nucleotide;
    }

    [[nodiscard]] bool is_standard_base() const {
        return base_type == BaseType::ADENINE ||
               base_type == BaseType::GUANINE ||
               base_type == BaseType::CYTOSINE ||
               base_type == BaseType::THYMINE ||
               base_type == BaseType::URACIL;
    }

    // === Query methods: Proteins ===

    [[nodiscard]] bool is_standard_amino_acid() const {
        return is_protein() && typing::is_standard_amino_acid(amino_acid_type);
    }

    [[nodiscard]] bool is_hydrophobic() const {
        return amino_acid_category == AminoAcidCategory::HYDROPHOBIC;
    }

    [[nodiscard]] bool is_polar() const {
        return amino_acid_category == AminoAcidCategory::POLAR;
    }

    [[nodiscard]] bool is_charged() const {
        return amino_acid_category == AminoAcidCategory::POSITIVE ||
               amino_acid_category == AminoAcidCategory::NEGATIVE;
    }

    [[nodiscard]] bool is_positive() const {
        return amino_acid_category == AminoAcidCategory::POSITIVE;
    }

    [[nodiscard]] bool is_negative() const {
        return amino_acid_category == AminoAcidCategory::NEGATIVE;
    }

    // === Query methods: Ions ===

    [[nodiscard]] bool is_cation() const {
        return is_ion() && typing::is_cation(ion_type);
    }

    [[nodiscard]] bool is_anion() const {
        return is_ion() && typing::is_anion(ion_type);
    }

    // === Legacy compatibility ===

    /**
     * @brief Convert to legacy ResidueType enum for backwards compatibility
     */
    [[nodiscard]] ResidueType to_legacy_type() const;

    /**
     * @brief Create classification from legacy ResidueType and residue name
     */
    [[nodiscard]] static ResidueClassification from_legacy(
        ResidueType type,
        const std::string& residue_name,
        bool is_purine_hint = false);
};

} // namespace typing

// Note: We don't expose ResidueClassification here to avoid conflicts
// with the existing core::ResidueClassification. The migration will
// update the existing class to inherit from typing::ResidueClassification.

} // namespace core
} // namespace x3dna
