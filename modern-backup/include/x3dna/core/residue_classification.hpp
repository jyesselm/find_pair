/**
 * @file residue_classification.hpp
 * @brief Hierarchical residue classification system
 */

#pragma once

#include <string>
#include <x3dna/core/molecule_type.hpp>
#include <x3dna/core/residue_type.hpp>

namespace x3dna {
namespace core {

/**
 * @struct ResidueClassification
 * @brief Complete classification of a residue with hierarchical type information
 *
 * Provides a unified way to query residue properties:
 * - Molecule type (nucleic acid, protein, water, ion, ligand)
 * - For nucleic acids: RNA vs DNA, canonical vs modified, base type
 * - Backwards compatible with legacy ResidueType enum
 */
struct ResidueClassification {
    // Level 1: Molecule type
    MoleculeType molecule_type = MoleculeType::UNKNOWN;

    // Level 2: For nucleic acids only
    NucleicAcidType nucleic_acid_type = NucleicAcidType::UNKNOWN;
    bool is_modified = false;

    // Level 3: For nucleotides only
    BaseType base_type = BaseType::UNKNOWN;
    BaseCategory base_category = BaseCategory::UNKNOWN;

    // Original residue name from PDB (e.g., "5MC", "PSU", "ATP", "HOH")
    std::string residue_name;

    // Canonical single-letter code ('A', 'C', 'G', 'T', 'U', or '?' for unknown)
    char canonical_code = '?';

    // === Query methods ===

    /// Is this a nucleotide (RNA or DNA)?
    [[nodiscard]] bool is_nucleotide() const;

    /// Is this RNA?
    [[nodiscard]] bool is_rna() const;

    /// Is this DNA?
    [[nodiscard]] bool is_dna() const;

    /// Is this a purine base (A, G, I)?
    [[nodiscard]] bool is_purine() const;

    /// Is this a pyrimidine base (C, T, U, PSU)?
    [[nodiscard]] bool is_pyrimidine() const;

    /// Is this a canonical (non-modified) nucleotide?
    [[nodiscard]] bool is_canonical() const;

    /// Is this a protein residue?
    [[nodiscard]] bool is_protein() const;

    /// Is this a water molecule?
    [[nodiscard]] bool is_water() const;

    /// Is this an ion?
    [[nodiscard]] bool is_ion() const;

    /// Is this a ligand?
    [[nodiscard]] bool is_ligand() const;

    // === Legacy compatibility ===

    /// Convert to legacy ResidueType enum for backwards compatibility
    [[nodiscard]] ResidueType to_legacy_type() const;

    /// Create classification from legacy ResidueType and residue name
    [[nodiscard]] static ResidueClassification from_legacy(ResidueType type, const std::string& residue_name,
                                                           bool is_purine_hint = false);
};

} // namespace core
} // namespace x3dna
