/**
 * @file residue.hpp
 * @brief Residue class representing a single residue in a PDB structure
 */

#pragma once

#include <string>
#include <vector>
#include <optional>
#include <algorithm>
#include <cctype>
#include <limits>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/residue_type.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>
#include <x3dna/core/typing/residue_classification.hpp>
#include <x3dna/core/typing/type_registry.hpp>
#include <x3dna/core/string_utils.hpp>

namespace x3dna {
namespace core {

/**
 * @class Residue
 * @brief Represents a single residue (nucleotide or amino acid) with atoms
 *
 * Residues can be constructed via:
 * 1. Simple constructor: Residue(name, seq_num, chain_id) - for basic creation
 * 2. Builder pattern: Residue::create(name, seq_num, chain_id).type(...).build() - for full control
 * 3. Residue::create_from_atoms() - recommended for proper property initialization
 */
class Residue {
public:
    class Builder; // Forward declaration

    /**
     * @brief Default constructor
     */
    Residue() = default;

    /**
     * @brief Constructor with name, sequence number, chain ID, and insertion code
     * @param name Residue name (e.g., "C", "G", "A", "T", "ADE", "GUA", etc.) - will be trimmed
     * @param seq_num Sequence number
     * @param chain_id Chain identifier
     * @param insertion Insertion code (PDB column 27, default "")
     */
    Residue(const std::string& name, int seq_num, const std::string& chain_id, const std::string& insertion = "")
        : name_(trim(name)), seq_num_(seq_num), chain_id_(chain_id), insertion_(insertion) {
        // Auto-initialize classification and one_letter_code from trimmed name
        classification_ = typing::TypeRegistry::instance().classify_residue(name_);
        one_letter_code_ = ModifiedNucleotideRegistry::get_one_letter_code(name_);
    }

    /**
     * @brief Create a Builder for fluent residue construction
     * @param name Residue name (required)
     * @param seq_num Sequence number (required)
     * @param chain_id Chain ID (required)
     * @return Builder instance for fluent construction
     */
    [[nodiscard]] static Builder create(const std::string& name, int seq_num, const std::string& chain_id);

    /**
     * @brief Create a fully-initialized residue from PDB/CIF data
     *
     * This is the recommended way to create residues. It uses ModifiedNucleotideRegistry
     * to determine one_letter_code, type, and is_purine properties.
     *
     * @param name Three-letter residue name (e.g., "A", "ATP", "PSU")
     * @param sequence_number PDB sequence number
     * @param chain_id Chain identifier
     * @param insertion_code Insertion code (default "")
     * @param atoms Vector of atoms in this residue
     * @return Residue with all properties initialized
     */
    [[nodiscard]] static Residue create_from_atoms(const std::string& name, int sequence_number,
                                                   const std::string& chain_id, const std::string& insertion_code,
                                                   const std::vector<Atom>& atoms);

    // Getters
    [[nodiscard]] const std::string& name() const {
        return name_;
    }
    [[nodiscard]] int seq_num() const {
        return seq_num_;
    }
    [[nodiscard]] const std::string& chain_id() const {
        return chain_id_;
    }
    [[nodiscard]] const std::string& insertion() const {
        return insertion_;
    }

    /**
     * @brief Get unique residue identifier
     * @return String in format "chain_id-res_name-res_num" or "chain_id-res_name-res_numX" if insertion code exists
     *
     * Examples:
     *   "A-G-5"    (chain A, guanine, position 5, no insertion)
     *   "A-C-10A"  (chain A, cytosine, position 10, insertion code A)
     *   "B-PSU-25" (chain B, pseudouridine, position 25)
     */
    [[nodiscard]] std::string res_id() const {
        std::string id = chain_id_ + "-" + name_ + "-" + std::to_string(seq_num_);
        if (!insertion_.empty()) {
            id += insertion_;
        }
        return id;
    }

    [[nodiscard]] const std::vector<Atom>& atoms() const {
        return atoms_;
    }
    std::vector<Atom>& atoms() {
        return atoms_;
    }
    [[nodiscard]] size_t num_atoms() const {
        return atoms_.size();
    }

    // Container-like interface for atoms
    [[nodiscard]] auto begin() const {
        return atoms_.begin();
    }
    [[nodiscard]] auto end() const {
        return atoms_.end();
    }
    [[nodiscard]] auto begin() {
        return atoms_.begin();
    }
    [[nodiscard]] auto end() {
        return atoms_.end();
    }
    [[nodiscard]] size_t size() const {
        return atoms_.size();
    }
    [[nodiscard]] bool empty() const {
        return atoms_.empty();
    }
    [[nodiscard]] const Atom& operator[](size_t idx) const {
        return atoms_[idx];
    }
    [[nodiscard]] Atom& operator[](size_t idx) {
        return atoms_[idx];
    }

    /**
     * @brief Get reference frame (if set)
     */
    [[nodiscard]] std::optional<ReferenceFrame> reference_frame() const {
        return reference_frame_;
    }

    // Modification methods (kept for essential post-construction updates)

    /**
     * @brief Set reference frame for this residue
     *
     * Reference frames are calculated after construction, so this setter
     * is retained for the frame calculation workflow.
     */
    void set_reference_frame(const ReferenceFrame& frame) {
        reference_frame_ = frame;
    }

    /**
     * @brief Add an atom to this residue
     *
     * Atoms are typically added one at a time during PDB parsing,
     * so this method is retained for the parsing workflow.
     */
    void add_atom(const Atom& atom) {
        atoms_.push_back(atom);
    }

    /**
     * @brief Find an atom by name
     * @param atom_name Atom name (can be trimmed or padded, e.g., "C1'" or " C1'")
     * @return Optional atom if found
     */
    [[nodiscard]] std::optional<Atom> find_atom(const std::string& atom_name) const {
        // Atom names are stored in PDB 4-character format, so compare directly
        // or normalize the search key to match
        for (const auto& atom : atoms_) {
            if (atom.name() == atom_name) {
                return atom;
            }
            // Also check trimmed version for convenience
            std::string trimmed = atom_name;
            trimmed.erase(0, trimmed.find_first_not_of(" \t\n\r"));
            trimmed.erase(trimmed.find_last_not_of(" \t\n\r") + 1);

            std::string atom_trimmed = atom.name();
            atom_trimmed.erase(0, atom_trimmed.find_first_not_of(" \t\n\r"));
            atom_trimmed.erase(atom_trimmed.find_last_not_of(" \t\n\r") + 1);

            if (atom_trimmed == trimmed) {
                return atom;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Get all ring atoms (base ring atoms for nucleotides)
     * @return Vector of ring atoms
     */
    [[nodiscard]] std::vector<Atom> ring_atoms() const {
        std::vector<Atom> ring;
        for (const auto& atom : atoms_) {
            if (atom.is_ring_atom()) {
                ring.push_back(atom);
            }
        }
        return ring;
    }

    /**
     * @brief Get one-letter code for this residue
     * @return One-letter code (A, C, G, T, U for nucleotides, lowercase for modified)
     *
     * Recognizes standard nucleotide names:
     * - Single letter: A, C, G, T, U
     * - Three letter: ADE, CYT, GUA, THY, URA
     * - DNA format: DA, DC, DG, DT (may have leading/trailing spaces)
     * - Modified nucleotides: PSU, 5MC, 5MU, H2U, A2M, etc.
     */
    /**
     * @brief Get one-letter code (stored value, not computed)
     * @return One-letter code ('A', 'C', 'G', 'U', 'a', 'c', etc.)
     */
    [[nodiscard]] char one_letter_code() const {
        return one_letter_code_;
    }

    /**
     * @brief Check if this is a nucleotide
     * @return True if nucleotide (A, C, G, T, U or modified a, c, g, t, u, P)
     */
    [[nodiscard]] bool is_nucleotide() const {
        return classification_.is_nucleotide();
    }

    /**
     * @brief Check if residue is a purine (A, G, I)
     * @return true if purine (A, G, I), false otherwise
     */
    [[nodiscard]] bool is_purine() const {
        return classification_.is_purine();
    }

    /**
     * @brief Get RY classification (Purine=1, Pyrimidine=0)
     * Uses classification system
     * @return 1 for purines, 0 for pyrimidines, -1 for non-nucleotides
     */
    [[nodiscard]] int ry_classification() const {
        if (classification_.is_purine()) {
            return 1; // Purine
        }
        // Check if it's a nucleotide
        char code = one_letter_code();
        if (code == 'C' || code == 'T' || code == 'U' || code == 'P' || code == 'c' || code == 't' || code == 'u') {
            return 0; // Pyrimidine
        }
        return -1; // Not a nucleotide
    }

    /**
     * @brief Get residue type
     */
    /**
     * @brief Get residue type (derived from classification)
     * @return ResidueType derived from classification
     */
    [[nodiscard]] ResidueType residue_type() const {
        return classification_.to_legacy_type();
    }

    /**
     * @brief Get full hierarchical classification of this residue
     * @return ResidueClassification with molecule type, base type, categories, etc.
     */
    [[nodiscard]] const typing::ResidueClassification& classification() const {
        return classification_;
    }

    /**
     * @brief Check if this is a pyrimidine (C, T, U, or pseudouridine)
     */
    [[nodiscard]] bool is_pyrimidine() const {
        return classification_.is_pyrimidine();
    }

    /**
     * @brief Check if this is a protein residue (amino acid)
     */
    [[nodiscard]] bool is_protein() const {
        return classification_.is_protein();
    }

    /**
     * @brief Check if this is a water molecule
     */
    [[nodiscard]] bool is_water() const {
        return classification_.is_water();
    }

    /**
     * @brief Check if this is an ion
     */
    [[nodiscard]] bool is_ion() const {
        return classification_.is_ion();
    }

    /**
     * @brief Check if this is an RNA nucleotide
     */
    [[nodiscard]] bool is_rna() const {
        return classification_.is_rna();
    }

    /**
     * @brief Check if this is a DNA nucleotide
     */
    [[nodiscard]] bool is_dna() const {
        return classification_.is_dna();
    }

    /**
     * @brief Get the base type (ADENINE, CYTOSINE, etc.)
     */
    [[nodiscard]] typing::BaseType base_type() const {
        return classification_.base_type;
    }

    /**
     * @brief Get the molecule type (NucleicAcid, Protein, Solvent, Unknown)
     */
    [[nodiscard]] typing::MoleculeType molecule_type() const {
        return classification_.molecule_type;
    }

    /**
     * @brief Get atom range for this residue (start_atom, end_atom)
     *
     * Returns the minimum and maximum legacy_atom_idx values from all atoms
     * in this residue. This is used for residue_indices JSON generation.
     *
     * @return Pair of (start_atom, end_atom) legacy atom indices, or (0, 0) if no atoms
     */
    [[nodiscard]] std::pair<int, int> atom_range() const {
        if (atoms_.empty()) {
            return {0, 0};
        }

        int min_idx = std::numeric_limits<int>::max();
        int max_idx = 0;

        for (const auto& atom : atoms_) {
            int legacy_idx = atom.legacy_atom_idx();
            if (legacy_idx > 0) {
                min_idx = std::min(min_idx, legacy_idx);
                max_idx = std::max(max_idx, legacy_idx);
            }
        }

        if (min_idx == std::numeric_limits<int>::max()) {
            return {0, 0};
        }

        return {min_idx, max_idx};
    }

    /**
     * @brief Get legacy residue index (1-based indexing from legacy code)
     * @return Legacy residue index for backward compatibility
     */
    [[nodiscard]] int legacy_residue_idx() const {
        return legacy_residue_idx_;
    }

    /**
     * @brief Set legacy residue index
     * @param idx Legacy residue index (1-based)
     */
    void set_legacy_residue_idx(int idx) {
        legacy_residue_idx_ = idx;
    }

private:
    friend class Builder;

    std::string name_;                              // Residue name (typically trimmed, e.g., "A", "ADE", "PSU")
    char one_letter_code_ = '?';                    // One-letter code (stored, not computed)
    int seq_num_ = 0;                               // Sequence number
    std::string chain_id_;                          // Chain identifier
    std::string insertion_;                         // Insertion code (PDB column 27)
    std::vector<Atom> atoms_;                       // Atoms in this residue
    std::optional<ReferenceFrame> reference_frame_; // Reference frame (if calculated)
    typing::ResidueClassification classification_;  // Full hierarchical classification
    int legacy_residue_idx_ = 0;                    // Legacy 1-based residue index for backward compatibility
};

/**
 * @class Residue::Builder
 * @brief Fluent builder for constructing Residue objects
 *
 * Usage:
 *   auto residue = Residue::create("  G", 42, 'A')
 *       .insertion(' ')
 *       .one_letter_code('G')
 *       .type(ResidueType::GUANINE)
 *       .is_purine(true)
 *       .atoms(atom_vector)
 *       .build();
 */
class Residue::Builder {
public:
    /**
     * @brief Constructor with required fields
     * @param name Residue name (will be trimmed, e.g., "  A" becomes "A")
     * @param seq_num Sequence number
     * @param chain_id Chain identifier
     */
    Builder(const std::string& name, int seq_num, const std::string& chain_id) {
        residue_.name_ = trim(name);
        residue_.seq_num_ = seq_num;
        residue_.chain_id_ = chain_id;
    }

    Builder& insertion(const std::string& ins) {
        residue_.insertion_ = ins;
        return *this;
    }

    Builder& one_letter_code(char code) {
        residue_.one_letter_code_ = code;
        return *this;
    }

    Builder& classification(const typing::ResidueClassification& c) {
        residue_.classification_ = c;
        return *this;
    }

    Builder& atoms(const std::vector<Atom>& atom_list) {
        residue_.atoms_ = atom_list;
        return *this;
    }

    Builder& atoms(std::vector<Atom>&& atom_list) {
        residue_.atoms_ = std::move(atom_list);
        return *this;
    }

    Builder& add_atom(const Atom& atom) {
        residue_.atoms_.push_back(atom);
        return *this;
    }

    Builder& legacy_residue_idx(int idx) {
        residue_.legacy_residue_idx_ = idx;
        return *this;
    }

    /**
     * @brief Build and return the constructed Residue
     * @return Constructed Residue object
     */
    [[nodiscard]] Residue build() const {
        return residue_;
    }

private:
    Residue residue_;
};

// Inline implementation of create()
inline Residue::Builder Residue::create(const std::string& name, int seq_num, const std::string& chain_id) {
    return Builder(name, seq_num, chain_id);
}

// Inline implementation of create_from_atoms()
inline Residue Residue::create_from_atoms(const std::string& name, int sequence_number, const std::string& chain_id,
                                          const std::string& insertion_code, const std::vector<Atom>& atoms) {
    // Helper to trim whitespace
    auto trim = [](const std::string& str) -> std::string {
        size_t start = str.find_first_not_of(" \t");
        if (start == std::string::npos)
            return "";
        size_t end = str.find_last_not_of(" \t");
        return str.substr(start, end - start + 1);
    };

    std::string trimmed = trim(name);

    // Get full classification from TypeRegistry
    auto classification = typing::TypeRegistry::instance().classify_residue(trimmed);

    // Get one-letter code from registry
    char one_letter = ModifiedNucleotideRegistry::get_one_letter_code(trimmed);

    // Build and return residue with trimmed name
    // Note: legacy_residue_idx is not set here - it should be set by the caller
    // after construction using set_legacy_residue_idx() if needed
    return Residue::create(trimmed, sequence_number, chain_id)
        .insertion(insertion_code)
        .one_letter_code(one_letter)
        .classification(classification)
        .atoms(atoms)
        .build();
}

} // namespace core
} // namespace x3dna
