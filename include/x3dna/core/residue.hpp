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
#include <nlohmann/json.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/residue_type.hpp>
#include <x3dna/core/modified_nucleotide_registry.hpp>

namespace x3dna {
namespace core {

/**
 * @class Residue
 * @brief Represents a single residue (nucleotide or amino acid) with atoms
 */
class Residue {
public:
    /**
     * @brief Default constructor
     */
    Residue() = default;

    /**
     * @brief Constructor with name, sequence number, chain ID, and insertion code
     * @param name Residue name (e.g., "  C", "  G", "  A", "  T")
     * @param seq_num Sequence number
     * @param chain_id Chain identifier
     * @param insertion Insertion code (PDB column 27, default ' ')
     */
    Residue(const std::string& name, int seq_num, char chain_id, char insertion = ' ')
        : name_(name), one_letter_code_('?'), type_(ResidueType::UNKNOWN), is_purine_(false), seq_num_(seq_num),
          chain_id_(chain_id), insertion_(insertion) {}

    /**
     * @brief Full constructor with all properties (used by ResidueFactory)
     * @param name Residue name
     * @param one_letter_code One-letter code ('A', 'C', 'G', 'U', etc.)
     * @param type Residue type
     * @param is_purine Whether this is a purine
     * @param seq_num Sequence number
     * @param chain_id Chain identifier
     * @param insertion Insertion code
     * @param atoms Vector of atoms
     */
    Residue(const std::string& name, char one_letter_code, ResidueType type, bool is_purine, int seq_num, char chain_id,
            char insertion, const std::vector<Atom>& atoms)
        : name_(name), one_letter_code_(one_letter_code), type_(type), is_purine_(is_purine), seq_num_(seq_num),
          chain_id_(chain_id), insertion_(insertion), atoms_(atoms) {}

    // Getters
    [[nodiscard]] const std::string& name() const {
        return name_;
    }
    [[nodiscard]] int seq_num() const {
        return seq_num_;
    }
    [[nodiscard]] char chain_id() const {
        return chain_id_;
    }
    [[nodiscard]] char insertion() const {
        return insertion_;
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

    /**
     * @brief Get reference frame (if set)
     */
    [[nodiscard]] std::optional<ReferenceFrame> reference_frame() const {
        return reference_frame_;
    }

    // Setters
    void set_name(const std::string& name) {
        name_ = name;
    }
    void set_seq_num(int seq_num) {
        seq_num_ = seq_num;
    }
    void set_chain_id(char chain_id) {
        chain_id_ = chain_id;
    }
    void set_insertion(char insertion) {
        insertion_ = insertion;
    }

    /**
     * @brief Set reference frame for this residue
     */
    void set_reference_frame(const ReferenceFrame& frame) {
        reference_frame_ = frame;
    }

    /**
     * @brief Add an atom to this residue
     */
    void add_atom(const Atom& atom) {
        atoms_.push_back(atom);
    }

    /**
     * @brief Find an atom by name
     * @param atom_name Atom name (e.g., " C1'", " N3 ")
     * @return Optional atom if found
     */
    [[nodiscard]] std::optional<Atom> find_atom(const std::string& atom_name) const {
        for (const auto& atom : atoms_) {
            if (atom.name() == atom_name) {
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
        char code = one_letter_code();
        // Standard nucleotides
        if (code == 'A' || code == 'C' || code == 'G' || code == 'T' || code == 'U' || code == 'I')
            return true;
        // Modified nucleotides (lowercase) and Pseudouridine (P)
        if (code == 'a' || code == 'c' || code == 'g' || code == 't' || code == 'u' || code == 'P')
            return true;
        return false;
    }

    /**
     * @brief Check if residue is a purine (stored value, not computed)
     * @return true if purine (A, G, I), false otherwise
     */
    [[nodiscard]] bool is_purine() const {
        return is_purine_;
    }

    /**
     * @brief Get RY classification (Purine=1, Pyrimidine=0)
     * Uses stored is_purine_ value
     * @return 1 for purines, 0 for pyrimidines, -1 for non-nucleotides
     */
    [[nodiscard]] int ry_classification() const {
        if (is_purine_) {
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
     * @brief Get residue type (stored value, not computed)
     * @return ResidueType set by ResidueFactory at creation
     */
    [[nodiscard]] ResidueType residue_type() const {
        return type_;
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
     * @brief Convert to legacy JSON format
     * Note: Residues are typically represented as part of pdb_atoms records
     * This creates a minimal representation for testing
     */
    [[nodiscard]] nlohmann::json to_json_legacy() const {
        nlohmann::json j;
        j["residue_name"] = name_;
        j["residue_seq"] = seq_num_;
        j["chain_id"] = std::string(1, chain_id_);
        j["atoms"] = nlohmann::json::array();
        for (const auto& atom : atoms_) {
            j["atoms"].push_back(atom.to_json_legacy());
        }
        if (reference_frame_.has_value()) {
            j["reference_frame"] = reference_frame_->to_json_legacy();
        }
        return j;
    }

    /**
     * @brief Create Residue from legacy JSON format
     */
    [[nodiscard]] static Residue from_json_legacy(const nlohmann::json& j) {
        std::string name = j.value("residue_name", "");
        int seq_num = j.value("residue_seq", 0);
        std::string chain_str = j.value("chain_id", "");
        char chain_id = chain_str.empty() ? '\0' : chain_str[0];

        Residue residue(name, seq_num, chain_id);

        if (j.contains("atoms") && j["atoms"].is_array()) {
            for (const auto& atom_json : j["atoms"]) {
                residue.add_atom(Atom::from_json_legacy(atom_json));
            }
        }

        if (j.contains("reference_frame")) {
            residue.set_reference_frame(ReferenceFrame::from_json_legacy(j["reference_frame"]));
        }

        return residue;
    }

    /**
     * @brief Convert to modern JSON format
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;
        j["name"] = name_;
        j["seq_num"] = seq_num_;
        j["chain_id"] = std::string(1, chain_id_);
        j["atoms"] = nlohmann::json::array();
        for (const auto& atom : atoms_) {
            j["atoms"].push_back(atom.to_json());
        }
        if (reference_frame_.has_value()) {
            j["reference_frame"] = reference_frame_->to_json();
        }
        return j;
    }

    /**
     * @brief Create Residue from modern JSON format
     */
    [[nodiscard]] static Residue from_json(const nlohmann::json& j) {
        std::string name = j.value("name", "");
        int seq_num = j.value("seq_num", 0);
        std::string chain_str = j.value("chain_id", "");
        char chain_id = chain_str.empty() ? '\0' : chain_str[0];

        Residue residue(name, seq_num, chain_id);

        if (j.contains("atoms") && j["atoms"].is_array()) {
            for (const auto& atom_json : j["atoms"]) {
                residue.add_atom(Atom::from_json(atom_json));
            }
        }

        if (j.contains("reference_frame")) {
            residue.set_reference_frame(ReferenceFrame::from_json(j["reference_frame"]));
        }

        return residue;
    }

private:
    std::string name_;                              // Residue name (e.g., "  C", "  G")
    char one_letter_code_ = '?';                    // One-letter code (stored, not computed)
    ResidueType type_ = ResidueType::UNKNOWN;       // Residue type (stored, not computed)
    bool is_purine_ = false;                        // Is purine flag (stored, not computed)
    int seq_num_ = 0;                               // Sequence number
    char chain_id_ = '\0';                          // Chain identifier
    char insertion_ = ' ';                          // Insertion code (PDB column 27)
    std::vector<Atom> atoms_;                       // Atoms in this residue
    std::optional<ReferenceFrame> reference_frame_; // Reference frame (if calculated)

    /**
     * @brief Trim whitespace from residue name
     */
    [[nodiscard]] std::string trim_name() const {
        std::string trimmed = name_;
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
        return trimmed;
    }
};

} // namespace core
} // namespace x3dna
