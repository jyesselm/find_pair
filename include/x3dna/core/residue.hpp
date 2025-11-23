/**
 * @file residue.hpp
 * @brief Residue class representing a single residue in a PDB structure
 */

#pragma once

#include <string>
#include <vector>
#include <optional>
#include <nlohmann/json.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/reference_frame.hpp>

namespace x3dna {
namespace core {

/**
 * @enum ResidueType
 * @brief Type of residue (nucleotide, amino acid, etc.)
 */
enum class ResidueType {
    UNKNOWN = -2,
    AMINO_ACID = -1,
    NUCLEOTIDE = 0,
    ADENINE = 1,
    CYTOSINE = 2,
    GUANINE = 3,
    THYMINE = 4,
    URACIL = 5
};

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
        : name_(name), seq_num_(seq_num), chain_id_(chain_id), insertion_(insertion) {}

    // Getters
    const std::string& name() const {
        return name_;
    }
    int seq_num() const {
        return seq_num_;
    }
    char chain_id() const {
        return chain_id_;
    }
    char insertion() const {
        return insertion_;
    }
    const std::vector<Atom>& atoms() const {
        return atoms_;
    }
    size_t num_atoms() const {
        return atoms_.size();
    }

    /**
     * @brief Get reference frame (if set)
     */
    std::optional<ReferenceFrame> reference_frame() const {
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
    std::optional<Atom> find_atom(const std::string& atom_name) const {
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
    std::vector<Atom> ring_atoms() const {
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
     * @return One-letter code (A, C, G, T, U for nucleotides)
     *
     * Recognizes standard nucleotide names:
     * - Single letter: A, C, G, T, U
     * - Three letter: ADE, CYT, GUA, THY, URA
     * - DNA format: DA, DC, DG, DT (may have leading/trailing spaces)
     */
    char one_letter_code() const {
        std::string trimmed = trim_name();

        // Single letter codes
        if (trimmed == "A" || trimmed == "ADE" || trimmed == "DA")
            return 'A';
        if (trimmed == "C" || trimmed == "CYT" || trimmed == "DC")
            return 'C';
        if (trimmed == "G" || trimmed == "GUA" || trimmed == "DG")
            return 'G';
        if (trimmed == "T" || trimmed == "THY" || trimmed == "DT")
            return 'T';
        if (trimmed == "U" || trimmed == "URA" || trimmed == "DU")
            return 'U';

        return '?';
    }

    /**
     * @brief Check if this is a nucleotide
     * @return True if nucleotide (A, C, G, T, U)
     */
    bool is_nucleotide() const {
        char code = one_letter_code();
        return (code == 'A' || code == 'C' || code == 'G' || code == 'T' || code == 'U');
    }

    /**
     * @brief Get RY classification (Purine=1, Pyrimidine=0)
     * @return 1 for purines (A, G), 0 for pyrimidines (C, T, U), -1 for non-nucleotides
     */
    int ry_classification() const {
        char code = one_letter_code();
        if (code == 'A' || code == 'G') {
            return 1; // Purine
        }
        if (code == 'C' || code == 'T' || code == 'U') {
            return 0; // Pyrimidine
        }
        return -1; // Not a nucleotide
    }

    /**
     * @brief Get residue type
     */
    ResidueType residue_type() const {
        char code = one_letter_code();
        if (code == 'A')
            return ResidueType::ADENINE;
        if (code == 'C')
            return ResidueType::CYTOSINE;
        if (code == 'G')
            return ResidueType::GUANINE;
        if (code == 'T')
            return ResidueType::THYMINE;
        if (code == 'U')
            return ResidueType::URACIL;
        if (code == '?') {
            // Could be amino acid or unknown
            return ResidueType::UNKNOWN;
        }
        return ResidueType::UNKNOWN;
    }

    /**
     * @brief Convert to legacy JSON format
     * Note: Residues are typically represented as part of pdb_atoms records
     * This creates a minimal representation for testing
     */
    nlohmann::json to_json_legacy() const {
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
    static Residue from_json_legacy(const nlohmann::json& j) {
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
    nlohmann::json to_json() const {
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
    static Residue from_json(const nlohmann::json& j) {
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
    int seq_num_ = 0;                               // Sequence number
    char chain_id_ = '\0';                          // Chain identifier
    char insertion_ = ' ';                          // Insertion code (PDB column 27)
    std::vector<Atom> atoms_;                       // Atoms in this residue
    std::optional<ReferenceFrame> reference_frame_; // Reference frame (if calculated)

    /**
     * @brief Trim whitespace from residue name
     */
    std::string trim_name() const {
        std::string trimmed = name_;
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
        return trimmed;
    }
};

} // namespace core
} // namespace x3dna
