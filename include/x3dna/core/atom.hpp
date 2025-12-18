/**
 * @file atom.hpp
 * @brief Atom class representing a single atom in a PDB structure
 */

#pragma once

#include <string>
#include <nlohmann/json.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/core/constants.hpp>

namespace x3dna {
namespace core {

/**
 * @class Atom
 * @brief Represents a single atom with name, position, and metadata
 *
 * Atoms are primarily constructed via the Builder pattern for clean initialization.
 * After construction, atoms are largely immutable except for legacy index metadata
 * which may be set during post-processing.
 */
class Atom {
public:
    class Builder; // Forward declaration

    /**
     * @brief Default constructor (for container compatibility)
     */
    Atom() = default;

    /**
     * @brief Constructor with name and position
     * @param name Atom name (e.g., " C1'", " N3 ")
     * @param position 3D position vector
     */
    Atom(const std::string& name, const geometry::Vector3D& position) : name_(name), position_(position) {}

    /**
     * @brief Constructor with full metadata
     * @param name Atom name
     * @param position 3D position vector
     * @param residue_name Residue name (e.g., "  C", "  G")
     * @param chain_id Chain identifier (e.g., 'A', 'B')
     * @param residue_seq Residue sequence number
     * @param record_type PDB record type (e.g., 'A' for ATOM, 'H' for HETATM)
     */
    Atom(const std::string& name, const geometry::Vector3D& position, const std::string& residue_name, char chain_id,
         int residue_seq, char record_type = 'A')
        : name_(name), position_(position), residue_name_(residue_name), chain_id_(chain_id), residue_seq_(residue_seq),
          record_type_(record_type) {}

    /**
     * @brief Create a Builder for fluent atom construction
     * @param name Atom name (required)
     * @param position 3D position (required)
     * @return Builder instance
     */
    [[nodiscard]] static Builder create(const std::string& name, const geometry::Vector3D& position);

    // Getters (all const, returning const references where appropriate)
    [[nodiscard]] const std::string& name() const {
        return name_;
    }
    [[nodiscard]] const geometry::Vector3D& position() const {
        return position_;
    }
    [[nodiscard]] const std::string& residue_name() const {
        return residue_name_;
    }
    [[nodiscard]] char chain_id() const {
        return chain_id_;
    }
    [[nodiscard]] int residue_seq() const {
        return residue_seq_;
    }
    [[nodiscard]] char record_type() const {
        return record_type_;
    }
    [[nodiscard]] char alt_loc() const {
        return alt_loc_;
    }
    [[nodiscard]] char insertion() const {
        return insertion_;
    }
    [[nodiscard]] double occupancy() const {
        return occupancy_;
    }
    [[nodiscard]] int atom_serial() const {
        return atom_serial_;
    }
    [[nodiscard]] int model_number() const {
        return model_number_;
    }
    [[nodiscard]] size_t line_number() const {
        return line_number_;
    }
    [[nodiscard]] double b_factor() const {
        return b_factor_;
    }
    [[nodiscard]] const std::string& element() const {
        return element_;
    }
    [[nodiscard]] const std::string& original_atom_name() const {
        return original_atom_name_;
    }
    [[nodiscard]] const std::string& original_residue_name() const {
        return original_residue_name_;
    }
    [[nodiscard]] int legacy_atom_idx() const {
        return legacy_atom_idx_;
    }
    [[nodiscard]] int legacy_residue_idx() const {
        return legacy_residue_idx_;
    }

    // Post-construction setters for parsing workflow
    // These are retained because the values depend on filtering/ordering decisions
    // made after initial atom construction from PDB line data.

    /**
     * @brief Set model number (from MODEL record, set after parsing atom line)
     */
    void set_model_number(int model_number) {
        model_number_ = model_number;
    }

    /**
     * @brief Set legacy atom index (set after filtering decision)
     */
    void set_legacy_atom_idx(int legacy_atom_idx) {
        legacy_atom_idx_ = legacy_atom_idx;
    }

    /**
     * @brief Set legacy residue index (for post-construction index fixing)
     */
    void set_legacy_residue_idx(int legacy_residue_idx) {
        legacy_residue_idx_ = legacy_residue_idx;
    }

    /**
     * @brief Calculate distance to another atom
     * @param other Another atom
     * @return Distance in Angstroms
     */
    [[nodiscard]] double distance_to(const Atom& other) const {
        return position_.distance_to(other.position_);
    }

    /**
     * @brief Check if this atom is a ring atom (part of base ring)
     * @return True if ring atom
     */
    [[nodiscard]] bool is_ring_atom() const {
        return constants::nucleotides::is_ring_atom(name_);
    }

    /**
     * @brief Check if this atom is a hydrogen bond donor
     * @return True if potential H-bond donor
     */
    [[nodiscard]] bool is_hydrogen_bond_donor() const {
        const std::string trimmed = trim_name();
        // Common H-bond donors: N with H (N1, N2, N3, N4, N6, N7, N9)
        return (trimmed.find("N") == 0 && trimmed.length() <= 2);
    }

    /**
     * @brief Check if this atom is a hydrogen bond acceptor
     * @return True if potential H-bond acceptor
     */
    [[nodiscard]] bool is_hydrogen_bond_acceptor() const {
        const std::string trimmed = trim_name();
        // Common H-bond acceptors: O (O2, O4, O6), N (N3, N7)
        return (trimmed.find("O") == 0 || trimmed == "N3" || trimmed == "N7");
    }

    /**
     * @brief Convert to legacy JSON format (for pdb_atoms record)
     * @return JSON object matching legacy format
     */
    [[nodiscard]] nlohmann::json to_json_legacy() const {
        nlohmann::json j;
        // Always include all fields to match legacy format exactly
        j["atom_name"] = name_;
        j["xyz"] = {position_.x(), position_.y(), position_.z()};
        j["residue_name"] = residue_name_.empty() ? "" : residue_name_;
        j["chain_id"] = (chain_id_ == '\0' || chain_id_ == ' ') ? "" : std::string(1, chain_id_);
        j["residue_seq"] = residue_seq_;
        j["record_type"] = (record_type_ == '\0') ? "" : std::string(1, record_type_);
        if (alt_loc_ != ' ' && alt_loc_ != '\0') {
            j["alt_loc"] = std::string(1, alt_loc_);
        }
        if (insertion_ != ' ' && insertion_ != '\0') {
            j["insertion"] = std::string(1, insertion_);
        }

        // Add debug information (always include for debugging)
        j["occupancy"] = occupancy_;
        if (atom_serial_ > 0) {
            j["atom_serial"] = atom_serial_;
        }
        if (model_number_ > 0) {
            j["model_number"] = model_number_;
        }
        if (line_number_ > 0) {
            j["line_number"] = line_number_;
        }
        if (b_factor_ != 0.0) { // Include if non-zero
            j["b_factor"] = b_factor_;
        }
        if (!element_.empty()) {
            j["element"] = element_;
        }
        if (!original_atom_name_.empty() && original_atom_name_ != name_) {
            j["original_atom_name"] = original_atom_name_;
        }
        if (!original_residue_name_.empty() && original_residue_name_ != residue_name_) {
            j["original_residue_name"] = original_residue_name_;
        }

        return j;
    }

    /**
     * @brief Create Atom from legacy JSON format
     * @param j JSON object from pdb_atoms record
     * @return Atom object
     */
    [[nodiscard]] static Atom from_json_legacy(const nlohmann::json& j) {
        std::string name = j.value("atom_name", "");
        std::vector<double> xyz = j.value("xyz", std::vector<double>{0.0, 0.0, 0.0});
        geometry::Vector3D position(xyz[0], xyz[1], xyz[2]);

        std::string residue_name = j.value("residue_name", "");
        std::string chain_str = j.value("chain_id", "");
        char chain_id;
        if (chain_str.empty()) {
            chain_id = '\0';
        } else {
            chain_id = chain_str[0];
        }
        int residue_seq = j.value("residue_seq", 0);
        std::string record_str = j.value("record_type", "A");
        char record_type;
        if (record_str.empty()) {
            record_type = 'A';
        } else {
            record_type = record_str[0];
        }

        return Atom(name, position, residue_name, chain_id, residue_seq, record_type);
    }

    /**
     * @brief Convert to JSON format for pdb_atoms record
     * @return JSON object matching pdb_atoms format (atom_idx, atom_name, xyz, etc.)
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;

        // Use legacy_atom_idx for atom_idx to match legacy exactly (1-based)
        // Always include atom_idx (atoms should have legacy_atom_idx set during parsing)
        int legacy_atom_idx = legacy_atom_idx_;
        j["atom_idx"] = legacy_atom_idx > 0 ? legacy_atom_idx : 0;

        // Core fields (match legacy exactly)
        j["atom_name"] = name_;
        j["residue_name"] = residue_name_;
        j["chain_id"] = (chain_id_ == '\0' || chain_id_ == ' ') ? "" : std::string(1, chain_id_);
        j["residue_seq"] = residue_seq_;

        // Insertion code (only if non-space, matches legacy format)
        if (insertion_ != ' ' && insertion_ != '\0') {
            j["insertion"] = std::string(1, insertion_);
        }

        // Coordinates
        j["xyz"] = nlohmann::json::array({position_.x(), position_.y(), position_.z()});

        // Record type (A for ATOM, H for HETATM)
        j["record_type"] = std::string(1, record_type_);

        // Optional metadata (only if present)
        if (line_number_ > 0) {
            j["line_number"] = line_number_;
        }

        return j;
    }

    /**
     * @brief Create Atom from JSON format (matches to_json() output)
     */
    [[nodiscard]] static Atom from_json(const nlohmann::json& j) {
        std::string name = j.value("atom_name", "");
        std::vector<double> xyz = j.value("xyz", std::vector<double>{0.0, 0.0, 0.0});
        geometry::Vector3D position(xyz[0], xyz[1], xyz[2]);

        std::string residue_name = j.value("residue_name", "");
        std::string chain_str = j.value("chain_id", "");
        char chain_id;
        if (chain_str.empty()) {
            chain_id = '\0';
        } else {
            chain_id = chain_str[0];
        }
        int residue_seq = j.value("residue_seq", 0);
        std::string record_str = j.value("record_type", "A");
        char record_type;
        if (record_str.empty()) {
            record_type = 'A';
        } else {
            record_type = record_str[0];
        }

        return Atom(name, position, residue_name, chain_id, residue_seq, record_type);
    }

private:
    friend class Builder;

    std::string name_;                  // Atom name (e.g., " C1'", " N3 ")
    geometry::Vector3D position_;       // 3D coordinates
    std::string residue_name_;          // Residue name (e.g., "  C", "  G")
    char chain_id_ = '\0';              // Chain identifier
    int residue_seq_ = 0;               // Residue sequence number
    char record_type_ = 'A';            // PDB record type ('A' = ATOM, 'H' = HETATM)
    char alt_loc_ = ' ';                // Alternate location indicator (PDB column 17)
    char insertion_ = ' ';              // Insertion code (PDB column 27)
    double occupancy_ = 1.0;            // Occupancy (PDB columns 55-60, default 1.0)
    int atom_serial_ = 0;               // Atom serial number (PDB column 7-11)
    int model_number_ = 0;              // Model number (from MODEL record, 0 if none)
    size_t line_number_ = 0;            // Line number in PDB file (for debugging)
    double b_factor_ = 0.0;             // B-factor/temperature factor (PDB column 61-66)
    std::string element_;               // Element symbol (PDB column 77-78)
    std::string original_atom_name_;    // Original atom name before normalization
    std::string original_residue_name_; // Original residue name before normalization
    int legacy_atom_idx_ = 0;           // Legacy atom index for direct comparison (0 if not set)
    int legacy_residue_idx_ = 0;        // Legacy residue index for direct comparison (0 if not set)

    /**
     * @brief Trim whitespace from atom name for comparison
     */
    [[nodiscard]] std::string trim_name() const {
        std::string trimmed = name_;
        // Remove leading/trailing spaces
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
        return trimmed;
    }
};

/**
 * @class Atom::Builder
 * @brief Fluent builder for constructing Atom objects
 *
 * Usage:
 *   auto atom = Atom::create("CA", position)
 *       .residue_name("ALA")
 *       .chain_id('A')
 *       .residue_seq(42)
 *       .build();
 */
class Atom::Builder {
public:
    /**
     * @brief Constructor with required fields
     * @param name Atom name
     * @param position 3D position
     */
    Builder(const std::string& name, const geometry::Vector3D& position) {
        atom_.name_ = name;
        atom_.position_ = position;
    }

    Builder& residue_name(const std::string& name) {
        atom_.residue_name_ = name;
        return *this;
    }

    Builder& chain_id(char id) {
        atom_.chain_id_ = id;
        return *this;
    }

    Builder& residue_seq(int seq) {
        atom_.residue_seq_ = seq;
        return *this;
    }

    Builder& record_type(char type) {
        atom_.record_type_ = type;
        return *this;
    }

    Builder& alt_loc(char loc) {
        atom_.alt_loc_ = loc;
        return *this;
    }

    Builder& insertion(char ins) {
        atom_.insertion_ = ins;
        return *this;
    }

    Builder& occupancy(double occ) {
        atom_.occupancy_ = occ;
        return *this;
    }

    Builder& atom_serial(int serial) {
        atom_.atom_serial_ = serial;
        return *this;
    }

    Builder& model_number(int num) {
        atom_.model_number_ = num;
        return *this;
    }

    Builder& line_number(size_t num) {
        atom_.line_number_ = num;
        return *this;
    }

    Builder& b_factor(double bf) {
        atom_.b_factor_ = bf;
        return *this;
    }

    Builder& element(const std::string& elem) {
        atom_.element_ = elem;
        return *this;
    }

    Builder& original_atom_name(const std::string& name) {
        atom_.original_atom_name_ = name;
        return *this;
    }

    Builder& original_residue_name(const std::string& name) {
        atom_.original_residue_name_ = name;
        return *this;
    }

    Builder& legacy_atom_idx(int idx) {
        atom_.legacy_atom_idx_ = idx;
        return *this;
    }

    Builder& legacy_residue_idx(int idx) {
        atom_.legacy_residue_idx_ = idx;
        return *this;
    }

    /**
     * @brief Build and return the constructed Atom
     * @return Constructed Atom object
     */
    [[nodiscard]] Atom build() const {
        return atom_;
    }

private:
    Atom atom_;
};

// Inline implementation of create()
inline Atom::Builder Atom::create(const std::string& name, const geometry::Vector3D& position) {
    return Builder(name, position);
}

} // namespace core
} // namespace x3dna
