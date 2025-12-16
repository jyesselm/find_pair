/**
 * @file atom.hpp
 * @brief Atom class representing a single atom in a PDB structure
 */

#pragma once

#include <string>
#include <nlohmann/json.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace core {

/**
 * @class Atom
 * @brief Represents a single atom with name, position, and metadata
 */
class Atom {
public:
    /**
     * @brief Default constructor
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

    // Getters
    const std::string& name() const {
        return name_;
    }
    const geometry::Vector3D& position() const {
        return position_;
    }
    const std::string& residue_name() const {
        return residue_name_;
    }
    char chain_id() const {
        return chain_id_;
    }
    int residue_seq() const {
        return residue_seq_;
    }
    char record_type() const {
        return record_type_;
    }
    char alt_loc() const {
        return alt_loc_;
    }
    char insertion() const {
        return insertion_;
    }
    double occupancy() const {
        return occupancy_;
    }
    int atom_serial() const {
        return atom_serial_;
    }
    int model_number() const {
        return model_number_;
    }
    size_t line_number() const {
        return line_number_;
    }
    double b_factor() const {
        return b_factor_;
    }
    const std::string& element() const {
        return element_;
    }
    const std::string& original_atom_name() const {
        return original_atom_name_;
    }
    const std::string& original_residue_name() const {
        return original_residue_name_;
    }
    int legacy_atom_idx() const {
        return legacy_atom_idx_;
    }
    int legacy_residue_idx() const {
        return legacy_residue_idx_;
    }

    // Setters
    void set_name(const std::string& name) {
        name_ = name;
    }
    void set_position(const geometry::Vector3D& position) {
        position_ = position;
    }
    void set_residue_name(const std::string& residue_name) {
        residue_name_ = residue_name;
    }
    void set_chain_id(char chain_id) {
        chain_id_ = chain_id;
    }
    void set_residue_seq(int residue_seq) {
        residue_seq_ = residue_seq;
    }
    void set_record_type(char record_type) {
        record_type_ = record_type;
    }
    void set_alt_loc(char alt_loc) {
        alt_loc_ = alt_loc;
    }
    void set_insertion(char insertion) {
        insertion_ = insertion;
    }
    void set_occupancy(double occupancy) {
        occupancy_ = occupancy;
    }
    void set_atom_serial(int atom_serial) {
        atom_serial_ = atom_serial;
    }
    void set_model_number(int model_number) {
        model_number_ = model_number;
    }
    void set_line_number(size_t line_number) {
        line_number_ = line_number;
    }
    void set_b_factor(double b_factor) {
        b_factor_ = b_factor;
    }
    void set_element(const std::string& element) {
        element_ = element;
    }
    void set_original_atom_name(const std::string& name) {
        original_atom_name_ = name;
    }
    void set_original_residue_name(const std::string& name) {
        original_residue_name_ = name;
    }
    void set_legacy_atom_idx(int legacy_atom_idx) {
        legacy_atom_idx_ = legacy_atom_idx;
    }
    void set_legacy_residue_idx(int legacy_residue_idx) {
        legacy_residue_idx_ = legacy_residue_idx;
    }

    /**
     * @brief Calculate distance to another atom
     * @param other Another atom
     * @return Distance in Angstroms
     */
    double distance_to(const Atom& other) const {
        return position_.distance_to(other.position_);
    }

    /**
     * @brief Check if this atom is a ring atom (part of base ring)
     * @return True if ring atom
     */
    bool is_ring_atom() const {
        // Common ring atoms in nucleic acids
        const std::string trimmed = trim_name();
        return (trimmed == "N1" || trimmed == "C2" || trimmed == "N3" || trimmed == "C4" || trimmed == "C5" ||
                trimmed == "C6" || trimmed == "N7" || trimmed == "C8" || trimmed == "N9");
    }

    /**
     * @brief Check if this atom is a hydrogen bond donor
     * @return True if potential H-bond donor
     */
    bool is_hydrogen_bond_donor() const {
        const std::string trimmed = trim_name();
        // Common H-bond donors: N with H (N1, N2, N3, N4, N6, N7, N9)
        return (trimmed.find("N") == 0 && trimmed.length() <= 2);
    }

    /**
     * @brief Check if this atom is a hydrogen bond acceptor
     * @return True if potential H-bond acceptor
     */
    bool is_hydrogen_bond_acceptor() const {
        const std::string trimmed = trim_name();
        // Common H-bond acceptors: O (O2, O4, O6), N (N3, N7)
        return (trimmed.find("O") == 0 || trimmed == "N3" || trimmed == "N7");
    }

    /**
     * @brief Convert to legacy JSON format (for pdb_atoms record)
     * @return JSON object matching legacy format
     */
    nlohmann::json to_json_legacy() const {
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
    static Atom from_json_legacy(const nlohmann::json& j) {
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
    nlohmann::json to_json() const {
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
     * @brief Create Atom from modern JSON format
     */
    static Atom from_json(const nlohmann::json& j) {
        std::string name = j.value("name", "");
        geometry::Vector3D position = geometry::Vector3D::from_json(j["position"]);

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
    std::string trim_name() const {
        std::string trimmed = name_;
        // Remove leading/trailing spaces
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
        return trimmed;
    }
};

} // namespace core
} // namespace x3dna
