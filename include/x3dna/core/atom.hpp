/**
 * @file atom.hpp
 * @brief Atom class representing a single atom in a PDB structure
 */

#pragma once

#include <string>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/core/constants.hpp>
#include <x3dna/core/string_utils.hpp>
#include <x3dna/core/typing/atom_classification.hpp>

namespace x3dna {
namespace core {

/**
 * @class Atom
 * @brief Represents a single atom with name, position, and metadata
 *
 * Atoms are primarily constructed via the Builder pattern for clean initialization.
 * After construction, atoms are largely immutable except for legacy index metadata
 * which may be set during post-processing.
 *
 * Note: Atom names are stored TRIMMED for consistent comparisons. The original
 * padded name is stored in original_atom_name_ for JSON output compatibility.
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
     * @param name Atom name (e.g., " C1'", " N3 ") - will be trimmed
     * @param position 3D position vector
     */
    Atom(const std::string& name, const geometry::Vector3D& position)
        : name_(trim(name)), original_atom_name_(name), position_(position) {}

    /**
     * @brief Constructor with full metadata
     * @param name Atom name - will be trimmed
     * @param position 3D position vector
     * @param residue_name Residue name (e.g., "  C", "  G") - will be trimmed
     * @param chain_id Chain identifier (e.g., "A", "B", "AA" for CIF) - will be trimmed
     * @param residue_seq Residue sequence number
     * @param record_type PDB record type (e.g., 'A' for ATOM, 'H' for HETATM)
     */
    Atom(const std::string& name, const geometry::Vector3D& position, const std::string& residue_name,
         const std::string& chain_id, int residue_seq, char record_type = 'A')
        : name_(trim(name)), original_atom_name_(name), position_(position),
          residue_name_(trim(residue_name)), original_residue_name_(residue_name),
          chain_id_(trim(chain_id)), residue_seq_(residue_seq), record_type_(record_type) {}

    /**
     * @brief Create a Builder for fluent atom construction
     * @param name Atom name (required)
     * @param position 3D position (required)
     * @return Builder instance
     */
    [[nodiscard]] static Builder create(const std::string& name, const geometry::Vector3D& position);

    // Getters (all const, returning const references where appropriate)

    /**
     * @brief Get trimmed atom name (for comparisons)
     */
    [[nodiscard]] const std::string& name() const {
        return name_;
    }

    /**
     * @brief Get original padded atom name (for JSON output)
     */
    [[nodiscard]] const std::string& original_atom_name() const {
        return original_atom_name_;
    }

    [[nodiscard]] const geometry::Vector3D& position() const {
        return position_;
    }

    /**
     * @brief Get trimmed residue name (for comparisons)
     */
    [[nodiscard]] const std::string& residue_name() const {
        return residue_name_;
    }

    /**
     * @brief Get original padded residue name (for JSON output)
     */
    [[nodiscard]] const std::string& original_residue_name() const {
        return original_residue_name_;
    }

    /**
     * @brief Get trimmed chain ID
     */
    [[nodiscard]] const std::string& chain_id() const {
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

    /**
     * @brief Get trimmed insertion code
     */
    [[nodiscard]] const std::string& insertion() const {
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
    [[nodiscard]] double b_factor() const {
        return b_factor_;
    }
    [[nodiscard]] const std::string& element() const {
        return element_;
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
        return typing::AtomClassifier::is_ring_atom(name_);
    }

    /**
     * @brief Check if this atom is a hydrogen bond donor
     * @return True if potential H-bond donor
     */
    [[nodiscard]] bool is_hydrogen_bond_donor() const {
        // Common H-bond donors: N with H (N1, N2, N3, N4, N6, N7, N9)
        return (name_.find("N") == 0 && name_.length() <= 2);
    }

    /**
     * @brief Check if this atom is a hydrogen bond acceptor
     * @return True if potential H-bond acceptor
     */
    [[nodiscard]] bool is_hydrogen_bond_acceptor() const {
        // Common H-bond acceptors: O (O2, O4, O6), N (N3, N7)
        return (name_.find("O") == 0 || name_ == "N3" || name_ == "N7");
    }

    /**
     * @brief Get full classification of this atom
     * @return AtomClassification with element, location, and H-bond role
     */
    [[nodiscard]] typing::AtomClassification classification() const {
        return typing::AtomClassifier::classify_nucleotide_atom(name_);
    }

    /**
     * @brief Check if this is a backbone atom (P, OP1, OP2, O5', O3', etc.)
     */
    [[nodiscard]] bool is_backbone_atom() const {
        return typing::AtomClassifier::is_backbone_atom(name_);
    }

    /**
     * @brief Check if this is a sugar atom (C1', C2', C3', C4', C5', O4', etc.)
     */
    [[nodiscard]] bool is_sugar_atom() const {
        return typing::AtomClassifier::is_sugar_atom(name_);
    }

    /**
     * @brief Check if this is a nucleobase atom (N1, C2, N3, C4, C5, C6, etc.)
     */
    [[nodiscard]] bool is_nucleobase_atom() const {
        return typing::AtomClassifier::is_nucleobase_atom(name_);
    }

private:
    friend class Builder;

    std::string name_;                  // Trimmed atom name (for comparisons)
    std::string original_atom_name_;    // Original atom name with padding (for JSON output)
    geometry::Vector3D position_;       // 3D coordinates
    std::string residue_name_;          // Trimmed residue name (for comparisons)
    std::string original_residue_name_; // Original residue name with padding (for JSON output)
    std::string chain_id_;              // Chain identifier (trimmed)
    int residue_seq_ = 0;               // Residue sequence number
    char record_type_ = 'A';            // PDB record type ('A' = ATOM, 'H' = HETATM)
    char alt_loc_ = ' ';                // Alternate location indicator (PDB column 17)
    std::string insertion_;             // Insertion code (trimmed)
    double occupancy_ = 1.0;            // Occupancy (PDB columns 55-60, default 1.0)
    int atom_serial_ = 0;               // Atom serial number (PDB column 7-11)
    int model_number_ = 0;              // Model number (from MODEL record, 0 if none)
    double b_factor_ = 0.0;             // B-factor/temperature factor (PDB column 61-66)
    std::string element_;               // Element symbol (PDB column 77-78)
    int legacy_atom_idx_ = 0;           // Legacy atom index for direct comparison (0 if not set)
    int legacy_residue_idx_ = 0;        // Legacy residue index for direct comparison (0 if not set)
};

/**
 * @class Atom::Builder
 * @brief Fluent builder for constructing Atom objects
 *
 * Usage:
 *   auto atom = Atom::create("CA", position)
 *       .residue_name("ALA")
 *       .chain_id("A")
 *       .residue_seq(42)
 *       .build();
 *
 * Note: All string values are automatically trimmed. Original values are preserved
 * for JSON output compatibility.
 */
class Atom::Builder {
public:
    /**
     * @brief Constructor with required fields
     * @param name Atom name (will be trimmed, original stored)
     * @param position 3D position
     */
    Builder(const std::string& name, const geometry::Vector3D& position) {
        atom_.name_ = trim(name);
        atom_.original_atom_name_ = name;
        atom_.position_ = position;
    }

    Builder& residue_name(const std::string& name) {
        atom_.residue_name_ = trim(name);
        atom_.original_residue_name_ = name;
        return *this;
    }

    Builder& chain_id(const std::string& id) {
        atom_.chain_id_ = trim(id);
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

    Builder& insertion(const std::string& ins) {
        atom_.insertion_ = trim(ins);
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

    Builder& b_factor(double bf) {
        atom_.b_factor_ = bf;
        return *this;
    }

    Builder& element(const std::string& elem) {
        atom_.element_ = elem;
        return *this;
    }

    /**
     * @brief Set original atom name explicitly (for parser normalization)
     *
     * This is used when the parser normalizes atom names (e.g., O1' → O4').
     * The normalized name should be set via the constructor, and the original
     * pre-normalization name can be set here.
     */
    Builder& original_atom_name(const std::string& name) {
        atom_.original_atom_name_ = name;
        return *this;
    }

    /**
     * @brief Set original residue name explicitly (for parser normalization)
     */
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
