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
 * Note: Atom names are trimmed on construction and stored without padding.
 * The original PDB 4-character format is not preserved (no longer needed after removing original_atom_name_).
 */
class Atom {
public:
    class Builder; // Forward declaration

    /**
     * @brief Default constructor (for container compatibility)
     */
    Atom() = default;

    Atom(const Atom& other) = default;
    Atom(Atom&& other) noexcept = default;
    Atom& operator=(const Atom& other) = default;
    Atom& operator=(Atom&& other) noexcept = default;

    /**
     * @brief Constructor with name and position
     * @param name Atom name (will be trimmed, e.g., " C1'" becomes "C1'", " N3 " becomes "N3")
     * @param position 3D position vector
     *
     * Also classifies the atom type at construction time for O(1) lookup later.
     */
    Atom(const std::string& name, const geometry::Vector3D& position)
        : name_(trim(name)),
          position_(position),
          standard_atom_(typing::AtomClassifier::get_atom_type(name_)) {}

    /**
     * @brief Create a Builder for fluent atom construction
     * @param name Atom name (required)
     * @param position 3D position (required)
     * @return Builder instance
     */
    [[nodiscard]] static Builder create(const std::string& name, const geometry::Vector3D& position);

    // Getters (all const, returning const references where appropriate)

    /**
     * @brief Get atom name (trimmed, without padding)
     */
    [[nodiscard]] const std::string& name() const {
        return name_;
    }

    /**
     * @brief Check if atom name matches (case-sensitive, handles both trimmed and padded input)
     * @param name_to_match Name to match against (can be trimmed like "C1'" or padded like " C1'")
     * @return True if names match after trimming
     */
    [[nodiscard]] bool name_matches(const std::string& name_to_match) const {
        return name_ == trim(name_to_match);
    }

    [[nodiscard]] const geometry::Vector3D& position() const {
        return position_;
    }

    [[nodiscard]] char alt_loc() const {
        return alt_loc_;
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
     * @brief Calculate distance to another atom
     * @param other Another atom
     * @return Distance in Angstroms
     */
    [[nodiscard]] double distance_to(const Atom& other) const {
        return position_.distance_to(other.position_);
    }

    /**
     * @brief Get the standard atom type (for fast enum comparison)
     * @return AtomType enum value (UNKNOWN for non-standard atoms)
     */
    [[nodiscard]] AtomType atom_type() const {
        return standard_atom_;
    }

    /**
     * @brief Alias for atom_type() - deprecated, use atom_type() instead
     */
    [[nodiscard]] AtomType standard_atom() const {
        return standard_atom_;
    }

    /**
     * @brief Check if this atom matches a specific AtomType (O(1) comparison)
     * @param type The AtomType to check against
     * @return True if atom type matches
     */
    [[nodiscard]] bool is(AtomType type) const {
        return standard_atom_ == type;
    }

    // === Specific atom type checks (O(1) comparison) ===

    [[nodiscard]] bool is_c1_prime() const { return standard_atom_ == AtomType::C1_PRIME; }
    [[nodiscard]] bool is_c2_prime() const { return standard_atom_ == AtomType::C2_PRIME; }
    [[nodiscard]] bool is_c3_prime() const { return standard_atom_ == AtomType::C3_PRIME; }
    [[nodiscard]] bool is_c4_prime() const { return standard_atom_ == AtomType::C4_PRIME; }
    [[nodiscard]] bool is_c5_prime() const { return standard_atom_ == AtomType::C5_PRIME; }
    [[nodiscard]] bool is_o2_prime() const { return standard_atom_ == AtomType::O2_PRIME; }
    [[nodiscard]] bool is_o3_prime() const { return standard_atom_ == AtomType::O3_PRIME; }
    [[nodiscard]] bool is_o4_prime() const { return standard_atom_ == AtomType::O4_PRIME; }
    [[nodiscard]] bool is_o5_prime() const { return standard_atom_ == AtomType::O5_PRIME; }

    // Ring atoms
    [[nodiscard]] bool is_n1() const { return standard_atom_ == AtomType::N1; }
    [[nodiscard]] bool is_n3() const { return standard_atom_ == AtomType::N3; }
    [[nodiscard]] bool is_n7() const { return standard_atom_ == AtomType::N7; }
    [[nodiscard]] bool is_n9() const { return standard_atom_ == AtomType::N9; }
    [[nodiscard]] bool is_c2() const { return standard_atom_ == AtomType::C2; }
    [[nodiscard]] bool is_c4() const { return standard_atom_ == AtomType::C4; }
    [[nodiscard]] bool is_c5() const { return standard_atom_ == AtomType::C5; }
    [[nodiscard]] bool is_c6() const { return standard_atom_ == AtomType::C6; }
    [[nodiscard]] bool is_c8() const { return standard_atom_ == AtomType::C8; }

    // Backbone atoms
    [[nodiscard]] bool is_phosphorus() const { return standard_atom_ == AtomType::P; }

    // === Atom type update ===

    /**
     * @brief Update atom type based on molecule context
     * @param molecule_type The type of molecule this atom belongs to
     *
     * This should be called after residue classification is known to ensure
     * atom types are correctly assigned based on context. For example, an atom
     * named "N7" will only get AtomType::N7 if the molecule is a nucleic acid.
     */
    void update_atom_type(typing::MoleculeType molecule_type) {
        standard_atom_ = typing::AtomClassifier::get_atom_type(name_, molecule_type);
    }

    /**
     * @brief Check if this atom is a ring atom (part of base ring)
     * @return True if ring atom
     *
     * Uses the cached standard_atom_ for O(1) lookup instead of string comparison.
     */
    [[nodiscard]] bool is_ring_atom() const {
        return typing::is_ring_atom(standard_atom_);
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

    std::string name_;            // Atom name (trimmed, without padding)
    geometry::Vector3D position_; // 3D coordinates
    AtomType standard_atom_ = AtomType::UNKNOWN; // Cached atom type for fast comparison
    char alt_loc_ = ' ';          // Alternate location indicator (PDB column 17)
    double occupancy_ = 1.0;      // Occupancy (PDB columns 55-60, default 1.0)
    int atom_serial_ = 0;         // Atom serial number (PDB column 7-11)
    int model_number_ = 0;        // Model number (from MODEL record, 0 if none)
    double b_factor_ = 0.0;       // B-factor/temperature factor (PDB column 61-66)
    std::string element_;         // Element symbol (PDB column 77-78)
    int legacy_atom_idx_ = 0;     // Legacy atom index for direct comparison (0 if not set)
};

/**
 * @class Atom::Builder
 * @brief Fluent builder for constructing Atom objects
 *
 * Usage:
 *   auto atom = Atom::create(" CA ", position)
 *       .alt_loc('A')
 *       .occupancy(1.0)
 *       .build();
 *
 * Note: Atom names are trimmed on construction, so " CA " becomes "CA".
 */
class Atom::Builder {
public:
    /**
     * @brief Constructor with required fields
     * @param name Atom name (will be trimmed)
     * @param position 3D position
     *
     * Also classifies the atom type at construction time for O(1) lookup later.
     */
    Builder(const std::string& name, const geometry::Vector3D& position) {
        atom_.name_ = trim(name);
        atom_.position_ = position;
        atom_.standard_atom_ = typing::AtomClassifier::get_atom_type(atom_.name_);
    }

    Builder& alt_loc(char loc) {
        atom_.alt_loc_ = loc;
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

    Builder& legacy_atom_idx(int idx) {
        atom_.legacy_atom_idx_ = idx;
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
