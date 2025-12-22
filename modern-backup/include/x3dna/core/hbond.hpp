/**
 * @file hbond.hpp
 * @brief Hydrogen bond representation
 */

#pragma once

#include <optional>
#include <string>
#include <x3dna/core/hbond_types.hpp>

namespace x3dna {
namespace core {

/**
 * @brief Represents a hydrogen bond between two atoms
 *
 * Note on naming: donor_atom_name and acceptor_atom_name are PROVISIONAL names
 * based on detection order. The actual donor/acceptor roles are determined
 * during classification and reflected in the 'classification' field.
 * These names are kept for JSON compatibility with legacy output.
 */
class HBond {
public:
    // === Atom identification (names match JSON output) ===
    std::string donor_atom_name;    // First atom found (JSON: "donor_atom")
    std::string acceptor_atom_name; // Second atom found (JSON: "acceptor_atom")

    // === Core geometry ===
    double distance = 0.0; // D...A distance in Angstroms

    // === Geometric angles (always calculated, heavy atoms only) ===
    // X-D...A angle where X is heavy atom bonded to D
    double donor_angle = 0.0;
    std::string donor_neighbor_atom; // X atom used for donor_angle

    // D...A-Y angle where Y is heavy atom bonded to A
    double acceptor_angle = 0.0;
    std::string acceptor_neighbor_atom; // Y atom used for acceptor_angle

    // X-D...A-Y dihedral (0.0 if neighbors not found)
    double dihedral_angle = 0.0;
    bool dihedral_valid = false; // True if both neighbors found

    // === Classification ===
    HBondClassification classification = HBondClassification::UNKNOWN;
    HBondContext context = HBondContext::UNKNOWN;

    // === Conflict resolution state ===
    ConflictState conflict_state = ConflictState::NO_CONFLICT;

    // === Indices ===
    std::optional<size_t> detection_index; // Order detected (JSON: "hbond_idx")
    size_t donor_residue_index = 0;        // 0-based residue index
    size_t acceptor_residue_index = 0;     // 0-based residue index

    // === Legacy compatibility ===
    [[nodiscard]] char legacy_type_char() const;
    [[nodiscard]] int legacy_linkage_type() const;
};

} // namespace core
} // namespace x3dna
