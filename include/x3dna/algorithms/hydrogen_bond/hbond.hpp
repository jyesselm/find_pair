/**
 * @file hbond.hpp
 * @brief Hydrogen bond representation
 */

#pragma once

#include <optional>
#include <string>
#include <x3dna/algorithms/hydrogen_bond/hbond_types.hpp>
#include <x3dna/algorithms/hydrogen_bond/hbond_quality.hpp>

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
    // === Atom identification ===
    std::string donor_atom_name;    // First atom (provisional name)
    std::string acceptor_atom_name; // Second atom (provisional name)

    // === Residue information ===
    size_t donor_residue_idx = 0;
    size_t acceptor_residue_idx = 0;
    std::string donor_res_id;    // e.g., "A-G-1"
    std::string acceptor_res_id; // e.g., "A-C-12"

    // === Core geometry ===
    double distance = 0.0;

    // === Angles (heavy atoms only) ===
    double donor_angle = 0.0; // X-D...A angle
    std::string donor_neighbor_atom;
    double acceptor_angle = 0.0; // D...A-Y angle
    std::string acceptor_neighbor_atom;
    double dihedral_angle = 0.0;
    bool dihedral_valid = false;

    // === Classification ===
    HBondClassification classification = HBondClassification::UNKNOWN;
    HBondContext context = HBondContext::UNKNOWN;
    ConflictState conflict_state = ConflictState::NO_CONFLICT;

    // === Leontis-Westhof edge classification ===
    BaseEdge donor_edge = BaseEdge::UNKNOWN;    // Which edge the donor atom is on
    BaseEdge acceptor_edge = BaseEdge::UNKNOWN; // Which edge the acceptor atom is on

    // === Detection metadata ===
    std::optional<size_t> detection_index;

    // === Quality scoring (optional - populated by HBondQualityScorer) ===
    std::optional<HBondQualityScore> quality_score;

    // === Legacy compatibility ===
    [[nodiscard]] char legacy_type_char() const {
        return to_legacy_char(classification);
    }

    [[nodiscard]] int legacy_linkage_type() const;

    [[nodiscard]] bool is_valid() const {
        return classification == HBondClassification::STANDARD || classification == HBondClassification::NON_STANDARD;
    }

    [[nodiscard]] bool is_standard() const {
        return classification == HBondClassification::STANDARD;
    }
};

} // namespace core
} // namespace x3dna
