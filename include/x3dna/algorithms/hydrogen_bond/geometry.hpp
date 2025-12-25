/**
 * @file geometry.hpp
 * @brief Geometric calculations for H-bonds using heavy atoms only
 */

#pragma once

#include <optional>
#include <string>
#include <x3dna/algorithms/hydrogen_bond/hbond_types.hpp>
#include <x3dna/core/typing/molecule_type.hpp>

// Forward declarations
namespace x3dna {
namespace geometry {
class Vector3D;
}
namespace core {
class Residue;
}
} // namespace x3dna

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {

/**
 * @brief Geometric calculations for H-bonds using heavy atoms only
 *
 * Provides static utilities for angle calculations, neighbor atom lookup,
 * and context determination. Works with trimmed atom names.
 */
class HBondGeometry {
public:
    // === Angle Calculations ===

    /**
     * @brief Calculate angle at vertex B for points A-B-C
     * @return Angle in degrees [0, 180]
     */
    [[nodiscard]] static double calculate_angle(const geometry::Vector3D& a, const geometry::Vector3D& b,
                                                const geometry::Vector3D& c);

    /**
     * @brief Calculate dihedral angle for points A-B-C-D
     * @return Dihedral in degrees [-180, 180]
     */
    [[nodiscard]] static double calculate_dihedral(const geometry::Vector3D& a, const geometry::Vector3D& b,
                                                   const geometry::Vector3D& c, const geometry::Vector3D& d);

    // === Neighbor Lookup ===

    /**
     * @brief Get reference neighbor atom for angle calculation
     * @param hbond_atom_name H-bond capable atom (trimmed, e.g., "N6", "O2'")
     * @return Neighbor atom name (trimmed), or empty if unknown
     */
    [[nodiscard]] static std::string get_neighbor_atom_name(const std::string& hbond_atom_name);

    /**
     * @brief Find neighbor atom position from residue
     * @return Position if found, nullopt otherwise
     */
    [[nodiscard]] static std::optional<geometry::Vector3D> find_neighbor_position(const std::string& hbond_atom_name,
                                                                                  const core::Residue& residue);

    // === Context Determination ===

    /**
     * @brief Determine H-bond context from two nucleotide atoms
     *
     * Uses AtomClassifier to classify atoms and determine context.
     */
    [[nodiscard]] static core::HBondContext determine_nucleotide_context(const std::string& atom1_name,
                                                                          const std::string& atom2_name);

    /**
     * @brief Determine H-bond context for any two atoms with molecule types
     *
     * Extended version that handles nucleic acid, protein, and ligand interactions.
     */
    [[nodiscard]] static core::HBondContext determine_context(const std::string& atom1_name,
                                                               const std::string& atom2_name,
                                                               core::typing::MoleculeType mol1_type,
                                                               core::typing::MoleculeType mol2_type);
};

} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
