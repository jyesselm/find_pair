/**
 * @file geometry.hpp
 * @brief Geometric calculations for H-bonds using heavy atoms only
 */

#pragma once

#include <optional>
#include <string>
#include <x3dna/core/hbond_types.hpp>

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

/**
 * @brief Geometric calculations for H-bonds using heavy atoms only
 *
 * Provides static utilities for angle calculations, neighbor atom lookup,
 * and atom classification (base/backbone/sugar).
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
     * @param hbond_atom_name H-bond capable atom (e.g., " N6 ")
     * @return Neighbor atom name, or empty if unknown
     */
    [[nodiscard]] static std::string get_neighbor_atom_name(const std::string& hbond_atom_name);

    /**
     * @brief Find neighbor atom position from residue
     * @return Position if found, nullopt otherwise
     */
    [[nodiscard]] static std::optional<geometry::Vector3D> find_neighbor_position(const std::string& hbond_atom_name,
                                                                                  const core::Residue& residue);

    // === Atom Classification ===

    /**
     * @brief Check if atoms can form H-bond based on elements
     */
    [[nodiscard]] static bool are_elements_hbond_compatible(const std::string& atom1_name,
                                                            const std::string& atom2_name,
                                                            const std::string& allowed_elements = ".O.N");

    /**
     * @brief Check if atom is a nucleobase atom (not backbone/sugar)
     */
    [[nodiscard]] static bool is_nucleobase_atom(const std::string& atom_name);

    /**
     * @brief Check if atom is a backbone atom (phosphate)
     */
    [[nodiscard]] static bool is_backbone_atom(const std::string& atom_name);

    /**
     * @brief Check if atom is a sugar atom (ribose)
     */
    [[nodiscard]] static bool is_sugar_atom(const std::string& atom_name);

    /**
     * @brief Determine H-bond context from atom names
     */
    [[nodiscard]] static core::HBondContext determine_context(const std::string& atom1_name,
                                                              const std::string& atom2_name);
};

} // namespace algorithms
} // namespace x3dna
