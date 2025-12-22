/**
 * @file overlap_calculator.hpp
 * @brief Calculator for base pair overlap area using polygon intersection
 *
 * Calculates the overlap area between two nucleotide bases by projecting
 * their ring atoms onto a common plane and computing polygon intersection.
 * This is used to determine how much two bases stack or overlap when
 * forming a base pair.
 */

#pragma once

#include <vector>
#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>

namespace x3dna {
namespace algorithms {
namespace validation {

// Forward declaration
class RingDataCache;

/**
 * @struct Point2D
 * @brief 2D point for polygon operations
 */
struct Point2D {
    double x;
    double y;
};

/**
 * @class OverlapCalculator
 * @brief Calculates overlap area between two nucleotide bases
 *
 * Projects ring atoms (with exocyclic substituents) onto a plane perpendicular
 * to the average z-axis and computes polygon intersection area using a
 * Sutherland-Hodgman style algorithm with integer arithmetic for precision.
 *
 * The overlap area helps identify stacking interactions vs. Watson-Crick pairing.
 * True base pairs typically have low overlap (bases are coplanar but not stacked).
 */
class OverlapCalculator {
public:
    /**
     * @brief Calculate overlap area between two residues
     * @param res1 First nucleotide residue
     * @param res2 Second nucleotide residue
     * @param average_origin Average origin of the two reference frames
     * @param average_z_axis Average z-axis of the two reference frames
     * @return Overlap area in square Angstroms
     *
     * Projects ring atoms (substituted with exocyclic atoms where bonded)
     * onto a plane perpendicular to the average z-axis, then computes
     * the intersection area of the two resulting polygons.
     */
    [[nodiscard]] static double calculate(const core::Residue& res1, const core::Residue& res2,
                                          const geometry::Vector3D& average_origin,
                                          const geometry::Vector3D& average_z_axis);

    /**
     * @brief Calculate overlap area using cached ring data (faster for batch processing)
     * @param res1 First nucleotide residue
     * @param res2 Second nucleotide residue
     * @param average_origin Average origin of the two reference frames
     * @param average_z_axis Average z-axis of the two reference frames
     * @param cache Ring data cache for O(1) ring coordinate lookup
     * @return Overlap area in square Angstroms
     *
     * Uses pre-computed ring atom indices and exocyclic mappings from the cache,
     * avoiding repeated O(n) lookups when the same residue appears in multiple pairs.
     */
    [[nodiscard]] static double calculate(const core::Residue& res1, const core::Residue& res2,
                                          const geometry::Vector3D& average_origin,
                                          const geometry::Vector3D& average_z_axis,
                                          RingDataCache& cache);

    /**
     * @brief Extract ring coordinates with exocyclic substituents for a residue
     * @param residue The nucleotide residue
     * @param average_origin Origin point (coordinates are returned relative to this)
     * @return Vector of 3D coordinates for ring positions (using exocyclic atoms where bonded)
     *
     * For each ring atom (C4, N3, C2, N1, C6, C5, N7, C8, N9), finds the closest
     * non-ring, non-hydrogen atom within bond distance. If found, uses that atom's
     * position instead of the ring atom. This expands the polygon to include
     * substituent groups like amino groups (N6 on adenine, N2 on guanine, etc.).
     */
    [[nodiscard]] static std::vector<geometry::Vector3D> get_ring_coordinates_with_exocyclic(
        const core::Residue& residue, const geometry::Vector3D& average_origin);

    /**
     * @brief Calculate intersection area of two 2D polygons
     * @param polygon_a First polygon vertices (in order)
     * @param polygon_b Second polygon vertices (in order)
     * @return Intersection area in square units
     *
     * Uses scaled integer arithmetic for robust edge intersection detection,
     * avoiding floating-point precision issues. The algorithm:
     * 1. Scales coordinates to a large integer range (GAMUT)
     * 2. Detects all edge-edge crossings between polygons
     * 3. Tracks winding numbers to identify interior regions
     * 4. Accumulates area using the shoelace formula
     */
    [[nodiscard]] static double calculate_polygon_intersection(const std::vector<Point2D>& polygon_a,
                                                               const std::vector<Point2D>& polygon_b);

private:
    /// Maximum polygon vertices supported (matches legacy MNPOLY constant)
    static constexpr long MAX_POLYGON_VERTICES = 1000;
};

}  // namespace validation
}  // namespace algorithms
}  // namespace x3dna
