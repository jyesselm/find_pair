/**
 * @file overlap_calculator.hpp
 * @brief Calculator for base pair overlap area using polygon intersection
 *
 * Calculates the overlap area between two nucleotide bases by projecting
 * their ring atoms onto a common plane and computing polygon intersection.
 */

#pragma once

#include <x3dna/core/residue.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <vector>

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
 * Projects ring atoms (with exocyclic atoms) onto a plane perpendicular
 * to the average z-axis and computes polygon intersection area.
 *
 * Matches legacy get_oarea() and pia_inter() functions.
 */
class OverlapCalculator {
public:
    /**
     * @brief Calculate overlap area between two residues
     * @param res1 First residue
     * @param res2 Second residue
     * @param oave Average origin of the two reference frames
     * @param zave Average z-axis of the two reference frames
     * @return Overlap area in square Angstroms
     */
    [[nodiscard]] static double calculate(const core::Residue& res1, const core::Residue& res2,
                                          const geometry::Vector3D& oave, const geometry::Vector3D& zave);

    /**
     * @brief Calculate overlap area using cached ring data (faster)
     * @param res1 First residue
     * @param res2 Second residue
     * @param oave Average origin of the two reference frames
     * @param zave Average z-axis of the two reference frames
     * @param cache Ring data cache for O(1) ring coordinate lookup
     * @return Overlap area in square Angstroms
     */
    [[nodiscard]] static double calculate(const core::Residue& res1, const core::Residue& res2,
                                          const geometry::Vector3D& oave, const geometry::Vector3D& zave,
                                          RingDataCache& cache);

    /**
     * @brief Get ring coordinates with exocyclic atoms for a residue
     * @param residue The nucleotide residue
     * @param oave Average origin (coordinates are relative to this)
     * @return Vector of 3D coordinates for ring + exocyclic atoms
     */
    [[nodiscard]] static std::vector<geometry::Vector3D> get_ring_coordinates_with_exocyclic(
        const core::Residue& residue, const geometry::Vector3D& oave);

    /**
     * @brief Calculate polygon intersection area
     * @param poly1 First polygon vertices
     * @param poly2 Second polygon vertices
     * @return Intersection area
     *
     * Uses the pia_inter algorithm (Polygon Intersection Area) for
     * accurate computation of overlapping polygon regions.
     */
    [[nodiscard]] static double calculate_polygon_intersection(const std::vector<Point2D>& poly1,
                                                               const std::vector<Point2D>& poly2);

private:
    // Maximum polygon vertices (matches legacy MNPOLY)
    static constexpr long MAX_POLYGON_VERTICES = 1000;
};

} // namespace validation
} // namespace algorithms
} // namespace x3dna
