/**
 * @file overlap_calculator.cpp
 * @brief Calculates overlap area between nucleotide bases using polygon intersection
 *
 * This module computes how much two nucleotide bases overlap when projected onto
 * a common plane perpendicular to their average z-axis. The algorithm:
 *
 * 1. Extracts ring atoms (with exocyclic substituents) for each base
 * 2. Projects 3D coordinates onto a 2D plane perpendicular to the average z-axis
 * 3. Computes polygon intersection area using integer arithmetic for precision
 *
 * The polygon intersection algorithm is adapted from the legacy X3DNA pia_inter
 * (Polygon Intersection Area) implementation in ana_fncs.c.
 */

#include <x3dna/algorithms/validation/overlap_calculator.hpp>
#include <x3dna/algorithms/validation/ring_data_cache.hpp>
#include <x3dna/algorithms/validation_constants.hpp>
#include <x3dna/core/atom.hpp>
#include <algorithm>
#include <cmath>
#include <set>

namespace x3dna {
namespace algorithms {
namespace validation {

// ============================================================================
// Polygon Intersection Data Structures
// ============================================================================

/**
 * @brief Vertex data for polygon intersection computation
 *
 * Each vertex stores both its scaled integer coordinates (for precise
 * intersection detection) and range information for its outgoing edge.
 */
struct PolygonVertex {
    Point2D scaled_position;  // Coordinates scaled to integer range for precision
    Point2D edge_x_range;     // Min/max x values along edge to next vertex
    Point2D edge_y_range;     // Min/max y values along edge to next vertex
    long winding_contribution;  // Contribution to winding number calculation
};

// ============================================================================
// Polygon Intersection Algorithm Helpers
// ============================================================================
//
// These functions implement the Sutherland-Hodgman style polygon intersection
// algorithm. The approach uses scaled integer coordinates to avoid floating
// point precision issues when detecting edge crossings.
//
// Algorithm overview:
// 1. Scale polygon coordinates to a large integer range (GAMUT)
// 2. Find all edge-edge intersections between the two polygons
// 3. Track winding numbers to determine which regions are inside both polygons
// 4. Accumulate signed area contributions from intersection points
// ============================================================================

namespace {

/**
 * @brief Compute signed area of triangle formed by point and edge
 * @param point Reference point
 * @param edge_start Start of edge
 * @param edge_end End of edge
 * @return Signed area (positive if point is left of edge, negative if right)
 *
 * Used to determine which side of an edge a point lies on, and to detect
 * edge crossings. The sign indicates orientation.
 */
double signed_triangle_area(const Point2D& point, const Point2D& edge_start, const Point2D& edge_end) {
    return edge_start.x * edge_end.y - edge_start.y * edge_end.x +
           point.x * (edge_start.y - edge_end.y) +
           point.y * (edge_end.x - edge_start.x);
}

/**
 * @brief Add trapezoid area contribution to running sum
 * @param area_sum Running area sum (modified in place)
 * @param from Start point of edge
 * @param to End point of edge
 * @param weight Winding number weight (+1 or -1)
 *
 * Computes the signed area of a trapezoid formed by projecting an edge
 * down to the x-axis. Used in the shoelace formula for polygon area.
 */
void add_trapezoid_area(double& area_sum, const Point2D& from, const Point2D& to, long weight) {
    area_sum += weight * (to.x - from.x) * (to.y + from.y) * 0.5;
}

/**
 * @brief Check if two 1D ranges overlap
 * @param range_a First range (x = min, y = max)
 * @param range_b Second range (x = min, y = max)
 * @return True if ranges overlap (excluding exact boundary touch)
 *
 * Used as a fast rejection test before checking for actual edge intersection.
 */
bool ranges_overlap(const Point2D& range_a, const Point2D& range_b) {
    return range_a.x < range_b.y && range_b.x < range_a.y;
}

/**
 * @brief Handle edge crossing between two polygon edges
 * @param area_sum Running area sum (modified)
 * @param edge_a_start First edge start vertex
 * @param edge_a_end First edge end vertex
 * @param edge_b_start Second edge start vertex
 * @param edge_b_end Second edge end vertex
 * @param signed_area_a1 Signed area: edge_a_start relative to edge_b
 * @param signed_area_a2 Signed area: edge_a_end relative to edge_b
 * @param signed_area_b1 Signed area: edge_b_start relative to edge_a
 * @param signed_area_b2 Signed area: edge_b_end relative to edge_a
 *
 * When two edges cross, this computes the intersection point and adds
 * the appropriate area contributions. Also updates winding contributions.
 */
void handle_edge_crossing(double& area_sum,
                          PolygonVertex* edge_a_start, PolygonVertex* edge_a_end,
                          PolygonVertex* edge_b_start, PolygonVertex* edge_b_end,
                          double signed_area_a1, double signed_area_a2,
                          double signed_area_b1, double signed_area_b2) {
    // Compute intersection parameters using linear interpolation
    double t_on_edge_a = signed_area_a1 / (signed_area_a1 + signed_area_a2);
    double t_on_edge_b = signed_area_b1 / (signed_area_b1 + signed_area_b2);

    // Compute intersection point on edge A
    Point2D intersection_on_a;
    intersection_on_a.x = edge_a_start->scaled_position.x +
                          t_on_edge_a * (edge_a_end->scaled_position.x - edge_a_start->scaled_position.x);
    intersection_on_a.y = edge_a_start->scaled_position.y +
                          t_on_edge_a * (edge_a_end->scaled_position.y - edge_a_start->scaled_position.y);

    // Add area contribution from edge A segment
    add_trapezoid_area(area_sum, intersection_on_a, edge_a_end->scaled_position, 1);

    // Compute intersection point on edge B
    Point2D intersection_on_b;
    intersection_on_b.x = edge_b_start->scaled_position.x +
                          t_on_edge_b * (edge_b_end->scaled_position.x - edge_b_start->scaled_position.x);
    intersection_on_b.y = edge_b_start->scaled_position.y +
                          t_on_edge_b * (edge_b_end->scaled_position.y - edge_b_start->scaled_position.y);

    // Add area contribution from edge B segment
    add_trapezoid_area(area_sum, edge_b_end->scaled_position, intersection_on_b, 1);

    // Update winding contributions
    ++edge_a_start->winding_contribution;
    --edge_b_start->winding_contribution;
}

/**
 * @brief Compute area contribution from vertices inside the other polygon
 * @param area_sum Running area sum (modified)
 * @param polygon_p First polygon vertices
 * @param vertex_count_p Number of vertices in first polygon
 * @param polygon_q Second polygon vertices
 * @param vertex_count_q Number of vertices in second polygon
 *
 * Determines which vertices of polygon P are inside polygon Q using
 * a ray casting approach, then adds their area contributions weighted
 * by the winding number.
 */
void add_interior_vertex_contributions(double& area_sum,
                                        PolygonVertex* polygon_p, long vertex_count_p,
                                        PolygonVertex* polygon_q, long vertex_count_q) {
    // Determine if first vertex of P is inside Q using ray casting
    long winding_number = 0;
    Point2D test_point = polygon_p[0].scaled_position;

    for (long i = vertex_count_q - 1; i >= 0; --i) {
        // Check if horizontal ray from test_point crosses this edge
        if (polygon_q[i].edge_x_range.x < test_point.x &&
            test_point.x < polygon_q[i].edge_x_range.y) {
            // Determine crossing direction using signed area
            bool point_left_of_edge = signed_triangle_area(
                test_point, polygon_q[i].scaled_position, polygon_q[i + 1].scaled_position) > 0;
            bool edge_goes_right = polygon_q[i].scaled_position.x < polygon_q[i + 1].scaled_position.x;

            if (point_left_of_edge == edge_goes_right) {
                winding_number += point_left_of_edge ? -1 : 1;
            }
        }
    }

    // Add area contributions for each edge of P, weighted by winding number
    for (long j = 0; j < vertex_count_p; ++j) {
        if (winding_number != 0) {
            add_trapezoid_area(area_sum,
                              polygon_p[j].scaled_position,
                              polygon_p[j + 1].scaled_position,
                              winding_number);
        }
        winding_number += polygon_p[j].winding_contribution;
    }
}

/**
 * @brief Scale polygon coordinates to integer range and compute edge bounds
 * @param min_x Minimum x coordinate in original space
 * @param min_y Minimum y coordinate in original space
 * @param midpoint Midpoint of scaled range (GAMUT/2)
 * @param scale_x X scaling factor
 * @param scale_y Y scaling factor
 * @param input_vertices Original polygon vertices
 * @param vertex_count Number of vertices
 * @param output_vertices Scaled vertex data (output)
 * @param fudge_bits Bits to OR into coordinates (0 or 2) for numerical stability
 *
 * Transforms floating-point polygon coordinates to a large integer range
 * to enable precise edge intersection detection. Also precomputes the
 * x and y ranges for each edge (used for fast rejection testing).
 */
void scale_polygon_to_integer_coords(double min_x, double min_y,
                                     double midpoint, double scale_x, double scale_y,
                                     const Point2D* input_vertices, long vertex_count,
                                     PolygonVertex* output_vertices, long fudge_bits) {
    // Scale each vertex to integer coordinates
    for (long i = vertex_count - 1; i >= 0; --i) {
        long scaled_x = static_cast<long>((input_vertices[i].x - min_x) * scale_x - midpoint);
        long scaled_y = static_cast<long>((input_vertices[i].y - min_y) * scale_y - midpoint);

        // Apply bit manipulation for numerical stability:
        // - Clear lowest 3 bits (& ~7)
        // - OR in fudge bits and vertex parity bit
        output_vertices[i].scaled_position.x = static_cast<double>((scaled_x & ~7) | fudge_bits | (i & 1));
        output_vertices[i].scaled_position.y = static_cast<double>((scaled_y & ~7) | fudge_bits);
    }

    // Add parity adjustment to first vertex's y coordinate
    output_vertices[0].scaled_position.y += vertex_count & 1;

    // Copy first vertex to end for wraparound
    output_vertices[vertex_count] = output_vertices[0];

    // Compute edge ranges for fast bounding box rejection
    for (long i = vertex_count - 1; i >= 0; --i) {
        double x1 = output_vertices[i].scaled_position.x;
        double x2 = output_vertices[i + 1].scaled_position.x;
        double y1 = output_vertices[i].scaled_position.y;
        double y2 = output_vertices[i + 1].scaled_position.y;

        // Store ranges as (min, max) pairs
        output_vertices[i].edge_x_range = (x1 < x2) ? Point2D{x1, x2} : Point2D{x2, x1};
        output_vertices[i].edge_y_range = (y1 < y2) ? Point2D{y1, y2} : Point2D{y2, y1};
        output_vertices[i].winding_contribution = 0;
    }
}

// ============================================================================
// Ring Atom Extraction
// ============================================================================

/// Names of atoms that form the nucleotide ring system (purines have all 9, pyrimidines have first 6)
static const char* RING_ATOM_NAMES[] = {"C4", "N3", "C2", "N1", "C6", "C5", "N7", "C8", "N9"};
static constexpr size_t NUM_RING_ATOMS = 9;

/**
 * @brief Extract ring coordinates with exocyclic substituents for overlap calculation
 * @tparam ResidueType Type supporting find_atom_ptr() and atoms() methods
 * @param residue The nucleotide residue
 * @param average_origin Origin point to subtract (coordinates returned relative to this)
 * @return Vector of 3D coordinates for ring atoms (using exocyclic atoms where bonded)
 *
 * For each ring atom, finds the closest non-ring, non-hydrogen atom within bond
 * distance. If found, uses that exocyclic atom's position instead of the ring atom.
 * This expands the effective polygon to include substituent groups like amino groups.
 */
template <typename ResidueType>
std::vector<geometry::Vector3D> extract_ring_coordinates(const ResidueType& residue,
                                                          const geometry::Vector3D& average_origin) {
    std::vector<geometry::Vector3D> coordinates;
    coordinates.reserve(NUM_RING_ATOMS);

    // Find all ring atoms present in this residue
    std::vector<const core::Atom*> ring_atoms;
    ring_atoms.reserve(NUM_RING_ATOMS);

    for (size_t i = 0; i < NUM_RING_ATOMS; ++i) {
        const core::Atom* atom = residue.find_atom_ptr(RING_ATOM_NAMES[i]);
        if (atom != nullptr) {
            ring_atoms.push_back(atom);
        }
    }

    // Build set of ring atom names for exclusion during exocyclic search
    std::set<std::string> ring_atom_name_set;
    for (const auto* atom : ring_atoms) {
        ring_atom_name_set.insert(atom->name());
    }

    // For each ring atom, find closest exocyclic substituent
    for (const auto* ring_atom : ring_atoms) {
        const core::Atom* exocyclic_atom = nullptr;
        double min_distance = validation_constants::BOND_DISTANCE;

        for (const auto& candidate : residue.atoms()) {
            // Skip ring atoms themselves
            if (ring_atom_name_set.count(candidate.name()) > 0) {
                continue;
            }
            // Skip hydrogen atoms (first character is 'H')
            if (!candidate.name().empty() && candidate.name()[0] == 'H') {
                continue;
            }

            double distance = (candidate.position() - ring_atom->position()).length();
            if (distance < min_distance && distance > validation_constants::MIN_ATOM_DISTANCE) {
                min_distance = distance;
                exocyclic_atom = &candidate;
            }
        }

        // Use exocyclic atom if found, otherwise the ring atom itself
        const core::Atom* atom_to_use = (exocyclic_atom != nullptr) ? exocyclic_atom : ring_atom;
        coordinates.push_back(atom_to_use->position() - average_origin);
    }

    return coordinates;
}

// ============================================================================
// 3D to 2D Projection
// ============================================================================

/**
 * @brief Project 3D ring coordinates onto 2D plane perpendicular to average z-axis
 * @param ring_coords_3d Vector of 3D ring coordinates
 * @param average_z_axis The average z-axis of the two base pair reference frames
 * @return Vector of 2D points representing the projected polygon
 *
 * Constructs an orthonormal basis where the z-axis is the given average_z_axis,
 * then projects all 3D coordinates onto the xy-plane of this basis.
 */
std::vector<Point2D> project_to_base_plane(const std::vector<geometry::Vector3D>& ring_coords_3d,
                                            const geometry::Vector3D& average_z_axis) {
    std::vector<Point2D> projected_points;
    projected_points.reserve(ring_coords_3d.size());

    // Normalize z-axis
    double z_length = average_z_axis.length();
    if (z_length < 1e-10) {
        return projected_points;  // Invalid z-axis
    }
    geometry::Vector3D z_unit = average_z_axis / z_length;

    // Check if z-axis is already aligned with global z
    geometry::Vector3D global_z(0.0, 0.0, 1.0);
    double alignment = z_unit.dot(global_z);
    double alignment_angle = std::acos(std::max(-1.0, std::min(1.0, alignment)));

    if (alignment_angle < 1e-6) {
        // Already aligned with global z - use x,y coordinates directly
        for (const auto& coord : ring_coords_3d) {
            projected_points.push_back({coord.x(), coord.y()});
        }
    } else {
        // Build orthonormal basis for projection
        geometry::Vector3D x_axis;

        // Choose initial x-axis not parallel to z
        if (std::abs(z_unit.x()) < 0.9) {
            x_axis = geometry::Vector3D(1.0, 0.0, 0.0);
        } else {
            x_axis = geometry::Vector3D(0.0, 1.0, 0.0);
        }

        // Gram-Schmidt: make x_axis orthogonal to z_unit
        x_axis = x_axis - z_unit * x_axis.dot(z_unit);
        double x_length = x_axis.length();
        if (x_length > 1e-10) {
            x_axis = x_axis / x_length;
        } else {
            x_axis = geometry::Vector3D(1.0, 0.0, 0.0);
        }

        // y_axis completes the right-handed basis
        geometry::Vector3D y_axis = z_unit.cross(x_axis);
        double y_length = y_axis.length();
        if (y_length > 1e-10) {
            y_axis = y_axis / y_length;
        }

        // Project each coordinate onto the xy-plane
        for (const auto& coord : ring_coords_3d) {
            double x_projected = coord.dot(x_axis);
            double y_projected = coord.dot(y_axis);
            projected_points.push_back({x_projected, y_projected});
        }
    }

    return projected_points;
}

/**
 * @brief Core overlap calculation given 3D ring coordinates
 * @tparam ResidueType Type supporting find_atom_ptr() and atoms() methods
 * @param res1 First residue
 * @param res2 Second residue
 * @param average_origin Average of the two reference frame origins
 * @param average_z_axis Average of the two reference frame z-axes
 * @return Overlap area in square Angstroms
 */
template <typename ResidueType>
double compute_overlap_area(const ResidueType& res1, const ResidueType& res2,
                            const geometry::Vector3D& average_origin,
                            const geometry::Vector3D& average_z_axis) {
    // Step 1: Extract ring coordinates with exocyclic atoms
    std::vector<geometry::Vector3D> ring_coords_1 = extract_ring_coordinates(res1, average_origin);
    std::vector<geometry::Vector3D> ring_coords_2 = extract_ring_coordinates(res2, average_origin);

    // Need at least 3 points to form a polygon
    if (ring_coords_1.size() < 3 || ring_coords_2.size() < 3) {
        return 0.0;
    }

    // Step 2: Project onto plane perpendicular to average z-axis
    std::vector<Point2D> polygon_1 = project_to_base_plane(ring_coords_1, average_z_axis);
    std::vector<Point2D> polygon_2 = project_to_base_plane(ring_coords_2, average_z_axis);

    if (polygon_1.size() < 3 || polygon_2.size() < 3) {
        return 0.0;
    }

    // Step 3: Calculate polygon intersection area
    return OverlapCalculator::calculate_polygon_intersection(polygon_1, polygon_2);
}

}  // anonymous namespace

// ============================================================================
// Public API Implementation
// ============================================================================

double OverlapCalculator::calculate(const core::Residue& res1, const core::Residue& res2,
                                    const geometry::Vector3D& average_origin,
                                    const geometry::Vector3D& average_z_axis) {
    return compute_overlap_area(res1, res2, average_origin, average_z_axis);
}

double OverlapCalculator::calculate(const core::Residue& res1, const core::Residue& res2,
                                    const geometry::Vector3D& average_origin,
                                    const geometry::Vector3D& average_z_axis,
                                    RingDataCache& cache) {
    // Use cached ring data for faster repeated lookups
    std::vector<geometry::Vector3D> ring_coords_1 = cache.get_ring_coords(res1, average_origin);
    std::vector<geometry::Vector3D> ring_coords_2 = cache.get_ring_coords(res2, average_origin);

    if (ring_coords_1.size() < 3 || ring_coords_2.size() < 3) {
        return 0.0;
    }

    std::vector<Point2D> polygon_1 = project_to_base_plane(ring_coords_1, average_z_axis);
    std::vector<Point2D> polygon_2 = project_to_base_plane(ring_coords_2, average_z_axis);

    if (polygon_1.size() < 3 || polygon_2.size() < 3) {
        return 0.0;
    }

    return calculate_polygon_intersection(polygon_1, polygon_2);
}

std::vector<geometry::Vector3D> OverlapCalculator::get_ring_coordinates_with_exocyclic(
    const core::Residue& residue, const geometry::Vector3D& average_origin) {
    return extract_ring_coordinates(residue, average_origin);
}

double OverlapCalculator::calculate_polygon_intersection(const std::vector<Point2D>& polygon_a,
                                                         const std::vector<Point2D>& polygon_b) {
    if (polygon_a.size() < 3 || polygon_b.size() < 3) {
        return 0.0;
    }

    const long num_vertices_a = static_cast<long>(polygon_a.size());
    const long num_vertices_b = static_cast<long>(polygon_b.size());

    // Compute bounding box of both polygons
    double min_x = validation_constants::XBIG;
    double min_y = validation_constants::XBIG;
    double max_x = -validation_constants::XBIG;
    double max_y = -validation_constants::XBIG;

    for (long i = 0; i < num_vertices_a; ++i) {
        min_x = std::min(min_x, polygon_a[i].x);
        min_y = std::min(min_y, polygon_a[i].y);
        max_x = std::max(max_x, polygon_a[i].x);
        max_y = std::max(max_y, polygon_a[i].y);
    }

    for (long i = 0; i < num_vertices_b; ++i) {
        min_x = std::min(min_x, polygon_b[i].x);
        min_y = std::min(min_y, polygon_b[i].y);
        max_x = std::max(max_x, polygon_b[i].x);
        max_y = std::max(max_y, polygon_b[i].y);
    }

    // Check for degenerate bounding box
    if (max_x <= min_x || max_y <= min_y) {
        return 0.0;
    }

    // Compute scaling factors to map coordinates to integer range
    const double scale_midpoint = 0.5 * validation_constants::GAMUT;
    const double scale_x = validation_constants::GAMUT / (max_x - min_x);
    const double scale_y = validation_constants::GAMUT / (max_y - min_y);
    const double inverse_scale = scale_x * scale_y;  // For converting area back

    // Check vertex count limits
    if (num_vertices_a > MAX_POLYGON_VERTICES || num_vertices_b > MAX_POLYGON_VERTICES) {
        return 0.0;
    }

    // Allocate scaled vertex arrays (+1 for wraparound to first vertex)
    std::vector<PolygonVertex> scaled_polygon_a(num_vertices_a + 1);
    std::vector<PolygonVertex> scaled_polygon_b(num_vertices_b + 1);

    // Scale polygons to integer coordinates
    // Fudge bits (0 vs 2) provide numerical separation between the two polygons
    scale_polygon_to_integer_coords(min_x, min_y, scale_midpoint, scale_x, scale_y,
                                    polygon_a.data(), num_vertices_a,
                                    scaled_polygon_a.data(), 0);
    scale_polygon_to_integer_coords(min_x, min_y, scale_midpoint, scale_x, scale_y,
                                    polygon_b.data(), num_vertices_b,
                                    scaled_polygon_b.data(), 2);

    // Find all edge-edge intersections and accumulate area
    double intersection_area = 0.0;

    for (long i = 0; i < num_vertices_a; ++i) {
        for (long j = 0; j < num_vertices_b; ++j) {
            // Fast rejection: check if edge bounding boxes overlap
            if (!ranges_overlap(scaled_polygon_a[i].edge_x_range, scaled_polygon_b[j].edge_x_range) ||
                !ranges_overlap(scaled_polygon_a[i].edge_y_range, scaled_polygon_b[j].edge_y_range)) {
                continue;
            }

            // Compute signed areas to check for edge crossing
            double area_a1 = -signed_triangle_area(scaled_polygon_a[i].scaled_position,
                                                   scaled_polygon_b[j].scaled_position,
                                                   scaled_polygon_b[j + 1].scaled_position);
            double area_a2 = signed_triangle_area(scaled_polygon_a[i + 1].scaled_position,
                                                  scaled_polygon_b[j].scaled_position,
                                                  scaled_polygon_b[j + 1].scaled_position);

            // Check if edge A endpoints are on opposite sides of edge B
            bool a1_negative = (area_a1 < 0);
            if (a1_negative != (area_a2 < 0)) {
                continue;  // No crossing possible
            }

            double area_b1 = signed_triangle_area(scaled_polygon_b[j].scaled_position,
                                                  scaled_polygon_a[i].scaled_position,
                                                  scaled_polygon_a[i + 1].scaled_position);
            double area_b2 = -signed_triangle_area(scaled_polygon_b[j + 1].scaled_position,
                                                   scaled_polygon_a[i].scaled_position,
                                                   scaled_polygon_a[i + 1].scaled_position);

            // Check if edge B endpoints are on opposite sides of edge A
            if ((area_b1 < 0) != (area_b2 < 0)) {
                continue;  // No crossing possible
            }

            // Edges cross - handle the intersection
            if (a1_negative) {
                handle_edge_crossing(intersection_area,
                                     &scaled_polygon_a[i], &scaled_polygon_a[i + 1],
                                     &scaled_polygon_b[j], &scaled_polygon_b[j + 1],
                                     area_a1, area_a2, area_b1, area_b2);
            } else {
                handle_edge_crossing(intersection_area,
                                     &scaled_polygon_b[j], &scaled_polygon_b[j + 1],
                                     &scaled_polygon_a[i], &scaled_polygon_a[i + 1],
                                     area_b1, area_b2, area_a1, area_a2);
            }
        }
    }

    // Add contributions from vertices inside the other polygon
    add_interior_vertex_contributions(intersection_area,
                                      scaled_polygon_a.data(), num_vertices_a,
                                      scaled_polygon_b.data(), num_vertices_b);
    add_interior_vertex_contributions(intersection_area,
                                      scaled_polygon_b.data(), num_vertices_b,
                                      scaled_polygon_a.data(), num_vertices_a);

    // Convert area back from scaled coordinates
    if (std::isnan(inverse_scale) || std::isinf(inverse_scale) || inverse_scale == 0.0) {
        return 0.0;
    }

    double result = std::fabs(intersection_area) / inverse_scale;

    if (std::isnan(result) || std::isinf(result)) {
        return 0.0;
    }

    return result;
}

}  // namespace validation
}  // namespace algorithms
}  // namespace x3dna
