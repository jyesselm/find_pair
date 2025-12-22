/**
 * @file overlap_calculator.cpp
 * @brief Implementation of overlap area calculation using polygon intersection
 *
 * Matches legacy pia_inter (Polygon Intersection Area) algorithm from ana_fncs.c
 */

#include <x3dna/algorithms/validation/overlap_calculator.hpp>
#include <x3dna/algorithms/validation/ring_data_cache.hpp>
#include <x3dna/algorithms/validation_constants.hpp>
#include <x3dna/core/atom.hpp>
#include <cmath>
#include <algorithm>
#include <set>

namespace x3dna {
namespace algorithms {
namespace validation {

// Internal vertex struct for integer arithmetic polygon intersection
struct Vertex {
    Point2D ip;  // Integer point (after scaling)
    Point2D rx;  // X range for this edge
    Point2D ry;  // Y range for this edge
    long inside; // Inside count
};

// Helper functions for pia_inter (polygon intersection area)
// Matches legacy implementation exactly for precision

namespace {

double pia_area(const Point2D& a, const Point2D& p, const Point2D& q) {
    // Matches legacy pia_area (line 3218-3221)
    return p.x * q.y - p.y * q.x + a.x * (p.y - q.y) + a.y * (q.x - p.x);
}

void pia_cntrib(double* s, const Point2D& f, const Point2D& t, long w) {
    // Matches legacy pia_cntrib (line 3223-3226)
    *s += w * (t.x - f.x) * (t.y + f.y) * 0.5;
}

bool pia_ovl(const Point2D& p, const Point2D& q) {
    // Matches legacy pia_ovl (line 3228-3231)
    return p.x < q.y && q.x < p.y;
}

void pia_cross(double* out_s, Vertex* a, Vertex* b, Vertex* c, Vertex* d, double a1, double a2, double a3, double a4) {
    // Matches legacy pia_cross (line 3233-3248)
    double r1, r2;
    Point2D dp;
    r1 = a1 / (a1 + a2);
    r2 = a3 / (a3 + a4);
    dp.x = a->ip.x + r1 * (b->ip.x - a->ip.x);
    dp.y = a->ip.y + r1 * (b->ip.y - a->ip.y);
    pia_cntrib(out_s, dp, b->ip, 1);
    dp.x = c->ip.x + r2 * (d->ip.x - c->ip.x);
    dp.y = c->ip.y + r2 * (d->ip.y - c->ip.y);
    pia_cntrib(out_s, d->ip, dp, 1);
    ++a->inside;
    --c->inside;
}

void pia_inness(double* out_s, Vertex* P, long cP, Vertex* Q, long cQ) {
    // Matches legacy pia_inness (line 3250-3264)
    long j, sgn, s = 0, c = cQ;
    Point2D p = P[0].ip;
    while (c--) {
        if (Q[c].rx.x < p.x && p.x < Q[c].rx.y) {
            sgn = (0 < pia_area(p, Q[c].ip, Q[c + 1].ip)) ? 1 : 0;
            s += sgn != (Q[c].ip.x < Q[c + 1].ip.x) ? 0 : (sgn ? -1 : 1);
        }
    }
    for (j = 0; j < cP; ++j) {
        if (s) {
            pia_cntrib(out_s, P[j].ip, P[j + 1].ip, s);
        }
        s += P[j].inside;
    }
}

void pia_fit(double minx, double miny, double mid, double sclx, double scly, const Point2D* x, long cx, Vertex* ix,
             long fudge) {
    // Matches legacy pia_fit (line 3188-3216)
    // Converts floating point coordinates to integer coordinates for precision
    long c, t;
    Point2D t1, t2;
    c = cx;
    while (c--) {
        t = static_cast<long>((x[c].x - minx) * sclx - mid);
        ix[c].ip.x = static_cast<double>((t & ~7) | fudge | (c & 1));
        t = static_cast<long>((x[c].y - miny) * scly - mid);
        ix[c].ip.y = static_cast<double>((t & ~7) | fudge);
    }
    ix[0].ip.y += cx & 1;
    ix[cx] = ix[0];
    c = cx;
    while (c--) {
        t1.x = ix[c].ip.x;
        t1.y = ix[c + 1].ip.x;
        t2.x = ix[c + 1].ip.x;
        t2.y = ix[c].ip.x;
        ix[c].rx = (ix[c].ip.x < ix[c + 1].ip.x) ? t1 : t2;
        t1.x = ix[c].ip.y;
        t1.y = ix[c + 1].ip.y;
        t2.x = ix[c + 1].ip.y;
        t2.y = ix[c].ip.y;
        ix[c].ry = (ix[c].ip.y < ix[c + 1].ip.y) ? t1 : t2;
        ix[c].inside = 0;
    }
}

} // anonymous namespace

// Helper function template to get ring coordinates - works with both Residue and IResidue
namespace {

// Static ring atom names (created once)
static const char* RING_ATOM_NAMES[] = {"C4", "N3", "C2", "N1", "C6", "C5", "N7", "C8", "N9"};
static constexpr size_t NUM_RING_ATOMS = 9;

template<typename ResidueType>
std::vector<geometry::Vector3D> get_ring_coords_impl(const ResidueType& residue, const geometry::Vector3D& oave) {
    std::vector<geometry::Vector3D> ring_coords;
    ring_coords.reserve(NUM_RING_ATOMS);  // Pre-allocate

    // Find ring atoms using O(1) lookup
    std::vector<const core::Atom*> ring_atoms;
    ring_atoms.reserve(NUM_RING_ATOMS);

    for (size_t i = 0; i < NUM_RING_ATOMS; ++i) {
        const core::Atom* atom_ptr = residue.find_atom_ptr(RING_ATOM_NAMES[i]);
        if (atom_ptr != nullptr) {
            ring_atoms.push_back(atom_ptr);
        }
    }

    // Build set of ring atom names for exclusion check
    std::set<std::string> ring_atom_names;
    for (const auto* atom : ring_atoms) {
        ring_atom_names.insert(atom->name());
    }

    // For each ring atom, find ONE exocyclic atom (connected non-ring, non-hydrogen atom)
    for (const auto* ring_atom : ring_atoms) {
        const core::Atom* exocyclic_atom = nullptr;
        double min_dist = validation_constants::BOND_DISTANCE;

        for (const auto& atom : residue.atoms()) {
            // Skip ring atoms
            if (ring_atom_names.count(atom.name()) > 0) {
                continue;
            }
            // Skip hydrogen atoms (matches legacy get_cntatom which skips idx==3)
            // With trimmed atom names, check first character for 'H'
            if (!atom.name().empty() && atom.name()[0] == 'H') {
                continue;
            }
            double dist = (atom.position() - ring_atom->position()).length();
            if (dist < min_dist && dist > validation_constants::MIN_ATOM_DISTANCE) {
                min_dist = dist;
                exocyclic_atom = &atom;
            }
        }

        // Use exocyclic atom if found, otherwise use ring atom itself (matches legacy)
        const core::Atom* atom_to_use = (exocyclic_atom != nullptr) ? exocyclic_atom : ring_atom;
        ring_coords.push_back(atom_to_use->position() - oave);
    }

    return ring_coords;
}

template<typename ResidueType>
double calculate_impl(const ResidueType& res1, const ResidueType& res2,
                      const geometry::Vector3D& oave, const geometry::Vector3D& zave) {
    // Match legacy get_oarea() logic (org/src/ana_fncs.c lines 3327-3358)
    // Steps:
    // 1. Get ring atoms for both residues (with exocyclic atoms)
    // 2. Align to z-axis (project to plane perpendicular to zave)
    // 3. Calculate polygon intersection area

    // Step 1: Get ring coordinates with exocyclic atoms for both residues
    std::vector<geometry::Vector3D> ring_coords1 = get_ring_coords_impl(res1, oave);
    std::vector<geometry::Vector3D> ring_coords2 = get_ring_coords_impl(res2, oave);

    // Need at least 3 points to form a polygon
    if (ring_coords1.size() < 3 || ring_coords2.size() < 3) {
        return 0.0;
    }

    // Step 2: Align to z-axis (project to plane perpendicular to zave)
    std::vector<Point2D> poly1, poly2;

    // Normalize z-axis
    geometry::Vector3D z_normalized = zave;
    double z_len = z_normalized.length();
    if (z_len < 1e-10) {
        return 0.0; // Invalid z-axis
    }
    z_normalized = z_normalized / z_len;

    // Build orthonormal basis: zave is z-axis, find x and y axes
    geometry::Vector3D z_target(0.0, 0.0, 1.0);
    geometry::Vector3D rot_axis = z_normalized.cross(z_target);
    double rot_angle = std::acos(std::max(-1.0, std::min(1.0, z_normalized.dot(z_target))));

    // If zave is already aligned with z-axis, no rotation needed
    if (rot_angle < 1e-6 || rot_axis.length() < 1e-6) {
        // Already aligned - just use x,y coordinates
        for (const auto& coord : ring_coords1) {
            poly1.push_back({coord.x(), coord.y()});
        }
        for (const auto& coord : ring_coords2) {
            poly2.push_back({coord.x(), coord.y()});
        }
    } else {
        // Build orthonormal basis for projection
        geometry::Vector3D x_axis, y_axis;

        // Choose a vector not parallel to zave for x-axis
        if (std::abs(z_normalized.x()) < 0.9) {
            x_axis = geometry::Vector3D(1.0, 0.0, 0.0);
        } else {
            x_axis = geometry::Vector3D(0.0, 1.0, 0.0);
        }

        // Make x_axis orthogonal to zave
        x_axis = x_axis - z_normalized * x_axis.dot(z_normalized);
        double x_len = x_axis.length();
        if (x_len > 1e-10) {
            x_axis = x_axis / x_len;
        } else {
            x_axis = geometry::Vector3D(1.0, 0.0, 0.0);
        }

        // y_axis = zave x x_axis
        y_axis = z_normalized.cross(x_axis);
        double y_len = y_axis.length();
        if (y_len > 1e-10) {
            y_axis = y_axis / y_len;
        }

        // Project coordinates onto x,y plane
        for (const auto& coord : ring_coords1) {
            double x_proj = coord.dot(x_axis);
            double y_proj = coord.dot(y_axis);
            poly1.push_back({x_proj, y_proj});
        }

        for (const auto& coord : ring_coords2) {
            double x_proj = coord.dot(x_axis);
            double y_proj = coord.dot(y_axis);
            poly2.push_back({x_proj, y_proj});
        }
    }

    // Step 3: Calculate polygon intersection area
    return OverlapCalculator::calculate_polygon_intersection(poly1, poly2);
}
} // anonymous namespace

double OverlapCalculator::calculate(const core::Residue& res1, const core::Residue& res2,
                                    const geometry::Vector3D& oave, const geometry::Vector3D& zave) {
    return calculate_impl(res1, res2, oave, zave);
}

double OverlapCalculator::calculate(const core::Residue& res1, const core::Residue& res2,
                                    const geometry::Vector3D& oave, const geometry::Vector3D& zave,
                                    RingDataCache& cache) {
    // Cache-aware version: reuse pre-computed ring atom indices and exocyclic mappings
    // Step 1: Get ring coordinates using cache
    std::vector<geometry::Vector3D> ring_coords1 = cache.get_ring_coords(res1, oave);
    std::vector<geometry::Vector3D> ring_coords2 = cache.get_ring_coords(res2, oave);

    // Need at least 3 points to form a polygon
    if (ring_coords1.size() < 3 || ring_coords2.size() < 3) {
        return 0.0;
    }

    // Step 2: Align to z-axis (project to plane perpendicular to zave)
    std::vector<Point2D> poly1, poly2;

    // Normalize z-axis
    geometry::Vector3D z_normalized = zave;
    double z_len = z_normalized.length();
    if (z_len < 1e-10) {
        return 0.0; // Invalid z-axis
    }
    z_normalized = z_normalized / z_len;

    // Build orthonormal basis: zave is z-axis, find x and y axes
    geometry::Vector3D z_target(0.0, 0.0, 1.0);
    geometry::Vector3D rot_axis = z_normalized.cross(z_target);
    double rot_angle = std::acos(std::max(-1.0, std::min(1.0, z_normalized.dot(z_target))));

    // If zave is already aligned with z-axis, no rotation needed
    if (rot_angle < 1e-6 || rot_axis.length() < 1e-6) {
        // Already aligned - just use x,y coordinates
        for (const auto& coord : ring_coords1) {
            poly1.push_back({coord.x(), coord.y()});
        }
        for (const auto& coord : ring_coords2) {
            poly2.push_back({coord.x(), coord.y()});
        }
    } else {
        // Build orthonormal basis for projection
        geometry::Vector3D x_axis, y_axis;

        // Choose a vector not parallel to zave for x-axis
        if (std::abs(z_normalized.x()) < 0.9) {
            x_axis = geometry::Vector3D(1.0, 0.0, 0.0);
        } else {
            x_axis = geometry::Vector3D(0.0, 1.0, 0.0);
        }

        // Make x_axis orthogonal to zave
        x_axis = x_axis - z_normalized * x_axis.dot(z_normalized);
        double x_len = x_axis.length();
        if (x_len > 1e-10) {
            x_axis = x_axis / x_len;
        } else {
            x_axis = geometry::Vector3D(1.0, 0.0, 0.0);
        }

        // y_axis = zave x x_axis
        y_axis = z_normalized.cross(x_axis);
        double y_len = y_axis.length();
        if (y_len > 1e-10) {
            y_axis = y_axis / y_len;
        }

        // Project coordinates onto x,y plane
        for (const auto& coord : ring_coords1) {
            double x_proj = coord.dot(x_axis);
            double y_proj = coord.dot(y_axis);
            poly1.push_back({x_proj, y_proj});
        }

        for (const auto& coord : ring_coords2) {
            double x_proj = coord.dot(x_axis);
            double y_proj = coord.dot(y_axis);
            poly2.push_back({x_proj, y_proj});
        }
    }

    // Step 3: Calculate polygon intersection area
    return calculate_polygon_intersection(poly1, poly2);
}

std::vector<geometry::Vector3D> OverlapCalculator::get_ring_coordinates_with_exocyclic(
    const core::Residue& residue, const geometry::Vector3D& oave) {
    return get_ring_coords_impl(residue, oave);
}

double OverlapCalculator::calculate_polygon_intersection(const std::vector<Point2D>& poly1,
                                                         const std::vector<Point2D>& poly2) {

    // Matches legacy pia_inter (line 3266-3325)
    if (poly1.size() < 3 || poly2.size() < 3) {
        return 0.0;
    }

    long na = static_cast<long>(poly1.size());
    long nb = static_cast<long>(poly2.size());

    // Find bounding box
    double minx = validation_constants::XBIG, miny = validation_constants::XBIG;
    double maxx = -validation_constants::XBIG, maxy = -validation_constants::XBIG;

    for (long j = 0; j < na; j++) {
        if (minx > poly1[j].x)
            minx = poly1[j].x;
        if (miny > poly1[j].y)
            miny = poly1[j].y;
        if (maxx < poly1[j].x)
            maxx = poly1[j].x;
        if (maxy < poly1[j].y)
            maxy = poly1[j].y;
    }

    for (long j = 0; j < nb; j++) {
        if (minx > poly2[j].x)
            minx = poly2[j].x;
        if (miny > poly2[j].y)
            miny = poly2[j].y;
        if (maxx < poly2[j].x)
            maxx = poly2[j].x;
        if (maxy < poly2[j].y)
            maxy = poly2[j].y;
    }

    // Check for degenerate case (zero area bounding box)
    if (maxx <= minx || maxy <= miny) {
        return 0.0;
    }

    // Scale to integer coordinates
    const double mid = 0.5 * validation_constants::GAMUT;
    double sclx = validation_constants::GAMUT / (maxx - minx);
    double scly = validation_constants::GAMUT / (maxy - miny);
    double ascale = sclx * scly;

    // Allocate vertex arrays (with +1 for wraparound)
    if (na > MAX_POLYGON_VERTICES || nb > MAX_POLYGON_VERTICES) {
        // Too many vertices - fallback to zero
        return 0.0;
    }

    std::vector<Vertex> ipa(na + 1);
    std::vector<Vertex> ipb(nb + 1);

    // Convert to integer coordinates
    pia_fit(minx, miny, mid, sclx, scly, poly1.data(), na, ipa.data(), 0);
    pia_fit(minx, miny, mid, sclx, scly, poly2.data(), nb, ipb.data(), 2);

    // Calculate intersection area
    double out_s = 0.0;
    double a1, a2, a3, a4;
    bool o;

    for (long j = 0; j < na; ++j) {
        for (long k = 0; k < nb; ++k) {
            if (pia_ovl(ipa[j].rx, ipb[k].rx) && pia_ovl(ipa[j].ry, ipb[k].ry)) {
                a1 = -pia_area(ipa[j].ip, ipb[k].ip, ipb[k + 1].ip);
                a2 = pia_area(ipa[j + 1].ip, ipb[k].ip, ipb[k + 1].ip);
                o = (a1 < 0);
                if (o == (a2 < 0)) {
                    a3 = pia_area(ipb[k].ip, ipa[j].ip, ipa[j + 1].ip);
                    a4 = -pia_area(ipb[k + 1].ip, ipa[j].ip, ipa[j + 1].ip);
                    if ((a3 < 0) == (a4 < 0)) {
                        if (o) {
                            pia_cross(&out_s, &ipa[j], &ipa[j + 1], &ipb[k], &ipb[k + 1], a1, a2, a3, a4);
                        } else {
                            pia_cross(&out_s, &ipb[k], &ipb[k + 1], &ipa[j], &ipa[j + 1], a3, a4, a1, a2);
                        }
                    }
                }
            }
        }
    }

    pia_inness(&out_s, ipa.data(), na, ipb.data(), nb);
    pia_inness(&out_s, ipb.data(), nb, ipa.data(), na);

    // Check for invalid ascale (division by zero protection)
    if (std::isnan(ascale) || std::isinf(ascale) || ascale == 0.0) {
        return 0.0;
    }

    double result = std::fabs(out_s) / ascale;

    // Ensure result is valid (not NaN or Inf)
    if (std::isnan(result) || std::isinf(result)) {
        return 0.0;
    }

    return result;
}

} // namespace validation
} // namespace algorithms
} // namespace x3dna
