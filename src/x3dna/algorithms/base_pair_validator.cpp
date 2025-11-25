/**
 * @file base_pair_validator.cpp
 * @brief Implementation of base pair validation (matches legacy check_pair)
 */

#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cctype>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>

// Helper struct for 2D points (matches legacy point struct)
struct Point2D {
    double x;
    double y;
};

// Helper struct for vertex (matches legacy vertex struct)
// Used in pia_inter for integer arithmetic polygon intersection
struct Vertex {
    Point2D ip;      // Integer point (after scaling)
    Point2D rx;      // X range for this edge
    Point2D ry;      // Y range for this edge
    long inside;     // Inside count
};

// Constants matching legacy
constexpr double XBIG = 1.0e+18;
constexpr long MNPOLY = 1000;
constexpr double GAMUT = 5.0e8;

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;
using namespace x3dna::geometry;

// Static members for atom list
std::map<std::string, std::string> BasePairValidator::atom_list_;
bool BasePairValidator::atom_list_loaded_ = false;

// Debug flag - set via environment variable DEBUG_GOOD_HB_ATOMS=1
static bool should_debug_good_hb_atoms() {
    static bool checked = false;
    static bool enabled = false;
    if (!checked) {
        const char* env = std::getenv("DEBUG_GOOD_HB_ATOMS");
        enabled = (env && std::string(env) == "1");
        checked = true;
    }
    return enabled;
}

// Debug flag for specific pairs - can be added later if needed
// Set via DEBUG_PAIRS="1VBY:20,21;3AVY:1204,1223"

ValidationResult BasePairValidator::validate(const Residue& res1, const Residue& res2) const {
    ValidationResult result;

    // Skip if same residue
    if (&res1 == &res2) {
        return result;
    }

    // Both residues must have reference frames
    auto frame1_opt = res1.reference_frame();
    auto frame2_opt = res2.reference_frame();
    if (!frame1_opt.has_value() || !frame2_opt.has_value()) {
        return result;
    }

    const ReferenceFrame& frame1 = frame1_opt.value();
    const ReferenceFrame& frame2 = frame2_opt.value();

    // Calculate average z-axis and origin
    Vector3D oave, zave;
    get_bp_zoave(frame1, frame2, oave, zave);

    // Calculate dorg (distance between origins)
    // Legacy: ddxyz(org[i], org[j], dorg) means dorg = org[i] - org[j]
    // Note: Legacy uses 1-based indexing, so org[i] is first residue, org[j] is second
    // Modern: frame1 is first (i), frame2 is second (j)
    Vector3D dorg = frame1.origin() - frame2.origin();
    result.dorg = dorg.length();

    // Calculate direction vectors (dot products of frame axes)
    calculate_direction_vectors(frame1, frame2, result.dir_x, result.dir_y, result.dir_z);

    // Calculate d_v (vertical distance)
    result.d_v = std::abs(dorg.dot(zave));

    // Calculate plane angle (angle between z-axes, 0-90 degrees)
    result.plane_angle = z1_z2_angle_in_0_to_90(frame1.z_axis(), frame2.z_axis());

    // Calculate dNN (distance between N1/N9 atoms)
    auto n1_n9_1 = find_n1_n9_position(res1);
    auto n1_n9_2 = find_n1_n9_position(res2);
    if (n1_n9_1.has_value() && n1_n9_2.has_value()) {
        Vector3D dNN_vec = n1_n9_1.value() - n1_n9_2.value();
        result.dNN = dNN_vec.length();
    } else {
        result.dNN = 1e10; // Large value if N1/N9 not found
    }

    // Calculate quality score (matches rtn_val[5])
    result.quality_score = result.dorg + 2.0 * result.d_v + result.plane_angle / 20.0;

    // Perform validation checks
    result.distance_check = in_range(result.dorg, params_.min_dorg, params_.max_dorg);
    result.d_v_check = in_range(result.d_v, params_.min_dv, params_.max_dv);
    result.plane_angle_check =
        in_range(result.plane_angle, params_.min_plane_angle, params_.max_plane_angle);
    result.dNN_check = in_range(result.dNN, params_.min_dNN, params_.max_dNN);

    // Check overlap area
    result.overlap_area = calculate_overlap_area(res1, res2, oave, zave);
    result.overlap_check = (result.overlap_area < params_.overlap_threshold);

    // If all distance/angle checks pass and overlap is acceptable, check hydrogen bonds
    bool cdns = result.distance_check && result.d_v_check && result.plane_angle_check &&
                result.dNN_check && result.overlap_check;

    if (cdns) {
        // Count H-bonds simply (BEFORE validation) - matches legacy check_pair behavior
        // This is the key fix: legacy counts H-bonds before validation for pair validation
        count_hydrogen_bonds_simple(res1, res2, result.num_base_hb, result.num_o2_hb);

        // Check H-bond requirement (matches legacy lines 4616-4617)
        if (params_.min_base_hb > 0) {
            result.hbond_check = (result.num_base_hb >= params_.min_base_hb);
        } else {
            result.hbond_check = (result.num_o2_hb > 0 || result.num_base_hb > 0);
        }

        // Pair is valid if all checks pass
        result.is_valid = result.hbond_check;

        // Find validated H-bonds (AFTER validation) - used for adjust_pairQuality
        // This matches legacy hb_numlist behavior which uses validated H-bonds
        if (result.is_valid) {
            result.hbonds = find_hydrogen_bonds(res1, res2);
        }

        // Determine base pair type (simplified - would need calculate_more_bppars)
        if (result.is_valid) {
            // For now, set to UNKNOWN - will be determined by calculate_more_bppars equivalent
            result.bp_type = BasePairType::UNKNOWN;
        }
    }

    return result;
}

void BasePairValidator::calculate_direction_vectors(const ReferenceFrame& frame1,
                                                    const ReferenceFrame& frame2, double& dir_x,
                                                    double& dir_y, double& dir_z) {

    // Dot products of corresponding axes
    dir_x = frame1.x_axis().dot(frame2.x_axis());
    dir_y = frame1.y_axis().dot(frame2.y_axis());
    dir_z = frame1.z_axis().dot(frame2.z_axis());
}

void BasePairValidator::get_bp_zoave(const ReferenceFrame& frame1, const ReferenceFrame& frame2,
                                     Vector3D& oave, Vector3D& zave) {

    // Average origin (matches legacy avexyz)
    oave = (frame1.origin() + frame2.origin()) * 0.5;

    // Average z-axis (matches legacy get_bp_zoave exactly)
    // Legacy: d = dot(z1, z2)
    //         if d > 0: zave = z2 + z1 (sumxyz)
    //         else: zave = z2 - z1 (ddxyz)
    //         then normalize
    Vector3D z1 = frame1.z_axis();
    Vector3D z2 = frame2.z_axis();
    double d = z1.dot(z2);

    if (d > 0.0) {
        // z-axes point in same direction: add them
        zave = z2 + z1;
    } else {
        // z-axes point in opposite directions: subtract them
        zave = z2 - z1;
    }

    // Normalize (matches legacy vec_norm)
    double len = zave.length();
    if (len > 1e-10) {
        zave = zave / len;
    } else {
        // Fallback if vectors are opposite and cancel out
        zave = z1;
    }
}

double BasePairValidator::z1_z2_angle_in_0_to_90(const Vector3D& z1, const Vector3D& z2) {

    double dot_product = z1.dot(z2);
    // Clamp to [-1, 1] to avoid numerical issues
    dot_product = std::max(-1.0, std::min(1.0, dot_product));

    double angle_rad = std::acos(dot_product);
    double angle_deg = angle_rad * 180.0 / M_PI;

    // Return angle in 0-90 degree range
    if (angle_deg > 90.0) {
        return 180.0 - angle_deg;
    }
    return angle_deg;
}

std::optional<Vector3D> BasePairValidator::find_n1_n9_position(const Residue& residue) {
    // Match legacy glyco_N logic exactly
    // RY[i] == 1 (purine): find N9
    // RY[i] == 0 (pyrimidine): find N1 (or C5 for P/p bases, but for dNN we use N1)

    ResidueType res_type = residue.residue_type();
    bool is_purine = (res_type == ResidueType::ADENINE || res_type == ResidueType::GUANINE);

    if (is_purine) {
        // Purine: find N9
        auto n9 = residue.find_atom(" N9 ");
        if (n9.has_value()) {
            return n9->position();
        }
    } else {
        // Pyrimidine: find N1
        // Special case: for P/p bases, legacy uses C5, but for dNN calculation
        // we should use N1 if available
        char base = residue.one_letter_code();
        if (base == 'P' || base == 'p') {
            // For P/p, try C5 first (matches legacy glyco_N for P/p)
            auto c5 = residue.find_atom(" C5 ");
            if (c5.has_value()) {
                return c5->position();
            }
        }
        // Otherwise, use N1
        auto n1 = residue.find_atom(" N1 ");
        if (n1.has_value()) {
            return n1->position();
        }
    }

    return std::nullopt;
}

// Forward declaration for helper function
static double calculate_polygon_intersection_area(const std::vector<Point2D>& poly1,
                                                  const std::vector<Point2D>& poly2);

double BasePairValidator::calculate_overlap_area(const Residue& res1, const Residue& res2,
                                                 const Vector3D& oave,
                                                 const Vector3D& zave) const {
    // Match legacy get_oarea() logic (org/src/ana_fncs.c lines 3327-3358)
    // Steps:
    // 1. Get ring atoms for both residues
    // 2. Translate coordinates relative to oave
    // 3. Align to z-axis (project to plane perpendicular to zave)
    // 4. Calculate polygon intersection area

    // Step 1: Get ring atoms (only_ring = 0 means include exocyclic atoms)
    // Legacy: ratom_xyz(ring_atom[r1], only_ring, xyz, oave, oxyz1)
    // We need to identify ring atoms from residues
    std::vector<Vector3D> ring_coords1, ring_coords2;
    
    // Get ring atoms using is_base_atom check (matches legacy is_baseatom)
    for (const auto& atom : res1.atoms()) {
        if (is_base_atom(atom.name())) {
            // Translate relative to oave (matches legacy ddxyz(oave, xyz[k], oxyz[i]))
            ring_coords1.push_back(atom.position() - oave);
        }
    }
    
    for (const auto& atom : res2.atoms()) {
        if (is_base_atom(atom.name())) {
            ring_coords2.push_back(atom.position() - oave);
        }
    }
    
    // Need at least 3 points to form a polygon
    if (ring_coords1.size() < 3 || ring_coords2.size() < 3) {
        return 0.0;
    }
    
    // Step 2: Align to z-axis (project to plane perpendicular to zave)
    // Legacy: align2zaxis(n1, zave, rotmat, oxyz1, oxyz1Z)
    // This rotates coordinates so zave becomes z-axis, then we use x,y as 2D coordinates
    std::vector<Point2D> poly1, poly2;
    
    // Build rotation matrix to align zave with z-axis
    // zave should become (0, 0, 1) after rotation
    Vector3D z_target(0.0, 0.0, 1.0);
    Vector3D z_normalized = zave;
    double z_len = z_normalized.length();
    if (z_len < 1e-10) {
        return 0.0; // Invalid z-axis
    }
    z_normalized = z_normalized / z_len;
    
    // Find rotation axis (cross product of zave and z_target)
    Vector3D rot_axis = z_normalized.cross(z_target);
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
        // Apply rotation using Rodrigues' rotation formula
        // For now, use a simpler approach: project onto plane perpendicular to zave
        // The plane is defined by x and y axes after rotating zave to z-axis
        
        // Build orthonormal basis: zave is z-axis, find x and y axes
        Vector3D x_axis, y_axis;
        
        // Choose a vector not parallel to zave for x-axis
        if (std::abs(z_normalized.x()) < 0.9) {
            x_axis = Vector3D(1.0, 0.0, 0.0);
        } else {
            x_axis = Vector3D(0.0, 1.0, 0.0);
        }
        
        // Make x_axis orthogonal to zave
        x_axis = x_axis - z_normalized * x_axis.dot(z_normalized);
        double x_len = x_axis.length();
        if (x_len > 1e-10) {
            x_axis = x_axis / x_len;
        } else {
            x_axis = Vector3D(1.0, 0.0, 0.0);
        }
        
        // y_axis = zave × x_axis
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
    // Legacy: pia_inter(a, n1, b, n2)
    return calculate_polygon_intersection_area(poly1, poly2);
}

// Helper functions for pia_inter (polygon intersection area)
// Matches legacy implementation exactly for precision

static double pia_area(const Point2D& a, const Point2D& p, const Point2D& q) {
    // Matches legacy pia_area (line 3218-3221)
    return p.x * q.y - p.y * q.x + a.x * (p.y - q.y) + a.y * (q.x - p.x);
}

static void pia_cntrib(double* s, const Point2D& f, const Point2D& t, long w) {
    // Matches legacy pia_cntrib (line 3223-3226)
    *s += w * (t.x - f.x) * (t.y + f.y) * 0.5;
}

static bool pia_ovl(const Point2D& p, const Point2D& q) {
    // Matches legacy pia_ovl (line 3228-3231)
    return p.x < q.y && q.x < p.y;
}

static void pia_cross(double* out_s, Vertex* a, Vertex* b, Vertex* c, Vertex* d,
                      double a1, double a2, double a3, double a4) {
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

static void pia_inness(double* out_s, Vertex* P, long cP, Vertex* Q, long cQ) {
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

static void pia_fit(double minx, double miny, const double mid, double sclx, double scly,
                    const Point2D* x, long cx, Vertex* ix, long fudge) {
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

// Main polygon intersection function (matches legacy pia_inter exactly)
static double calculate_polygon_intersection_area(const std::vector<Point2D>& poly1,
                                                  const std::vector<Point2D>& poly2) {
    // Matches legacy pia_inter (line 3266-3325)
    if (poly1.size() < 3 || poly2.size() < 3) {
        return 0.0;
    }
    
    long na = static_cast<long>(poly1.size());
    long nb = static_cast<long>(poly2.size());
    
    // Find bounding box
    double minx = XBIG, miny = XBIG;
    double maxx = -XBIG, maxy = -XBIG;
    
    for (long j = 0; j < na; j++) {
        if (minx > poly1[j].x) minx = poly1[j].x;
        if (miny > poly1[j].y) miny = poly1[j].y;
        if (maxx < poly1[j].x) maxx = poly1[j].x;
        if (maxy < poly1[j].y) maxy = poly1[j].y;
    }
    
    for (long j = 0; j < nb; j++) {
        if (minx > poly2[j].x) minx = poly2[j].x;
        if (miny > poly2[j].y) miny = poly2[j].y;
        if (maxx < poly2[j].x) maxx = poly2[j].x;
        if (maxy < poly2[j].y) maxy = poly2[j].y;
    }
    
    // Check for degenerate case (zero area bounding box)
    if (maxx <= minx || maxy <= miny) {
        return 0.0;
    }
    
    // Scale to integer coordinates
    const double mid = 0.5 * GAMUT;
    double sclx = GAMUT / (maxx - minx);
    double scly = GAMUT / (maxy - miny);
    double ascale = sclx * scly;
    
    // Allocate vertex arrays (with +1 for wraparound)
    if (na > MNPOLY || nb > MNPOLY) {
        // Too many vertices - fallback to bounding box approximation
        // This shouldn't happen for base pairs (typically < 10 vertices)
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

void BasePairValidator::count_hydrogen_bonds_simple(const Residue& res1, const Residue& res2,
                                                    int& num_base_hb, int& num_o2_hb) const {
    // Matches legacy check_pair H-bond counting (lines 4605-4614 in cmn_fncs.c)
    // Counts H-bonds BEFORE validation - this is the key difference

    num_base_hb = 0;
    num_o2_hb = 0;
    bool debug_enabled = should_debug_good_hb_atoms();

    // Loop through all atom pairs (matches legacy nested loops)
    for (const auto& atom1 : res1.atoms()) {
        for (const auto& atom2 : res2.atoms()) {
            // Check distance using within_limits equivalent
            Vector3D diff = atom1.position() - atom2.position();
            double dist = diff.length();

            // Check if distance is in range [hb_lower, hb_dist1]
            if (dist < params_.hb_lower || dist > params_.hb_dist1) {
                continue;
            }

            // Check if both are base atoms and can form H-bond
            // Legacy: O2' is counted separately as num_o2_hb, not num_base_hb
            // So exclude O2' from base H-bond counting
            bool atom1_is_base = is_base_atom(atom1.name());
            bool atom2_is_base = is_base_atom(atom2.name());
            bool both_base = atom1_is_base && atom2_is_base;
            bool not_o2prime = (atom1.name() != " O2'" && atom2.name() != " O2'");
            bool good_hb = false;

            if (both_base && not_o2prime) {
                good_hb = good_hb_atoms(atom1.name(), atom2.name());
                if (good_hb) {
                    num_base_hb++;
                }
            }

            // Debug output for H-bond counting
            if (debug_enabled && both_base && dist <= params_.hb_dist1 &&
                dist >= params_.hb_lower) {
                std::cerr << "[DEBUG count_hb] " << atom1.name() << " - " << atom2.name()
                          << " (dist=" << std::fixed << std::setprecision(3) << dist << "): "
                          << "atom1_is_base=" << atom1_is_base
                          << ", atom2_is_base=" << atom2_is_base << ", not_o2prime=" << not_o2prime
                          << ", good_hb_atoms=" << good_hb << ", num_base_hb=" << num_base_hb
                          << std::endl;
            }

            // Check if either atom is O2'
            if (atom1.name() == " O2'" || atom2.name() == " O2'") {
                num_o2_hb++;
            }
        }
    }

    if (debug_enabled) {
        std::cerr << "[DEBUG count_hb] Final counts: num_base_hb=" << num_base_hb
                  << ", num_o2_hb=" << num_o2_hb << std::endl;
    }
}

std::vector<core::hydrogen_bond> BasePairValidator::find_hydrogen_bonds(const Residue& res1,
                                                                        const Residue& res2) const {
    // Use HydrogenBondFinder to get ALL H-bonds including invalid ones (type=' ')
    // This matches legacy behavior where ALL H-bonds are recorded to JSON
    // Legacy records ALL H-bonds in get_hbond_ij (including type=' ')
    using namespace x3dna::algorithms;
    
    double hb_lower = params_.hb_lower;
    double hb_dist1 = params_.hb_dist1;
    double hb_dist2 = 4.5; // Default value (matches legacy)
    
    DetailedHBondResult detailed = HydrogenBondFinder::find_hydrogen_bonds_detailed(
        res1, res2, hb_lower, hb_dist1, hb_dist2);
    
    // Return after_validation which includes ALL H-bonds (including type=' ')
    // This matches legacy behavior where all H-bonds are recorded to JSON
    // Note: Legacy records ALL H-bonds to JSON, not just validated ones
    // Legacy uses fabs() for distance in json_writer_record_hbond_list (line 1165)
    std::vector<core::hydrogen_bond> all_hbonds;
    for (const auto& hbond_result : detailed.after_validation) {
        core::hydrogen_bond hbond;
        hbond.donor_atom = hbond_result.donor_atom;
        hbond.acceptor_atom = hbond_result.acceptor_atom;
        // Use absolute distance (matches legacy fabs(hb_dist[i]) in json_writer.c:1165)
        hbond.distance = std::abs(hbond_result.distance);
        hbond.type = hbond_result.type;
        all_hbonds.push_back(hbond);
    }
    
    return all_hbonds;
}

bool BasePairValidator::is_base_atom(const std::string& atom_name) {
    // Matches legacy is_baseatom (line 4652 in cmn_fncs.c)
    // Base atoms: C5M or atoms matching pattern " HP" where H is not H or P
    if (atom_name == " C5M") {
        return true;
    }

    // Pattern: space, character (not H or P), digit, space
    // Legacy: atomname[0] == ' ' && strchr("HP", atomname[1]) == NULL
    //         && isdigit(atomname[2]) && atomname[3] == ' '
    if (atom_name.length() >= 4 && atom_name[0] == ' ' && atom_name[1] != 'H' &&
        atom_name[1] != 'P' && std::isdigit(static_cast<unsigned char>(atom_name[2])) &&
        atom_name[3] == ' ') {
        return true;
    }

    return false;
}

bool BasePairValidator::good_hb_atoms(const std::string& atom1, const std::string& atom2) const {
    // Match legacy good_hbatoms() logic EXACTLY (lines 3864-3877 in cmn_fncs.c)
    // Legacy signature: long good_hbatoms(miscPars *misc_pars, char *atom1, char *atom2, long idx1,
    // long idx2)

    bool debug_enabled = should_debug_good_hb_atoms();

    if (debug_enabled) {
        std::cerr << "\n[DEBUG good_hb_atoms] Checking: " << atom1 << " - " << atom2 << std::endl;
    }

    // Step 1: PO list check (matches legacy lines 3866-3870)
    // PO list: atoms that should be rejected if BOTH are in this list
    // Legacy: static char *PO[] = { " O1P", " O2P", " O3'", " O4'", " O5'", " N7 " };
    // Legacy: static long numPO = sizeof PO / sizeof PO[0] - 1;  // = 5 (6 elements - 1)
    static const std::vector<std::string> PO = {" O1P", " O2P", " O3'", " O4'", " O5'", " N7 "};

    // Check if both atoms match PO list using exact string comparison (matches num_strmatch)
    // Legacy: if (num_strmatch(atom1, PO, 0, numPO) && num_strmatch(atom2, PO, 0, numPO)) return
    // FALSE; num_strmatch does exact strcmp comparison, std::find does exact string comparison
    bool atom1_in_po = std::find(PO.begin(), PO.end(), atom1) != PO.end();
    bool atom2_in_po = std::find(PO.begin(), PO.end(), atom2) != PO.end();
    if (debug_enabled) {
        std::cerr << "  PO check: atom1_in_po=" << atom1_in_po << ", atom2_in_po=" << atom2_in_po
                  << std::endl;
    }
    if (atom1_in_po && atom2_in_po) {
        if (debug_enabled) {
            std::cerr << "  -> REJECTED (both in PO list)" << std::endl;
        }
        return false; // Matches legacy: return FALSE
    }

    // Step 2: idx-based check (matches legacy lines 3871-3874)
    // Legacy: if ((idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4) &&
    //             (lval_in_set(idx1, 1, natom, misc_pars->hb_idx) &&
    //              lval_in_set(idx2, 1, natom, misc_pars->hb_idx)))
    //             return TRUE;
    // NOTE: idx1 and idx2 are atom type indices from atom_idx() function
    //       Based on atoms_list: idx 0=UNK, 1=C, 2=O, 3=H, 4=N, 5=S, 6=P, ...
    //       hb_idx array is built from hb_atoms string (default ".O.N") via
    //       reset_alist_symbol_idx()

    // Build hb_idx array from hb_atoms string (matches reset_alist_symbol_idx)
    // For default ".O.N": hb_idx = [2, 4] (O and N)
    std::vector<int> hb_idx;
    std::string hb_atoms_upper = params_.hb_atoms;
    std::transform(hb_atoms_upper.begin(), hb_atoms_upper.end(), hb_atoms_upper.begin(), ::toupper);

    // Replace '.' with ' ' (matches legacy reset_alist_symbol_idx)
    std::replace(hb_atoms_upper.begin(), hb_atoms_upper.end(), '.', ' ');

    // Parse pairs of characters as atomic symbols
    // atoms_list mapping: " C"=1, " O"=2, " H"=3, " N"=4, " S"=5, " P"=6, ...
    // Use function-local static to avoid capture issues
    static const std::map<std::string, int> atoms_list_idx = {{" C", 1}, {" O", 2}, {" H", 3},
                                                              {" N", 4}, {" S", 5}, {" P", 6}};

    for (size_t i = 0; i + 1 < hb_atoms_upper.length(); i += 2) {
        std::string sym = hb_atoms_upper.substr(i, 2);
        auto it = atoms_list_idx.find(sym);
        if (it != atoms_list_idx.end()) {
            hb_idx.push_back(it->second);
        }
    }

    // Extract atomic symbol from atom name and map to idx (matches aname2asym + asym_idx)
    // Legacy aname2asym: replaces non-alphabetic with '.', matches against atomlist, extracts
    // symbol Load atomlist if not already loaded
    if (!atom_list_loaded_) {
        load_atom_list();
    }

    auto get_atom_idx = [debug_enabled](const std::string& atom_name) -> int {
        if (atom_name.length() < 4) {
            return 0; // UNK - need at least 4 chars for PDB format
        }

        // Match legacy aname2asym logic EXACTLY (lines 4160-4193 in cmn_fncs.c)
        // Step 1: Create pattern by replacing non-alphabetic with '.'
        std::string aname_pattern = atom_name;
        for (size_t i = 0; i < 4 && i < aname_pattern.length(); i++) {
            if (!std::isalpha(static_cast<unsigned char>(aname_pattern[i]))) {
                aname_pattern[i] = '.';
            }
        }

        if (debug_enabled) {
            std::cerr << "    Pattern for '" << atom_name << "': '" << aname_pattern << "'"
                      << std::endl;
        }

        // Step 2: Try to match against atomlist (matches legacy lines 4168-4170)
        // Legacy uses str_pmatch(atomlist[i], aname) which does:
        // strncmp(atomlist[i], aname, strlen(aname))
        // This checks if atomlist[i] (like ".O4. O") starts with aname (like ".O4.")
        // So we check if pattern starts with aname_pattern
        std::string my_asym;
        bool found_match = false;
        std::string matched_pattern;
        for (const auto& [pattern, sym] : atom_list_) {
            // str_pmatch(atomlist[i], aname) checks if atomlist[i] starts with aname
            // So check if pattern starts with aname_pattern
            if (pattern.length() >= aname_pattern.length() &&
                pattern.substr(0, aname_pattern.length()) == aname_pattern) {
                my_asym = sym;
                matched_pattern = pattern;
                found_match = true;
                if (debug_enabled) {
                    std::cerr << "    Atomlist match: pattern='" << pattern << "' -> sym='" << sym
                              << "'" << std::endl;
                }
                break;
            }
        }

        // Step 3: If no match, use fallback logic (matches legacy lines 4171-4184)
        if (!found_match) {
            bool unknown = (aname_pattern == ".UNK");

            if (debug_enabled) {
                std::cerr << "    No atomlist match, using fallback (unknown=" << unknown << ")"
                          << std::endl;
            }

            // Fallback case 1: aname[0]!='.' && aname[1]!='.' && aname[2]=='.' && aname[3]=='.'
            // Extract aname[0] + aname[1] as atomic symbol
            if (aname_pattern.length() >= 4 && aname_pattern[0] != '.' && aname_pattern[1] != '.' &&
                aname_pattern[2] == '.' && aname_pattern[3] == '.') {
                my_asym = std::string(" ") + aname_pattern[1]; // Space + second char
                if (debug_enabled) {
                    std::cerr << "    Fallback case 1: sym='" << my_asym << "'" << std::endl;
                }
            }
            // Fallback case 2: aname[0]=='.' && aname[1]!='.' && !unknown
            // Extract " " + aname[1] as atomic symbol (matches legacy line 4178-4180)
            else if (aname_pattern.length() >= 2 && aname_pattern[0] == '.' &&
                     aname_pattern[1] != '.' && !unknown) {
                my_asym = std::string(" ") + aname_pattern[1];
                if (debug_enabled) {
                    std::cerr << "    Fallback case 2: sym='" << my_asym << "'" << std::endl;
                }
            }
            // Fallback case 3: aname[0]=='H' -> " H" (matches legacy line 4181-4182)
            else if (aname_pattern.length() >= 1 && aname_pattern[0] == 'H') {
                my_asym = " H";
                if (debug_enabled) {
                    std::cerr << "    Fallback case 3: sym='" << my_asym << "'" << std::endl;
                }
            }
            // No match - use UNK
            else {
                my_asym = "XX"; // UNKATM
                if (debug_enabled) {
                    std::cerr << "    Fallback: no match, using UNK" << std::endl;
                }
            }
        }

        // Step 4: Map atomic symbol to idx (matches legacy asym_idx)
        auto it = atoms_list_idx.find(my_asym);
        int idx = 0;
        if (it != atoms_list_idx.end()) {
            idx = it->second;
        }

        if (debug_enabled) {
            std::cerr << "    Final: sym='" << my_asym << "' -> idx=" << idx << std::endl;
        }

        return idx;
    };

    if (debug_enabled) {
        std::cerr << "  hb_atoms string: '" << params_.hb_atoms << "'" << std::endl;
        std::cerr << "  hb_idx array: [";
        for (size_t i = 0; i < hb_idx.size(); i++) {
            if (i > 0)
                std::cerr << ", ";
            std::cerr << hb_idx[i];
        }
        std::cerr << "]" << std::endl;
    }

    int idx1 = get_atom_idx(atom1);
    int idx2 = get_atom_idx(atom2);

    if (debug_enabled) {
        std::cerr << "  idx1=" << idx1 << ", idx2=" << idx2 << std::endl;
    }

    // Legacy check: (idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4)
    //               AND both in hb_idx
    bool at_least_one_on = (idx1 == 2 || idx1 == 4 || idx2 == 2 || idx2 == 4);
    if (debug_enabled) {
        std::cerr << "  at_least_one_on (idx==2||idx==4): " << at_least_one_on << std::endl;
    }

    if (at_least_one_on) {
        // Check if both idx are in hb_idx (matches lval_in_set)
        bool idx1_in_hb = std::find(hb_idx.begin(), hb_idx.end(), idx1) != hb_idx.end();
        bool idx2_in_hb = std::find(hb_idx.begin(), hb_idx.end(), idx2) != hb_idx.end();
        if (debug_enabled) {
            std::cerr << "  idx1_in_hb: " << idx1_in_hb << ", idx2_in_hb: " << idx2_in_hb
                      << std::endl;
        }
        if (idx1_in_hb && idx2_in_hb) {
            if (debug_enabled) {
                std::cerr << "  -> ACCEPTED (both idx in hb_idx)" << std::endl;
            }
            return true; // Matches legacy: return TRUE
        } else {
            if (debug_enabled) {
                std::cerr << "  -> REJECTED (not both idx in hb_idx)" << std::endl;
            }
        }
    } else {
        if (debug_enabled) {
            std::cerr << "  -> REJECTED (no idx==2||idx==4)" << std::endl;
        }
    }

    // Legacy: else return FALSE;
    return false;
}

char BasePairValidator::donor_acceptor(char base1, char base2, const std::string& atom1,
                                       const std::string& atom2) {
    // Matches legacy donor_acceptor function
    // CB_LIST = "ACGITU" (A=0, C=1, G=2, I=3, T=4, U=5)
    static const char* CB_LIST = "ACGITU";
    static const char* da_types[] = {"AD", "AX", "XD", "XX", "DA", "DX", "XA"};
    static const char* bb_da[] = {" O1P", " O2P", " O5'", " O4'", " O3'", " O2'"};
    static const char bb_da_type[] = {'A', 'A', 'A', 'A', 'A', 'X'};

    // Base-specific atom patterns: [base_index][atom_index] = "atom_name_type"
    // base_da[base_index][atom_index] format: " N9 _?" means N9 with type '?' (either)
    static const char* base_da[6][6] = {// A (Adenine, index 0)
                                        {" N9 ", " N7 ", " N6 ", " N1 ", " N3 ", nullptr},
                                        // C (Cytosine, index 1)
                                        {" N1 ", " O2 ", " N3 ", " N4 ", nullptr, nullptr},
                                        // G (Guanine, index 2)
                                        {" N9 ", " N7 ", " O6 ", " N1 ", " N2 ", " N3 "},
                                        // I (Inosine, index 3)
                                        {" N9 ", " N7 ", " O6 ", " N1 ", " N3 ", nullptr},
                                        // T (Thymine, index 4)
                                        {" N1 ", " O2 ", " N3 ", " O4 ", nullptr, nullptr},
                                        // U (Uracil, index 5)
                                        {" N1 ", " O2 ", " N3 ", " O4 ", nullptr, nullptr}};

    // Base-specific atom types: [base_index][atom_index] = type
    static const char base_da_type[6][6] = {// A
                                            {'?', 'A', 'D', 'A', 'A', '\0'},
                                            // C
                                            {'?', 'A', 'A', 'D', '\0', '\0'},
                                            // G
                                            {'?', 'A', 'A', 'D', 'D', 'A'},
                                            // I
                                            {'?', 'A', 'A', 'D', 'A', '\0'},
                                            // T
                                            {'?', 'A', 'D', 'A', '\0', '\0'},
                                            // U
                                            {'?', 'A', 'D', 'A', '\0', '\0'}};

    char hbatom_type = '*';
    char ia = '\0', ja = '\0';

    // Find base indices in CB_LIST
    int inum = -1, jnum = -1;
    const char* pchar = strchr(CB_LIST, std::toupper(static_cast<unsigned char>(base1)));
    if (pchar)
        inum = pchar - CB_LIST;
    pchar = strchr(CB_LIST, std::toupper(static_cast<unsigned char>(base2)));
    if (pchar)
        jnum = pchar - CB_LIST;

    if (inum < 0 || jnum < 0) {
        return hbatom_type; // Invalid base
    }

    // Check backbone atoms first
    for (int i = 0; i < 6; ++i) {
        if (atom1.length() >= 4 && atom1.substr(0, 4) == bb_da[i]) {
            ia = bb_da_type[i];
        }
        if (atom2.length() >= 4 && atom2.substr(0, 4) == bb_da[i]) {
            ja = bb_da_type[i];
        }
    }

    // Check base-specific atoms if not found in backbone
    if (!ia && inum >= 0 && inum < 6) {
        for (int i = 0; i < 6; ++i) {
            if (base_da[inum][i] == nullptr)
                break;
            if (atom1.length() >= 4 && atom1.substr(0, 4) == base_da[inum][i]) {
                ia = base_da_type[inum][i];
                break;
            }
        }
    }

    if (!ja && jnum >= 0 && jnum < 6) {
        for (int i = 0; i < 6; ++i) {
            if (base_da[jnum][i] == nullptr)
                break;
            if (atom2.length() >= 4 && atom2.substr(0, 4) == base_da[jnum][i]) {
                ja = base_da_type[jnum][i];
                break;
            }
        }
    }

    // Check if donor-acceptor combination is valid
    if (ia && ja) {
        char da[3] = {ia, ja, '\0'};
        // Check against valid da_types: "AD", "AX", "XD", "XX", "DA", "DX", "XA"
        for (int i = 0; i < 7; ++i) {
            if (std::strcmp(da, da_types[i]) == 0) {
                hbatom_type = '-'; // Valid standard H-bond
                break;
            }
        }
    }

    return hbatom_type;
}

bool BasePairValidator::pattern_match(const std::string& str, const std::string& pattern) {
    // Matches legacy str_pmatch: pattern matching where '.' matches any character
    if (str.length() != pattern.length()) {
        return false;
    }
    for (size_t i = 0; i < str.length() && i < pattern.length(); i++) {
        if (pattern[i] != '.' && pattern[i] != str[i]) {
            return false;
        }
    }
    return true;
}

void BasePairValidator::load_atom_list(const std::string& x3dna_home) {
    if (atom_list_loaded_) {
        return; // Already loaded
    }

    // Get X3DNA home directory
    std::string home_dir = x3dna_home;
    if (home_dir.empty()) {
        const char* env_home = std::getenv("X3DNA");
        if (env_home) {
            home_dir = env_home;
        } else {
            // Default location
            home_dir = "/Users/jyesselman2/installs/x3dna";
        }
    }

    std::string atomlist_file = home_dir + "/config/atomlist.dat";
    std::ifstream file(atomlist_file);
    if (!file.is_open()) {
        // If file not found, use fallback logic only
        atom_list_loaded_ = true;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse line: atom_name_pattern atomic_symbol
        std::istringstream iss(line);
        std::string aname4, asym;
        if (!(iss >> aname4 >> asym)) {
            continue;
        }

        // Skip if starts with #
        if (aname4[0] == '#' || asym[0] == '#') {
            continue;
        }

        // Validate format (matches legacy get_atomlist)
        if (aname4.length() != 4) {
            continue;
        }
        if (asym.length() != 1 && asym.length() != 2) {
            continue;
        }

        // Normalize: convert to uppercase and pad atomic symbol to 2 chars
        std::transform(aname4.begin(), aname4.end(), aname4.begin(), ::toupper);
        std::transform(asym.begin(), asym.end(), asym.begin(), ::toupper);

        // Pad atomic symbol to 2 chars with space if needed
        if (asym.length() == 1) {
            asym = " " + asym;
        }

        // Store in map: pattern -> atomic symbol
        atom_list_[aname4] = asym;
    }

    file.close();
    atom_list_loaded_ = true;
}

} // namespace algorithms
} // namespace x3dna
