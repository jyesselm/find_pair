/**
 * @file parameter_calculator.cpp
 * @brief Implementation of ParameterCalculator class
 */

#include <x3dna/algorithms/parameter_calculator.hpp>
#include <cmath>
#include <algorithm>

namespace x3dna {
namespace algorithms {

// Geometry utility functions (matching legacy implementations)

double ParameterCalculator::deg2rad(double ang) {
    return ang * PI / 180.0;
}

double ParameterCalculator::rad2deg(double ang) {
    return ang * 180.0 / PI;
}

double ParameterCalculator::magang(const geometry::Vector3D& va, const geometry::Vector3D& vb) {
    double vlen_a = va.length();
    double vlen_b = vb.length();
    if (vlen_a < XEPS || vlen_b < XEPS) {
        return 0.0;
    }
    geometry::Vector3D va_norm = va / vlen_a;
    geometry::Vector3D vb_norm = vb / vlen_b;
    double dot_val = va_norm.dot(vb_norm);
    // Clamp to [-1, 1] to avoid numerical issues
    dot_val = std::max(-1.0, std::min(1.0, dot_val));
    return rad2deg(std::acos(dot_val));
}

geometry::Matrix3D ParameterCalculator::arb_rotation(const geometry::Vector3D& axis, double angle_deg) {
    double vlen = axis.length();
    if (vlen < XEPS) {
        return geometry::Matrix3D::identity();
    }
    geometry::Vector3D va = axis / vlen;
    double angle_rad = deg2rad(angle_deg);
    double c = std::cos(angle_rad);
    double s = std::sin(angle_rad);
    double dc = 1.0 - c;

    // Build rotation matrix (Rodrigues' rotation formula)
    return geometry::Matrix3D(c + dc * va.x() * va.x(), va.x() * va.y() * dc - va.z() * s,
                              va.x() * va.z() * dc + va.y() * s, va.x() * va.y() * dc + va.z() * s,
                              c + dc * va.y() * va.y(), va.y() * va.z() * dc - va.x() * s,
                              va.x() * va.z() * dc - va.y() * s, va.y() * va.z() * dc + va.x() * s,
                              c + dc * va.z() * va.z());
}

double ParameterCalculator::vec_ang(const geometry::Vector3D& va, const geometry::Vector3D& vb,
                                    const geometry::Vector3D& vref) {
    // Project va and vb onto plane perpendicular to vref
    geometry::Vector3D vref_norm = vref / vref.length();
    geometry::Vector3D va_proj = va - vref_norm * va.dot(vref_norm);
    geometry::Vector3D vb_proj = vb - vref_norm * vb.dot(vref_norm);

    double va_len = va_proj.length();
    double vb_len = vb_proj.length();

    if (va_len < XEPS || vb_len < XEPS) {
        return 0.0;
    }

    va_proj = va_proj / va_len;
    vb_proj = vb_proj / vb_len;

    double ang_deg = magang(va_proj, vb_proj);

    // Determine sign using cross product
    geometry::Vector3D cross_prod = va_proj.cross(vb_proj);
    if (cross_prod.dot(vref_norm) < 0.0) {
        ang_deg = -ang_deg;
    }

    return ang_deg;
}

geometry::Vector3D ParameterCalculator::get_vector(const geometry::Vector3D& va, const geometry::Vector3D& vref,
                                                   double deg_ang) {
    geometry::Vector3D va_proj = va;
    geometry::Vector3D vref_norm = vref / vref.length();

    // Project va onto plane perpendicular to vref
    if (va_proj.dot(vref_norm) > XEPS) {
        va_proj = va_proj - vref_norm * va_proj.dot(vref_norm);
    }

    // Rotate va_proj around vref by deg_ang
    geometry::Matrix3D rot = arb_rotation(vref_norm, deg_ang);
    geometry::Vector3D vo = rot * va_proj;
    vo = vo / vo.length();

    return vo;
}

geometry::Matrix3D ParameterCalculator::x_y_z_2_mtx(const geometry::Vector3D& x, const geometry::Vector3D& y,
                                                    const geometry::Vector3D& z) {
    // Build matrix with x, y, z as columns
    return geometry::Matrix3D(x.x(), y.x(), z.x(), x.y(), y.y(), z.y(), x.z(), y.z(), z.z());
}

// Core bpstep_par implementation (matching legacy code exactly)
void ParameterCalculator::bpstep_par_impl(const geometry::Matrix3D& r1, const geometry::Vector3D& o1,
                                          const geometry::Matrix3D& r2, const geometry::Vector3D& o2,
                                          core::BasePairStepParameters& params,
                                          core::ReferenceFrame& midstep_frame) const {
    // Legacy uses 1-based indexing, we use 0-based
    // r1[i][j] in legacy = r1.at(i-1, j-1) in modern

    // Get z-axes (third column of rotation matrices)
    geometry::Vector3D t1 = r1.column(2); // z-axis of frame 1
    geometry::Vector3D t2 = r2.column(2); // z-axis of frame 2

    // Calculate hinge vector (cross product of z-axes)
    geometry::Vector3D hinge = t1.cross(t2);
    double rolltilt = magang(t1, t2);

    // Handle degenerate case (parallel or anti-parallel z-axes)
    if (hinge.length() < XEPS && (std::abs(rolltilt - 180.0) < XEPS || rolltilt < XEPS)) {
        // Use sum of x and y axes as hinge
        geometry::Vector3D x1 = r1.column(0);
        geometry::Vector3D y1 = r1.column(1);
        geometry::Vector3D x2 = r2.column(0);
        geometry::Vector3D y2 = r2.column(1);
        hinge = x1 + y1 + x2 + y2;
    }

    // Calculate para_bp1 and para_bp2 (rotated frames)
    // Legacy: multi_matrix(temp, rot2, para_bp2) = temp * rot2
    geometry::Matrix3D temp_rot = arb_rotation(hinge, -0.5 * rolltilt);
    geometry::Matrix3D para_bp2 = temp_rot * r2;

    temp_rot = arb_rotation(hinge, 0.5 * rolltilt);
    geometry::Matrix3D para_bp1 = temp_rot * r1;

    // Calculate mstz (midstep z-axis)
    geometry::Vector3D mstz = para_bp2.column(2);

    // Calculate twist (pars[6])
    geometry::Vector3D y1_para = para_bp1.column(1);
    geometry::Vector3D y2_para = para_bp2.column(1);
    params.twist = vec_ang(y1_para, y2_para, mstz);

    // Calculate msty (midstep y-axis)
    geometry::Vector3D msty = get_vector(y1_para, mstz, 0.5 * params.twist);

    // Calculate mstx (midstep x-axis)
    geometry::Vector3D mstx = msty.cross(mstz);

    // Calculate midstep origin (average of o1 and o2)
    geometry::Vector3D mst_org = (o1 + o2) * 0.5;

    // Calculate midstep orientation matrix
    geometry::Matrix3D mst_orien = x_y_z_2_mtx(mstx, msty, mstz);
    midstep_frame = core::ReferenceFrame(mst_orien, mst_org);

    // Calculate displacement vector
    geometry::Vector3D t1_displacement = o2 - o1;

    // Calculate Shift, Slide, Rise (pars[1], pars[2], pars[3])
    // These are the components of displacement in the midstep frame
    params.shift = t1_displacement.dot(mstx);
    params.slide = t1_displacement.dot(msty);
    params.rise = t1_displacement.dot(mstz);

    // Calculate Tilt and Roll (pars[4], pars[5])
    double phi = deg2rad(vec_ang(hinge, msty, mstz));
    params.roll = rolltilt * std::cos(phi);
    params.tilt = rolltilt * std::sin(phi);
}

core::BasePairStepParameters ParameterCalculator::calculate_step_parameters(const core::ReferenceFrame& frame1,
                                                                            const core::ReferenceFrame& frame2) const {
    core::BasePairStepParameters params;
    core::ReferenceFrame midstep_frame;

    bpstep_par_impl(frame1.rotation(), frame1.origin(), frame2.rotation(), frame2.origin(), params, midstep_frame);

    params.midstep_frame = midstep_frame;
    return params;
}

core::BasePairStepParameters ParameterCalculator::calculate_step_parameters(const core::BasePair& pair1,
                                                                            const core::BasePair& pair2) const {
    // For step parameters between consecutive pairs:
    // Use frame1 from each pair (the first residue's frame)
    return calculate_step_parameters(pair1.frame1(), pair2.frame1());
}

core::BasePairStepParameters ParameterCalculator::calculate_step_parameters_for_pair(const core::BasePair& pair) const {
    // For a single base pair (for bp_type_id calculation):
    // Legacy: bpstep_par(r2, org[j], r1, org[i], ...)
    // This means: frame2 (residue j) is first, frame1 (residue i) is second
    // Note: Legacy order is reversed (frame2 first, frame1 second)
    return calculate_step_parameters(pair.frame2(), pair.frame1());
}

core::HelicalParameters ParameterCalculator::calculate_helical_parameters(const core::BasePair& pair1,
                                                                          const core::BasePair& pair2) const {
    // Use frame1 from each pair (matching legacy refs_i_j behavior)
    return calculate_helical_parameters_impl(pair1.frame1(), pair2.frame1());
}

core::HelicalParameters ParameterCalculator::calculate_helical_parameters_impl(
    const core::ReferenceFrame& frame1, const core::ReferenceFrame& frame2) const {
    core::HelicalParameters params;

    // Get rotation matrices and origins
    const geometry::Matrix3D& rot1 = frame1.rotation();
    const geometry::Vector3D& org1 = frame1.origin();
    const geometry::Matrix3D& rot2 = frame2.rotation();
    const geometry::Vector3D& org2 = frame2.origin();

    // Legacy uses 1-based indexing, we use 0-based
    // rot[i][j] in legacy = rot.at(i-1, j-1) in modern
    // Columns: 0=x, 1=y, 2=z

    // Calculate axis_h from cross product of (rot2 x-axis - rot1 x-axis) and (rot2 y-axis - rot1
    // y-axis)
    geometry::Vector3D x1 = rot1.column(0); // rot1 x-axis
    geometry::Vector3D y1 = rot1.column(1); // rot1 y-axis
    geometry::Vector3D z1 = rot1.column(2); // rot1 z-axis
    geometry::Vector3D x2 = rot2.column(0); // rot2 x-axis
    geometry::Vector3D y2 = rot2.column(1); // rot2 y-axis
    geometry::Vector3D z2 = rot2.column(2); // rot2 z-axis

    geometry::Vector3D t1 = x2 - x1; // rot2[i][1] - rot1[i][1]
    geometry::Vector3D t2 = y2 - y1; // rot2[i][2] - rot1[i][2]
    geometry::Vector3D axis_h = t1.cross(t2);

    // Normalize axis_h
    double vlen = axis_h.length();
    if (vlen < XEPS) {
        // Default to z-axis if degenerate
        axis_h = geometry::Vector3D(0.0, 0.0, 1.0);
    } else {
        axis_h = axis_h / vlen;
    }

    // Calculate TipInc1: angle between axis_h and rot1 z-axis
    double TipInc1 = magang(axis_h, z1);
    geometry::Vector3D hinge1 = axis_h.cross(z1);
    if (hinge1.length() < XEPS) {
        hinge1 = geometry::Vector3D(1.0, 0.0, 0.0); // Default hinge
    } else {
        hinge1 = hinge1 / hinge1.length();
    }

    // Rotate rot1 to align z-axis with axis_h
    geometry::Matrix3D temp_rot = arb_rotation(hinge1, -TipInc1);
    geometry::Matrix3D rot1_h = temp_rot * rot1;

    // Calculate TipInc2: angle between axis_h and rot2 z-axis
    double TipInc2 = magang(axis_h, z2);
    geometry::Vector3D hinge2 = axis_h.cross(z2);
    if (hinge2.length() < XEPS) {
        hinge2 = geometry::Vector3D(1.0, 0.0, 0.0); // Default hinge
    } else {
        hinge2 = hinge2 / hinge2.length();
    }

    // Rotate rot2 to align z-axis with axis_h
    temp_rot = arb_rotation(hinge2, -TipInc2);
    geometry::Matrix3D rot2_h = temp_rot * rot2;

    // Calculate helical midstep frame orientation
    // Sum of x and y axes from rotated frames
    geometry::Vector3D h_x = rot1_h.column(0) + rot2_h.column(0);
    geometry::Vector3D h_y = rot1_h.column(1) + rot2_h.column(1);

    // Normalize
    double h_x_len = h_x.length();
    double h_y_len = h_y.length();
    if (h_x_len > XEPS)
        h_x = h_x / h_x_len;
    if (h_y_len > XEPS)
        h_y = h_y / h_y_len;

    // Build helical midstep orientation matrix
    geometry::Matrix3D mst_orienH = x_y_z_2_mtx(h_x, h_y, axis_h);

    // Calculate h-Twist (pars[6]): angle between y-axes of rotated frames around axis_h
    geometry::Vector3D y1_h = rot1_h.column(1); // rot1_h y-axis (column 2 in 1-based)
    geometry::Vector3D y2_h = rot2_h.column(1); // rot2_h y-axis (column 2 in 1-based)
    params.twist = vec_ang(y1_h, y2_h, axis_h);

    // Calculate h-Rise (pars[3]): projection of org2-org1 onto axis_h
    geometry::Vector3D org_diff = org2 - org1;
    params.rise = org_diff.dot(axis_h);

    // Calculate Tip and Inclination (pars[5] and pars[4])
    // phi: angle between hinge1 and y1_h (rot1_h y-axis) around axis_h
    // Matching legacy: phi = vec_ang(hinge1, t1, axis_h) where t1 = rot1_h y-axis
    double phi = deg2rad(vec_ang(hinge1, y1_h, axis_h));
    params.tip = TipInc1 * std::cos(phi);
    params.inclination = TipInc1 * std::sin(phi);

    // Calculate X-disp and Y-disp (pars[1] and pars[2])
    // Project displacement onto helical frame
    geometry::Vector3D t1_proj = org_diff - (params.rise * axis_h);

    // Calculate org1_h and org2_h on helical axis
    constexpr double HTWIST0 = 0.05; // Threshold for small twist angle (matches legacy)
    geometry::Vector3D org1_h, org2_h;

    if (std::abs(params.twist) < HTWIST0) {
        // Small twist: use midpoint
        org1_h = org1 + 0.5 * t1_proj;
    } else {
        // Large twist: calculate axis displacement
        geometry::Vector3D AD_axis = get_vector(t1_proj, axis_h, 90.0 - params.twist / 2.0);
        double AD_mag = 0.5 * t1_proj.length() / std::sin(deg2rad(params.twist / 2.0));
        org1_h = org1 + AD_mag * AD_axis;
    }

    org2_h = org1_h + params.rise * axis_h;

    // Update helical midstep origin now that we have org1_h and org2_h
    geometry::Vector3D mst_orgH = (org1_h + org2_h) * 0.5;
    params.midstep_frame = core::ReferenceFrame(mst_orienH, mst_orgH);

    // Calculate X-disp and Y-disp from org1 to org1_h
    // Legacy: ddxyz(org1_h, org1, t1) which calculates t1 = org1 - org1_h
    geometry::Vector3D disp = org1 - org1_h;

    // Project onto rotated frame axes
    geometry::Vector3D rot1_h_x = rot1_h.column(0);
    geometry::Vector3D rot1_h_y = rot1_h.column(1);

    params.x_displacement = disp.dot(rot1_h_x);
    params.y_displacement = disp.dot(rot1_h_y);

    return params;
}

std::vector<core::BasePairStepParameters> ParameterCalculator::calculate_all_step_parameters(
    const std::vector<core::BasePair>& pairs) const {
    std::vector<core::BasePairStepParameters> result;
    if (pairs.size() < 2) {
        return result;
    }

    for (size_t i = 0; i < pairs.size() - 1; ++i) {
        result.push_back(calculate_step_parameters(pairs[i], pairs[i + 1]));
    }

    return result;
}

core::ReferenceFrame ParameterCalculator::calculate_midstep_frame(const core::ReferenceFrame& frame1,
                                                                  const core::ReferenceFrame& frame2) const {
    core::BasePairStepParameters params;
    core::ReferenceFrame midstep_frame;

    bpstep_par_impl(frame1.rotation(), frame1.origin(), frame2.rotation(), frame2.origin(), params, midstep_frame);

    return midstep_frame;
}

core::ReferenceFrame ParameterCalculator::calculate_pair_frame(const core::ReferenceFrame& frame1,
                                                               const core::ReferenceFrame& frame2) const {
    // This matches legacy cehs_average behavior for a 2-base pair
    // Legacy code:
    //   1. Start with mst = frame1
    //   2. For frame2: if z-axes anti-parallel, reverse y and z columns
    //   3. Call bpstep_par(frame2_modified, frame1)

    // Get z-axes
    geometry::Vector3D z1 = frame1.z_axis();
    geometry::Vector3D z2 = frame2.z_axis();

    // Check if z-axes are anti-parallel (negative dot product)
    double z_dot = z1.dot(z2);

    geometry::Matrix3D r2_modified = frame2.rotation();

    if (z_dot < 0.0) {
        // Reverse y and z columns (legacy: reverse_y_z_columns)
        // Column 1 (y-axis) and Column 2 (z-axis) are negated
        r2_modified = geometry::Matrix3D(r2_modified.at(0, 0), -r2_modified.at(0, 1), -r2_modified.at(0, 2),
                                         r2_modified.at(1, 0), -r2_modified.at(1, 1), -r2_modified.at(1, 2),
                                         r2_modified.at(2, 0), -r2_modified.at(2, 1), -r2_modified.at(2, 2));
    }

    // Calculate midstep frame using bpstep_par
    // Legacy order: bpstep_par(bi, org[ik], mst, morg, ...)
    // where bi=frame2_modified, org[ik]=frame2.origin, mst=frame1, morg=frame1.origin
    core::BasePairStepParameters params;
    core::ReferenceFrame midstep_frame;

    bpstep_par_impl(r2_modified, frame2.origin(), frame1.rotation(), frame1.origin(), params, midstep_frame);

    return midstep_frame;
}

} // namespace algorithms
} // namespace x3dna
