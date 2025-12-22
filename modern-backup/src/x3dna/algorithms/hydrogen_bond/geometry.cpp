/**
 * @file geometry.cpp
 * @brief Implementation of H-bond geometric calculations
 */

#include <x3dna/algorithms/hydrogen_bond/geometry.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/atom_classification.hpp>
#include <x3dna/algorithms/hydrogen_bond/atom_utils.hpp>
#include <unordered_map>
#include <cmath>
#include <cctype>
#include <algorithm>

namespace x3dna {
namespace algorithms {

using namespace x3dna::geometry;
using namespace x3dna::core;

namespace {

/**
 * @brief Neighbor atom lookup table
 *
 * Maps H-bond capable atoms to their reference neighbor atoms for angle calculations.
 * Organized by base type for clarity.
 */
const std::unordered_map<std::string, std::string> NEIGHBOR_LOOKUP = {
    // === ADENINE ===
    {" N6 ", " C6 "}, // Amino N - attached to C6
    {" N1 ", " C2 "}, // Ring N1 - use C2 as reference
    {" N3 ", " C2 "}, // Ring N3 - use C2 as reference
    {" N7 ", " C5 "}, // Ring N7 - use C5 as reference
    {" N9 ", " C4 "}, // Glycosidic N - use C4

    // === GUANINE ===
    {" O6 ", " C6 "}, // Carbonyl O
    {" N2 ", " C2 "}, // Amino N
    // N1, N3, N7, N9 same as adenine

    // === CYTOSINE ===
    {" N4 ", " C4 "}, // Amino N
    {" O2 ", " C2 "}, // Carbonyl O

    // === URACIL / THYMINE ===
    {" O4 ", " C4 "}, // Carbonyl O
    // O2, N3 patterns covered above

    // === BACKBONE (PHOSPHATE) ===
    {" O1P", " P  "},
    {" O2P", " P  "},
    {" OP1", " P  "}, // PDBv3 name
    {" OP2", " P  "}, // PDBv3 name
    {" O5'", " C5'"},
    {" O3'", " C3'"},

    // === SUGAR (RIBOSE) ===
    {" O4'", " C4'"}, // Ring O
    {" O2'", " C2'"}, // 2'-OH (RNA only)
};

} // anonymous namespace

double HBondGeometry::calculate_angle(const Vector3D& a, const Vector3D& b, const Vector3D& c) {
    Vector3D ba = a - b;
    Vector3D bc = c - b;

    const double dot = ba.dot(bc);
    const double mag_ba = ba.length();
    const double mag_bc = bc.length();

    if (mag_ba < 1e-10 || mag_bc < 1e-10) {
        return 0.0;
    }

    double cos_angle = dot / (mag_ba * mag_bc);
    cos_angle = std::clamp(cos_angle, -1.0, 1.0);

    return std::acos(cos_angle) * 180.0 / M_PI;
}

double HBondGeometry::calculate_dihedral(const Vector3D& a, const Vector3D& b, const Vector3D& c, const Vector3D& d) {
    // Compute vectors along bonds
    const Vector3D b1 = b - a;
    const Vector3D b2 = c - b;
    const Vector3D b3 = d - c;

    // Compute normals to planes
    const Vector3D n1 = b1.cross(b2);
    const Vector3D n2 = b2.cross(b3);

    // Compute dihedral angle
    const double n1_len = n1.length();
    const double n2_len = n2.length();

    if (n1_len < 1e-10 || n2_len < 1e-10) {
        return 0.0;
    }

    // Angle between normals
    double cos_angle = n1.dot(n2) / (n1_len * n2_len);
    cos_angle = std::clamp(cos_angle, -1.0, 1.0);

    double angle = std::acos(cos_angle) * 180.0 / M_PI;

    // Determine sign using triple product
    const double sign = b2.dot(n1.cross(n2));
    if (sign < 0.0) {
        angle = -angle;
    }

    return angle;
}

std::string HBondGeometry::get_neighbor_atom_name(const std::string& hbond_atom_name) {
    const auto it = NEIGHBOR_LOOKUP.find(hbond_atom_name);
    if (it != NEIGHBOR_LOOKUP.end()) {
        return it->second;
    }
    return "";
}

std::optional<Vector3D> HBondGeometry::find_neighbor_position(const std::string& hbond_atom_name,
                                                              const Residue& residue) {
    const std::string neighbor_name = get_neighbor_atom_name(hbond_atom_name);
    if (neighbor_name.empty()) {
        return std::nullopt;
    }

    const auto atom_opt = residue.find_atom(neighbor_name);
    if (!atom_opt.has_value()) {
        return std::nullopt;
    }

    return atom_opt.value().position();
}

bool HBondGeometry::are_elements_hbond_compatible(const std::string& atom1_name, const std::string& atom2_name,
                                                  const std::string& allowed_elements) {
    return hydrogen_bond::good_hb_atoms(atom1_name, atom2_name, allowed_elements);
}

bool HBondGeometry::is_nucleobase_atom(const std::string& atom_name) {
    return core::atom_classification::is_nucleobase_atom(atom_name);
}

bool HBondGeometry::is_backbone_atom(const std::string& atom_name) {
    return core::atom_classification::is_backbone_atom(atom_name);
}

bool HBondGeometry::is_sugar_atom(const std::string& atom_name) {
    return core::atom_classification::is_sugar_atom(atom_name);
}

HBondContext HBondGeometry::determine_context(const std::string& atom1_name, const std::string& atom2_name) {
    const bool a1_base = is_nucleobase_atom(atom1_name);
    const bool a2_base = is_nucleobase_atom(atom2_name);
    const bool a1_backbone = is_backbone_atom(atom1_name);
    const bool a2_backbone = is_backbone_atom(atom2_name);
    const bool a1_sugar = is_sugar_atom(atom1_name);
    const bool a2_sugar = is_sugar_atom(atom2_name);

    if (a1_base && a2_base) {
        return HBondContext::BASE_BASE;
    }
    if ((a1_base && a2_backbone) || (a1_backbone && a2_base)) {
        return HBondContext::BASE_BACKBONE;
    }
    if (a1_backbone && a2_backbone) {
        return HBondContext::BACKBONE_BACKBONE;
    }
    if ((a1_base && a2_sugar) || (a1_sugar && a2_base)) {
        return HBondContext::BASE_SUGAR;
    }
    if (a1_sugar && a2_sugar) {
        return HBondContext::SUGAR_SUGAR;
    }

    return HBondContext::UNKNOWN;
}

} // namespace algorithms
} // namespace x3dna
