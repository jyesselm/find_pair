/**
 * @file slot_predictor.cpp
 * @brief Implementation of geometry-based slot prediction
 */

#include <x3dna/algorithms/hydrogen_bond/slot/slot_predictor.hpp>
#include <x3dna/algorithms/hydrogen_bond/slot/atom_capacity.hpp>
#include <cmath>
#include <algorithm>

namespace x3dna {
namespace algorithms {
namespace hydrogen_bond {
namespace slot {

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD = PI / 180.0;

// Base connectivity tables (which atoms are bonded to which)
struct ConnectivityKey {
    char base_type;
    std::string atom_name;
    bool operator==(const ConnectivityKey& o) const {
        return base_type == o.base_type && atom_name == o.atom_name;
    }
};

struct ConnectivityHash {
    size_t operator()(const ConnectivityKey& k) const {
        return std::hash<char>()(k.base_type) ^
               (std::hash<std::string>()(k.atom_name) << 1);
    }
};

using ConnectivityMap = std::unordered_map<ConnectivityKey, std::vector<std::string>, ConnectivityHash>;

const ConnectivityMap& get_connectivity_map() {
    static const ConnectivityMap map = {
        // Adenine
        {{'A', "N6"}, {"C6"}},
        {{'A', "N1"}, {"C2", "C6"}},
        {{'A', "N3"}, {"C2", "C4"}},
        {{'A', "N7"}, {"C5", "C8"}},
        {{'A', "N9"}, {"C4", "C8"}},
        // Guanine
        {{'G', "N1"}, {"C2", "C6"}},
        {{'G', "N2"}, {"C2"}},
        {{'G', "O6"}, {"C6"}},
        {{'G', "N3"}, {"C2", "C4"}},
        {{'G', "N7"}, {"C5", "C8"}},
        {{'G', "N9"}, {"C4", "C8"}},
        // Cytosine
        {{'C', "N4"}, {"C4"}},
        {{'C', "O2"}, {"C2"}},
        {{'C', "N3"}, {"C2", "C4"}},
        {{'C', "N1"}, {"C2", "C6"}},
        // Uracil
        {{'U', "N3"}, {"C2", "C4"}},
        {{'U', "O2"}, {"C2"}},
        {{'U', "O4"}, {"C4"}},
        {{'U', "N1"}, {"C2", "C6"}},
        // Thymine
        {{'T', "N3"}, {"C2", "C4"}},
        {{'T', "O2"}, {"C2"}},
        {{'T', "O4"}, {"C4"}},
        {{'T', "N1"}, {"C2", "C6"}},
        // Ribose atoms (common to all)
        {{'A', "O2'"}, {"C2'"}},
        {{'G', "O2'"}, {"C2'"}},
        {{'C', "O2'"}, {"C2'"}},
        {{'U', "O2'"}, {"C2'"}},
        {{'T', "O2'"}, {"C2'"}},
        {{'A', "O4'"}, {"C1'", "C4'"}},
        {{'G', "O4'"}, {"C1'", "C4'"}},
        {{'C', "O4'"}, {"C1'", "C4'"}},
        {{'U', "O4'"}, {"C1'", "C4'"}},
        {{'T', "O4'"}, {"C1'", "C4'"}},
    };
    return map;
}

// Amino group atoms (NH2) - donate 2 H
const std::vector<std::pair<char, std::string>> AMINO_ATOMS = {
    {'A', "N6"}, {'C', "N4"}, {'G', "N2"}
};

// Imino group atoms (NH) - donate 1 H
const std::vector<std::pair<char, std::string>> IMINO_ATOMS = {
    {'G', "N1"}, {'U', "N3"}, {'T', "N3"}
};

// Carbonyl oxygens - accept 2 (sp2)
const std::vector<std::pair<char, std::string>> CARBONYL_ATOMS = {
    {'G', "O6"}, {'U', "O2"}, {'U', "O4"}, {'C', "O2"}, {'T', "O2"}, {'T', "O4"}
};

// Ring nitrogens (sp2 acceptors) - accept 1
const std::vector<std::pair<char, std::string>> RING_N_ACCEPTORS = {
    {'A', "N1"}, {'A', "N3"}, {'A', "N7"}, {'G', "N3"}, {'G', "N7"}, {'C', "N3"}
};

bool is_amino(char base, const std::string& atom) {
    return std::find(AMINO_ATOMS.begin(), AMINO_ATOMS.end(),
                     std::make_pair(base, atom)) != AMINO_ATOMS.end();
}

bool is_imino(char base, const std::string& atom) {
    return std::find(IMINO_ATOMS.begin(), IMINO_ATOMS.end(),
                     std::make_pair(base, atom)) != IMINO_ATOMS.end();
}

bool is_carbonyl(char base, const std::string& atom) {
    return std::find(CARBONYL_ATOMS.begin(), CARBONYL_ATOMS.end(),
                     std::make_pair(base, atom)) != CARBONYL_ATOMS.end();
}

bool is_ring_n_acceptor(char base, const std::string& atom) {
    return std::find(RING_N_ACCEPTORS.begin(), RING_N_ACCEPTORS.end(),
                     std::make_pair(base, atom)) != RING_N_ACCEPTORS.end();
}

bool is_ribose_oxygen(const std::string& atom) {
    return atom == "O2'" || atom == "O3'" || atom == "O4'" || atom == "O5'";
}

} // anonymous namespace

geometry::Vector3D SlotPredictor::rotate_around_axis(
    const geometry::Vector3D& v,
    const geometry::Vector3D& axis,
    double angle_deg) {

    double angle_rad = angle_deg * DEG_TO_RAD;
    double cos_a = std::cos(angle_rad);
    double sin_a = std::sin(angle_rad);

    // Rodrigues' rotation formula
    geometry::Vector3D v_rot = v * cos_a +
                               axis.cross(v) * sin_a +
                               axis * (axis.dot(v) * (1.0 - cos_a));
    return v_rot;
}

std::vector<std::string> SlotPredictor::get_connectivity(char base_type, const std::string& atom_name) {
    const auto& map = get_connectivity_map();
    auto it = map.find({base_type, atom_name});
    if (it != map.end()) {
        return it->second;
    }
    return {};
}

geometry::Vector3D SlotPredictor::compute_base_normal(const core::Residue& residue) {
    // Use ring atoms to compute plane normal
    std::vector<std::string> ring_atoms = {"C2", "C4", "C6", "N1", "N3"};
    std::vector<geometry::Vector3D> positions;

    for (const auto& name : ring_atoms) {
        auto atom = residue.find_atom(name);
        if (atom) {
            positions.push_back(atom->position());
        }
    }

    if (positions.size() < 3) {
        return geometry::Vector3D(0, 0, 1); // Default
    }

    // Cross product of two edge vectors
    geometry::Vector3D v1 = positions[1] - positions[0];
    geometry::Vector3D v2 = positions[2] - positions[0];
    return v1.cross(v2).normalized();
}

std::vector<HSlot> SlotPredictor::predict_sp2_amino_slots(
    const geometry::Vector3D& donor_pos,
    const geometry::Vector3D& bonded_pos,
    const geometry::Vector3D& base_normal) {

    // NH2: two hydrogens at ~120° from C-N bond, in the base plane
    geometry::Vector3D bond_dir = (donor_pos - bonded_pos).normalized();

    // First H: rotate 120° around base normal
    geometry::Vector3D h1_dir = rotate_around_axis(bond_dir, base_normal, 120.0);

    // Second H: rotate -120° around base normal
    geometry::Vector3D h2_dir = rotate_around_axis(bond_dir, base_normal, -120.0);

    return {HSlot(h1_dir, 1), HSlot(h2_dir, 1)}; // Each H can only donate once
}

std::vector<HSlot> SlotPredictor::predict_sp2_imino_slots(
    const geometry::Vector3D& donor_pos,
    const std::vector<geometry::Vector3D>& bonded_positions,
    const geometry::Vector3D& base_normal) {

    if (bonded_positions.size() < 2) {
        return {};
    }

    // NH in ring: single H pointing opposite to average of ring bonds
    geometry::Vector3D avg_bond(0, 0, 0);
    for (const auto& pos : bonded_positions) {
        avg_bond += (pos - donor_pos).normalized();
    }
    avg_bond = avg_bond.normalized();

    // H points opposite to average bond direction, in the base plane
    geometry::Vector3D h_dir = -avg_bond;

    // Project onto base plane
    double dot = h_dir.dot(base_normal);
    h_dir = (h_dir - base_normal * dot).normalized();

    return {HSlot(h_dir, 2)}; // Single H, but can bifurcate
}

std::vector<LPSlot> SlotPredictor::predict_sp2_carbonyl_slots(
    const geometry::Vector3D& acceptor_pos,
    const geometry::Vector3D& bonded_pos,
    const geometry::Vector3D& base_normal) {

    // C=O: two lone pairs at ~120° from C=O bond, in the base plane
    geometry::Vector3D bond_dir = (acceptor_pos - bonded_pos).normalized();

    // LP1: rotate 120° around base normal
    geometry::Vector3D lp1_dir = rotate_around_axis(bond_dir, base_normal, 120.0);

    // LP2: rotate -120° around base normal
    geometry::Vector3D lp2_dir = rotate_around_axis(bond_dir, base_normal, -120.0);

    return {LPSlot(lp1_dir, 1), LPSlot(lp2_dir, 1)};
}

std::vector<LPSlot> SlotPredictor::predict_sp2_ring_nitrogen_slots(
    const geometry::Vector3D& acceptor_pos,
    const std::vector<geometry::Vector3D>& bonded_positions) {

    if (bonded_positions.size() < 2) {
        return {};
    }

    // Ring N: single LP pointing opposite to average of ring bonds (in plane)
    geometry::Vector3D avg_bond(0, 0, 0);
    for (const auto& pos : bonded_positions) {
        avg_bond += (pos - acceptor_pos).normalized();
    }
    avg_bond = avg_bond.normalized();

    geometry::Vector3D lp_dir = -avg_bond;

    return {LPSlot(lp_dir, 1)}; // Single LP, cannot bifurcate (only 1 acceptor)
}

std::vector<HSlot> SlotPredictor::predict_sp3_hydroxyl_h_slots(
    const geometry::Vector3D& oxygen_pos,
    const geometry::Vector3D& bonded_carbon_pos) {

    // O-H: single H, tetrahedral geometry
    // Point away from C, allow rotation
    geometry::Vector3D bond_dir = (oxygen_pos - bonded_carbon_pos).normalized();

    // H is approximately along C-O extension (simplified)
    return {HSlot(bond_dir, 2)};
}

std::vector<LPSlot> SlotPredictor::predict_sp3_hydroxyl_lp_slots(
    const geometry::Vector3D& oxygen_pos,
    const geometry::Vector3D& bonded_carbon_pos) {

    // sp3 O: two lone pairs in tetrahedral arrangement
    geometry::Vector3D bond_dir = (bonded_carbon_pos - oxygen_pos).normalized();

    // Approximate LP directions (opposite to H, spread ~109°)
    geometry::Vector3D perp = geometry::Vector3D(1, 0, 0);
    if (std::abs(bond_dir.dot(perp)) > 0.9) {
        perp = geometry::Vector3D(0, 1, 0);
    }
    perp = bond_dir.cross(perp).normalized();

    geometry::Vector3D lp1 = rotate_around_axis(-bond_dir, perp, 54.75);
    geometry::Vector3D lp2 = rotate_around_axis(-bond_dir, perp, -54.75);

    return {LPSlot(lp1, 1), LPSlot(lp2, 1)};
}

std::vector<HSlot> SlotPredictor::predict_h_slots(
    char base_type,
    const std::string& atom_name,
    const core::Residue& residue,
    const geometry::Vector3D& base_normal) {

    int capacity = AtomCapacity::get_donor_capacity(std::string(1, base_type), atom_name);
    if (capacity == 0) {
        return {};
    }

    auto donor_atom = residue.find_atom(atom_name);
    if (!donor_atom) {
        return {};
    }
    geometry::Vector3D donor_pos = donor_atom->position();

    // Get connectivity
    auto connected = get_connectivity(base_type, atom_name);
    std::vector<geometry::Vector3D> bonded_positions;
    for (const auto& name : connected) {
        auto atom = residue.find_atom(name);
        if (atom) {
            bonded_positions.push_back(atom->position());
        }
    }

    // Amino group (NH2)
    if (is_amino(base_type, atom_name) && !bonded_positions.empty()) {
        return predict_sp2_amino_slots(donor_pos, bonded_positions[0], base_normal);
    }

    // Imino group (NH in ring)
    if (is_imino(base_type, atom_name) && bonded_positions.size() >= 2) {
        return predict_sp2_imino_slots(donor_pos, bonded_positions, base_normal);
    }

    // Ribose hydroxyl (O2', O3', O5')
    if (is_ribose_oxygen(atom_name) && !bonded_positions.empty()) {
        return predict_sp3_hydroxyl_h_slots(donor_pos, bonded_positions[0]);
    }

    return {};
}

std::vector<LPSlot> SlotPredictor::predict_lp_slots(
    char base_type,
    const std::string& atom_name,
    const core::Residue& residue,
    const geometry::Vector3D& base_normal) {

    int capacity = AtomCapacity::get_acceptor_capacity(std::string(1, base_type), atom_name);
    if (capacity == 0) {
        return {};
    }

    auto acceptor_atom = residue.find_atom(atom_name);
    if (!acceptor_atom) {
        return {};
    }
    geometry::Vector3D acceptor_pos = acceptor_atom->position();

    // Get connectivity
    auto connected = get_connectivity(base_type, atom_name);
    std::vector<geometry::Vector3D> bonded_positions;
    for (const auto& name : connected) {
        auto atom = residue.find_atom(name);
        if (atom) {
            bonded_positions.push_back(atom->position());
        }
    }

    // Carbonyl oxygen (C=O)
    if (is_carbonyl(base_type, atom_name) && !bonded_positions.empty()) {
        return predict_sp2_carbonyl_slots(acceptor_pos, bonded_positions[0], base_normal);
    }

    // Ring nitrogen acceptor
    if (is_ring_n_acceptor(base_type, atom_name) && bonded_positions.size() >= 2) {
        return predict_sp2_ring_nitrogen_slots(acceptor_pos, bonded_positions);
    }

    // Ribose oxygens
    if (is_ribose_oxygen(atom_name) && !bonded_positions.empty()) {
        return predict_sp3_hydroxyl_lp_slots(acceptor_pos, bonded_positions[0]);
    }

    return {};
}

} // namespace slot
} // namespace hydrogen_bond
} // namespace algorithms
} // namespace x3dna
