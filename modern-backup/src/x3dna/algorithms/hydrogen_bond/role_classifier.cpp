/**
 * @file role_classifier.cpp
 * @brief Implementation of H-bond validation and donor/acceptor classification
 */

#include <x3dna/algorithms/hydrogen_bond/role_classifier.hpp>
#include <cctype>
#include <cstring>

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;

namespace {

// Legacy BASE_LIST = "ACGITU" (A=0, C=1, G=2, I=3, T=4, U=5)
constexpr const char* BASE_LIST = "ACGITU";

// Valid donor-acceptor combinations
// Legacy: "AD", "AX", "XD", "XX", "DA", "DX", "XA"
constexpr const char* DA_TYPES[] = {"AD", "AX", "XD", "XX", "DA", "DX", "XA"};

// Backbone atom patterns with roles
// Legacy bb_da array format: " O1P_A" means O1P is an Acceptor
struct BackboneAtom {
    const char* name;
    char role; // 'A' = Acceptor, 'X' = Either
};

constexpr BackboneAtom BACKBONE_ATOMS[] = {{" O1P", 'A'}, {" O2P", 'A'}, {" O5'", 'A'},
                                           {" O4'", 'A'}, {" O3'", 'A'}, {" O2'", 'X'}};

// Base-specific atom patterns
// Legacy base_da array format: [base_index][atom_index] = " N9 _?"
// where '?' is the role (D=Donor, A=Acceptor, X=Either)
struct BaseAtom {
    const char* name;
    char role;
};

// Adenine (index 0 in BASE_LIST)
constexpr BaseAtom ADENINE_ATOMS[] = {{" N9 ", '?'}, // Glycosidic - Either
                                      {" N7 ", 'A'}, // Ring N - Acceptor
                                      {" N6 ", 'D'}, // Amino group - Donor
                                      {" N1 ", 'A'}, // Ring N - Acceptor
                                      {" N3 ", 'A'}, // Ring N - Acceptor
                                      {nullptr, '\0'}};

// Cytosine (index 1)
constexpr BaseAtom CYTOSINE_ATOMS[] = {{" N1 ", '?'}, // Glycosidic - Either
                                       {" O2 ", 'A'}, // Carbonyl - Acceptor
                                       {" N3 ", 'A'}, // Ring N - Acceptor
                                       {" N4 ", 'D'}, // Amino group - Donor
                                       {nullptr, '\0'}};

// Guanine (index 2)
constexpr BaseAtom GUANINE_ATOMS[] = {{" N9 ", '?'}, // Glycosidic - Either
                                      {" N7 ", 'A'}, // Ring N - Acceptor
                                      {" O6 ", 'A'}, // Carbonyl - Acceptor
                                      {" N1 ", 'D'}, // NH in ring - Donor
                                      {" N2 ", 'D'}, // Amino group - Donor
                                      {" N3 ", 'A'}, // Ring N - Acceptor
                                      {nullptr, '\0'}};

// Inosine (index 3) - like Guanine but no N2
constexpr BaseAtom INOSINE_ATOMS[] = {{" N9 ", '?'}, // Glycosidic - Either
                                      {" N7 ", 'A'}, // Ring N - Acceptor
                                      {" O6 ", 'A'}, // Carbonyl - Acceptor
                                      {" N1 ", 'D'}, // NH in ring - Donor
                                      {" N3 ", 'A'}, // Ring N - Acceptor
                                      {nullptr, '\0'}};

// Thymine (index 4) and Uracil (index 5) - same pattern
constexpr BaseAtom THYMINE_URACIL_ATOMS[] = {{" N1 ", '?'}, // Glycosidic - Either
                                             {" O2 ", 'A'}, // Carbonyl - Acceptor
                                             {" N3 ", 'D'}, // NH in ring - Donor
                                             {" O4 ", 'A'}, // Carbonyl - Acceptor
                                             {nullptr, '\0'}};

// Map base index to atom table
const BaseAtom* get_base_atoms(int base_index) {
    switch (base_index) {
        case 0:
            return ADENINE_ATOMS;
        case 1:
            return CYTOSINE_ATOMS;
        case 2:
            return GUANINE_ATOMS;
        case 3:
            return INOSINE_ATOMS;
        case 4:
        case 5:
            return THYMINE_URACIL_ATOMS;
        default:
            return nullptr;
    }
}

// Convert role character to HBondAtomRole
HBondAtomRole char_to_role(char role_char) {
    switch (role_char) {
        case 'D':
            return HBondAtomRole::DONOR;
        case 'A':
            return HBondAtomRole::ACCEPTOR;
        case '?':
        case 'X':
            return HBondAtomRole::EITHER;
        default:
            return HBondAtomRole::UNKNOWN;
    }
}

} // anonymous namespace

HBondAtomRole HBondRoleClassifier::get_atom_role(char base, const std::string& atom_name) {
    // Ensure we have 4-character atom name
    if (atom_name.length() < 4) {
        return HBondAtomRole::UNKNOWN;
    }

    // Check backbone atoms first
    for (const auto& bb_atom : BACKBONE_ATOMS) {
        if (atom_name.substr(0, 4) == bb_atom.name) {
            return char_to_role(bb_atom.role);
        }
    }

    // Find base index
    const char* base_ptr = std::strchr(BASE_LIST, std::toupper(static_cast<unsigned char>(base)));
    if (base_ptr == nullptr) {
        return HBondAtomRole::UNKNOWN;
    }

    int base_index = base_ptr - BASE_LIST;
    const BaseAtom* base_atoms = get_base_atoms(base_index);
    if (base_atoms == nullptr) {
        return HBondAtomRole::UNKNOWN;
    }

    // Check base-specific atoms
    for (int i = 0; base_atoms[i].name != nullptr; ++i) {
        if (atom_name.substr(0, 4) == base_atoms[i].name) {
            return char_to_role(base_atoms[i].role);
        }
    }

    return HBondAtomRole::UNKNOWN;
}

// Get the raw role character for an atom (preserves 'X' vs '?' distinction)
// Legacy behavior: If base is unknown (not in BASE_LIST), return '\0' for ALL atoms
// including backbone atoms. This ensures bonds with unknown bases return NON_STANDARD.
char get_atom_role_char(char base, const std::string& atom_name) {
    // Ensure we have 4-character atom name
    if (atom_name.length() < 4) {
        return '\0';
    }

    // Find base index FIRST (legacy behavior)
    // If base is unknown, we can't classify ANY atoms from that residue
    const char* base_ptr = std::strchr(BASE_LIST, std::toupper(static_cast<unsigned char>(base)));
    if (base_ptr == nullptr) {
        // Unknown base (e.g., PSU='P', not in "ACGITU")
        // Legacy returns '*' for all atoms from unknown bases
        return '\0';
    }

    int base_index = base_ptr - BASE_LIST;

    // Check backbone atoms - these use 'X' for EITHER
    for (const auto& bb_atom : BACKBONE_ATOMS) {
        if (atom_name.substr(0, 4) == bb_atom.name) {
            return bb_atom.role; // Returns 'X' for O2'
        }
    }

    // Check base-specific atoms - these use '?' for EITHER
    const BaseAtom* base_atoms = get_base_atoms(base_index);
    if (base_atoms == nullptr) {
        return '\0';
    }

    for (int i = 0; base_atoms[i].name != nullptr; ++i) {
        if (atom_name.substr(0, 4) == base_atoms[i].name) {
            return base_atoms[i].role; // Returns '?' for N9, N1 (glycosidic)
        }
    }

    return '\0';
}

HBondClassification HBondRoleClassifier::classify_bond(char base1, char base2, const std::string& atom1,
                                                       const std::string& atom2) {
    // Get raw role characters (preserves 'X' vs '?' distinction)
    // Legacy code uses 'X' for backbone EITHER (O2') but '?' for base EITHER (N9)
    // DA_TYPES has 'X', so backbone 'X' matches but base '?' fails strcmp
    char role1_char = get_atom_role_char(base1, atom1);
    char role2_char = get_atom_role_char(base2, atom2);

    // Unknown atoms → NON_STANDARD (matches legacy '*')
    if (role1_char == '\0' || role2_char == '\0') {
        return HBondClassification::NON_STANDARD;
    }

    // Build donor-acceptor string
    char da[3] = {role1_char, role2_char, '\0'};

    // Check against valid combinations
    // Valid: "AD", "AX", "XD", "XX", "DA", "DX", "XA"
    for (const auto& valid_type : DA_TYPES) {
        if (std::strcmp(da, valid_type) == 0) {
            return HBondClassification::STANDARD;
        }
    }

    // Not in valid list → NON_STANDARD (matches legacy '*')
    // Invalid combinations like "AA", "DD" return '*' in legacy
    return HBondClassification::NON_STANDARD;
}

bool HBondRoleClassifier::is_good_hbond_distance(double distance, double min_dist, double max_dist) {
    return distance >= min_dist && distance <= max_dist;
}

int HBondRoleClassifier::count_good_hbonds(const std::vector<HBond>& bonds, double min_dist, double max_dist) {
    int count = 0;
    for (const auto& bond : bonds) {
        // Count STANDARD H-bonds with good distance
        if (bond.classification == HBondClassification::STANDARD &&
            is_good_hbond_distance(bond.distance, min_dist, max_dist)) {
            ++count;
        }
    }
    return count;
}

} // namespace algorithms
} // namespace x3dna
