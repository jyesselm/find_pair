/**
 * @file atom_classification.cpp
 * @brief Implementation of AtomClassifier
 */

#include "x3dna/core/typing/atom_classification.hpp"

#include <cctype>
#include <set>
#include <unordered_map>
#include <algorithm>

namespace x3dna {
namespace core {
namespace typing {

ElementType AtomClassifier::get_element(const std::string& atom_name) {
    if (atom_name.empty()) {
        return ElementType::UNKNOWN;
    }

    // For trimmed names, element is at position 0 (e.g., "C1'", "N3", "O2'")
    const char element_char = atom_name[0];

    switch (element_char) {
        case 'C': return ElementType::CARBON;
        case 'N': return ElementType::NITROGEN;
        case 'O': return ElementType::OXYGEN;
        case 'H': return ElementType::HYDROGEN;
        case 'P': return ElementType::PHOSPHORUS;
        case 'S': return ElementType::SULFUR;
        default: return ElementType::UNKNOWN;
    }
}

int AtomClassifier::get_legacy_element_index(const std::string& atom_name) {
    if (atom_name.empty()) {
        return 0;
    }

    // For trimmed names, element is at position 0
    const char element_char = atom_name[0];

    switch (element_char) {
        case 'C': return 1;
        case 'O': return 2;
        case 'H': return 3;
        case 'N': return 4;
        case 'S': return 5;
        case 'P': return 6;
        default: return 0;
    }
}

bool AtomClassifier::is_backbone_atom(const std::string& atom_name) {
    static const std::set<std::string> backbone_atoms = {
        "P", "OP1", "OP2", "O1P", "O2P", "O5'", "O3'"
    };
    return backbone_atoms.count(atom_name) > 0;
}

bool AtomClassifier::is_sugar_atom(const std::string& atom_name) {
    static const std::set<std::string> sugar_atoms = {
        "C1'", "C2'", "C3'", "C4'", "C5'", "O4'", "O2'"
    };
    return sugar_atoms.count(atom_name) > 0;
}

bool AtomClassifier::is_nucleobase_atom(const std::string& atom_name) {
    return !is_backbone_atom(atom_name) && !is_sugar_atom(atom_name);
}

bool AtomClassifier::is_ring_atom(const std::string& atom_name) {
    // Atom names are stored trimmed - direct lookup
    static const std::set<std::string> ring_atoms = {
        "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"
    };

    return ring_atoms.count(atom_name) > 0;
}

bool AtomClassifier::is_mainchain_atom(const std::string& atom_name) {
    static const std::set<std::string> mainchain_atoms = {
        "N", "CA", "C", "O", "OXT"
    };
    return mainchain_atoms.count(atom_name) > 0;
}

bool AtomClassifier::is_sidechain_atom(const std::string& atom_name) {
    if (is_mainchain_atom(atom_name)) {
        return false;
    }
    // Exclude hydrogens (trimmed names start with H)
    if (!atom_name.empty() && atom_name[0] == 'H') {
        return false;
    }
    return true;
}

bool AtomClassifier::can_form_hbond(const std::string& atom_name,
                                     const std::string& allowed_elements) {
    if (atom_name.empty()) {
        return false;
    }

    // For trimmed names, element is at position 0
    std::string pattern(1, '.');
    pattern += atom_name[0];
    pattern += '.';

    return allowed_elements.find(pattern) != std::string::npos;
}

bool AtomClassifier::can_form_hbond_pair(const std::string& atom1,
                                          const std::string& atom2,
                                          const std::string& allowed_elements) {
    return can_form_hbond(atom1, allowed_elements) &&
           can_form_hbond(atom2, allowed_elements);
}

bool AtomClassifier::is_base_atom_for_hbond(const std::string& atom_name) {
    // Atom names are stored trimmed - direct comparison

    if (atom_name == "C5M") {
        return true;
    }

    if (atom_name.size() < 2) {
        return false;
    }

    // For trimmed names like "N1", "C2", etc.
    const char c0 = atom_name[0];
    const char c1 = atom_name[1];

    if (c0 == 'H' || c0 == 'P') {
        return false;
    }

    return std::isdigit(static_cast<unsigned char>(c1));
}

AtomClassification AtomClassifier::classify(const std::string& atom_name, MoleculeType molecule_type) {
    switch (molecule_type) {
        case MoleculeType::NUCLEIC_ACID:
            return classify_nucleotide_atom(atom_name);
        case MoleculeType::PROTEIN:
            return classify_protein_atom(atom_name);
        default:
            // For unknown/other types, return basic element info only
            AtomClassification result;
            result.element = get_element(atom_name);
            result.legacy_element_index = get_legacy_element_index(atom_name);
            result.location = AtomLocation::UNKNOWN;
            result.hbond_role = HBondRole::UNKNOWN;
            return result;
    }
}

AtomClassification AtomClassifier::classify_nucleotide_atom(const std::string& atom_name) {
    AtomClassification result;

    result.element = get_element(atom_name);
    result.legacy_element_index = get_legacy_element_index(atom_name);
    result.is_ring_atom = is_ring_atom(atom_name);

    if (is_backbone_atom(atom_name)) {
        result.location = AtomLocation::BACKBONE;
    } else if (is_sugar_atom(atom_name)) {
        result.location = AtomLocation::SUGAR;
    } else {
        result.location = AtomLocation::NUCLEOBASE;
    }

    // Determine H-bond role based on element
    if (result.element == ElementType::OXYGEN) {
        result.hbond_role = HBondRole::ACCEPTOR;
    } else if (result.element == ElementType::NITROGEN) {
        result.hbond_role = HBondRole::BOTH;  // N can be donor or acceptor
    } else {
        result.hbond_role = HBondRole::NONE;
    }

    return result;
}

AtomClassification AtomClassifier::classify_protein_atom(const std::string& atom_name) {
    AtomClassification result;

    result.element = get_element(atom_name);
    result.legacy_element_index = get_legacy_element_index(atom_name);
    result.is_ring_atom = false;  // Proteins don't have nucleobase rings

    if (is_mainchain_atom(atom_name)) {
        result.location = AtomLocation::MAINCHAIN;
    } else {
        result.location = AtomLocation::SIDECHAIN;
    }

    // Determine H-bond role based on element
    if (result.element == ElementType::OXYGEN) {
        result.hbond_role = HBondRole::ACCEPTOR;
    } else if (result.element == ElementType::NITROGEN) {
        result.hbond_role = HBondRole::BOTH;
    } else {
        result.hbond_role = HBondRole::NONE;
    }

    return result;
}

namespace {
// Nucleotide-specific atoms (ring, exocyclic, sugar, backbone)
const std::unordered_map<std::string, AtomType> NUCLEOTIDE_ATOM_MAP = {
    // Ring atoms
    {"C4", AtomType::C4},
    {"N3", AtomType::N3},
    {"C2", AtomType::C2},
    {"N1", AtomType::N1},
    {"C6", AtomType::C6},
    {"C5", AtomType::C5},
    {"N7", AtomType::N7},
    {"C8", AtomType::C8},
    {"N9", AtomType::N9},
    // Exocyclic atoms
    {"O6", AtomType::O6},
    {"N6", AtomType::N6},
    {"O2", AtomType::O2},
    {"N2", AtomType::N2},
    {"O4", AtomType::O4},
    {"N4", AtomType::N4},
    {"C5M", AtomType::C5M},
    {"C7", AtomType::C7},
    // Sugar atoms
    {"C1'", AtomType::C1_PRIME},
    {"C2'", AtomType::C2_PRIME},
    {"C3'", AtomType::C3_PRIME},
    {"C4'", AtomType::C4_PRIME},
    {"C5'", AtomType::C5_PRIME},
    {"O2'", AtomType::O2_PRIME},
    {"O3'", AtomType::O3_PRIME},
    {"O4'", AtomType::O4_PRIME},
    {"O5'", AtomType::O5_PRIME},
    // Backbone atoms
    {"P", AtomType::P},
    {"OP1", AtomType::OP1},
    {"OP2", AtomType::OP2},
    {"OP3", AtomType::OP3},
    {"O1P", AtomType::OP1},
    {"O2P", AtomType::OP2},
};

// Protein-specific atoms (backbone and side chain)
const std::unordered_map<std::string, AtomType> PROTEIN_ATOM_MAP = {
    // Backbone atoms
    {"N", AtomType::N},
    {"CA", AtomType::CA},
    {"C", AtomType::C},
    {"O", AtomType::O},
    {"OXT", AtomType::OXT},
    // Side chain atoms
    {"CB", AtomType::CB},
    {"CG", AtomType::CG},
    {"CG1", AtomType::CG1},
    {"CG2", AtomType::CG2},
    {"CD", AtomType::CD},
    {"CD1", AtomType::CD1},
    {"CD2", AtomType::CD2},
    {"CE", AtomType::CE},
    {"CE1", AtomType::CE1},
    {"CE2", AtomType::CE2},
    {"CE3", AtomType::CE3},
    {"CZ", AtomType::CZ},
    {"CZ2", AtomType::CZ2},
    {"CZ3", AtomType::CZ3},
    {"CH2", AtomType::CH2},
    {"OG", AtomType::OG},
    {"OG1", AtomType::OG1},
    {"OD1", AtomType::OD1},
    {"OD2", AtomType::OD2},
    {"OE1", AtomType::OE1},
    {"OE2", AtomType::OE2},
    {"OH", AtomType::OH},
    {"ND1", AtomType::ND1},
    {"ND2", AtomType::ND2},
    {"NE", AtomType::NE},
    {"NE1", AtomType::NE1},
    {"NE2", AtomType::NE2},
    {"NH1", AtomType::NH1},
    {"NH2", AtomType::NH2},
    {"NZ", AtomType::NZ},
    {"SD", AtomType::SD},
    {"SG", AtomType::SG},
};

// Water atoms
const std::unordered_map<std::string, AtomType> WATER_ATOM_MAP = {
    {"OW", AtomType::OW},
    {"O", AtomType::OW},  // Water oxygen is often just "O"
};
} // anonymous namespace

AtomType AtomClassifier::get_atom_type(const std::string& atom_name, MoleculeType molecule_type) {
    switch (molecule_type) {
        case MoleculeType::NUCLEIC_ACID: {
            auto it = NUCLEOTIDE_ATOM_MAP.find(atom_name);
            return (it != NUCLEOTIDE_ATOM_MAP.end()) ? it->second : AtomType::UNKNOWN;
        }
        case MoleculeType::PROTEIN: {
            auto it = PROTEIN_ATOM_MAP.find(atom_name);
            return (it != PROTEIN_ATOM_MAP.end()) ? it->second : AtomType::UNKNOWN;
        }
        case MoleculeType::WATER: {
            auto it = WATER_ATOM_MAP.find(atom_name);
            return (it != WATER_ATOM_MAP.end()) ? it->second : AtomType::UNKNOWN;
        }
        default:
            return AtomType::UNKNOWN;
    }
}

AtomType AtomClassifier::get_atom_type(const std::string& atom_name) {
    // Legacy behavior: check all maps (used when molecule type is unknown at construction)
    // First check nucleotide atoms (most common use case)
    auto it = NUCLEOTIDE_ATOM_MAP.find(atom_name);
    if (it != NUCLEOTIDE_ATOM_MAP.end()) {
        return it->second;
    }
    // Then check protein atoms
    it = PROTEIN_ATOM_MAP.find(atom_name);
    if (it != PROTEIN_ATOM_MAP.end()) {
        return it->second;
    }
    // Finally check water
    auto wit = WATER_ATOM_MAP.find(atom_name);
    if (wit != WATER_ATOM_MAP.end()) {
        return wit->second;
    }
    return AtomType::UNKNOWN;
}

} // namespace typing
} // namespace core
} // namespace x3dna
