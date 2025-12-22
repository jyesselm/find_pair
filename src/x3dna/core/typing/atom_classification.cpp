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

StandardAtom AtomClassifier::get_standard_atom(const std::string& atom_name) {
    // Static map for O(1) lookup - initialized once
    static const std::unordered_map<std::string, StandardAtom> atom_map = {
        // Nucleotide ring atoms
        {"C4", StandardAtom::C4},
        {"N3", StandardAtom::N3},
        {"C2", StandardAtom::C2},
        {"N1", StandardAtom::N1},
        {"C6", StandardAtom::C6},
        {"C5", StandardAtom::C5},
        {"N7", StandardAtom::N7},
        {"C8", StandardAtom::C8},
        {"N9", StandardAtom::N9},

        // Nucleotide exocyclic atoms
        {"O6", StandardAtom::O6},
        {"N6", StandardAtom::N6},
        {"O2", StandardAtom::O2},
        {"N2", StandardAtom::N2},
        {"O4", StandardAtom::O4},
        {"N4", StandardAtom::N4},
        {"C5M", StandardAtom::C5M},
        {"C7", StandardAtom::C7},

        // Nucleotide sugar atoms
        {"C1'", StandardAtom::C1_PRIME},
        {"C2'", StandardAtom::C2_PRIME},
        {"C3'", StandardAtom::C3_PRIME},
        {"C4'", StandardAtom::C4_PRIME},
        {"C5'", StandardAtom::C5_PRIME},
        {"O2'", StandardAtom::O2_PRIME},
        {"O3'", StandardAtom::O3_PRIME},
        {"O4'", StandardAtom::O4_PRIME},
        {"O5'", StandardAtom::O5_PRIME},

        // Nucleotide backbone atoms
        {"P", StandardAtom::P},
        {"OP1", StandardAtom::OP1},
        {"OP2", StandardAtom::OP2},
        {"OP3", StandardAtom::OP3},
        // Legacy phosphate oxygen names
        {"O1P", StandardAtom::OP1},
        {"O2P", StandardAtom::OP2},

        // Amino acid backbone atoms
        {"N", StandardAtom::N},
        {"CA", StandardAtom::CA},
        {"C", StandardAtom::C},
        {"O", StandardAtom::O},
        {"OXT", StandardAtom::OXT},

        // Common amino acid side chain atoms
        {"CB", StandardAtom::CB},
        {"CG", StandardAtom::CG},
        {"CG1", StandardAtom::CG1},
        {"CG2", StandardAtom::CG2},
        {"CD", StandardAtom::CD},
        {"CD1", StandardAtom::CD1},
        {"CD2", StandardAtom::CD2},
        {"CE", StandardAtom::CE},
        {"CE1", StandardAtom::CE1},
        {"CE2", StandardAtom::CE2},
        {"CE3", StandardAtom::CE3},
        {"CZ", StandardAtom::CZ},
        {"CZ2", StandardAtom::CZ2},
        {"CZ3", StandardAtom::CZ3},
        {"CH2", StandardAtom::CH2},
        {"OG", StandardAtom::OG},
        {"OG1", StandardAtom::OG1},
        {"OD1", StandardAtom::OD1},
        {"OD2", StandardAtom::OD2},
        {"OE1", StandardAtom::OE1},
        {"OE2", StandardAtom::OE2},
        {"OH", StandardAtom::OH},
        {"ND1", StandardAtom::ND1},
        {"ND2", StandardAtom::ND2},
        {"NE", StandardAtom::NE},
        {"NE1", StandardAtom::NE1},
        {"NE2", StandardAtom::NE2},
        {"NH1", StandardAtom::NH1},
        {"NH2", StandardAtom::NH2},
        {"NZ", StandardAtom::NZ},
        {"SD", StandardAtom::SD},
        {"SG", StandardAtom::SG},

        // Water
        {"OW", StandardAtom::OW},
    };

    auto it = atom_map.find(atom_name);
    if (it != atom_map.end()) {
        return it->second;
    }
    return StandardAtom::UNKNOWN;
}

} // namespace typing
} // namespace core
} // namespace x3dna
