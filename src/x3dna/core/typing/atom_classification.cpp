/**
 * @file atom_classification.cpp
 * @brief Implementation of AtomClassifier
 */

#include "x3dna/core/typing/atom_classification.hpp"

#include <cctype>
#include <set>
#include <algorithm>

namespace x3dna {
namespace core {
namespace typing {

namespace {

// Trim whitespace from string_view
std::string_view trim(std::string_view sv) {
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front()))) {
        sv.remove_prefix(1);
    }
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back()))) {
        sv.remove_suffix(1);
    }
    return sv;
}

} // anonymous namespace

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
    // Trim whitespace
    std::string_view sv = trim(atom_name);

    // Ring atoms: N1, C2, N3, C4, C5, C6, N7, C8, N9
    static const std::set<std::string> ring_atoms = {
        "N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"
    };

    return ring_atoms.count(std::string(sv)) > 0;
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
    // Trim the name for comparison
    std::string name = atom_name;
    name.erase(0, name.find_first_not_of(" \t\n\r"));
    name.erase(name.find_last_not_of(" \t\n\r") + 1);

    if (name == "C5M") {
        return true;
    }

    if (name.size() < 2) {
        return false;
    }

    // For trimmed names like "N1", "C2", etc.
    const char c0 = name[0];
    const char c1 = name[1];

    if (c0 == 'H' || c0 == 'P') {
        return false;
    }

    return std::isdigit(static_cast<unsigned char>(c1));
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

} // namespace typing
} // namespace core
} // namespace x3dna
