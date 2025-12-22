#include "x3dna/core/atom_classification.hpp"

#include <cctype>
#include <set>

namespace x3dna::core::atom_classification {

bool is_backbone_atom(const std::string& atom_name) {
    static const std::set<std::string> backbone_atoms = {" P  ", " OP1", " OP2", " O1P", " O2P", " O5'", " O3'"};
    return backbone_atoms.count(atom_name) > 0;
}

bool is_sugar_atom(const std::string& atom_name) {
    static const std::set<std::string> sugar_atoms = {" C1'", " C2'", " C3'", " C4'", " C5'", " O4'", " O2'"};
    return sugar_atoms.count(atom_name) > 0;
}

bool is_nucleobase_atom(const std::string& atom_name) {
    return !is_backbone_atom(atom_name) && !is_sugar_atom(atom_name);
}

bool is_base_atom_for_hbond(const std::string& atom_name) {
    if (atom_name == " C5M") {
        return true;
    }

    if (atom_name.size() != 4) {
        return false;
    }

    const char c1 = atom_name[1];
    const char c2 = atom_name[2];

    if (c1 == 'H' || c1 == 'P') {
        return false;
    }

    return std::isdigit(static_cast<unsigned char>(c2));
}

bool can_form_hbond(const std::string& atom1, const std::string& atom2, const std::string& allowed_elements) {
    if (atom1.size() < 2 || atom2.size() < 2) {
        return false;
    }

    std::string s1(1, '.');
    s1 += atom1[1];
    s1 += '.';

    std::string s2(1, '.');
    s2 += atom2[1];
    s2 += '.';

    const bool atom1_allowed = allowed_elements.find(s1) != std::string::npos;
    const bool atom2_allowed = allowed_elements.find(s2) != std::string::npos;

    return atom1_allowed && atom2_allowed;
}

int get_element_index(const std::string& atom_name) {
    if (atom_name.length() < 2) {
        return 0;
    }

    const char element_char = atom_name[1];

    switch (element_char) {
        case 'C':
            return 1;
        case 'O':
            return 2;
        case 'H':
            return 3;
        case 'N':
            return 4;
        case 'S':
            return 5;
        case 'P':
            return 6;
        default:
            return 0;
    }
}

} // namespace x3dna::core::atom_classification
