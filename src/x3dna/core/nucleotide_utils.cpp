/**
 * @file nucleotide_utils.cpp
 * @brief Implementation of nucleotide utility functions
 */

#include <x3dna/core/nucleotide_utils.hpp>
#include <x3dna/core/residue.hpp>

namespace x3dna {
namespace core {

std::vector<Atom> ring_atoms(const Residue& residue) {
    std::vector<Atom> ring;
    for (const auto& atom : residue.atoms()) {
        if (atom.is_ring_atom()) {
            ring.push_back(atom);
        }
    }
    return ring;
}

char one_letter_code(const Residue& residue) {
    return TypeRegistry::instance().get_one_letter_code(residue.name());
}

int ry_classification(const Residue& residue) {
    if (residue.classification().is_purine()) {
        return 1;
    }
    char code = one_letter_code(residue);
    if (code == 'C' || code == 'T' || code == 'U' || code == 'P' ||
        code == 'c' || code == 't' || code == 'u') {
        return 0;
    }
    return -1;
}

} // namespace core
} // namespace x3dna
