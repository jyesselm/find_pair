/**
 * @file ring_atom_matcher.cpp
 * @brief Implementation of ring atom matcher
 */

#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/core/constants.hpp>
#include <algorithm>

namespace x3dna {
namespace algorithms {

MatchedAtoms RingAtomMatcher::match(const core::Residue& residue, const core::Structure& standard_template,
                                    std::optional<core::typing::BaseType> detected_type) {
    MatchedAtoms result;

    // Determine base type and get appropriate ring atom list
    // Use detected type if provided, otherwise use residue's base_type
    core::typing::BaseType base_type = detected_type.has_value() ? detected_type.value() : residue.base_type();
    std::vector<std::string> ring_atom_names = get_ring_atom_names(base_type);

    // Match atoms by name
    for (const auto& atom_name : ring_atom_names) {
        auto exp_atom = find_atom_by_name(residue, atom_name);
        auto std_atom = find_atom_by_name(standard_template, atom_name);

        // Only include matched atoms (both experimental and standard must be found)
        bool is_matched = exp_atom.has_value() && std_atom.has_value();
        if (is_matched) {
            result.experimental.push_back(exp_atom.value());
            result.standard.push_back(std_atom.value());
            result.atom_names.push_back(atom_name);
            result.num_matched++;
        }
    }

    return result;
}

std::vector<std::string> RingAtomMatcher::get_ring_atom_names(core::typing::BaseType base_type) {
    // Use constants as single source of truth for ring atom names
    // Return trimmed names directly - no padding needed since atoms are stored trimmed
    return constants::nucleotides::ring_atoms_for_type(base_type);
}

std::optional<core::Atom> RingAtomMatcher::find_atom_by_name(const core::Residue& residue,
                                                             const std::string& atom_name) {
    // Atom names are stored trimmed - direct comparison
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == atom_name) {
            return atom;
        }
    }
    return std::nullopt;
}

std::optional<core::Atom> RingAtomMatcher::find_atom_by_name(const core::Structure& structure,
                                                             const std::string& atom_name) {
    // Search through all chains and residues
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            auto atom = find_atom_by_name(residue, atom_name);
            if (atom.has_value()) {
                return atom;
            }
        }
    }
    return std::nullopt;
}

bool RingAtomMatcher::is_purine(core::typing::BaseType type) {
    return core::typing::is_purine(type);
}

} // namespace algorithms
} // namespace x3dna
