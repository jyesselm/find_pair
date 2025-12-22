/**
 * @file ring_atom_matcher.cpp
 * @brief Implementation of ring atom matcher
 */

#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/core/constants.hpp>
#include <algorithm>

namespace x3dna {
namespace algorithms {

namespace {

/**
 * @brief Pad an atom name to 4 characters (PDB format)
 *
 * PDB atom names are 4 characters with specific positioning:
 * - 2-character elements right-aligned in first 2 chars
 * - 1-character elements in column 2
 *
 * For simplicity, we pad to " XX " format (space + name + space)
 */
[[nodiscard]] std::string pad_atom_name(const std::string& name) {
    if (name.length() == 2) {
        return " " + name + " ";
    }
    // Already padded or longer - return as-is
    return name;
}

/**
 * @brief Convert registry atom names to PDB-padded format
 */
[[nodiscard]] std::vector<std::string> get_padded_names(const std::vector<std::string>& names) {
    std::vector<std::string> padded;
    padded.reserve(names.size());
    for (const auto& name : names) {
        padded.push_back(pad_atom_name(name));
    }
    return padded;
}

} // anonymous namespace

MatchedAtoms RingAtomMatcher::match(const core::Residue& residue, const core::Structure& standard_template,
                                    std::optional<core::ResidueType> detected_type) {
    MatchedAtoms result;

    // Determine residue type and get appropriate ring atom list
    // Use detected type if provided, otherwise use residue's type
    core::ResidueType residue_type = detected_type.has_value() ? detected_type.value() : residue.residue_type();
    std::vector<std::string> ring_atom_names = get_ring_atom_names(residue_type);

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

std::vector<std::string> RingAtomMatcher::get_ring_atom_names(core::ResidueType residue_type) {
    // Use constants as single source of truth for ring atom names
    // Pad names to PDB format for atom matching
    return get_padded_names(constants::nucleotides::ring_atoms_for_type(residue_type));
}

std::optional<core::Atom> RingAtomMatcher::find_atom_by_name(const core::Residue& residue,
                                                             const std::string& atom_name) {
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

bool RingAtomMatcher::is_purine(core::ResidueType type) {
    return constants::nucleotides::is_purine(type);
}

} // namespace algorithms
} // namespace x3dna
