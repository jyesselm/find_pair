/**
 * @file ring_atom_matcher.cpp
 * @brief Implementation of ring atom matcher
 */

#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <algorithm>
#include <iostream>

namespace x3dna {
namespace algorithms {

// Ring atom names from RA_LIST definition
static const std::vector<std::string> RING_ATOMS_PURINE = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ",
                                                           " C5 ", " N7 ", " C8 ", " N9 "};

static const std::vector<std::string> RING_ATOMS_PYRIMIDINE = {" C4 ", " N3 ", " C2 ",
                                                               " N1 ", " C6 ", " C5 "};

// C1' is a sugar atom, not a ring atom, so it should NOT be included in ring atom matching
// (even for RNA). The legacy code confirms this - it never includes C1' in base_frame_calc.

MatchedAtoms RingAtomMatcher::match(const core::Residue& residue,
                                    const core::Structure& standard_template,
                                    std::optional<core::ResidueType> detected_type) {
    MatchedAtoms result;

    // Determine residue type and get appropriate ring atom list
    // Use detected type if provided, otherwise use residue's type
    core::ResidueType residue_type =
        detected_type.has_value() ? detected_type.value() : residue.residue_type();
    std::vector<std::string> ring_atom_names = get_ring_atom_names(residue_type);

#ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: RingAtomMatcher - residue: " << residue.name()
              << " type: " << static_cast<int>(residue_type) << "\n";
    std::cerr << "DEBUG: Looking for " << ring_atom_names.size() << " ring atoms\n";
    std::cerr << "DEBUG: Residue has " << residue.num_atoms() << " atoms\n";
#endif

    // Match atoms by name
    for (const auto& atom_name : ring_atom_names) {
        auto exp_atom = find_atom_by_name(residue, atom_name);
        auto std_atom = find_atom_by_name(standard_template, atom_name);

#ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Atom " << atom_name << " (repr: " << std::hex << std::showbase;
        for (char c : atom_name) {
            std::cerr << static_cast<int>(c) << " ";
        }
        std::cerr << std::dec << "): " << (exp_atom.has_value() ? "FOUND" : "NOT FOUND")
                  << " in residue, " << (std_atom.has_value() ? "FOUND" : "NOT FOUND")
                  << " in template\n";
        if (!exp_atom.has_value()) {
            // Show available atoms
            std::cerr << "DEBUG: Available atoms in residue: ";
            for (const auto& atom : residue.atoms()) {
                std::cerr << atom.name() << " ";
            }
            std::cerr << "\n";
        }
#endif

        // Only include matched atoms (both experimental and standard must be found)
        bool is_matched = exp_atom.has_value() && std_atom.has_value();
        if (is_matched) {
            result.experimental.push_back(exp_atom.value());
            result.standard.push_back(std_atom.value());
            result.atom_names.push_back(atom_name);
            result.num_matched++;
        }
    }

#ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: Total matched: " << result.num_matched << "\n";
#endif

    return result;
}

std::vector<std::string> RingAtomMatcher::get_ring_atom_names(core::ResidueType residue_type) {
    // Use correct ring atom lists (includes C4, excludes C1' and H)
    // NOTE: C1' is NOT a ring atom (it's a sugar atom), so it should NOT be included
    // in ring atom matching, even for RNA. The legacy code confirms this - it never
    // includes C1' in base_frame_calc matched_atoms.

    if (is_purine(residue_type)) {
        return RING_ATOMS_PURINE;
    } else {
        return RING_ATOMS_PYRIMIDINE;
    }
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
    return type == core::ResidueType::ADENINE || type == core::ResidueType::GUANINE;
}

} // namespace algorithms
} // namespace x3dna
