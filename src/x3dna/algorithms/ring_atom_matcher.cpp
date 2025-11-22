/**
 * @file ring_atom_matcher.cpp
 * @brief Implementation of ring atom matcher
 */

#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <algorithm>

namespace x3dna {
namespace algorithms {

// Ring atom names from RA_LIST definition
static const std::vector<std::string> RING_ATOMS_PURINE = {
    " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
};

static const std::vector<std::string> RING_ATOMS_PYRIMIDINE = {
    " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "
};

// RNA adds C1' to the beginning
static const std::string C1_PRIME = " C1'";

MatchedAtoms RingAtomMatcher::match(const core::Residue& residue,
                                     const core::Structure& standard_template,
                                     bool is_rna) {
    MatchedAtoms result;
    
    // Determine residue type and get appropriate ring atom list
    core::ResidueType residue_type = residue.residue_type();
    std::vector<std::string> ring_atom_names = get_ring_atom_names(residue_type, is_rna);
    
    // Match atoms by name
    for (const auto& atom_name : ring_atom_names) {
        auto exp_atom = find_atom_by_name(residue, atom_name);
        auto std_atom = find_atom_by_name(standard_template, atom_name);
        
        if (exp_atom.has_value() && std_atom.has_value()) {
            result.experimental.push_back(exp_atom.value());
            result.standard.push_back(std_atom.value());
            result.atom_names.push_back(atom_name);
            result.num_matched++;
        }
    }
    
    return result;
}

std::vector<std::string> RingAtomMatcher::get_ring_atom_names(core::ResidueType residue_type,
                                                               bool is_rna) {
    std::vector<std::string> atom_names;
    
    // For RNA, add C1' at the beginning
    if (is_rna) {
        atom_names.push_back(C1_PRIME);
    }
    
    // Add ring atoms based on purine vs pyrimidine
    if (is_purine(residue_type)) {
        atom_names.insert(atom_names.end(),
                         RING_ATOMS_PURINE.begin(),
                         RING_ATOMS_PURINE.end());
    } else {
        atom_names.insert(atom_names.end(),
                         RING_ATOMS_PYRIMIDINE.begin(),
                         RING_ATOMS_PYRIMIDINE.end());
    }
    
    return atom_names;
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

