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
                                     bool is_rna,
                                     bool exclude_c4) {
    MatchedAtoms result;
    
    // Determine residue type and get appropriate ring atom list
    core::ResidueType residue_type = residue.residue_type();
    std::vector<std::string> ring_atom_names = get_ring_atom_names(residue_type, is_rna, exclude_c4);
    
    #ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: RingAtomMatcher - residue: " << residue.name()
              << " type: " << static_cast<int>(residue_type)
              << " is_rna: " << is_rna << "\n";
    std::cerr << "DEBUG: Looking for " << ring_atom_names.size() << " ring atoms\n";
    std::cerr << "DEBUG: Residue has " << residue.num_atoms() << " atoms\n";
    #endif
    
    // Match atoms by name
    for (const auto& atom_name : ring_atom_names) {
        auto exp_atom = find_atom_by_name(residue, atom_name);
        auto std_atom = find_atom_by_name(standard_template, atom_name);
        
        #ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Atom " << atom_name << " (repr: " 
                  << std::hex << std::showbase;
        for (char c : atom_name) {
            std::cerr << static_cast<int>(c) << " ";
        }
        std::cerr << std::dec << "): "
                  << (exp_atom.has_value() ? "FOUND" : "NOT FOUND") << " in residue, "
                  << (std_atom.has_value() ? "FOUND" : "NOT FOUND") << " in template\n";
        if (!exp_atom.has_value()) {
            // Show available atoms
            std::cerr << "DEBUG: Available atoms in residue: ";
            for (const auto& atom : residue.atoms()) {
                std::cerr << atom.name() << " ";
            }
            std::cerr << "\n";
        }
        #endif
        
        if (exp_atom.has_value() && std_atom.has_value()) {
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

std::vector<std::string> RingAtomMatcher::get_ring_atom_names(core::ResidueType residue_type,
                                                               bool is_rna,
                                                               bool exclude_c4) {
    std::vector<std::string> atom_names;
    
    // For RNA, add C1' at the beginning
    if (is_rna) {
        atom_names.push_back(C1_PRIME);
    }
    
    // Add ring atoms based on purine vs pyrimidine
    std::vector<std::string> ring_atoms;
    if (is_purine(residue_type)) {
        ring_atoms = RING_ATOMS_PURINE;
    } else {
        ring_atoms = RING_ATOMS_PYRIMIDINE;
    }
    
    // Add ring atoms, optionally excluding C4 (legacy compatibility mode)
    for (const auto& atom_name : ring_atoms) {
        if (exclude_c4 && atom_name == " C4 ") {
            // Skip C4 atom in legacy mode
            continue;
        }
        atom_names.push_back(atom_name);
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

