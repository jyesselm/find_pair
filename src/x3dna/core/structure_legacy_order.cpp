/**
 * @file structure_legacy_order.cpp
 * @brief Implementation of legacy order utilities
 */

#include <x3dna/core/structure_legacy_order.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/atom.hpp>
#include <algorithm>

namespace x3dna {
namespace core {

std::vector<const Residue*> get_residues_in_legacy_order(const Structure& structure) {
    // Collect all atoms with their line numbers and residue info
    struct AtomWithResidue {
        const Atom* atom;
        const Residue* residue;
        size_t line_number;
    };
    
    std::vector<AtomWithResidue> atoms_with_residues;
    
    // Collect all atoms from all residues
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                atoms_with_residues.push_back({&atom, &residue, atom.line_number()});
            }
        }
    }
    
    // Sort by line number (PDB file order)
    std::sort(atoms_with_residues.begin(), atoms_with_residues.end(),
              [](const AtomWithResidue& a, const AtomWithResidue& b) {
                  return a.line_number < b.line_number;
              });
    
    // Group by (ResName, ChainID, ResSeq, insertion) and collect unique residues in order
    std::set<std::tuple<std::string, char, int, char>> seen_residues;
    std::vector<const Residue*> residues_in_order;
    
    for (const auto& awr : atoms_with_residues) {
        // Legacy groups by (ResName, ChainID, ResSeq, insertion)
        auto legacy_key = std::make_tuple(
            awr.atom->residue_name(),
            awr.atom->chain_id(),
            awr.atom->residue_seq(),
            awr.atom->insertion()
        );
        
        // Only add residue when we first see it
        if (seen_residues.find(legacy_key) == seen_residues.end()) {
            seen_residues.insert(legacy_key);
            residues_in_order.push_back(awr.residue);
        }
    }
    
    return residues_in_order;
}

const Residue* get_residue_by_legacy_idx(const Structure& structure, int legacy_idx) {
    if (legacy_idx < 1) {
        return nullptr;
    }
    
    auto residues = get_residues_in_legacy_order(structure);
    size_t idx = static_cast<size_t>(legacy_idx - 1); // Convert to 0-based
    
    if (idx >= residues.size()) {
        return nullptr;
    }
    
    return residues[idx];
}

int get_legacy_idx_for_residue(const Structure& structure, const Residue* residue) {
    if (!residue) {
        return 0;
    }
    
    auto residues = get_residues_in_legacy_order(structure);
    
    for (size_t i = 0; i < residues.size(); i++) {
        if (residues[i] == residue) {
            return static_cast<int>(i + 1); // Convert to 1-based
        }
    }
    
    return 0;
}

} // namespace core
} // namespace x3dna

