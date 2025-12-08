/**
 * @file check_residue_indices.cpp
 * @brief Check for duplicate legacy_residue_idx values in parsed structure
 */

#include <iostream>
#include <map>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>

using namespace x3dna;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> [target_idx]\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int target_idx = (argc > 2) ? std::stoi(argv[2]) : 1102;

    std::cout << "Checking residue indices in: " << pdb_file << "\n";
    std::cout << "Target index: " << target_idx << "\n\n";

    // Parse PDB
    io::PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);

    // Build map and check for duplicates
    std::map<int, std::vector<core::Residue*>> residues_by_legacy_idx;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residues_by_legacy_idx[legacy_idx].push_back(&residue);
                }
            }
        }
    }

    // Check for duplicates
    std::cout << "Checking for duplicate legacy_residue_idx values...\n";
    bool found_duplicates = false;
    for (const auto& [idx, residues] : residues_by_legacy_idx) {
        if (residues.size() > 1) {
            found_duplicates = true;
            std::cout << "\n⚠️  Duplicate legacy_residue_idx " << idx << " found in " << residues.size()
                      << " residues:\n";
            for (auto* res : residues) {
                std::cout << "  - " << res->name() << " Chain " << res->chain_id() << " Seq " << res->seq_num() << "\n";
            }
        }
    }

    if (!found_duplicates) {
        std::cout << "✓ No duplicates found\n";
    }

    // Check target index
    std::cout << "\nResidues with legacy_residue_idx = " << target_idx << ":\n";
    auto it = residues_by_legacy_idx.find(target_idx);
    if (it != residues_by_legacy_idx.end()) {
        for (auto* res : it->second) {
            std::cout << "  - " << res->name() << " Chain " << res->chain_id() << " Seq " << res->seq_num() << "\n";
        }
    } else {
        std::cout << "  Not found\n";
    }

    return 0;
}
