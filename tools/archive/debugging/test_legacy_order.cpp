/**
 * @file test_legacy_order.cpp
 * @brief Test the legacy order function
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/structure.hpp>
#include <iostream>
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> [residue_idx]\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    int target_idx = -1;
    if (argc >= 3) {
        target_idx = std::stoi(argv[2]);
    }

    // Parse PDB - include all residues to match legacy (HETATMs, waters, etc.)
    PdbParser parser;
    parser.set_include_hetatm(true); // Include HETATM records (legacy includes all)
    parser.set_include_waters(true); // Include water molecules (legacy includes all)
    auto structure = parser.parse_file(pdb_file);

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Legacy Order Test\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "PDB: " << pdb_file << "\n\n";

    // Get residues in legacy order (using Structure's built-in method)
    auto residues = structure.residues_in_legacy_order();

    std::cout << "Total residues in legacy order: " << residues.size() << "\n\n";

    if (target_idx > 0) {
        std::cout << "Residue at legacy index " << target_idx << ":\n";
        const Residue* res = structure.get_residue_by_legacy_idx(target_idx);
        if (res) {
            std::cout << "  " << res->name() << " (chain " << res->chain_id() << ", seq " << res->seq_num() << ")\n";
        } else {
            std::cout << "  Not found!\n";
        }
        std::cout << "\n";
    }

    // Show first 10 and around target
    std::cout << "First 10 residues in legacy order:\n";
    for (size_t i = 0; i < std::min<size_t>(10, residues.size()); i++) {
        const Residue* res = residues[i];
        std::cout << "  " << std::setw(4) << (i + 1) << ". " << std::setw(3) << res->name() << " (chain "
                  << res->chain_id() << ", seq " << std::setw(4) << res->seq_num() << ")\n";
    }

    if (target_idx > 0 && target_idx <= static_cast<int>(residues.size())) {
        std::cout << "\nResidues around index " << target_idx << ":\n";
        int start = std::max(1, target_idx - 2);
        int end = std::min(static_cast<int>(residues.size()), target_idx + 2);
        for (int i = start; i <= end; i++) {
            const Residue* res = residues[i - 1];
            const char* marker = (i == target_idx) ? " <--" : "";
            std::cout << "  " << std::setw(4) << i << ". " << std::setw(3) << res->name() << " (chain "
                      << res->chain_id() << ", seq " << std::setw(4) << res->seq_num() << ")" << marker << "\n";
        }
    }

    return 0;
}
