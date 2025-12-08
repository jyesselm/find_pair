/**
 * @file compare_atom_selection.cpp
 * @brief Compare which atoms legacy and modern check for H-bond detection
 *
 * This tool helps debug why modern finds different H-bonds by comparing:
 * - Which atoms legacy checks (seidx range)
 * - Which atoms modern checks (all atoms in residue)
 * - Atom-by-atom comparison
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <set>
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

void print_atom_list(const std::string& label, const std::vector<std::string>& atoms) {
    std::cout << "\n" << label << " (" << atoms.size() << " atoms):\n";
    for (size_t i = 0; i < atoms.size(); i++) {
        std::cout << "  " << std::setw(3) << (i + 1) << ". " << atoms[i] << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3G8T.pdb 946 947\n";
        return 1;
    }

    std::filesystem::path pdb_file = argv[1];
    int residue1_idx = std::stoi(argv[2]);
    int residue2_idx = std::stoi(argv[3]);

    // Parse PDB
    PdbParser parser;
    Structure structure = parser.parse_file(pdb_file);

    // Find residues (legacy uses 1-based indexing)
    Residue* res1 = nullptr;
    Residue* res2 = nullptr;
    size_t legacy_idx = 1;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (legacy_idx == static_cast<size_t>(residue1_idx)) {
                res1 = &residue;
            }
            if (legacy_idx == static_cast<size_t>(residue2_idx)) {
                res2 = &residue;
            }
            legacy_idx++;
        }
    }

    if (!res1 || !res2) {
        std::cerr << "Error: Could not find residues " << residue1_idx << " and " << residue2_idx << "\n";
        return 1;
    }

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Atom Selection Comparison\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Pair: (" << residue1_idx << ", " << residue2_idx << ")\n";
    std::cout << "Residue 1: " << res1->name() << " (chain " << res1->chain_id() << ", seq " << res1->seq_num()
              << ")\n";
    std::cout << "Residue 2: " << res2->name() << " (chain " << res2->chain_id() << ", seq " << res2->seq_num()
              << ")\n";

    // Get modern atom lists
    std::vector<std::string> modern_atoms1;
    std::vector<std::string> modern_atoms2;

    for (const auto& atom : res1->atoms()) {
        modern_atoms1.push_back(atom.name());
    }

    for (const auto& atom : res2->atoms()) {
        modern_atoms2.push_back(atom.name());
    }

    print_atom_list("Modern Residue 1 Atoms", modern_atoms1);
    print_atom_list("Modern Residue 2 Atoms", modern_atoms2);

    // Legacy seidx information
    // Note: We can't directly get seidx from modern code, but we can infer
    // Legacy typically includes: base atoms + backbone atoms (P, O1P, O2P, O5', C5', C4', O4', C3',
    // O3', C2', O2', C1') For modified nucleotides, it might include additional atoms

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Legacy seidx Information\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Legacy uses seidx[i][1] to seidx[i][2] to define atom range.\n";
    std::cout << "This is typically all atoms in the residue from PDB file.\n";
    std::cout << "However, legacy might filter or exclude certain atoms.\n";
    std::cout << "\nTo get exact seidx range, need to:\n";
    std::cout << "  1. Add debug output to legacy code\n";
    std::cout << "  2. Or parse legacy's atom selection logic\n";

    // Check for potential H-bond atoms
    std::set<std::string> hbond_atoms = {"N1", "N2", "N3",  "N4",  "N6",  "N7",  "N9",  "O2",
                                         "O4", "O6", "O4'", "O2'", "O3'", "O5'", "O1P", "O2P"};

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Potential H-bond Atoms\n";
    std::cout << std::string(60, '=') << "\n";

    std::vector<std::string> res1_hbond_atoms;
    std::vector<std::string> res2_hbond_atoms;

    for (const auto& atom_name : modern_atoms1) {
        // Check if atom name contains H-bond atoms
        for (const auto& hb_atom : hbond_atoms) {
            if (atom_name.find(hb_atom) != std::string::npos) {
                res1_hbond_atoms.push_back(atom_name);
                break;
            }
        }
    }

    for (const auto& atom_name : modern_atoms2) {
        for (const auto& hb_atom : hbond_atoms) {
            if (atom_name.find(hb_atom) != std::string::npos) {
                res2_hbond_atoms.push_back(atom_name);
                break;
            }
        }
    }

    print_atom_list("Residue 1 Potential H-bond Atoms", res1_hbond_atoms);
    print_atom_list("Residue 2 Potential H-bond Atoms", res2_hbond_atoms);

    // Check for specific atoms that differ
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Key Observations\n";
    std::cout << std::string(60, '=') << "\n";

    bool has_o4 = false, has_o4prime = false, has_n6 = false;
    for (const auto& atom : modern_atoms1) {
        if (atom.find(" O4 ") != std::string::npos)
            has_o4 = true;
        if (atom.find("O4'") != std::string::npos)
            has_o4prime = true;
        if (atom.find(" N6 ") != std::string::npos)
            has_n6 = true;
    }
    for (const auto& atom : modern_atoms2) {
        if (atom.find(" O4 ") != std::string::npos)
            has_o4 = true;
        if (atom.find("O4'") != std::string::npos)
            has_o4prime = true;
        if (atom.find(" N6 ") != std::string::npos)
            has_n6 = true;
    }

    std::cout << "Has O4 (base oxygen): " << (has_o4 ? "YES" : "NO") << "\n";
    std::cout << "Has O4' (backbone oxygen): " << (has_o4prime ? "YES" : "NO") << "\n";
    std::cout << "Has N6: " << (has_n6 ? "YES" : "NO") << "\n";
    std::cout << "\nNote: Legacy finds N3->O4, but modern finds N3->O4'\n";
    std::cout << "      This suggests legacy checks O4 but modern checks O4'\n";

    return 0;
}
