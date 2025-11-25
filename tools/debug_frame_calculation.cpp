/**
 * @file debug_frame_calculation.cpp
 * @brief Debug tool for frame calculation on specific residues
 * 
 * This tool helps debug why frame calculation fails for specific residues.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <iomanip>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

void print_residue_info(const Residue& residue) {
    std::cout << "\n========================================\n";
    std::cout << "Residue Information\n";
    std::cout << "========================================\n";
    std::cout << "  Name: " << residue.name() << "\n";
    std::cout << "  Chain: " << residue.chain_id() << "\n";
    std::cout << "  Sequence: " << residue.seq_num() << "\n";
    std::cout << "  Insertion: " << (residue.insertion() == ' ' ? "(none)" : std::string(1, residue.insertion())) << "\n";
    std::cout << "  Residue Type: " << static_cast<int>(residue.residue_type()) << "\n";
    std::cout << "  One Letter: " << residue.one_letter_code() << "\n";
    std::cout << "  Is Nucleotide: " << (residue.is_nucleotide() ? "yes" : "no") << "\n";
    std::cout << "  Number of Atoms: " << residue.num_atoms() << "\n";
    
    std::cout << "\n  Atoms:\n";
    for (const auto& atom : residue.atoms()) {
        std::cout << "    " << atom.name() << " at (" 
                  << std::fixed << std::setprecision(3)
                  << atom.position().x() << ", "
                  << atom.position().y() << ", "
                  << atom.position().z() << ")\n";
    }
}

void check_ring_atoms(const Residue& residue) {
    std::cout << "\n========================================\n";
    std::cout << "Ring Atom Detection\n";
    std::cout << "========================================\n";
    
    static const std::vector<std::string> common_ring_atoms = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "};
    static const std::vector<std::string> purine_ring_atoms = {" N7 ", " C8 ", " N9 "};
    
    int ring_count = 0;
    std::cout << "  Common ring atoms:\n";
    for (const auto& atom_name : common_ring_atoms) {
        auto atom = residue.find_atom(atom_name);
        if (atom.has_value()) {
            std::cout << "    ✓ " << atom_name << " found\n";
            ring_count++;
        } else {
            std::cout << "    ✗ " << atom_name << " missing\n";
        }
    }
    
    std::cout << "\n  Purine-specific atoms:\n";
    bool has_purine = false;
    for (const auto& atom_name : purine_ring_atoms) {
        auto atom = residue.find_atom(atom_name);
        if (atom.has_value()) {
            std::cout << "    ✓ " << atom_name << " found\n";
            has_purine = true;
        } else {
            std::cout << "    ✗ " << atom_name << " missing\n";
        }
    }
    
    std::cout << "\n  Summary:\n";
    std::cout << "    Common ring atoms found: " << ring_count << " / " << common_ring_atoms.size() << "\n";
    std::cout << "    Has purine atoms: " << (has_purine ? "yes" : "no") << "\n";
    std::cout << "    Sufficient for nucleotide: " << (ring_count >= 3 ? "yes" : "no") << "\n";
}

void debug_frame_calculation(const Residue& residue, BaseFrameCalculator& calculator) {
    std::cout << "\n========================================\n";
    std::cout << "Frame Calculation Debug\n";
    std::cout << "========================================\n";
    
    // Make a mutable copy for calculation
    Residue mutable_residue = residue;
    
    FrameCalculationResult result = calculator.calculate_frame(mutable_residue);
    
    std::cout << "  Is Valid: " << (result.is_valid ? "yes" : "no") << "\n";
    std::cout << "  Template File: " << (result.template_file.empty() ? "(none)" : result.template_file) << "\n";
    std::cout << "  RMS Fit: " << std::fixed << std::setprecision(6) << result.rms_fit << "\n";
    std::cout << "  Matched Atoms: " << result.num_matched << "\n";
    
    if (!result.matched_atoms.empty()) {
        std::cout << "  Matched Atom Names:\n";
        for (const auto& atom_name : result.matched_atoms) {
            std::cout << "    - " << atom_name << "\n";
        }
    }
    
    if (result.is_valid) {
        std::cout << "\n  Frame Origin: (" 
                  << std::fixed << std::setprecision(6)
                  << result.frame.origin().x() << ", "
                  << result.frame.origin().y() << ", "
                  << result.frame.origin().z() << ")\n";
        
        auto rot = result.frame.rotation().as_array();
        std::cout << "  Rotation Matrix:\n";
        for (int i = 0; i < 3; i++) {
            std::cout << "    [";
            for (int j = 0; j < 3; j++) {
                std::cout << std::setw(10) << std::fixed << std::setprecision(6) << rot[i*3 + j];
                if (j < 2) std::cout << ", ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "\n  ❌ Frame calculation FAILED\n";
        
        if (result.template_file.empty()) {
            std::cout << "  Reason: Template file not found or couldn't be loaded\n";
        } else if (result.num_matched < 3) {
            std::cout << "  Reason: Insufficient atom matching (" << result.num_matched << " < 3 required)\n";
        } else {
            std::cout << "  Reason: Unknown (template loaded, atoms matched, but still failed)\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <chain_id> <seq_num> [insertion]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3KNC.pdb B 1\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/5UJ2.pdb T 2\n";
        return 1;
    }
    
    std::filesystem::path pdb_file = argv[1];
    char chain_id = argv[2][0];
    int seq_num = std::stoi(argv[3]);
    char insertion = (argc >= 5) ? argv[4][0] : ' ';
    
    if (!std::filesystem::exists(pdb_file)) {
        std::cerr << "Error: PDB file not found: " << pdb_file << "\n";
        return 1;
    }
    
    try {
        // Parse PDB
        std::cout << "Parsing PDB file: " << pdb_file << "\n";
        PdbParser parser;
        parser.set_include_hetatm(true);
        parser.set_include_waters(true);
        
        Structure structure = parser.parse_file(pdb_file);
        
        // Find the residue
        const Residue* target_residue = nullptr;
        for (const auto& chain : structure.chains()) {
            if (chain.chain_id() == chain_id) {
                for (const auto& residue : chain.residues()) {
                    if (residue.seq_num() == seq_num && residue.insertion() == insertion) {
                        target_residue = &residue;
                        break;
                    }
                }
            }
            if (target_residue) break;
        }
        
        if (!target_residue) {
            std::cerr << "Error: Residue not found: " << chain_id << ":" << seq_num;
            if (insertion != ' ') std::cerr << insertion;
            std::cerr << "\n";
            return 1;
        }
        
        // Print residue info
        print_residue_info(*target_residue);
        
        // Check ring atoms
        check_ring_atoms(*target_residue);
        
        // Debug frame calculation
        BaseFrameCalculator calculator("data/templates");
        debug_frame_calculation(*target_residue, calculator);
        
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

