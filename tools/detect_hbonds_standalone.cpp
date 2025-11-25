/**
 * @file detect_hbonds_standalone.cpp
 * @brief Standalone H-bond detection tool - completely separate from pair validation
 * 
 * This tool detects H-bonds between two residues WITHOUT any pair validation.
 * It uses the same H-bond detection logic but runs it independently, showing:
 * - All H-bonds found (before any validation filtering)
 * - H-bond distances
 * - H-bond types (standard '-' or non-standard '*')
 * - Which H-bonds are "good" for quality adjustment (type='-' and dist in [2.5, 3.5])
 * 
 * This is separate from pair validation which has additional checks and filtering.
 */

#include <iostream>
#include <fstream>
#include <filesystem>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <iomanip>
#include <vector>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

struct HBondInfo {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;  // '-' for standard, '*' for non-standard
    bool is_good;  // type='-' AND distance in [2.5, 3.5]
};

void print_hbond_info(const std::vector<HBondInfo>& hbonds, const std::string& label) {
    std::cout << "\n" << label << "\n";
    std::cout << "========================================\n";
    
    if (hbonds.empty()) {
        std::cout << "  (no H-bonds found)\n";
        return;
    }
    
    int good_count = 0;
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < hbonds.size(); i++) {
        const auto& hb = hbonds[i];
        std::cout << "  " << std::setw(3) << (i+1) << ". ";
        std::cout << std::setw(6) << hb.donor_atom << " -> " << std::setw(6) << hb.acceptor_atom;
        std::cout << "  dist=" << std::setw(10) << hb.distance;
        std::cout << "  type=" << hb.type;
        std::cout << "  good=" << (hb.is_good ? "YES" : "NO ");
        std::cout << "\n";
        
        if (hb.is_good) {
            good_count++;
        }
    }
    
    std::cout << "\n  Summary:\n";
    std::cout << "    Total H-bonds: " << hbonds.size() << "\n";
    std::cout << "    Good H-bonds (type='-' and dist in [2.5, 3.5]): " << good_count << "\n";
    std::cout << "    adjust_pairQuality: ";
    if (good_count >= 2) {
        std::cout << "-3.0 (2+ good H-bonds)\n";
    } else {
        std::cout << "-" << good_count << ".0 (" << good_count << " good H-bond";
        if (good_count != 1) std::cout << "s";
        std::cout << ")\n";
    }
}

std::vector<HBondInfo> detect_hbonds_standalone(const Residue& res1, const Residue& res2) {
    std::vector<HBondInfo> hbonds;
    
    // Use default validation parameters for H-bond detection
    ValidationParameters params = ValidationParameters::defaults();
    double hb_lower = params.hb_lower;
    double hb_dist1 = params.hb_dist1;
    double hb_dist2 = 4.5; // Default value (matches legacy)
    
    // Use HydrogenBondFinder directly - this is separate from pair validation
    DetailedHBondResult detailed = HydrogenBondFinder::find_hydrogen_bonds_detailed(
        res1, res2, hb_lower, hb_dist1, hb_dist2);
    
    // Get ALL H-bonds found (after validation, but this is H-bond validation, not pair validation)
    for (const auto& hbond_result : detailed.after_validation) {
        HBondInfo info;
        info.donor_atom = hbond_result.donor_atom;
        info.acceptor_atom = hbond_result.acceptor_atom;
        info.distance = std::abs(hbond_result.distance); // Use absolute distance
        info.type = hbond_result.type;
        
        // Good H-bond: type='-' AND distance in [2.5, 3.5]
        info.is_good = (info.type == '-' && info.distance >= 2.5 && info.distance <= 3.5);
        
        hbonds.push_back(info);
    }
    
    return hbonds;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <chain1> <seq1> <chain2> <seq2> [insertion1] [insertion2]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3G8T.pdb A 92 A 160\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb A 75 A 78\n";
        std::cerr << "\n";
        std::cerr << "This tool detects H-bonds between two residues INDEPENDENTLY.\n";
        std::cerr << "It does NOT use pair validation - only pure H-bond detection.\n";
        return 1;
    }
    
    std::filesystem::path pdb_file = argv[1];
    char chain1 = argv[2][0];
    int seq1 = std::stoi(argv[3]);
    char chain2 = argv[4][0];
    int seq2 = std::stoi(argv[5]);
    char insertion1 = (argc >= 7) ? argv[6][0] : ' ';
    char insertion2 = (argc >= 8) ? argv[7][0] : ' ';
    
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
        
        // Calculate frames for all residues
        BaseFrameCalculator calculator("data/templates");
        calculator.calculate_all_frames(structure);
        
        // Find the two residues
        const Residue* res1_ptr = nullptr;
        const Residue* res2_ptr = nullptr;
        
        for (const auto& chain : structure.chains()) {
            if (chain.chain_id() == chain1) {
                for (const auto& residue : chain.residues()) {
                    if (residue.seq_num() == seq1 && residue.insertion() == insertion1) {
                        res1_ptr = &residue;
                        break;
                    }
                }
            }
            if (chain.chain_id() == chain2) {
                for (const auto& residue : chain.residues()) {
                    if (residue.seq_num() == seq2 && residue.insertion() == insertion2) {
                        res2_ptr = &residue;
                        break;
                    }
                }
            }
            if (res1_ptr && res2_ptr) break;
        }
        
        if (!res1_ptr) {
            std::cerr << "Error: Residue not found: " << chain1 << ":" << seq1;
            if (insertion1 != ' ') std::cerr << insertion1;
            std::cerr << "\n";
            return 1;
        }
        
        if (!res2_ptr) {
            std::cerr << "Error: Residue not found: " << chain2 << ":" << seq2;
            if (insertion2 != ' ') std::cerr << insertion2;
            std::cerr << "\n";
            return 1;
        }
        
        // Print residue info
        std::cout << "\nResidue 1: " << res1_ptr->name() << " " << chain1 << ":" << seq1;
        if (insertion1 != ' ') std::cout << insertion1;
        std::cout << " (one_letter=" << res1_ptr->one_letter_code() << ")\n";
        
        std::cout << "Residue 2: " << res2_ptr->name() << " " << chain2 << ":" << seq2;
        if (insertion2 != ' ') std::cout << insertion2;
        std::cout << " (one_letter=" << res2_ptr->one_letter_code() << ")\n";
        
        // Check if residues have frames
        if (!res1_ptr->reference_frame().has_value()) {
            std::cerr << "Warning: Residue 1 does not have a reference frame\n";
        }
        if (!res2_ptr->reference_frame().has_value()) {
            std::cerr << "Warning: Residue 2 does not have a reference frame\n";
        }
        
        // Detect H-bonds independently (no pair validation)
        std::vector<HBondInfo> hbonds = detect_hbonds_standalone(*res1_ptr, *res2_ptr);
        
        // Print results
        print_hbond_info(hbonds, "Standalone H-bond Detection");
        
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

