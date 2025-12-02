/**
 * @file debug_frame_on_residue.cpp
 * @brief Debug tool to check frames on residue objects after calculation
 * 
 * This tool helps identify frame retrieval bugs by:
 * 1. Calculating frames on a PDB structure
 * 2. Retrieving residues by legacy index (like base_pair_finder does)
 * 3. Checking what frames those residue objects actually have
 * 4. Comparing to expected frames from frame_calc JSON
 * 
 * Usage:
 *   ./build/tools/debug_frame_on_residue <PDB_FILE> <legacy_idx1> [legacy_idx2]
 *   ./build/tools/debug_frame_on_residue data/pdb/9CF3.pdb 25 27
 */

#include <iostream>
#include <filesystem>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/core/structure.hpp>
#include <iomanip>
#include <cmath>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <PDB_FILE> <legacy_idx1> [legacy_idx2]\n";
        return 1;
    }
    
    std::filesystem::path pdb_file = argv[1];
    int legacy_idx1 = std::stoi(argv[2]);
    int legacy_idx2 = (argc > 3) ? std::stoi(argv[3]) : 0;
    
    // Parse PDB
    x3dna::io::PdbParser parser;
    x3dna::core::Structure structure = parser.parse_file(pdb_file);
    
    std::cout << "Parsed PDB: " << pdb_file << "\n";
    std::cout << "Looking for residues with legacy_idx: " << legacy_idx1;
    if (legacy_idx2 > 0) {
        std::cout << " and " << legacy_idx2;
    }
    std::cout << "\n\n";
    
    // Calculate frames (this should set frames on residue objects)
    x3dna::algorithms::BaseFrameCalculator calculator("data/templates");
    
    // Detect RNA
    bool is_rna = false;
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == " O2'") {
                    is_rna = true;
                    break;
                }
            }
            if (is_rna) break;
        }
        if (is_rna) break;
    }
    calculator.set_is_rna(is_rna);
    
    std::cout << "Calculating frames...\n";
    calculator.calculate_all_frames(structure);
    std::cout << "Frames calculated.\n\n";
    
    // Build residue_by_legacy_idx mapping (like base_pair_finder does)
    std::map<int, const x3dna::core::Residue*> residue_by_legacy_idx;
    
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residue_by_legacy_idx[legacy_idx] = &residue;
                }
            }
        }
    }
    
    std::cout << "Built residue_by_legacy_idx mapping with " 
              << residue_by_legacy_idx.size() << " residues\n\n";
    
    // Check frames on residues
    auto check_residue = [&](int legacy_idx) {
        auto it = residue_by_legacy_idx.find(legacy_idx);
        if (it == residue_by_legacy_idx.end()) {
            std::cout << "Residue " << legacy_idx << ": NOT FOUND in mapping\n";
            return;
        }
        
        const x3dna::core::Residue* residue = it->second;
        
        std::cout << "Residue " << legacy_idx << ":\n";
        std::cout << "  Name: " << residue->name() 
                  << ", Chain: " << residue->chain_id()
                  << ", Seq: " << residue->seq_num() << "\n";
        
        auto frame_opt = residue->reference_frame();
        if (!frame_opt.has_value()) {
            std::cout << "  ❌ NO FRAME SET!\n";
            return;
        }
        
        const auto& frame = frame_opt.value();
        auto origin = frame.origin();
        
        std::cout << "  ✅ Frame exists\n";
        std::cout << "  Origin: (" 
                  << std::fixed << std::setprecision(6)
                  << origin.x() << ", "
                  << origin.y() << ", "
                  << origin.z() << ")\n";
        
        // Also check atoms to verify this is the right residue
        if (!residue->atoms().empty()) {
            auto first_atom = residue->atoms()[0];
            std::cout << "  First atom: " << first_atom.name()
                      << " (legacy_idx=" << first_atom.legacy_residue_idx() << ")\n";
        }
        std::cout << "\n";
    };
    
    check_residue(legacy_idx1);
    if (legacy_idx2 > 0) {
        check_residue(legacy_idx2);
        
        // Calculate dorg if both residues found
        auto it1 = residue_by_legacy_idx.find(legacy_idx1);
        auto it2 = residue_by_legacy_idx.find(legacy_idx2);
        
        if (it1 != residue_by_legacy_idx.end() && it2 != residue_by_legacy_idx.end()) {
            const auto* res1 = it1->second;
            const auto* res2 = it2->second;
            
            auto frame1_opt = res1->reference_frame();
            auto frame2_opt = res2->reference_frame();
            
            if (frame1_opt.has_value() && frame2_opt.has_value()) {
                const auto& frame1 = frame1_opt.value();
                const auto& frame2 = frame2_opt.value();
                
                auto dorg_vec = frame1.origin() - frame2.origin();
                double dorg = std::sqrt(
                    dorg_vec.x() * dorg_vec.x() +
                    dorg_vec.y() * dorg_vec.y() +
                    dorg_vec.z() * dorg_vec.z()
                );
                
                std::cout << "Calculated dorg from residue frames: " 
                          << std::fixed << std::setprecision(6) << dorg << " Å\n";
                std::cout << "Expected from frame_calc JSON: 4.874563 Å\n";
                
                if (std::abs(dorg - 4.874563) < 0.01) {
                    std::cout << "✅ dorg matches!\n";
                } else {
                    std::cout << "❌ dorg DOES NOT MATCH! Difference: " 
                              << std::abs(dorg - 4.874563) << " Å\n";
                    std::cout << "   This suggests frames on residue objects are WRONG!\n";
                }
            }
        }
    }
    
    return 0;
}

