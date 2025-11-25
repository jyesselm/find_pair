/**
 * @file compare_hbond_detection.cpp
 * @brief Tool to compare modern and legacy H-bond detection for specific pairs
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/hydrogen_bond_finder.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/json_reader.hpp>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <iomanip>

using namespace x3dna;
using namespace x3dna::io;
using namespace x3dna::algorithms;

using HydrogenBondResult = HydrogenBondResult;
using DetailedHBondResult = DetailedHBondResult;

void print_hbond_result(const HydrogenBondResult& hb) {
    std::cout << "    " << hb.donor_atom << " -> " << hb.acceptor_atom 
              << ": " << std::fixed << std::setprecision(3) << hb.distance
              << " (type=" << hb.type << ")\n";
}

void compare_pair_hbonds(const std::string& pdb_file, int legacy_idx1, int legacy_idx2) {
    std::cout << "\n=== Comparing H-bond detection for pair (" 
              << legacy_idx1 << ", " << legacy_idx2 << ") ===\n";
    
    // Parse PDB
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    auto structure = parser.parse_file(pdb_file);
    
    // Find residues by legacy index
    const core::Residue* res1 = nullptr;
    const core::Residue* res2 = nullptr;
    
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx == legacy_idx1) {
                    res1 = &residue;
                }
                if (legacy_idx == legacy_idx2) {
                    res2 = &residue;
                }
            }
        }
    }
    
    if (!res1 || !res2) {
        std::cerr << "ERROR: Could not find residues " << legacy_idx1 
                  << " and/or " << legacy_idx2 << "\n";
        return;
    }
    
    std::cout << "Residue 1: " << res1->name() << " " << res1->chain_id() 
              << ":" << res1->seq_num() << " (type=" << res1->one_letter_code() << ")\n";
    std::cout << "Residue 2: " << res2->name() << " " << res2->chain_id() 
              << ":" << res2->seq_num() << " (type=" << res2->one_letter_code() << ")\n";
    
    // Calculate frames if needed
    BaseFrameCalculator calculator("data/templates");
    if (!res1->reference_frame().has_value()) {
        calculator.calculate_frame(const_cast<core::Residue&>(*res1));
    }
    if (!res2->reference_frame().has_value()) {
        calculator.calculate_frame(const_cast<core::Residue&>(*res2));
    }
    
    // Find H-bonds using modern code
    double hb_lower = 2.0;  // Default values - should match legacy
    double hb_dist1 = 4.0;
    
    auto detailed = HydrogenBondFinder::find_hydrogen_bonds_detailed(
        *res1, *res2, hb_lower, hb_dist1
    );
    
    std::cout << "\nModern H-bond detection:\n";
    std::cout << "  Initial H-bonds found: " << detailed.initial_hbonds.size() << "\n";
    for (const auto& hb : detailed.initial_hbonds) {
        print_hbond_result(hb);
    }
    
    std::cout << "\n  After conflict resolution: " << detailed.after_conflict_resolution.size() << "\n";
    for (const auto& hb : detailed.after_conflict_resolution) {
        print_hbond_result(hb);
    }
    
    std::cout << "\n  After validation: " << detailed.after_validation.size() << "\n";
    for (const auto& hb : detailed.after_validation) {
        print_hbond_result(hb);
    }
    
    std::cout << "\n  Final H-bonds: " << detailed.final_hbonds.size() << "\n";
    for (const auto& hb : detailed.final_hbonds) {
        print_hbond_result(hb);
    }
    
    std::cout << "\n  Good H-bonds (type='-' and 2.5-3.5): " << detailed.num_good_hb << "\n";
    
    // Load legacy H-bond data from JSON if available
    std::filesystem::path pdb_path(pdb_file);
    std::string pdb_id = pdb_path.stem().string();
    std::filesystem::path legacy_json = "data/json_legacy/" + pdb_id + "_hbond_list.json";
    
    if (std::filesystem::exists(legacy_json)) {
        std::cout << "\nLegacy H-bond data from JSON:\n";
        std::ifstream file(legacy_json);
        nlohmann::json data;
        file >> data;
        
        for (const auto& record : data) {
            if (record.contains("base_i") && record.contains("base_j")) {
                int i = record["base_i"];
                int j = record["base_j"];
                if ((i == legacy_idx1 && j == legacy_idx2) ||
                    (i == legacy_idx2 && j == legacy_idx1)) {
                    if (record.contains("hbonds")) {
                        std::cout << "  Found " << record["hbonds"].size() << " H-bonds:\n";
                        for (const auto& hb : record["hbonds"]) {
                            std::cout << "    " << hb.value("donor_atom", "") 
                                      << " -> " << hb.value("acceptor_atom", "")
                                      << ": " << hb.value("distance", 0.0)
                                      << " (type=" << hb.value("type", "?") << ")\n";
                        }
                    }
                    break;
                }
            }
        }
    } else {
        std::cout << "\nLegacy JSON not found: " << legacy_json << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <legacy_idx1> <legacy_idx2>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/3G8T.pdb 92 160\n";
        return 1;
    }
    
    std::string pdb_file = argv[1];
    int legacy_idx1 = std::stoi(argv[2]);
    int legacy_idx2 = std::stoi(argv[3]);
    
    compare_pair_hbonds(pdb_file, legacy_idx1, legacy_idx2);
    
    return 0;
}

