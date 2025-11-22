/**
 * @file debug_frame_calculation.cpp
 * @brief Debugging tool to investigate frame calculation differences
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"
#include <filesystem>
#include <vector>
#include <set>
#include <map>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::test;

void print_atom_info(const Atom& atom) {
    std::cout << "    " << std::setw(6) << atom.name() 
              << " (" << std::fixed << std::setprecision(3)
              << std::setw(8) << atom.position().x() << ", "
              << std::setw(8) << atom.position().y() << ", "
              << std::setw(8) << atom.position().z() << ")";
}

void debug_residue(const Residue& residue, const nlohmann::json& legacy_record, 
                   BaseFrameCalculator& calculator) {
    std::cout << "\n=== Residue " << residue.seq_num() << " (" << residue.name() << ") ===" << std::endl;
    
    // Get residue type
    ResidueType residue_type = residue.residue_type();
    std::cout << "Residue Type: ";
    switch (residue_type) {
        case ResidueType::ADENINE: std::cout << "ADENINE (A)"; break;
        case ResidueType::CYTOSINE: std::cout << "CYTOSINE (C)"; break;
        case ResidueType::GUANINE: std::cout << "GUANINE (G)"; break;
        case ResidueType::THYMINE: std::cout << "THYMINE (T)"; break;
        case ResidueType::URACIL: std::cout << "URACIL (U)"; break;
        default: std::cout << "UNKNOWN"; break;
    }
    std::cout << std::endl;
    
    // Legacy info
    if (legacy_record.contains("base_type")) {
        std::cout << "Legacy base_type: " << legacy_record["base_type"].get<std::string>() << std::endl;
    }
    if (legacy_record.contains("num_matched_atoms")) {
        std::cout << "Legacy num_matched_atoms: " << legacy_record["num_matched_atoms"].get<size_t>() << std::endl;
    }
    if (legacy_record.contains("matched_atoms")) {
        auto legacy_atoms = legacy_record["matched_atoms"].get<std::vector<std::string>>();
        std::cout << "Legacy matched_atoms: [";
        for (const auto& atom : legacy_atoms) {
            std::cout << atom << " ";
        }
        std::cout << "]" << std::endl;
    }
    
    // Get expected ring atom names
    std::vector<std::string> expected_ring_atoms = RingAtomMatcher::get_ring_atom_names(residue_type, false);
    std::cout << "\nExpected ring atoms (" << expected_ring_atoms.size() << "): [";
    for (const auto& atom_name : expected_ring_atoms) {
        std::cout << atom_name << " ";
    }
    std::cout << "]" << std::endl;
    
    // Check which atoms exist in residue
    std::cout << "\nAtoms in residue:" << std::endl;
    std::set<std::string> residue_atom_names;
    for (const auto& atom : residue.atoms()) {
        residue_atom_names.insert(atom.name());
        print_atom_info(atom);
        std::cout << std::endl;
    }
    
    // Check which expected ring atoms are present
    std::cout << "\nRing atom availability:" << std::endl;
    std::vector<std::string> available_ring_atoms;
    for (const auto& atom_name : expected_ring_atoms) {
        bool found = residue_atom_names.find(atom_name) != residue_atom_names.end();
        std::cout << "  " << atom_name << ": " << (found ? "✓" : "✗") << std::endl;
        if (found) {
            available_ring_atoms.push_back(atom_name);
        }
    }
    
    // Load standard template
    StandardBaseTemplates templates("data/templates");
    try {
        Structure standard_template = templates.load_template(residue_type);
        
        std::cout << "\nStandard template atoms:" << std::endl;
        std::set<std::string> template_atom_names;
        for (const auto& chain : standard_template.chains()) {
            for (const auto& template_residue : chain.residues()) {
                for (const auto& atom : template_residue.atoms()) {
                    template_atom_names.insert(atom.name());
                    print_atom_info(atom);
                    std::cout << std::endl;
                }
            }
        }
        
        // Try matching
        MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template, false);
        
        std::cout << "\n=== Matching Results ===" << std::endl;
        std::cout << "Number of matched atoms: " << matched.num_matched << std::endl;
        std::cout << "Matched atom names: [";
        for (const auto& name : matched.atom_names) {
            std::cout << name << " ";
        }
        std::cout << "]" << std::endl;
        
        std::cout << "\nMatched atom pairs:" << std::endl;
        for (size_t i = 0; i < matched.num_matched; ++i) {
            std::cout << "  " << matched.atom_names[i] << ":" << std::endl;
            std::cout << "    Experimental: ";
            print_atom_info(matched.experimental[i]);
            std::cout << std::endl;
            std::cout << "    Standard:     ";
            print_atom_info(matched.standard[i]);
            std::cout << std::endl;
        }
        
        // Calculate frame
        if (matched.is_valid()) {
            FrameCalculationResult result = calculator.calculate_frame_const(residue);
            
            std::cout << "\n=== Frame Calculation Result ===" << std::endl;
            std::cout << "Valid: " << (result.is_valid ? "Yes" : "No") << std::endl;
            std::cout << "RMS fit: " << std::fixed << std::setprecision(6) << result.rms_fit << std::endl;
            std::cout << "Num matched: " << result.num_matched << std::endl;
            
            if (legacy_record.contains("rms_fit") && !legacy_record["rms_fit"].is_null()) {
                double legacy_rms = legacy_record["rms_fit"].get<double>();
                std::cout << "Legacy RMS: " << legacy_rms << std::endl;
                std::cout << "Difference: " << std::abs(result.rms_fit - legacy_rms) << std::endl;
            }
            
            std::cout << "\nRotation matrix:" << std::endl;
            for (size_t i = 0; i < 3; ++i) {
                std::cout << "  [";
                for (size_t j = 0; j < 3; ++j) {
                    std::cout << std::setw(10) << std::fixed << std::setprecision(6) 
                              << result.rotation_matrix.at(i, j);
                    if (j < 2) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
            
            std::cout << "\nTranslation:" << std::endl;
            std::cout << "  [" << std::fixed << std::setprecision(6)
                      << result.translation.x() << ", "
                      << result.translation.y() << ", "
                      << result.translation.z() << "]" << std::endl;
            
            if (legacy_record.contains("translation")) {
                auto legacy_trans = legacy_record["translation"];
                std::cout << "Legacy translation: [" 
                          << legacy_trans[0].get<double>() << ", "
                          << legacy_trans[1].get<double>() << ", "
                          << legacy_trans[2].get<double>() << "]" << std::endl;
            }
        } else {
            std::cout << "\nMatching failed - not enough atoms matched" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "\nError loading template: " << e.what() << std::endl;
    }
}

int main(int /* argc */, char* /* argv */[]) {
    // Discover PDB/JSON pairs
    auto pairs = test_data_discovery::discover_pairs();
    
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }
    
    // Initialize calculator
    BaseFrameCalculator calculator("data/templates");
    
    // Debug first 5 PDBs
    size_t num_to_debug = std::min(pairs.size(), size_t(5));
    
    for (size_t i = 0; i < num_to_debug; ++i) {
        const auto& pair = pairs[i];
        
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "PDB: " << pair.pdb_name << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        
        // Load PDB file
        PdbParser parser;
        Structure structure = parser.parse_file(pair.pdb_file);
        
        // Load legacy JSON
        std::ifstream json_file(pair.json_file);
        nlohmann::json legacy_json;
        json_file >> legacy_json;
        
        // Find base_frame_calc records
        std::vector<nlohmann::json> base_frame_records;
        if (legacy_json.contains("calculations")) {
            for (const auto& calc : legacy_json["calculations"]) {
                if (calc.contains("type") && calc["type"] == "base_frame_calc") {
                    base_frame_records.push_back(calc);
                }
            }
        }
        
        // Get pdb_atoms to build ordered residue list (all residues, not just nucleotides)
        std::vector<std::tuple<char, int, std::string>> ordered_residues;
        if (legacy_json.contains("calculations")) {
            for (const auto& calc : legacy_json["calculations"]) {
                if (calc.contains("type") && calc["type"] == "pdb_atoms" &&
                    calc.contains("atoms") && calc["atoms"].is_array()) {
                    std::set<std::tuple<char, int, std::string>> seen;
                    for (const auto& atom : calc["atoms"]) {
                        std::string chain_str = atom.value("chain_id", "");
                        char chain_id = chain_str.empty() ? '\0' : chain_str[0];
                        int seq_num = atom.value("residue_seq", 0);
                        std::string res_name = atom.value("residue_name", "");
                        auto key = std::make_tuple(chain_id, seq_num, res_name);
                        if (seen.find(key) == seen.end()) {
                            ordered_residues.push_back(key);
                            seen.insert(key);
                        }
                    }
                    break;
                }
            }
        }
        
        // Debug first few base_frame_calc records
        size_t num_to_debug_per_pdb = std::min(base_frame_records.size(), size_t(3));
        
        for (size_t i = 0; i < num_to_debug_per_pdb; ++i) {
            const auto& legacy_record = base_frame_records[i];
            if (!legacy_record.contains("residue_idx")) continue;
            
            size_t legacy_idx = legacy_record["residue_idx"].get<size_t>();
            
            // Convert 1-based to 0-based index
            if (legacy_idx == 0 || legacy_idx > ordered_residues.size()) continue;
            
            auto [legacy_chain, legacy_seq, legacy_name] = ordered_residues[legacy_idx - 1];
            
            // Find matching residue in structure
            for (const auto& chain : structure.chains()) {
                if (chain.chain_id() != legacy_chain) continue;
                
                for (const auto& residue : chain.residues()) {
                    if (residue.seq_num() == legacy_seq) {
                        debug_residue(residue, legacy_record, calculator);
                        break;
                    }
                }
            }
        }
    }
    
    return 0;
}

