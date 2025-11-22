/**
 * @file debug_atom_matching.cpp
 * @brief Detailed debugging of atom matching logic vs legacy
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"
#include <cmath>
#include <algorithm>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::test;
using namespace x3dna::geometry;

// RA_LIST as in legacy: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
static const std::vector<std::string> RA_LIST = {
    " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
};

std::vector<std::tuple<char, int, std::string>> build_ordered_residue_list(const nlohmann::json& legacy_json) {
    std::vector<std::tuple<char, int, std::string>> ordered_residues;
    if (legacy_json.contains("calculations")) {
        for (const auto& calc : legacy_json["calculations"]) {
            if (calc.contains("type") && calc["type"] == "pdb_atoms" &&
                calc.contains("atoms") && calc["atoms"].is_array()) {
                std::vector<std::tuple<char, int, std::string>> seen;
                for (const auto& atom : calc["atoms"]) {
                    std::string chain_str = atom.value("chain_id", "");
                    char chain_id = chain_str.empty() ? '\0' : chain_str[0];
                    int seq_num = atom.value("residue_seq", 0);
                    std::string res_name = atom.value("residue_name", "");
                    auto key = std::make_tuple(chain_id, seq_num, res_name);
                    bool found = false;
                    for (const auto& s : seen) {
                        if (s == key) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        ordered_residues.push_back(key);
                        seen.push_back(key);
                    }
                }
                break;
            }
        }
    }
    return ordered_residues;
}

std::vector<nlohmann::json> find_records_by_type(const nlohmann::json& json, const std::string& type) {
    std::vector<nlohmann::json> results;
    if (json.contains("calculations")) {
        for (const auto& calc : json["calculations"]) {
            if (calc.contains("type") && calc["type"] == type) {
                results.push_back(calc);
            }
        }
    }
    return results;
}

std::optional<const Residue*> find_residue_by_legacy_idx(
    const Structure& structure, 
    size_t legacy_residue_idx,
    const std::vector<std::tuple<char, int, std::string>>& ordered_residues) {
    
    if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) {
        return std::nullopt;
    }
    
    auto [legacy_chain, legacy_seq, legacy_name] = ordered_residues[legacy_residue_idx - 1];
    
    for (const auto& chain : structure.chains()) {
        if (chain.chain_id() != legacy_chain) continue;
        for (const auto& residue : chain.residues()) {
            if (residue.seq_num() == legacy_seq) {
                return &residue;
            }
        }
    }
    
    return std::nullopt;
}

void debug_residue_matching(const Residue& residue, const nlohmann::json& frame_calc_record,
                            StandardBaseTemplates& templates) {
    std::cout << "\n=== Detailed Atom Matching Debug ===" << std::endl;
    std::cout << "Residue: " << residue.name() << " seq " << residue.seq_num() << std::endl;
    
    ResidueType residue_type = residue.residue_type();
    bool is_purine = (residue_type == ResidueType::ADENINE || residue_type == ResidueType::GUANINE);
    int RingAtom_num = is_purine ? 9 : 6;
    
    std::cout << "Residue type: " << (is_purine ? "Purine" : "Pyrimidine") << std::endl;
    std::cout << "RingAtom_num: " << RingAtom_num << std::endl;
    
    // Load standard template
    Structure standard_template = templates.load_template(residue_type);
    
    // Legacy matching order (as in RA_LIST)
    std::cout << "\n--- Legacy Matching Order (RA_LIST) ---" << std::endl;
    std::vector<std::string> our_matched;
    std::vector<Vector3D> our_exp_coords;
    std::vector<Vector3D> our_std_coords;
    
    int nmatch = 0;
    for (int j = 0; j < RingAtom_num; ++j) {
        const std::string& atom_name = RA_LIST[j];
        std::cout << "Checking [" << j << "] " << atom_name << ": ";
        
        // Find in experimental residue
        std::optional<Atom> exp_atom;
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                exp_atom = atom;
                break;
            }
        }
        
        // Find in standard template
        std::optional<Atom> std_atom;
        for (const auto& chain : standard_template.chains()) {
            for (const auto& res : chain.residues()) {
                for (const auto& atom : res.atoms()) {
                    if (atom.name() == atom_name) {
                        std_atom = atom;
                        break;
                    }
                }
                if (std_atom.has_value()) break;
            }
            if (std_atom.has_value()) break;
        }
        
        if (exp_atom.has_value() && std_atom.has_value()) {
            nmatch++;
            our_matched.push_back(atom_name);
            our_exp_coords.push_back(exp_atom.value().position());
            our_std_coords.push_back(std_atom.value().position());
            std::cout << "✓ MATCHED (nmatch=" << nmatch << ")" << std::endl;
            std::cout << "  Exp: (" << std::fixed << std::setprecision(3)
                      << exp_atom.value().position().x() << ", "
                      << exp_atom.value().position().y() << ", "
                      << exp_atom.value().position().z() << ")" << std::endl;
            std::cout << "  Std: (" << std_atom.value().position().x() << ", "
                      << std_atom.value().position().y() << ", "
                      << std_atom.value().position().z() << ")" << std::endl;
        } else {
            std::cout << "✗ NOT FOUND" << std::endl;
            if (!exp_atom.has_value()) std::cout << "  (missing in experimental)" << std::endl;
            if (!std_atom.has_value()) std::cout << "  (missing in standard)" << std::endl;
        }
    }
    
    std::cout << "\nTotal matched: " << nmatch << " atoms" << std::endl;
    std::cout << "Matched atoms: [";
    for (const auto& name : our_matched) {
        std::cout << name << " ";
    }
    std::cout << "]" << std::endl;
    
    // Compare with legacy frame_calc record if available
    if (frame_calc_record.contains("matched_coordinates")) {
        std::cout << "\n--- Legacy matched_coordinates ---" << std::endl;
        auto legacy_coords = frame_calc_record["matched_coordinates"];
        std::cout << "Legacy num_matched: " << legacy_coords.size() << std::endl;
        
        for (size_t i = 0; i < legacy_coords.size() && i < our_matched.size(); ++i) {
            auto legacy = legacy_coords[i];
            auto std_xyz = legacy.value("std_xyz", std::vector<double>{});
            auto exp_xyz = legacy.value("exp_xyz", std::vector<double>{});
            
            if (std_xyz.size() >= 3 && exp_xyz.size() >= 3) {
                std::cout << "  [" << i << "] " << our_matched[i] << ":" << std::endl;
                std::cout << "    Our Exp: (" << std::fixed << std::setprecision(3)
                          << our_exp_coords[i].x() << ", "
                          << our_exp_coords[i].y() << ", "
                          << our_exp_coords[i].z() << ")" << std::endl;
                std::cout << "    Legacy Exp: (" << exp_xyz[0] << ", " << exp_xyz[1] << ", " << exp_xyz[2] << ")" << std::endl;
                double diff_x = std::abs(our_exp_coords[i].x() - exp_xyz[0]);
                double diff_y = std::abs(our_exp_coords[i].y() - exp_xyz[1]);
                double diff_z = std::abs(our_exp_coords[i].z() - exp_xyz[2]);
                if (diff_x > 0.01 || diff_y > 0.01 || diff_z > 0.01) {
                    std::cout << "    ⚠ DIFFERENCE! max diff: " 
                              << std::max({diff_x, diff_y, diff_z}) << std::endl;
                }
                std::cout << "    Our Std: (" << our_std_coords[i].x() << ", "
                          << our_std_coords[i].y() << ", "
                          << our_std_coords[i].z() << ")" << std::endl;
                std::cout << "    Legacy Std: (" << std_xyz[0] << ", " << std_xyz[1] << ", " << std_xyz[2] << ")" << std::endl;
                diff_x = std::abs(our_std_coords[i].x() - std_xyz[0]);
                diff_y = std::abs(our_std_coords[i].y() - std_xyz[1]);
                diff_z = std::abs(our_std_coords[i].z() - std_xyz[2]);
                if (diff_x > 0.01 || diff_y > 0.01 || diff_z > 0.01) {
                    std::cout << "    ⚠ DIFFERENCE! max diff: " 
                              << std::max({diff_x, diff_y, diff_z}) << std::endl;
                }
            }
        }
    }
    
    // Now check what our RingAtomMatcher does
    std::cout << "\n--- Our RingAtomMatcher Result ---" << std::endl;
    MatchedAtoms our_result = RingAtomMatcher::match(residue, standard_template, false);
    std::cout << "Our num_matched: " << our_result.num_matched << std::endl;
    std::cout << "Our matched atoms: [";
    for (const auto& name : our_result.atom_names) {
        std::cout << name << " ";
    }
    std::cout << "]" << std::endl;
    
    // Compare
    if (our_result.atom_names != our_matched) {
        std::cout << "\n⚠⚠⚠ ATOM MATCHING ORDER/DIFFERENCE DETECTED! ⚠⚠⚠" << std::endl;
        std::cout << "Legacy order: [";
        for (const auto& name : our_matched) std::cout << name << " ";
        std::cout << "]" << std::endl;
        std::cout << "Our order: [";
        for (const auto& name : our_result.atom_names) std::cout << name << " ";
        std::cout << "]" << std::endl;
    } else {
        std::cout << "\n✓ Atom matching order matches!" << std::endl;
    }
}

int main() {
    auto pairs = test_data_discovery::discover_pairs();
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }
    
    BaseFrameCalculator calculator("data/templates");
    StandardBaseTemplates templates("data/templates");
    
    // Debug first PDB, first 3 residues that have frame_calc records
    const auto& pair = pairs[0];
    
    std::cout << "Debugging PDB: " << pair.pdb_name << std::endl;
    
    // Load PDB
    PdbParser parser;
    Structure structure = parser.parse_file(pair.pdb_file);
    
    // Load legacy JSON
    std::ifstream json_file(pair.json_file);
    nlohmann::json legacy_json;
    json_file >> legacy_json;
    
    // Get frame_calc records
    auto frame_calc_records = find_records_by_type(legacy_json, "frame_calc");
    auto ordered_residues = build_ordered_residue_list(legacy_json);
    
    // Debug first 3 residues with frame_calc records
    size_t debugged = 0;
    for (const auto& frame_calc_record : frame_calc_records) {
        if (debugged >= 3) break;
        
        if (!frame_calc_record.contains("residue_idx")) continue;
        
        size_t legacy_residue_idx = frame_calc_record["residue_idx"].get<size_t>();
        auto residue_opt = find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);
        
        if (!residue_opt.has_value()) continue;
        
        const Residue* residue_ptr = residue_opt.value();
        debug_residue_matching(*residue_ptr, frame_calc_record, templates);
        debugged++;
    }
    
    return 0;
}

