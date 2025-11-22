/**
 * @file compare_ls_fitting.cpp
 * @brief Directly compare our ls fitting with legacy results
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <x3dna/geometry/least_squares_fitter.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;
using namespace x3dna::test;

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

int main() {
    auto pairs = test_data_discovery::discover_pairs();
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }
    
    StandardBaseTemplates templates("data/templates");
    
    // Get first PDB with failures
    const auto& pair = pairs[0];
    
    // Load PDB
    PdbParser parser;
    Structure structure = parser.parse_file(pair.pdb_file);
    
    // Load legacy JSON
    std::ifstream json_file(pair.json_file);
    nlohmann::json legacy_json;
    json_file >> legacy_json;
    
    // Get ls_fitting and frame_calc records
    auto ls_records = find_records_by_type(legacy_json, "ls_fitting");
    auto frame_calc_records = find_records_by_type(legacy_json, "frame_calc");
    auto ordered_residues = build_ordered_residue_list(legacy_json);
    
    // Find a failing residue (one where num_points differs or RMS differs significantly)
    for (const auto& ls_record : ls_records) {
        if (!ls_record.contains("residue_idx")) continue;
        
        size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();
        auto residue_opt = find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);
        if (!residue_opt.has_value()) continue;
        
        const Residue* residue_ptr = residue_opt.value();
        
        // Get frame_calc record for coordinates
        const nlohmann::json* frame_calc = nullptr;
        for (const auto& fc : frame_calc_records) {
            if (fc.contains("residue_idx") && fc["residue_idx"].get<size_t>() == legacy_residue_idx) {
                frame_calc = &fc;
                break;
            }
        }
        
        if (!frame_calc || !frame_calc->contains("matched_coordinates")) continue;
        
        // Extract legacy coordinates
        auto legacy_coords = (*frame_calc)["matched_coordinates"];
        if (legacy_coords.size() < 3) continue;
        
        std::cout << "\n=== Comparing LS Fitting for Residue " << legacy_residue_idx << " ===" << std::endl;
        
        // Build coordinate arrays from legacy
        std::vector<Vector3D> legacy_standard;
        std::vector<Vector3D> legacy_experimental;
        
        for (const auto& coord : legacy_coords) {
            auto std_xyz = coord.value("std_xyz", std::vector<double>{});
            auto exp_xyz = coord.value("exp_xyz", std::vector<double>{});
            if (std_xyz.size() >= 3 && exp_xyz.size() >= 3) {
                legacy_standard.push_back(Vector3D(std_xyz[0], std_xyz[1], std_xyz[2]));
                legacy_experimental.push_back(Vector3D(exp_xyz[0], exp_xyz[1], exp_xyz[2]));
            }
        }
        
        std::cout << "Legacy num_points: " << legacy_standard.size() << std::endl;
        
        // Match atoms using our code
        Structure standard_template = templates.load_template(residue_ptr->residue_type());
        MatchedAtoms matched = RingAtomMatcher::match(*residue_ptr, standard_template, false);
        
        if (!matched.is_valid()) {
            std::cout << "Our matching failed!" << std::endl;
            continue;
        }
        
        std::cout << "Our num_matched: " << matched.num_matched << std::endl;
        
        // Extract our coordinates
        std::vector<Vector3D> our_standard;
        std::vector<Vector3D> our_experimental;
        
        for (size_t i = 0; i < matched.num_matched; ++i) {
            our_standard.push_back(matched.standard[i].position());
            our_experimental.push_back(matched.experimental[i].position());
        }
        
        // Compare coordinates
        std::cout << "\nCoordinate Comparison:" << std::endl;
        bool coords_match = true;
        for (size_t i = 0; i < std::min(legacy_standard.size(), our_standard.size()); ++i) {
            double diff_std = (legacy_standard[i] - our_standard[i]).length();
            double diff_exp = (legacy_experimental[i] - our_experimental[i]).length();
            if (diff_std > 0.01 || diff_exp > 0.01) {
                std::cout << "  [" << i << "] DIFFERENCE! std_diff=" << diff_std 
                          << ", exp_diff=" << diff_exp << std::endl;
                coords_match = false;
            }
        }
        
        if (coords_match) {
            std::cout << "  ✓ Coordinates match!" << std::endl;
        }
        
        // Perform least squares fitting with our code
        LeastSquaresFitter fitter;
        auto our_result = fitter.fit(our_standard, our_experimental);
        
        // Get legacy results
        auto legacy_rot = ls_record["rotation_matrix"];
        auto legacy_trans = ls_record["translation"];
        double legacy_rms = ls_record.value("rms_fit", 0.0);
        
        std::cout << "\nLS Fitting Comparison:" << std::endl;
        std::cout << "  Legacy RMS: " << std::fixed << std::setprecision(6) << legacy_rms << std::endl;
        std::cout << "  Our RMS:    " << our_result.rms << std::endl;
        std::cout << "  Difference: " << std::abs(our_result.rms - legacy_rms) << std::endl;
        
        std::cout << "\nRotation Matrix:" << std::endl;
        double max_rot_diff = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            std::cout << "  Row " << i << ": ";
            for (size_t j = 0; j < 3; ++j) {
                double our_val = our_result.rotation.at(i, j);
                double leg_val = legacy_rot[i][j].get<double>();
                double diff = std::abs(our_val - leg_val);
                if (diff > max_rot_diff) max_rot_diff = diff;
                std::cout << std::fixed << std::setprecision(6) << std::setw(10) << diff << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "  Max difference: " << max_rot_diff << std::endl;
        
        std::cout << "\nTranslation:" << std::endl;
        double max_trans_diff = 0.0;
        std::vector<double> our_trans = {our_result.translation.x(), 
                                         our_result.translation.y(), 
                                         our_result.translation.z()};
        const char* labels[] = {"X", "Y", "Z"};
        for (size_t i = 0; i < 3; ++i) {
            double our_val = our_trans[i];
            double leg_val = legacy_trans[i].get<double>();
            double diff = std::abs(our_val - leg_val);
            if (diff > max_trans_diff) max_trans_diff = diff;
            std::cout << "  " << labels[i] << " Our: " << std::fixed << std::setprecision(6) 
                      << std::setw(12) << our_val
                      << " Legacy: " << std::setw(12) << leg_val
                      << " Diff: " << std::setw(12) << diff << std::endl;
        }
        std::cout << "  Max difference: " << max_trans_diff << std::endl;
        
        // Only show first failing residue, then break
        if (max_rot_diff > 0.05 || max_trans_diff > 0.05 || std::abs(our_result.rms - legacy_rms) > 0.005) {
            std::cout << "\n*** FOUND FAILING RESIDUE ***" << std::endl;
            break;
        }
    }
    
    return 0;
}

