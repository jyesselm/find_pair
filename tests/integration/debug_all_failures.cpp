/**
 * @file debug_all_failures.cpp
 * @brief Detailed debugging tool to record all frame calculation failures
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <nlohmann/json.hpp>
#include "test_data_discovery.hpp"
#include <cmath>
#include <sstream>
#include <filesystem>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::test;
using namespace x3dna::geometry;

struct FailureDetail {
    std::string pdb_name;
    size_t legacy_residue_idx;
    char chain_id;
    int seq_num;
    std::string residue_name;
    std::string base_type;
    std::string failure_reason;

    // Our values
    size_t our_num_matched;
    double our_rms;
    std::vector<double> our_rotation;
    std::vector<double> our_translation;

    // Legacy values
    size_t legacy_num_matched;
    double legacy_rms;
    std::vector<double> legacy_rotation;
    std::vector<double> legacy_translation;

    // Differences
    double max_rot_diff;
    double max_trans_diff;
    double rms_diff;

    // Atom matching info
    std::vector<std::string> our_atoms;
    std::vector<std::string> legacy_atoms;
    bool atoms_differ;

    // Coordinate info (first matched atom)
    std::vector<double> first_exp_coord;
    std::vector<double> first_std_coord;
};

std::vector<std::tuple<char, int, std::string>>
build_ordered_residue_list(const nlohmann::json& legacy_json) {
    std::vector<std::tuple<char, int, std::string>> ordered_residues;
    if (legacy_json.contains("calculations")) {
        for (const auto& calc : legacy_json["calculations"]) {
            if (calc.contains("type") && calc["type"] == "pdb_atoms" && calc.contains("atoms") &&
                calc["atoms"].is_array()) {
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

std::vector<nlohmann::json> find_records_by_type(const nlohmann::json& json,
                                                 const std::string& type) {
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
    const Structure& structure, size_t legacy_residue_idx,
    const std::vector<std::tuple<char, int, std::string>>& ordered_residues) {

    if (legacy_residue_idx == 0 || legacy_residue_idx > ordered_residues.size()) {
        return std::nullopt;
    }

    auto [legacy_chain, legacy_seq, legacy_name] = ordered_residues[legacy_residue_idx - 1];

    for (const auto& chain : structure.chains()) {
        if (chain.chain_id() != legacy_chain)
            continue;
        for (const auto& residue : chain.residues()) {
            if (residue.seq_num() == legacy_seq) {
                return &residue;
            }
        }
    }

    return std::nullopt;
}

double max_rotation_diff(const Matrix3D& m1, const nlohmann::json& json_m2) {
    if (!json_m2.is_array() || json_m2.size() != 3)
        return -1.0;
    double max_diff = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        if (!json_m2[i].is_array() || json_m2[i].size() != 3)
            return -1.0;
        for (size_t j = 0; j < 3; ++j) {
            double v1 = m1.at(i, j);
            double v2 = json_m2[i][j].get<double>();
            double diff = std::abs(v1 - v2);
            if (diff > max_diff)
                max_diff = diff;
        }
    }
    return max_diff;
}

double max_translation_diff(const Vector3D& v1, const nlohmann::json& json_v2) {
    if (!json_v2.is_array() || json_v2.size() != 3)
        return -1.0;
    double v2_x = json_v2[0].get<double>();
    double v2_y = json_v2[1].get<double>();
    double v2_z = json_v2[2].get<double>();
    double diff_x = std::abs(v1.x() - v2_x);
    double diff_y = std::abs(v1.y() - v2_y);
    double diff_z = std::abs(v1.z() - v2_z);
    return std::max({diff_x, diff_y, diff_z});
}

void extract_rotation_matrix(const Matrix3D& m, std::vector<double>& out) {
    out.clear();
    out.reserve(9);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            out.push_back(m.at(i, j));
        }
    }
}

void extract_translation_vector(const Vector3D& v, std::vector<double>& out) {
    out.clear();
    out.reserve(3);
    out.push_back(v.x());
    out.push_back(v.y());
    out.push_back(v.z());
}

int main() {
    auto pairs = test_data_discovery::discover_pairs();
    if (pairs.empty()) {
        std::cerr << "No PDB/JSON pairs found" << std::endl;
        return 1;
    }

    BaseFrameCalculator calculator("data/templates");
    StandardBaseTemplates templates("data/templates");

    std::vector<FailureDetail> all_failures;
    std::vector<FailureDetail> real_failures; // Prioritize actual numerical failures
    std::map<std::string, size_t> failure_reasons;

    std::cout << "Analyzing all " << pairs.size() << " PDB files..." << std::endl;
    std::cout << "Progress will be shown every 100 PDBs..." << std::endl;

    for (size_t pdb_idx = 0; pdb_idx < pairs.size(); ++pdb_idx) {
        if (pdb_idx % 100 == 0) {
            std::cout << "Processing PDB " << pdb_idx << "/" << pairs.size() << " ("
                      << all_failures.size() << " failures so far)..." << std::endl;
        }

        const auto& pair = pairs[pdb_idx];

        try {
            // Load PDB
            PdbParser parser;
            Structure structure = parser.parse_file(pair.pdb_file);

            // Load legacy JSON
            std::ifstream json_file(pair.json_file);
            if (!json_file.is_open())
                continue;

            nlohmann::json legacy_json;
            json_file >> legacy_json;

            // Get records
            auto ls_records = find_records_by_type(legacy_json, "ls_fitting");
            auto base_frame_records = find_records_by_type(legacy_json, "base_frame_calc");
            auto frame_calc_records = find_records_by_type(legacy_json, "frame_calc");
            auto ordered_residues = build_ordered_residue_list(legacy_json);

            // Calculate frames
            calculator.calculate_all_frames(structure);

            // Analyze each residue
            for (const auto& ls_record : ls_records) {
                if (!ls_record.contains("residue_idx"))
                    continue;

                size_t legacy_residue_idx = ls_record["residue_idx"].get<size_t>();

                auto [legacy_chain, legacy_seq, legacy_name] =
                    ordered_residues[legacy_residue_idx - 1];

                FailureDetail detail;
                detail.pdb_name = pair.pdb_name;
                detail.legacy_residue_idx = legacy_residue_idx;
                detail.chain_id = legacy_chain;
                detail.seq_num = legacy_seq;
                detail.residue_name = legacy_name;

                auto residue_opt =
                    find_residue_by_legacy_idx(structure, legacy_residue_idx, ordered_residues);

                if (!residue_opt.has_value()) {
                    detail.failure_reason = "RESIDUE_NOT_FOUND";
                    all_failures.push_back(detail);
                    failure_reasons["RESIDUE_NOT_FOUND"]++;
                    continue;
                }

                const Residue* residue_ptr = residue_opt.value();

                if (!residue_ptr->reference_frame().has_value()) {
                    // Try to calculate to see why
                    FrameCalculationResult test_result =
                        calculator.calculate_frame_const(*residue_ptr);
                    if (!test_result.is_valid) {
                        detail.failure_reason = "NO_FRAME_INVALID";
                        detail.our_num_matched = test_result.num_matched;
                    } else {
                        detail.failure_reason = "NO_FRAME_NOT_STORED";
                    }
                    all_failures.push_back(detail);
                    failure_reasons[detail.failure_reason]++;
                    continue;
                }

                FrameCalculationResult result = calculator.calculate_frame_const(*residue_ptr);

                if (!result.is_valid) {
                    detail.failure_reason = "INVALID_CALCULATION";
                    detail.our_num_matched = result.num_matched;
                    all_failures.push_back(detail);
                    failure_reasons["INVALID_CALCULATION"]++;
                    continue;
                }

                // Get base_frame_calc record
                const nlohmann::json* base_frame_record = nullptr;
                for (const auto& bf : base_frame_records) {
                    if (bf.contains("residue_idx") &&
                        bf["residue_idx"].get<size_t>() == legacy_residue_idx) {
                        base_frame_record = &bf;
                        break;
                    }
                }

                if (base_frame_record) {
                    detail.base_type = base_frame_record->value("base_type", "");
                    if (base_frame_record->contains("matched_atoms")) {
                        detail.legacy_atoms =
                            base_frame_record->value("matched_atoms", std::vector<std::string>());
                    }
                    detail.legacy_num_matched = base_frame_record->value("num_matched_atoms", 0UL);
                }

                detail.our_num_matched = result.num_matched;
                detail.our_atoms = result.matched_atoms;
                detail.our_rms = result.rms_fit;

                extract_rotation_matrix(result.rotation_matrix, detail.our_rotation);
                extract_translation_vector(result.translation, detail.our_translation);

                // Compare with legacy
                bool failed = false;
                std::stringstream failure_details;

                if (ls_record.contains("rotation_matrix")) {
                    detail.legacy_rotation.clear();
                    detail.legacy_rotation.reserve(9);
                    auto leg_rot = ls_record["rotation_matrix"];
                    for (size_t i = 0; i < 3; ++i) {
                        for (size_t j = 0; j < 3; ++j) {
                            detail.legacy_rotation.push_back(leg_rot[i][j].get<double>());
                        }
                    }
                    detail.max_rot_diff =
                        max_rotation_diff(result.rotation_matrix, ls_record["rotation_matrix"]);
                    if (detail.max_rot_diff > 0.05) {
                        failed = true;
                        failure_details << "ROT_";
                    }
                }

                if (ls_record.contains("translation")) {
                    detail.legacy_translation.clear();
                    detail.legacy_translation.reserve(3);
                    auto leg_trans = ls_record["translation"];
                    for (size_t i = 0; i < 3; ++i) {
                        detail.legacy_translation.push_back(leg_trans[i].get<double>());
                    }
                    detail.max_trans_diff =
                        max_translation_diff(result.translation, ls_record["translation"]);
                    if (detail.max_trans_diff > 0.05) {
                        failed = true;
                        failure_details << "TRANS_";
                    }
                }

                if (ls_record.contains("rms_fit") && !ls_record["rms_fit"].is_null()) {
                    detail.legacy_rms = ls_record["rms_fit"].get<double>();
                    detail.rms_diff = std::abs(result.rms_fit - detail.legacy_rms);
                    if (detail.rms_diff > 0.005) {
                        failed = true;
                        failure_details << "RMS_";
                    }
                }

                if (ls_record.contains("num_points")) {
                    detail.legacy_num_matched = ls_record["num_points"].get<size_t>();
                    if (detail.our_num_matched != detail.legacy_num_matched) {
                        failed = true;
                        failure_details << "NUM_MATCHED_";
                    }
                }

                // Get frame_calc record for coordinates
                const nlohmann::json* frame_calc = nullptr;
                for (const auto& fc : frame_calc_records) {
                    if (fc.contains("residue_idx") &&
                        fc["residue_idx"].get<size_t>() == legacy_residue_idx) {
                        frame_calc = &fc;
                        break;
                    }
                }

                if (frame_calc && frame_calc->contains("matched_coordinates")) {
                    auto coords = (*frame_calc)["matched_coordinates"];
                    if (coords.size() > 0) {
                        auto first_coord = coords[0];
                        auto std_xyz = first_coord.value("std_xyz", std::vector<double>{});
                        auto exp_xyz = first_coord.value("exp_xyz", std::vector<double>{});
                        if (std_xyz.size() >= 3)
                            detail.first_std_coord = std_xyz;
                        if (exp_xyz.size() >= 3)
                            detail.first_exp_coord = exp_xyz;
                    }
                }

                // Check if atoms differ
                detail.atoms_differ = (detail.our_atoms != detail.legacy_atoms);
                if (detail.atoms_differ) {
                    failed = true;
                    failure_details << "ATOMS_";
                }

                if (failed) {
                    detail.failure_reason = failure_details.str();
                    // Remove trailing underscore
                    if (!detail.failure_reason.empty() && detail.failure_reason.back() == '_') {
                        detail.failure_reason.pop_back();
                    }

                    // Check if this is a real numerical failure or just atom ordering
                    bool is_real_failure = false;
                    if (detail.failure_reason != "ATOMS" || detail.max_rot_diff > 1e-6 ||
                        detail.max_trans_diff > 1e-6 || detail.rms_diff > 1e-6 ||
                        detail.our_num_matched != detail.legacy_num_matched) {
                        is_real_failure = true;
                    }

                    if (is_real_failure) {
                        real_failures.push_back(detail);
                    } else {
                        all_failures.push_back(detail);
                    }
                    failure_reasons[detail.failure_reason]++;
                }
            }
        } catch (const std::exception& e) {
            // Skip PDBs that fail to load
            continue;
        }
    }

    // Write detailed failure report
    std::ofstream report("docs/frame_calculation_failures.json");
    nlohmann::json failure_json;

    failure_json["summary"] = nlohmann::json::object();
    failure_json["summary"]["total_failures"] = all_failures.size();
    failure_json["summary"]["failure_breakdown"] = nlohmann::json::object();
    for (const auto& [reason, count] : failure_reasons) {
        failure_json["summary"]["failure_breakdown"][reason] = count;
    }

    failure_json["failures"] = nlohmann::json::array();

    // Prioritize real failures first, then add atom ordering differences
    // Write up to 2000 total (prioritize real failures)
    size_t num_real_to_write = std::min(real_failures.size(), size_t(1000));
    size_t num_atom_ordering_to_write = std::min(all_failures.size(), size_t(1000));

    // First add real failures
    for (size_t i = 0; i < num_real_to_write; ++i) {
        const auto& f = real_failures[i];
        nlohmann::json failure;
        failure["pdb_name"] = f.pdb_name;
        failure["legacy_residue_idx"] = f.legacy_residue_idx;
        failure["chain_id"] = std::string(1, f.chain_id);
        failure["seq_num"] = f.seq_num;
        failure["residue_name"] = f.residue_name;
        failure["base_type"] = f.base_type;
        failure["failure_reason"] = f.failure_reason;

        failure["our"] = nlohmann::json::object();
        failure["our"]["num_matched"] = f.our_num_matched;
        failure["our"]["rms"] = f.our_rms;
        failure["our"]["rotation"] = f.our_rotation;
        failure["our"]["translation"] = f.our_translation;
        failure["our"]["matched_atoms"] = f.our_atoms;

        failure["legacy"] = nlohmann::json::object();
        failure["legacy"]["num_matched"] = f.legacy_num_matched;
        failure["legacy"]["rms"] = f.legacy_rms;
        failure["legacy"]["rotation"] = f.legacy_rotation;
        failure["legacy"]["translation"] = f.legacy_translation;
        failure["legacy"]["matched_atoms"] = f.legacy_atoms;

        failure["differences"] = nlohmann::json::object();
        failure["differences"]["max_rot_diff"] = f.max_rot_diff;
        failure["differences"]["max_trans_diff"] = f.max_trans_diff;
        failure["differences"]["rms_diff"] = f.rms_diff;
        failure["differences"]["atoms_differ"] = f.atoms_differ;

        if (!f.first_exp_coord.empty()) {
            failure["first_coordinates"] = nlohmann::json::object();
            failure["first_coordinates"]["experimental"] = f.first_exp_coord;
            failure["first_coordinates"]["standard"] = f.first_std_coord;
        }

        failure_json["failures"].push_back(failure);
    }

    // Then add atom ordering differences (up to remaining space)
    for (size_t i = 0; i < num_atom_ordering_to_write && failure_json["failures"].size() < 2000;
         ++i) {
        const auto& f = all_failures[i];
        nlohmann::json failure;
        failure["pdb_name"] = f.pdb_name;
        failure["legacy_residue_idx"] = f.legacy_residue_idx;
        failure["chain_id"] = std::string(1, f.chain_id);
        failure["seq_num"] = f.seq_num;
        failure["residue_name"] = f.residue_name;
        failure["base_type"] = f.base_type;
        failure["failure_reason"] = f.failure_reason;

        failure["our"] = nlohmann::json::object();
        failure["our"]["num_matched"] = f.our_num_matched;
        failure["our"]["rms"] = f.our_rms;
        failure["our"]["rotation"] = f.our_rotation;
        failure["our"]["translation"] = f.our_translation;
        failure["our"]["matched_atoms"] = f.our_atoms;

        failure["legacy"] = nlohmann::json::object();
        failure["legacy"]["num_matched"] = f.legacy_num_matched;
        failure["legacy"]["rms"] = f.legacy_rms;
        failure["legacy"]["rotation"] = f.legacy_rotation;
        failure["legacy"]["translation"] = f.legacy_translation;
        failure["legacy"]["matched_atoms"] = f.legacy_atoms;

        failure["differences"] = nlohmann::json::object();
        failure["differences"]["max_rot_diff"] = f.max_rot_diff;
        failure["differences"]["max_trans_diff"] = f.max_trans_diff;
        failure["differences"]["rms_diff"] = f.rms_diff;
        failure["differences"]["atoms_differ"] = f.atoms_differ;

        if (!f.first_exp_coord.empty()) {
            failure["first_coordinates"] = nlohmann::json::object();
            failure["first_coordinates"]["experimental"] = f.first_exp_coord;
            failure["first_coordinates"]["standard"] = f.first_std_coord;
        }

        failure_json["failures"].push_back(failure);
    }

    report << failure_json.dump(2);
    report.close();

    // Print summary
    std::cout << "\n=== Failure Analysis Complete ===" << std::endl;
    std::cout << "Total failures: " << (all_failures.size() + real_failures.size()) << std::endl;
    std::cout << "  Real numerical failures: " << real_failures.size() << std::endl;
    std::cout << "  Atom ordering only: " << all_failures.size() << std::endl;
    std::cout << "\nFailure breakdown:" << std::endl;
    for (const auto& [reason, count] : failure_reasons) {
        std::cout << "  " << reason << ": " << count << std::endl;
    }

    std::cout << "\nDetailed failure report written to: docs/frame_calculation_failures.json"
              << std::endl;
    std::cout << "(" << num_real_to_write << " real failures + "
              << failure_json["failures"].size() - num_real_to_write
              << " atom ordering differences included in JSON)" << std::endl;

    // Print examples of different failure types (prioritize real failures)
    std::cout << "\n=== Sample Failures by Type ===" << std::endl;
    std::map<std::string, size_t> shown;

    // First show real failures
    for (const auto& f : real_failures) {
        if (shown[f.failure_reason] < 3) {
            std::cout << "\n"
                      << f.failure_reason << " - " << f.pdb_name << " residue_idx "
                      << f.legacy_residue_idx << " (" << f.chain_id << ":" << f.seq_num << " "
                      << f.residue_name << ")" << std::endl;

            if (f.failure_reason.find("COMPARISON") != std::string::npos ||
                f.failure_reason.find("ROT") != std::string::npos ||
                f.failure_reason.find("TRANS") != std::string::npos ||
                f.failure_reason.find("RMS") != std::string::npos) {
                std::cout << "  Our RMS: " << std::fixed << std::setprecision(6) << f.our_rms
                          << ", Legacy RMS: " << f.legacy_rms << std::endl;
                std::cout << "  Max rot diff: " << f.max_rot_diff << std::endl;
                std::cout << "  Max trans diff: " << f.max_trans_diff << std::endl;
                if (f.atoms_differ) {
                    std::cout << "  Atoms differ! Our: " << f.our_num_matched
                              << ", Legacy: " << f.legacy_num_matched << std::endl;
                }
            }

            shown[f.failure_reason]++;
            if (shown.size() >= 10)
                break; // Show up to 10 different failure types
        }
    }

    return 0;
}
