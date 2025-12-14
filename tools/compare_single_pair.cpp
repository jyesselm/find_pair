/**
 * @file compare_single_pair.cpp
 * @brief Tool to compare pair validation between legacy and modern code for a specific pair
 *
 * Usage:
 *   ./compare_single_pair <pdb_file> <base_i> <base_j> [--json-dir <dir>]
 *
 * Example:
 *   ./compare_single_pair data/pdb/1EHZ.pdb 1 72 --json-dir data/json_legacy
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <optional>

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/quality_score_calculator.hpp>
#include <x3dna/debug/pair_validation_debugger.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::debug;

void print_usage(const char* program) {
    std::cerr << "Usage: " << program << " <pdb_file> <base_i> <base_j> [options]\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  --json-dir <dir>   Legacy JSON directory (default: data/json_legacy)\n";
    std::cerr << "  --verbose          Show detailed output\n";
    std::cerr << "  --help             Show this help\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << program << " data/pdb/1EHZ.pdb 1 72\n";
}

struct LegacyPairData {
    bool found = false;
    int is_valid = 0;
    int bp_type_id = -1;
    double dir_x = 0.0, dir_y = 0.0, dir_z = 0.0;
    double dorg = 0.0, d_v = 0.0, plane_angle = 0.0, dNN = 0.0;
    double quality_score = 0.0;
    bool distance_check = false;
    bool d_v_check = false;
    bool plane_angle_check = false;
    bool dNN_check = false;
};

LegacyPairData load_legacy_pair_validation(const std::string& json_dir,
                                            const std::string& pdb_id,
                                            int base_i, int base_j) {
    LegacyPairData data;

    std::string path = json_dir + "/pair_validation/" + pdb_id + ".json";
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open " << path << "\n";
        return data;
    }

    try {
        nlohmann::json json_data = nlohmann::json::parse(file);

        // Normalize the pair key
        int norm_i = std::min(base_i, base_j);
        int norm_j = std::max(base_i, base_j);

        for (const auto& record : json_data) {
            int rec_i = record.value("base_i", 0);
            int rec_j = record.value("base_j", 0);

            int rec_norm_i = std::min(rec_i, rec_j);
            int rec_norm_j = std::max(rec_i, rec_j);

            if (rec_norm_i == norm_i && rec_norm_j == norm_j) {
                data.found = true;
                data.is_valid = record.value("is_valid", 0);
                data.bp_type_id = record.value("bp_type_id", -1);

                if (record.contains("direction_vectors")) {
                    auto& dir = record["direction_vectors"];
                    data.dir_x = dir.value("dir_x", 0.0);
                    data.dir_y = dir.value("dir_y", 0.0);
                    data.dir_z = dir.value("dir_z", 0.0);
                }

                if (record.contains("calculated_values")) {
                    auto& calc = record["calculated_values"];
                    data.dorg = calc.value("dorg", 0.0);
                    data.d_v = calc.value("d_v", 0.0);
                    data.plane_angle = calc.value("plane_angle", 0.0);
                    data.dNN = calc.value("dNN", 0.0);
                    data.quality_score = calc.value("quality_score", 0.0);
                }

                if (record.contains("validation_checks")) {
                    auto& checks = record["validation_checks"];
                    data.distance_check = checks.value("distance_check", false);
                    data.d_v_check = checks.value("d_v_check", false);
                    data.plane_angle_check = checks.value("plane_angle_check", false);
                    data.dNN_check = checks.value("dNN_check", false);
                }
                break;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << "\n";
    }

    return data;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string pdb_file = argv[1];
    int base_i = std::atoi(argv[2]);
    int base_j = std::atoi(argv[3]);
    std::string json_dir = "data/json_legacy";
    bool verbose = false;

    // Parse options
    for (int i = 4; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--json-dir" && i + 1 < argc) {
            json_dir = argv[++i];
        } else if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
    }

    // Extract PDB ID from filename
    std::string pdb_id;
    size_t last_slash = pdb_file.find_last_of("/\\");
    size_t last_dot = pdb_file.find_last_of('.');
    if (last_slash != std::string::npos) {
        pdb_id = pdb_file.substr(last_slash + 1, last_dot - last_slash - 1);
    } else {
        pdb_id = pdb_file.substr(0, last_dot);
    }

    std::cout << "=== Comparing Pair (" << base_i << ", " << base_j << ") in " << pdb_id << " ===\n\n";

    // Load and parse PDB file
    if (verbose) std::cout << "Loading PDB file: " << pdb_file << "\n";
    PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);

    // Calculate frames
    if (verbose) std::cout << "Calculating reference frames...\n";
    BaseFrameCalculator frame_calc;
    frame_calc.calculate_all_frames(structure);

    // Find residues by legacy index
    core::Residue* res1 = nullptr;
    core::Residue* res2 = nullptr;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx == base_i) res1 = &residue;
                if (legacy_idx == base_j) res2 = &residue;
            }
        }
    }

    if (!res1 || !res2) {
        std::cerr << "Error: Could not find residues with legacy indices "
                  << base_i << " and " << base_j << "\n";
        return 1;
    }

    if (verbose) {
        std::cout << "Found residue " << base_i << ": " << res1->name()
                  << " (chain " << res1->chain_id() << ", seq " << res1->seq_num() << ")\n";
        std::cout << "Found residue " << base_j << ": " << res2->name()
                  << " (chain " << res2->chain_id() << ", seq " << res2->seq_num() << ")\n\n";
        std::cout << "Running modern validation...\n";
    }
    BasePairValidator validator;
    auto modern_result = validator.validate(*res1, *res2);

    // Calculate quality score components
    QualityScoreCalculator score_calc;
    double hbond_adjustment = score_calc.adjust_pair_quality(modern_result.hbonds);
    int modern_bp_type_id = score_calc.calculate_bp_type_id(*res1, *res2, modern_result);
    double adjusted_quality = modern_result.quality_score + hbond_adjustment;
    if (modern_bp_type_id == 2) {
        adjusted_quality -= 2.0;  // WC bonus
    }

    // Load legacy data
    if (verbose) std::cout << "Loading legacy validation from: " << json_dir << "\n\n";
    auto legacy_data = load_legacy_pair_validation(json_dir, pdb_id, base_i, base_j);

    // Print comparison
    std::cout << std::fixed << std::setprecision(6);

    auto compare_float = [](const std::string& name, double leg, double mod, double tol = 1e-5) {
        double diff = std::abs(leg - mod);
        bool match = diff <= tol;
        std::cout << "  " << std::setw(20) << std::left << name << ": "
                  << std::setw(14) << leg << " vs " << std::setw(14) << mod;
        if (!match) {
            std::cout << " [DIFF: " << diff << "]";
        } else {
            std::cout << " [OK]";
        }
        std::cout << "\n";
        return match;
    };

    auto compare_bool = [](const std::string& name, bool leg, bool mod) {
        bool match = (leg == mod);
        std::cout << "  " << std::setw(20) << std::left << name << ": "
                  << std::setw(14) << (leg ? "true" : "false")
                  << " vs " << std::setw(14) << (mod ? "true" : "false");
        if (!match) {
            std::cout << " [MISMATCH]";
        } else {
            std::cout << " [OK]";
        }
        std::cout << "\n";
        return match;
    };

    auto compare_int = [](const std::string& name, int leg, int mod) {
        bool match = (leg == mod);
        std::cout << "  " << std::setw(20) << std::left << name << ": "
                  << std::setw(14) << leg << " vs " << std::setw(14) << mod;
        if (!match) {
            std::cout << " [MISMATCH]";
        } else {
            std::cout << " [OK]";
        }
        std::cout << "\n";
        return match;
    };

    bool all_match = true;

    if (!legacy_data.found) {
        std::cout << "WARNING: Pair not found in legacy JSON!\n";
        std::cout << "This might indicate different pair selection.\n\n";

        std::cout << "--- Modern Results Only ---\n";
        std::cout << "  is_valid: " << modern_result.is_valid << "\n";
        std::cout << "  dorg: " << modern_result.dorg << "\n";
        std::cout << "  d_v: " << modern_result.d_v << "\n";
        std::cout << "  plane_angle: " << modern_result.plane_angle << "\n";
        std::cout << "  dNN: " << modern_result.dNN << "\n";
        std::cout << "  overlap_area: " << modern_result.overlap_area << "\n";
        std::cout << "  num_base_hb: " << modern_result.num_base_hb << "\n";
        std::cout << "  num_o2_hb: " << modern_result.num_o2_hb << "\n";
        std::cout << "  quality_score: " << modern_result.quality_score << "\n";
        std::cout << "  hbond_adjustment: " << hbond_adjustment << "\n";
        std::cout << "  adjusted_quality: " << adjusted_quality << "\n";
        std::cout << "  Validation checks:\n";
        std::cout << "    distance_check: " << modern_result.distance_check << "\n";
        std::cout << "    d_v_check: " << modern_result.d_v_check << "\n";
        std::cout << "    plane_angle_check: " << modern_result.plane_angle_check << "\n";
        std::cout << "    dNN_check: " << modern_result.dNN_check << "\n";
        std::cout << "    overlap_check: " << modern_result.overlap_check << "\n";
        std::cout << "    hbond_check: " << modern_result.hbond_check << "\n";
        all_match = false;
    } else {
        std::cout << "--- Geometry ---\n";
        all_match &= compare_float("dorg", legacy_data.dorg, modern_result.dorg);
        all_match &= compare_float("d_v", legacy_data.d_v, modern_result.d_v);
        all_match &= compare_float("plane_angle", legacy_data.plane_angle, modern_result.plane_angle);
        all_match &= compare_float("dNN", legacy_data.dNN, modern_result.dNN);

        std::cout << "\n--- Direction Vectors ---\n";
        all_match &= compare_float("dir_x", legacy_data.dir_x, modern_result.dir_x);
        all_match &= compare_float("dir_y", legacy_data.dir_y, modern_result.dir_y);
        all_match &= compare_float("dir_z", legacy_data.dir_z, modern_result.dir_z);

        std::cout << "\n--- Validation Results ---\n";
        all_match &= compare_bool("is_valid", legacy_data.is_valid == 1, modern_result.is_valid);
        all_match &= compare_int("bp_type_id", legacy_data.bp_type_id, modern_bp_type_id);

        std::cout << "\n--- Quality Score ---\n";
        all_match &= compare_float("quality_score", legacy_data.quality_score, adjusted_quality);
    }

    std::cout << "\n=== SUMMARY ===\n";
    if (all_match) {
        std::cout << "RESULT: MATCH - Legacy and modern validation agree\n";
        return 0;
    } else {
        std::cout << "RESULT: MISMATCH - Differences found between legacy and modern\n";
        return 1;
    }
}
