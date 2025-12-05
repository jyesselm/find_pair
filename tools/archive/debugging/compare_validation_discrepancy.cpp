/**
 * @file compare_validation_discrepancy.cpp
 * @brief Compare validation results between legacy and modern for specific pairs
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;
using json = nlohmann::json;

struct LegacyValidationResult {
    bool is_valid;
    double dorg;
    double d_v;
    double plane_angle;
    double dNN;
    double quality_score;
    int bp_type_id;
    double dir_x, dir_y, dir_z;
    bool distance_check, d_v_check, plane_angle_check, dNN_check;
    int num_base_hb, num_o2_hb;
};

LegacyValidationResult load_legacy_validation(const std::string& json_file, int idx1, int idx2) {
    LegacyValidationResult result = {};
    result.is_valid = false;

    std::ifstream file(json_file);
    if (!file.is_open()) {
        return result;
    }

    json data;
    file >> data;

    // Find the pair (0-based indices in JSON)
    int json_idx1 = idx1 - 1;
    int json_idx2 = idx2 - 1;

    for (const auto& entry : data) {
        int r1 = entry.value("residue1_idx", -1);
        int r2 = entry.value("residue2_idx", -1);

        if ((r1 == json_idx1 && r2 == json_idx2) || (r1 == json_idx2 && r2 == json_idx1)) {
            result.is_valid = entry.value("is_valid", 0) != 0;
            result.dorg = entry.value("dorg", 0.0);
            result.d_v = entry.value("d_v", 0.0);
            result.plane_angle = entry.value("plane_angle", 0.0);
            result.dNN = entry.value("dNN", 0.0);
            result.quality_score = entry.value("quality_score", 0.0);
            result.bp_type_id = entry.value("bp_type_id", -1);
            result.dir_x = entry.value("dir_x", 0.0);
            result.dir_y = entry.value("dir_y", 0.0);
            result.dir_z = entry.value("dir_z", 0.0);
            result.distance_check = entry.value("distance_check", false);
            result.d_v_check = entry.value("d_v_check", false);
            result.plane_angle_check = entry.value("plane_angle_check", false);
            result.dNN_check = entry.value("dNN_check", false);
            result.num_base_hb = entry.value("num_base_hb", 0);
            result.num_o2_hb = entry.value("num_o2_hb", 0);
            break;
        }
    }

    return result;
}

void print_comparison(const LegacyValidationResult& legacy, const ValidationResult& modern,
                      const ReferenceFrame& frame1, const ReferenceFrame& frame2, int idx1,
                      int idx2) {
    std::cout << "\n============================================================\n";
    std::cout << "COMPARISON: Pair (" << idx1 << ", " << idx2 << ")\n";
    std::cout << "============================================================\n\n";

    std::cout << "FRAMES:\n";
    std::cout << "  Frame 1 origin: [" << frame1.origin().x() << ", " << frame1.origin().y() << ", "
              << frame1.origin().z() << "]\n";
    std::cout << "  Frame 2 origin: [" << frame2.origin().x() << ", " << frame2.origin().y() << ", "
              << frame2.origin().z() << "]\n";
    std::cout << "  Origin distance: " << (frame1.origin() - frame2.origin()).length() << " Å\n";

    std::cout << "\nDIRECTION VECTORS:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Legacy:  dir_x=" << legacy.dir_x << ", dir_y=" << legacy.dir_y
              << ", dir_z=" << legacy.dir_z << "\n";
    std::cout << "  Modern:  dir_x=" << modern.dir_x << ", dir_y=" << modern.dir_y
              << ", dir_z=" << modern.dir_z << "\n";
    std::cout << "  Diff:    dir_x=" << (modern.dir_x - legacy.dir_x)
              << ", dir_y=" << (modern.dir_y - legacy.dir_y)
              << ", dir_z=" << (modern.dir_z - legacy.dir_z) << "\n";

    std::cout << "\nGEOMETRIC PARAMETERS:\n";
    std::cout << "  dorg:        Legacy=" << legacy.dorg << ", Modern=" << modern.dorg
              << ", Diff=" << (modern.dorg - legacy.dorg) << "\n";
    std::cout << "  d_v:          Legacy=" << legacy.d_v << ", Modern=" << modern.d_v
              << ", Diff=" << (modern.d_v - legacy.d_v) << "\n";
    std::cout << "  plane_angle:  Legacy=" << legacy.plane_angle
              << ", Modern=" << modern.plane_angle
              << ", Diff=" << (modern.plane_angle - legacy.plane_angle) << "\n";
    std::cout << "  dNN:          Legacy=" << legacy.dNN << ", Modern=" << modern.dNN
              << ", Diff=" << (modern.dNN - legacy.dNN) << "\n";

    std::cout << "\nVALIDATION CHECKS:\n";
    std::cout << "  distance_check:    Legacy=" << (legacy.distance_check ? "PASS" : "FAIL")
              << ", Modern=" << (modern.distance_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  d_v_check:          Legacy=" << (legacy.d_v_check ? "PASS" : "FAIL")
              << ", Modern=" << (modern.d_v_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  plane_angle_check:  Legacy=" << (legacy.plane_angle_check ? "PASS" : "FAIL")
              << ", Modern=" << (modern.plane_angle_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  dNN_check:          Legacy=" << (legacy.dNN_check ? "PASS" : "FAIL")
              << ", Modern=" << (modern.dNN_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  overlap_check:      Legacy=N/A, Modern="
              << (modern.overlap_check ? "PASS" : "FAIL") << "\n";

    std::cout << "\nHYDROGEN BONDS:\n";
    std::cout << "  num_base_hb:  Legacy=" << legacy.num_base_hb
              << ", Modern=" << modern.num_base_hb << "\n";
    std::cout << "  num_o2_hb:     Legacy=" << legacy.num_o2_hb << ", Modern=" << modern.num_o2_hb
              << "\n";

    std::cout << "\nFINAL RESULT:\n";
    std::cout << "  Legacy:  is_valid=" << (legacy.is_valid ? "YES" : "NO")
              << ", quality=" << legacy.quality_score << ", bp_type_id=" << legacy.bp_type_id
              << "\n";
    std::cout << "  Modern:  is_valid=" << (modern.is_valid ? "YES" : "NO")
              << ", quality=" << modern.quality_score
              << ", bp_type=" << static_cast<int>(modern.bp_type) << "\n";

    if (legacy.is_valid != modern.is_valid) {
        std::cout << "\n*** DISCREPANCY: Validation results differ! ***\n";
        if (!legacy.is_valid && modern.is_valid) {
            std::cout << "  Legacy marks as INVALID, Modern marks as VALID\n";
            std::cout << "  Checking which validation step fails in legacy:\n";
            if (!legacy.distance_check)
                std::cout << "    - distance_check FAILED\n";
            if (!legacy.d_v_check)
                std::cout << "    - d_v_check FAILED\n";
            if (!legacy.plane_angle_check)
                std::cout << "    - plane_angle_check FAILED\n";
            if (!legacy.dNN_check)
                std::cout << "    - dNN_check FAILED\n";
        } else {
            std::cout << "  Legacy marks as VALID, Modern marks as INVALID\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <pdb_file> <legacy_json> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0]
                  << " data/pdb/6CAQ.pdb data/json_legacy/pair_validation/6CAQ.json 980 997\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    std::string legacy_json = argv[2];
    int idx1 = std::stoi(argv[3]);
    int idx2 = std::stoi(argv[4]);

    std::cout << "============================================================\n";
    std::cout << "Validation Discrepancy Comparison Tool\n";
    std::cout << "============================================================\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "Legacy JSON: " << legacy_json << "\n";
    std::cout << "Pair: (" << idx1 << ", " << idx2 << ")\n\n";

    // Load legacy validation result
    LegacyValidationResult legacy = load_legacy_validation(legacy_json, idx1, idx2);
    if (legacy.dorg == 0.0 && legacy.d_v == 0.0) {
        std::cerr << "Warning: Could not find pair in legacy JSON\n";
    }

    // Parse PDB and find residues
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);

    Residue* mutable_res1 = nullptr;
    Residue* mutable_res2 = nullptr;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx == idx1) {
                    mutable_res1 = &residue;
                }
                if (legacy_idx == idx2) {
                    mutable_res2 = &residue;
                }
            }
        }
    }

    if (!mutable_res1 || !mutable_res2) {
        std::cerr << "Error: Could not find residues " << idx1 << " and/or " << idx2 << "\n";
        return 1;
    }

    std::cout << "Residue 1 (legacy_idx=" << idx1 << "): " << mutable_res1->name() << " Chain "
              << mutable_res1->chain_id() << " Seq " << mutable_res1->seq_num() << "\n";
    std::cout << "Residue 2 (legacy_idx=" << idx2 << "): " << mutable_res2->name() << " Chain "
              << mutable_res2->chain_id() << " Seq " << mutable_res2->seq_num() << "\n";

    // Calculate frames
    BaseFrameCalculator calculator("data/templates");
    bool is_rna = false;
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == " O2'" || atom.name() == " O2*") {
                    is_rna = true;
                    break;
                }
            }
            if (is_rna)
                break;
        }
        if (is_rna)
            break;
    }

    FrameCalculationResult frame1_result = calculator.calculate_frame(*mutable_res1);
    FrameCalculationResult frame2_result = calculator.calculate_frame(*mutable_res2);

    if (!frame1_result.is_valid || !frame2_result.is_valid) {
        std::cerr << "Error: Frame calculation failed\n";
        return 1;
    }

    // Store frames on residues
    mutable_res1->set_reference_frame(frame1_result.frame);
    mutable_res2->set_reference_frame(frame2_result.frame);

    ReferenceFrame frame1 = frame1_result.frame;
    ReferenceFrame frame2 = frame2_result.frame;

    // Run modern validation
    BasePairValidator validator;
    ValidationResult modern = validator.validate(*mutable_res1, *mutable_res2);

    // Print comparison
    print_comparison(legacy, modern, frame1, frame2, idx1, idx2);

    return 0;
}
