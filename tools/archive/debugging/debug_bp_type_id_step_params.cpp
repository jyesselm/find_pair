/**
 * @file debug_bp_type_id_step_params.cpp
 * @brief Debug tool to investigate step parameters for pairs with bp_type_id differences
 *
 * This tool calculates step parameters for specific pairs and outputs them
 * to help identify why bp_type_id differs between legacy and modern.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/residue_index_fixer.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <nlohmann/json.hpp>
#include <filesystem>

using namespace x3dna;
using json = nlohmann::json;

void print_separator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <legacy_idx1> <legacy_idx2> [pdb_id] [legacy_json_file]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 1141 1151 6CAQ\n";
        std::cerr << "         " << argv[0]
                  << " data/pdb/6CAQ.pdb 1141 1151 6CAQ data/json_legacy/base_frame_calc/6CAQ.json\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int target_idx1 = std::stoi(argv[2]);
    int target_idx2 = std::stoi(argv[3]);
    std::string pdb_id = (argc > 4) ? argv[4] : "";
    std::string legacy_json_file = (argc > 5) ? argv[5] : "";

    std::cout << "Debug step parameters for bp_type_id calculation\n";
    std::cout << "Pair (" << target_idx1 << ", " << target_idx2 << ")\n";
    std::cout << "PDB file: " << pdb_file << "\n";

    // Step 1: Parse PDB
    print_separator("STEP 1: Parse PDB and build residue mapping");

    io::PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);

    // Optionally fix residue indices from legacy JSON
    if (!legacy_json_file.empty() && std::filesystem::exists(legacy_json_file)) {
        std::cout << "Fixing residue indices from: " << legacy_json_file << "\n";
        int fixed_count = io::fix_residue_indices_from_json(structure, legacy_json_file);
        std::cout << "Fixed " << fixed_count << " residue indices\n";
    } else if (legacy_json_file.empty()) {
        // Try to auto-detect legacy JSON file
        std::filesystem::path pdb_path(pdb_file);
        std::string pdb_id_from_file = pdb_path.stem().string();
        std::filesystem::path auto_json = std::filesystem::path("data/json_legacy/base_frame_calc") /
                                          (pdb_id_from_file + ".json");
        if (std::filesystem::exists(auto_json)) {
            std::cout << "Auto-detected legacy JSON: " << auto_json << "\n";
            int fixed_count = io::fix_residue_indices_from_json(structure, auto_json.string());
            std::cout << "Fixed " << fixed_count << " residue indices\n";
        }
    }

    // Build mapping from legacy_residue_idx to residue
    std::map<int, core::Residue*> residue_by_legacy_idx;
    int max_legacy_idx = 0;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residue_by_legacy_idx[legacy_idx] = &residue;
                    if (legacy_idx > max_legacy_idx) {
                        max_legacy_idx = legacy_idx;
                    }
                }
            }
        }
    }

    // Step 2: Find target residues
    print_separator("STEP 2: Find target residues");

    auto it1 = residue_by_legacy_idx.find(target_idx1);
    auto it2 = residue_by_legacy_idx.find(target_idx2);

    if (it1 == residue_by_legacy_idx.end()) {
        std::cerr << "ERROR: Residue at legacy_idx " << target_idx1 << " not found!\n";
        return 1;
    }
    if (it2 == residue_by_legacy_idx.end()) {
        std::cerr << "ERROR: Residue at legacy_idx " << target_idx2 << " not found!\n";
        return 1;
    }

    core::Residue* res1 = it1->second;
    core::Residue* res2 = it2->second;

    std::cout << "Residue 1 (legacy_idx=" << target_idx1 << "):\n";
    std::cout << "  Name: " << res1->name() << "\n";
    std::cout << "  Type: " << static_cast<int>(res1->residue_type()) << "\n";
    std::cout << "  Chain: " << res1->chain_id() << "\n";
    std::cout << "  Seq: " << res1->seq_num() << "\n";
    std::cout << "  Insertion: '" << res1->insertion() << "'\n";
    std::cout << "  One-letter: " << res1->one_letter_code() << "\n";
    std::cout << "  Num atoms: " << res1->num_atoms() << "\n";

    std::cout << "Residue 2 (legacy_idx=" << target_idx2 << "):\n";
    std::cout << "  Name: " << res2->name() << "\n";
    std::cout << "  Type: " << static_cast<int>(res2->residue_type()) << "\n";
    std::cout << "  Chain: " << res2->chain_id() << "\n";
    std::cout << "  Seq: " << res2->seq_num() << "\n";
    std::cout << "  Insertion: '" << res2->insertion() << "'\n";
    std::cout << "  One-letter: " << res2->one_letter_code() << "\n";
    std::cout << "  Num atoms: " << res2->num_atoms() << "\n";

    // Step 3: Calculate frames
    print_separator("STEP 3: Calculate frames");

    // Try resources/templates first, then data/templates
    std::filesystem::path template_path = "resources/templates";
    if (!std::filesystem::exists(template_path)) {
        template_path = "data/templates";
    }
    algorithms::BaseFrameCalculator calculator(template_path.string());
    auto frame_result1 = calculator.calculate_frame(*res1);
    auto frame_result2 = calculator.calculate_frame(*res2);

    if (!frame_result1.is_valid || !frame_result2.is_valid) {
        std::cerr << "ERROR: Frame calculation failed for one or both residues\n";
        return 1;
    }

    std::cout << "Frames calculated successfully\n";
    std::cout << "Frame 1 origin: [" << frame_result1.frame.origin().x() << ", " << frame_result1.frame.origin().y()
              << ", " << frame_result1.frame.origin().z() << "]\n";
    std::cout << "Frame 1 RMS fit: " << frame_result1.rms_fit << "\n";
    std::cout << "Frame 1 matched atoms (" << frame_result1.num_matched << "): ";
    for (const auto& name : frame_result1.matched_atoms) {
        std::cout << name << " ";
    }
    std::cout << "\n";
    std::cout << "Frame 2 origin: [" << frame_result2.frame.origin().x() << ", " << frame_result2.frame.origin().y()
              << ", " << frame_result2.frame.origin().z() << "]\n";
    std::cout << "Frame 2 RMS fit: " << frame_result2.rms_fit << "\n";
    std::cout << "Frame 2 matched atoms (" << frame_result2.num_matched << "): ";
    for (const auto& name : frame_result2.matched_atoms) {
        std::cout << name << " ";
    }
    std::cout << "\n";

    // Calculate dorg from frame origins
    auto dorg_vec = frame_result1.frame.origin() - frame_result2.frame.origin();
    double dorg_calc = std::sqrt(dorg_vec.x() * dorg_vec.x() + dorg_vec.y() * dorg_vec.y() +
                                 dorg_vec.z() * dorg_vec.z());
    std::cout << "Calculated dorg from frame origins: " << dorg_calc << "\n";
    std::cout << "\n";

    // Show some atom coordinates for debugging
    std::cout << "Residue 1 atoms (first 5):\n";
    for (size_t i = 0; i < std::min(size_t(5), res1->atoms().size()); ++i) {
        const auto& atom = res1->atoms()[i];
        auto pos = atom.position();
        std::cout << "  " << atom.name() << ": [" << pos.x() << ", " << pos.y() << ", " << pos.z() << "]\n";
    }
    std::cout << "Residue 2 atoms (first 5):\n";
    for (size_t i = 0; i < std::min(size_t(5), res2->atoms().size()); ++i) {
        const auto& atom = res2->atoms()[i];
        auto pos = atom.position();
        std::cout << "  " << atom.name() << ": [" << pos.x() << ", " << pos.y() << ", " << pos.z() << "]\n";
    }
    std::cout << "\n";

    // Step 4: Run validation
    print_separator("STEP 4: Run validation");

    algorithms::BasePairValidator validator;
    algorithms::ValidationResult validation_result = validator.validate(*res1, *res2);

    std::cout << "Validation result:\n";
    std::cout << "  dorg: " << validation_result.dorg << "\n";
    std::cout << "  d_v: " << validation_result.d_v << "\n";
    std::cout << "  plane_angle: " << validation_result.plane_angle << "\n";
    std::cout << "  dir_x: " << validation_result.dir_x << "\n";
    std::cout << "  dir_y: " << validation_result.dir_y << "\n";
    std::cout << "  dir_z: " << validation_result.dir_z << "\n";
    std::cout << "  is_valid: " << (validation_result.is_valid ? "YES" : "NO") << "\n";

    // Check direction vector condition
    bool condition_met = (validation_result.dir_x > 0.0 && validation_result.dir_y < 0.0 &&
                          validation_result.dir_z < 0.0);
    std::cout << "\nDirection Vector Condition (dir_x > 0 && dir_y < 0 && dir_z < 0): "
              << (condition_met ? "MET" : "NOT MET") << "\n";

    if (!condition_met) {
        std::cout << "  ⚠️  Condition not met - bp_type_id should remain -1\n";
        return 0;
    }

    // Step 5: Calculate step parameters
    print_separator("STEP 5: Calculate step parameters");

    algorithms::ParameterCalculator param_calculator;

    // Get frames
    auto frame1_opt = res1->reference_frame();
    auto frame2_opt = res2->reference_frame();

    if (!frame1_opt.has_value() || !frame2_opt.has_value()) {
        std::cerr << "ERROR: Frames not available on residues\n";
        return 1;
    }

    core::ReferenceFrame frame1 = frame1_opt.value();
    core::ReferenceFrame frame2 = frame2_opt.value();

    // Apply frame reversal if dir_z <= 0 (matches legacy logic)
    if (validation_result.dir_z <= 0.0) {
        std::cout << "Applying frame reversal (dir_z <= 0)\n";
        geometry::Matrix3D rot2 = frame2.rotation();
        geometry::Vector3D y_col = rot2.column(1);
        geometry::Vector3D z_col = rot2.column(2);
        rot2.set_column(1, -y_col);
        rot2.set_column(2, -z_col);
        frame2 = core::ReferenceFrame(rot2, frame2.origin());
    } else {
        std::cout << "No frame reversal needed (dir_z > 0)\n";
    }

    // Calculate step parameters (frame2 first, frame1 second - matches legacy order)
    core::BasePairStepParameters params = param_calculator.calculate_step_parameters(frame2, frame1);

    std::cout << "\nStep Parameters:\n";
    std::cout << "  Shift:  " << std::fixed << std::setprecision(6) << params.shift << "\n";
    std::cout << "  Slide:  " << params.slide << " (shear)\n";
    std::cout << "  Rise:   " << params.rise << " (stretch)\n";
    std::cout << "  Tilt:   " << params.tilt << "\n";
    std::cout << "  Roll:   " << params.roll << "\n";
    std::cout << "  Twist:  " << params.twist << " (opening, degrees)\n";

    // Step 6: Check thresholds
    print_separator("STEP 6: Check thresholds and classification");

    double shear = params.slide;
    double stretch = params.rise;
    double opening = params.twist;

    std::cout << "Geometric Parameters:\n";
    std::cout << "  shear (slide):   " << std::abs(shear) << "\n";
    std::cout << "  stretch (rise):  " << std::abs(stretch) << "\n";
    std::cout << "  opening (twist): " << std::abs(opening) << " degrees\n";

    std::cout << "\nThreshold Checks:\n";
    bool stretch_ok = (std::abs(stretch) <= 2.0);
    bool opening_ok = (std::abs(opening) <= 60.0);
    std::cout << "  stretch <= 2.0:  " << (stretch_ok ? "PASS" : "FAIL") << " (abs=" << std::abs(stretch) << ")\n";
    std::cout << "  opening <= 60.0: " << (opening_ok ? "PASS" : "FAIL") << " (abs=" << std::abs(opening) << ")\n";

    if (!stretch_ok || !opening_ok) {
        std::cout << "  ⚠️  Threshold check failed - bp_type_id should remain -1\n";
        return 0;
    }

    // Step 7: Check base pair type
    print_separator("STEP 7: Check base pair type (WC_LIST matching)");

    // Get base letters
    char base1_char = '?';
    char base2_char = '?';

    switch (res1->residue_type()) {
        case core::ResidueType::ADENINE:
            base1_char = 'A';
            break;
        case core::ResidueType::CYTOSINE:
            base1_char = 'C';
            break;
        case core::ResidueType::GUANINE:
            base1_char = 'G';
            break;
        case core::ResidueType::THYMINE:
            base1_char = 'T';
            break;
        case core::ResidueType::URACIL:
            base1_char = 'U';
            break;
        default:
            base1_char = '?';
            break;
    }

    switch (res2->residue_type()) {
        case core::ResidueType::ADENINE:
            base2_char = 'A';
            break;
        case core::ResidueType::CYTOSINE:
            base2_char = 'C';
            break;
        case core::ResidueType::GUANINE:
            base2_char = 'G';
            break;
        case core::ResidueType::THYMINE:
            base2_char = 'T';
            break;
        case core::ResidueType::URACIL:
            base2_char = 'U';
            break;
        default:
            base2_char = '?';
            break;
    }

    std::string bp_type = std::string(1, base1_char) + std::string(1, base2_char);

    std::cout << "Base Pair Type: \"" << bp_type << "\"\n";
    std::cout << "  Residue 1: " << res1->name() << " -> " << base1_char << "\n";
    std::cout << "  Residue 2: " << res2->name() << " -> " << base2_char << "\n";

    static const std::vector<std::string> WC_LIST = {"XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"};

    bool in_wc_list = false;
    for (const auto& wc : WC_LIST) {
        if (bp_type == wc) {
            in_wc_list = true;
            break;
        }
    }

    std::cout << "  In WC_LIST: " << (in_wc_list ? "YES" : "NO") << "\n";

    // Step 8: Determine bp_type_id
    print_separator("STEP 8: Determine bp_type_id");

    int bp_type_id = -1;

    // Check for wobble pair (fabs(shear) in [1.8, 2.8])
    if (std::abs(shear) >= 1.8 && std::abs(shear) <= 2.8) {
        bp_type_id = 1; // Wobble
        std::cout << "Wobble pair detected: abs(shear) = " << std::abs(shear) << " in [1.8, 2.8]\n";
    }

    // Check for Watson-Crick pair (fabs(shear) <= 1.8 AND in WC_LIST)
    if (std::abs(shear) <= 1.8) {
        std::cout << "Shear check: abs(shear) = " << std::abs(shear) << " <= 1.8\n";
        if (in_wc_list) {
            bp_type_id = 2; // Watson-Crick
            std::cout << "Watson-Crick pair detected: in WC_LIST\n";
        } else {
            std::cout << "Not in WC_LIST - keeping previous assignment\n";
        }
    }

    std::cout << "\nFinal bp_type_id: " << bp_type_id << "\n";
    if (bp_type_id == -1) {
        std::cout << "  (-1 = not classified)\n";
    } else if (bp_type_id == 1) {
        std::cout << "  (1 = wobble pair)\n";
    } else if (bp_type_id == 2) {
        std::cout << "  (2 = Watson-Crick pair)\n";
    }

    // Step 9: Compare with legacy if JSON available
    if (!pdb_id.empty()) {
        print_separator("STEP 9: Compare with legacy JSON");

        std::string legacy_file = "data/json_legacy/pair_validation/" + pdb_id + ".json";
        std::ifstream f(legacy_file);
        if (f.good()) {
            json legacy_data;
            f >> legacy_data;

            for (const auto& rec : legacy_data) {
                if (rec.contains("base_i") && rec.contains("base_j")) {
                    int base_i = rec["base_i"];
                    int base_j = rec["base_j"];
                    if ((base_i == target_idx1 && base_j == target_idx2) ||
                        (base_i == target_idx2 && base_j == target_idx1)) {
                        int legacy_bp_type_id = rec.value("bp_type_id", 0);
                        std::cout << "Legacy bp_type_id: " << legacy_bp_type_id << "\n";
                        std::cout << "Modern bp_type_id: " << bp_type_id << "\n";
                        if (legacy_bp_type_id != bp_type_id) {
                            std::cout << "  ⚠️  MISMATCH DETECTED!\n";
                        } else {
                            std::cout << "  ✅ Match\n";
                        }
                        break;
                    }
                }
            }
        } else {
            std::cout << "Legacy JSON file not found: " << legacy_file << "\n";
        }
    }

    print_separator("INVESTIGATION COMPLETE");

    return 0;
}
