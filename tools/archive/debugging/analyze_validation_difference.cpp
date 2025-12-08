/**
 * @file analyze_validation_difference.cpp
 * @brief Tool to analyze why a pair passes validation in modern but not in legacy
 *
 * Usage:
 *   build/analyze_validation_difference <pdb_file> <residue1_idx> <residue2_idx>
 *
 * Example:
 *   build/analyze_validation_difference data/pdb/1T0K.pdb 491 492
 */

#include <iostream>
#include <iomanip>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/geometry/vector3d.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;

void print_validation_result(const ValidationResult& result, const ValidationParameters& params) {
    std::cout << "\n=== Validation Result ===\n";
    std::cout << "is_valid: " << (result.is_valid ? "YES" : "NO") << "\n";
    std::cout << "\nGeometric Parameters:\n";
    std::cout << "  dorg: " << std::fixed << std::setprecision(6) << result.dorg << "\n";
    std::cout << "  d_v: " << result.d_v << "\n";
    std::cout << "  plane_angle: " << result.plane_angle << "\n";
    std::cout << "  dNN: " << result.dNN << "\n";
    std::cout << "  overlap_area: " << result.overlap_area << "\n";

    std::cout << "\nValidation Checks:\n";
    std::cout << "  distance_check (dorg): " << (result.distance_check ? "PASS" : "FAIL") << " [" << params.min_dorg
              << " <= " << result.dorg << " <= " << params.max_dorg << "]\n";
    std::cout << "  d_v_check: " << (result.d_v_check ? "PASS" : "FAIL") << " [" << params.min_dv
              << " <= " << result.d_v << " <= " << params.max_dv << "]\n";
    std::cout << "  plane_angle_check: " << (result.plane_angle_check ? "PASS" : "FAIL") << " ["
              << params.min_plane_angle << " <= " << result.plane_angle << " <= " << params.max_plane_angle << "]\n";
    std::cout << "  dNN_check: " << (result.dNN_check ? "PASS" : "FAIL") << " [" << params.min_dNN
              << " <= " << result.dNN << " <= " << params.max_dNN << "]\n";
    std::cout << "  overlap_check: " << (result.overlap_check ? "PASS" : "FAIL") << " [overlap_area < "
              << params.overlap_threshold << "]\n";
    std::cout << "  hbond_check: " << (result.hbond_check ? "PASS" : "FAIL") << " [num_base_hb=" << result.num_base_hb
              << ", min_base_hb=" << params.min_base_hb << "]\n";

    std::cout << "\nDirection Vectors:\n";
    std::cout << "  dir_x: " << result.dir_x << "\n";
    std::cout << "  dir_y: " << result.dir_y << "\n";
    std::cout << "  dir_z: " << result.dir_z << "\n";

    std::cout << "\nH-bonds:\n";
    std::cout << "  num_base_hb: " << result.num_base_hb << "\n";
    std::cout << "  num_o2_hb: " << result.num_o2_hb << "\n";
    std::cout << "  total hbonds: " << result.hbonds.size() << "\n";

    std::cout << "\nQuality Score:\n";
    std::cout << "  base_score: " << result.quality_score << "\n";
    std::cout << "  formula: dorg + 2.0 * d_v + plane_angle / 20.0\n";
    std::cout << "  = " << result.dorg << " + 2.0 * " << result.d_v << " + " << result.plane_angle << " / 20.0\n";
    std::cout << "  = " << result.dorg << " + " << (2.0 * result.d_v) << " + " << (result.plane_angle / 20.0) << "\n";
    std::cout << "  = " << result.quality_score << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/1T0K.pdb 491 492\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);

    std::cout << "Analyzing validation for pair (" << idx1 << ", " << idx2 << ") in " << pdb_file << "\n";
    std::cout << "=" << std::string(70, '=') << "\n";

    // Parse PDB
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);

    // Find residues (1-based legacy indices)
    Residue* mutable_res1 = nullptr;
    Residue* mutable_res2 = nullptr;

    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (residue.atoms().empty())
                continue;
            int legacy_idx = residue.atoms()[0].legacy_residue_idx();
            if (legacy_idx == idx1) {
                mutable_res1 = &residue;
            }
            if (legacy_idx == idx2) {
                mutable_res2 = &residue;
            }
        }
    }

    if (!mutable_res1 || !mutable_res2) {
        std::cerr << "Error: Could not find residues " << idx1 << " and/or " << idx2 << "\n";
        return 1;
    }

    std::cout << "Found residues:\n";
    std::cout << "  Residue " << idx1 << ": " << mutable_res1->name() << " (chain " << mutable_res1->chain_id()
              << ")\n";
    std::cout << "  Residue " << idx2 << ": " << mutable_res2->name() << " (chain " << mutable_res2->chain_id()
              << ")\n";

    // Calculate frames
    BaseFrameCalculator frame_calc("data/templates");
    std::cout << "\nCalculating frames...\n";
    auto frame1_result = frame_calc.calculate_frame(*mutable_res1);
    auto frame2_result = frame_calc.calculate_frame(*mutable_res2);

    if (!frame1_result.is_valid || !frame2_result.is_valid) {
        std::cerr << "Error: Failed to calculate frames\n";
        if (!frame1_result.is_valid)
            std::cerr << "  Frame 1 (residue " << idx1 << ") failed\n";
        if (!frame2_result.is_valid)
            std::cerr << "  Frame 2 (residue " << idx2 << ") failed\n";
        return 1;
    }
    std::cout << "Frames calculated successfully\n";

    // Validate pair
    BasePairValidator validator;
    auto result = validator.validate(*mutable_res1, *mutable_res2);

    print_validation_result(result, validator.parameters());

    // Determine why legacy might reject
    std::cout << "\n=== Analysis: Why Legacy Might Reject ===\n";

    bool passes_cdns = result.distance_check && result.d_v_check && result.plane_angle_check && result.dNN_check;

    if (!passes_cdns) {
        std::cout << "❌ Fails cdns (distance/angle checks)\n";
        if (!result.distance_check)
            std::cout << "  - dorg check failed\n";
        if (!result.d_v_check)
            std::cout << "  - d_v check failed\n";
        if (!result.plane_angle_check)
            std::cout << "  - plane_angle check failed\n";
        if (!result.dNN_check)
            std::cout << "  - dNN check failed\n";
    } else {
        std::cout << "✅ Passes cdns (distance/angle checks)\n";
    }

    if (!result.overlap_check) {
        std::cout << "❌ Fails overlap check\n";
    } else {
        std::cout << "✅ Passes overlap check\n";
    }

    if (!result.hbond_check) {
        std::cout << "❌ Fails H-bond check\n";
    } else {
        std::cout << "✅ Passes H-bond check\n";
    }

    if (result.is_valid) {
        std::cout << "\n✅ Modern validation: PASSES\n";
        std::cout << "If legacy rejects, possible causes:\n";
        std::cout << "  1. Different validation thresholds\n";
        std::cout << "  2. Different frame calculations\n";
        std::cout << "  3. Different overlap calculation\n";
        std::cout << "  4. Different H-bond detection\n";
        std::cout << "  5. Early rejection (before validation)\n";
    } else {
        std::cout << "\n❌ Modern validation: FAILS\n";
        std::cout << "Both legacy and modern reject this pair\n";
    }

    return 0;
}
