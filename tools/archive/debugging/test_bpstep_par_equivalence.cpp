/**
 * @file test_bpstep_par_equivalence.cpp
 * @brief Test to verify bpstep_par implementation matches legacy exactly
 *
 * This tool tests the bpstep_par implementation with known values
 * to ensure numerical precision matches legacy.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

void test_bpstep_par_with_frames(const ReferenceFrame& frame1, const ReferenceFrame& frame2, const std::string& label) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "Testing: " << label << "\n";
    std::cout << std::string(60, '=') << "\n";

    ParameterCalculator calc;
    BasePairStepParameters params = calc.calculate_step_parameters(frame1, frame2);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Step Parameters:\n";
    std::cout << "  Shift:  " << std::setw(12) << params.shift << "\n";
    std::cout << "  Slide:  " << std::setw(12) << params.slide << " (shear for bp_type_id)\n";
    std::cout << "  Rise:   " << std::setw(12) << params.rise << " (stretch for bp_type_id)\n";
    std::cout << "  Tilt:   " << std::setw(12) << params.tilt << "\n";
    std::cout << "  Roll:   " << std::setw(12) << params.roll << "\n";
    std::cout << "  Twist:  " << std::setw(12) << params.twist << " (opening for bp_type_id)\n";

    // Check bp_type_id thresholds
    std::cout << "\nbp_type_id Threshold Checks:\n";
    double abs_shear = std::abs(params.slide);
    double abs_stretch = std::abs(params.rise);
    double abs_opening = std::abs(params.twist);

    std::cout << "  fabs(shear) = " << abs_shear;
    if (abs_shear <= 1.8) {
        std::cout << " <= 1.8 ✅ (Watson-Crick candidate)\n";
    } else if (abs_shear >= 1.8 && abs_shear <= 2.8) {
        std::cout << " in [1.8, 2.8] ✅ (Wobble candidate)\n";
    } else {
        std::cout << " > 2.8 ❌ (Outside range)\n";
    }

    std::cout << "  fabs(stretch) = " << abs_stretch;
    if (abs_stretch <= 2.0) {
        std::cout << " <= 2.0 ✅\n";
    } else {
        std::cout << " > 2.0 ❌ (Exceeds threshold)\n";
    }

    std::cout << "  fabs(opening) = " << abs_opening;
    if (abs_opening <= 60.0) {
        std::cout << " <= 60.0 ✅\n";
    } else {
        std::cout << " > 60.0 ❌ (Exceeds threshold)\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <idx1> <idx2>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 1024 1188\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);

    // Load structure
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);

    // Find residues by legacy index
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
    BaseFrameCalculator frame_calc("data/templates");
    FrameCalculationResult frame1_result = frame_calc.calculate_frame(*mutable_res1);
    FrameCalculationResult frame2_result = frame_calc.calculate_frame(*mutable_res2);

    if (!frame1_result.is_valid || !frame2_result.is_valid) {
        std::cerr << "Error: Frame calculation failed for one or both residues\n";
        return 1;
    }

    ReferenceFrame frame1 = frame1_result.frame;
    ReferenceFrame frame2 = frame2_result.frame;

    // Test with original frames
    test_bpstep_par_with_frames(frame1, frame2, "Original Frames");

    // Test with reversed frame2 (if dir_z <= 0)
    // Calculate direction vectors
    Vector3D x1 = frame1.rotation().column(0);
    Vector3D y1 = frame1.rotation().column(1);
    Vector3D z1 = frame1.rotation().column(2);
    Vector3D x2 = frame2.rotation().column(0);
    Vector3D y2 = frame2.rotation().column(1);
    Vector3D z2 = frame2.rotation().column(2);

    double dir_x = x1.dot(x2);
    double dir_y = y1.dot(y2);
    double dir_z = z1.dot(z2);

    std::cout << "\nDirection Vectors:\n";
    std::cout << "  dir_x: " << dir_x << "\n";
    std::cout << "  dir_y: " << dir_y << "\n";
    std::cout << "  dir_z: " << dir_z << "\n";

    if (dir_z <= 0.0) {
        std::cout << "\nApplying frame reversal (dir_z <= 0)...\n";
        Matrix3D rot2 = frame2.rotation();
        Vector3D y_col = rot2.column(1);
        Vector3D z_col = rot2.column(2);
        rot2.set_column(1, -y_col);
        rot2.set_column(2, -z_col);
        ReferenceFrame frame2_reversed(rot2, frame2.origin());

        // Test with reversed frame2 (legacy order: r2, r1)
        test_bpstep_par_with_frames(frame2_reversed, frame1, "After Frame Reversal (r2, r1)");
    }

    return 0;
}
