/**
 * @file compare_frames_and_step_params.cpp
 * @brief Compare final frames and 6 base pair step parameters between legacy and modern
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <x3dna/io/json_reader.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

void print_frame(const std::string& label, const ReferenceFrame& frame) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n" << label << ":\n";
    std::cout << "  Origin: [" << frame.origin().x() << ", " 
              << frame.origin().y() << ", " << frame.origin().z() << "]\n";
    
    Matrix3D rot = frame.rotation();
    std::cout << "  Rotation matrix:\n";
    for (size_t i = 0; i < 3; i++) {
        std::cout << "    [";
        for (size_t j = 0; j < 3; j++) {
            std::cout << std::setw(10) << rot.at(i, j);
            if (j < 2) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    // Print as 9-element array (legacy format)
    std::cout << "  As 9-element array (legacy format):\n    [";
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            std::cout << rot.at(i, j);
            if (i < 2 || j < 2) std::cout << ", ";
        }
    }
    std::cout << "]\n";
}

void print_step_parameters(const BasePairStepParameters& params) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n6 BASE PAIR STEP PARAMETERS (from bpstep_par):\n";
    std::cout << "  1. Shift:          " << std::setw(10) << params.shift << "\n";
    std::cout << "  2. Slide (shear):  " << std::setw(10) << params.slide << "\n";
    std::cout << "  3. Rise (stretch): " << std::setw(10) << params.rise << "\n";
    std::cout << "  4. Tilt:          " << std::setw(10) << params.tilt << "\n";
    std::cout << "  5. Roll:          " << std::setw(10) << params.roll << "\n";
    std::cout << "  6. Twist (opening):" << std::setw(10) << params.twist << "\n";
    
    // Print as array (legacy format: pars[1..6])
    // Legacy uses: pars[1]=Slide, pars[2]=Rise, pars[6]=Twist for bp_type_id
    std::cout << "\n  As array [pars[1], pars[2], pars[3], pars[4], pars[5], pars[6]]:\n";
    std::cout << "    [" << params.shift << ", " << params.slide << ", " 
              << params.rise << ", " << params.tilt << ", " 
              << params.roll << ", " << params.twist << "]\n";
    std::cout << "\n  Parameters used for bp_type_id:\n";
    std::cout << "    pars[1] (Slide/Shear): " << params.slide << "\n";
    std::cout << "    pars[2] (Rise/Stretch): " << params.rise << "\n";
    std::cout << "    pars[6] (Twist/Opening): " << params.twist << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 1024 1188\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);

    std::cout << "============================================================\n";
    std::cout << "Frame and Step Parameter Comparison Tool\n";
    std::cout << "============================================================\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "Pair: (" << idx1 << ", " << idx2 << ")\n";

    // Parse PDB
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);

    // Find residues
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

    std::cout << "\nResidue 1 (legacy_idx=" << idx1 << "): " << mutable_res1->name() 
              << " Chain " << mutable_res1->chain_id() << " Seq " << mutable_res1->seq_num() << "\n";
    std::cout << "Residue 2 (legacy_idx=" << idx2 << "): " << mutable_res2->name() 
              << " Chain " << mutable_res2->chain_id() << " Seq " << mutable_res2->seq_num() << "\n";

    // Calculate and store frames
    BaseFrameCalculator calculator("data/templates");
    FrameCalculationResult frame1_result = calculator.calculate_frame(*mutable_res1);
    FrameCalculationResult frame2_result = calculator.calculate_frame(*mutable_res2);

    if (!frame1_result.is_valid || !frame2_result.is_valid) {
        std::cerr << "Error: Frame calculation failed\n";
        return 1;
    }

    mutable_res1->set_reference_frame(frame1_result.frame);
    mutable_res2->set_reference_frame(frame2_result.frame);

    // Print frames
    print_frame("FRAME 1 (Residue " + std::to_string(idx1) + ")", frame1_result.frame);
    print_frame("FRAME 2 (Residue " + std::to_string(idx2) + ")", frame2_result.frame);

    // Calculate direction vectors
    BasePairValidator validator;
    ValidationResult result = validator.validate(*mutable_res1, *mutable_res2);

    std::cout << "\n============================================================\n";
    std::cout << "DIRECTION VECTORS:\n";
    std::cout << "  dir_x: " << result.dir_x << "\n";
    std::cout << "  dir_y: " << result.dir_y << "\n";
    std::cout << "  dir_z: " << result.dir_z << "\n";

    // Apply frame reversal if dir_z <= 0 (matches legacy logic)
    ReferenceFrame frame1 = frame1_result.frame;
    ReferenceFrame frame2 = frame2_result.frame;
    
    if (result.dir_z <= 0.0) {
        std::cout << "\nApplying frame reversal (dir_z <= 0):\n";
        std::cout << "  Reversing y and z columns of frame2\n";
        Matrix3D rot2 = frame2.rotation();
        Vector3D y_col = rot2.column(1);
        Vector3D z_col = rot2.column(2);
        rot2.set_column(1, -y_col);
        rot2.set_column(2, -z_col);
        frame2 = ReferenceFrame(rot2, frame2.origin());
        print_frame("FRAME 2 (after reversal)", frame2);
    }

    // Calculate step parameters (using frame2 first, frame1 second - matches legacy)
    ParameterCalculator param_calc;
    BasePairStepParameters params = param_calc.calculate_step_parameters(frame2, frame1);

    print_step_parameters(params);

    // Calculate bp_type_id
    int bp_type_id = -1;
    if (result.is_valid && result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) {
        if (std::abs(params.rise) <= 2.0 && std::abs(params.twist) <= 60.0) {
            // Check wobble
            if (std::abs(params.slide) >= 1.8 && std::abs(params.slide) <= 2.8) {
                bp_type_id = 1;
            }
            // Check Watson-Crick
            if (std::abs(params.slide) <= 1.8) {
                char base1 = std::toupper(mutable_res1->one_letter_code());
                char base2 = std::toupper(mutable_res2->one_letter_code());
                std::string bp_type = std::string(1, base1) + std::string(1, base2);
                static const std::vector<std::string> WC_LIST = {"XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"};
                for (const auto& wc : WC_LIST) {
                    if (bp_type == wc) {
                        bp_type_id = 2;
                        break;
                    }
                }
            }
        }
    }

    std::cout << "\n============================================================\n";
    std::cout << "bp_type_id CALCULATION:\n";
    std::cout << "  Direction condition (dir_x>0 && dir_y<0 && dir_z<0): " 
              << ((result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) ? "PASS" : "FAIL") << "\n";
    std::cout << "  fabs(stretch) <= 2.0: " << (std::abs(params.rise) <= 2.0 ? "PASS" : "FAIL") 
              << " (value: " << std::abs(params.rise) << ")\n";
    std::cout << "  fabs(opening) <= 60.0: " << (std::abs(params.twist) <= 60.0 ? "PASS" : "FAIL") 
              << " (value: " << std::abs(params.twist) << ")\n";
    std::cout << "  fabs(shear) <= 1.8: " << (std::abs(params.slide) <= 1.8 ? "PASS" : "FAIL") 
              << " (value: " << std::abs(params.slide) << ")\n";
    std::cout << "  Base pair type: " << mutable_res1->one_letter_code() 
              << mutable_res2->one_letter_code() << "\n";
    std::cout << "  Final bp_type_id: " << bp_type_id << "\n";

    // Note: Legacy comparison can be done via Python script
    std::cout << "\n============================================================\n";
    std::cout << "NOTE: Use Python script to compare with legacy JSON data\n";

    return 0;
}

