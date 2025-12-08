/**
 * @file compare_bp_type_id_calculation.cpp
 * @brief Compare bp_type_id calculation between modern and legacy for specific pairs
 */

#include <iostream>
#include <iomanip>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/parameter_calculator.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

void analyze_bp_type_id(int idx1, int idx2, const Residue* res1, const Residue* res2, const ValidationResult& result) {
    std::cout << "\n============================================================\n";
    std::cout << "bp_type_id ANALYSIS: Pair (" << idx1 << ", " << idx2 << ")\n";
    std::cout << "============================================================\n\n";

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "DIRECTION VECTORS:\n";
    std::cout << "  dir_x: " << result.dir_x << "\n";
    std::cout << "  dir_y: " << result.dir_y << "\n";
    std::cout << "  dir_z: " << result.dir_z << "\n";

    // Check direction vector condition
    bool dir_condition = (result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0);
    std::cout << "  Condition (dir_x>0 && dir_y<0 && dir_z<0): " << (dir_condition ? "PASS" : "FAIL") << "\n";

    if (!dir_condition) {
        std::cout << "\nbp_type_id = -1 (direction vector condition not met)\n";
        return;
    }

    // Get frames
    if (!res1->reference_frame().has_value() || !res2->reference_frame().has_value()) {
        std::cout << "\nbp_type_id = -1 (frames not available)\n";
        return;
    }

    ReferenceFrame frame1 = res1->reference_frame().value();
    ReferenceFrame frame2 = res2->reference_frame().value();

    // Apply frame reversal if dir_z <= 0
    if (result.dir_z <= 0.0) {
        Matrix3D rot2 = frame2.rotation();
        Vector3D y_col = rot2.column(1);
        Vector3D z_col = rot2.column(2);
        rot2.set_column(1, -y_col);
        rot2.set_column(2, -z_col);
        frame2 = ReferenceFrame(rot2, frame2.origin());
    }

    // Calculate step parameters
    ParameterCalculator param_calc;
    BasePairStepParameters params = param_calc.calculate_step_parameters(frame2, frame1);

    double shear = params.slide;
    double stretch = params.rise;
    double opening = params.twist;

    std::cout << "\nSTEP PARAMETERS:\n";
    std::cout << "  shear (slide): " << shear << "\n";
    std::cout << "  stretch (rise): " << stretch << "\n";
    std::cout << "  opening (twist): " << opening << " degrees\n";

    std::cout << "\nTHRESHOLD CHECKS:\n";
    bool stretch_ok = (std::abs(stretch) <= 2.0);
    bool opening_ok = (std::abs(opening) <= 60.0);
    std::cout << "  fabs(stretch) <= 2.0: " << (stretch_ok ? "PASS" : "FAIL") << " (value: " << std::abs(stretch)
              << ")\n";
    std::cout << "  fabs(opening) <= 60.0: " << (opening_ok ? "PASS" : "FAIL") << " (value: " << std::abs(opening)
              << ")\n";

    if (!stretch_ok || !opening_ok) {
        std::cout << "\nbp_type_id = -1 (stretch or opening threshold exceeded)\n";
        return;
    }

    // Get base pair type
    char base1 = res1->one_letter_code();
    char base2 = res2->one_letter_code();
    std::string bp_type = std::string(1, base1) + std::string(1, base2);

    std::cout << "\nBASE PAIR TYPE:\n";
    std::cout << "  bp_type: " << bp_type << "\n";

    static const std::vector<std::string> WC_LIST = {"XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"};
    bool in_wc_list = false;
    for (const auto& wc : WC_LIST) {
        if (bp_type == wc) {
            in_wc_list = true;
            break;
        }
    }
    std::cout << "  In WC_LIST: " << (in_wc_list ? "YES" : "NO") << "\n";

    // Check wobble
    bool is_wobble = (std::abs(shear) >= 1.8 && std::abs(shear) <= 2.8);
    std::cout << "\nWOBBLE CHECK:\n";
    std::cout << "  fabs(shear) in [1.8, 2.8]: " << (is_wobble ? "YES" : "NO") << " (value: " << std::abs(shear)
              << ")\n";

    // Check Watson-Crick
    bool is_wc = (std::abs(shear) <= 1.8 && in_wc_list);
    std::cout << "\nWATSON-CRICK CHECK:\n";
    std::cout << "  fabs(shear) <= 1.8: " << (std::abs(shear) <= 1.8 ? "YES" : "NO") << " (value: " << std::abs(shear)
              << ")\n";
    std::cout << "  In WC_LIST: " << (in_wc_list ? "YES" : "NO") << "\n";
    std::cout << "  Both conditions met: " << (is_wc ? "YES" : "NO") << "\n";

    std::cout << "\nFINAL bp_type_id:\n";
    int bp_type_id = -1;
    if (is_wobble) {
        bp_type_id = 1;
        std::cout << "  bp_type_id = 1 (Wobble)\n";
    }
    if (is_wc) {
        bp_type_id = 2;
        std::cout << "  bp_type_id = 2 (Watson-Crick) - OVERWRITES wobble\n";
    }
    if (bp_type_id == -1) {
        std::cout << "  bp_type_id = -1 (Not classified)\n";
    }

    std::cout << "\nQUALITY SCORE ADJUSTMENT:\n";
    double base_quality = result.quality_score;
    std::cout << "  Base quality: " << base_quality << "\n";
    if (bp_type_id == 2) {
        double adjusted = base_quality - 2.0;
        std::cout << "  After bp_type_id=2 adjustment (-2.0): " << adjusted << "\n";
    } else {
        std::cout << "  No adjustment (bp_type_id != 2)\n";
    }
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
    std::cout << "bp_type_id Calculation Analysis Tool\n";
    std::cout << "============================================================\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "Pair: (" << idx1 << ", " << idx2 << ")\n\n";

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

    std::cout << "Residue 1 (legacy_idx=" << idx1 << "): " << mutable_res1->name() << " Chain "
              << mutable_res1->chain_id() << " Seq " << mutable_res1->seq_num() << "\n";
    std::cout << "Residue 2 (legacy_idx=" << idx2 << "): " << mutable_res2->name() << " Chain "
              << mutable_res2->chain_id() << " Seq " << mutable_res2->seq_num() << "\n";

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

    // Run validation
    BasePairValidator validator;
    ValidationResult result = validator.validate(*mutable_res1, *mutable_res2);

    // Analyze bp_type_id calculation
    analyze_bp_type_id(idx1, idx2, mutable_res1, mutable_res2, result);

    return 0;
}
