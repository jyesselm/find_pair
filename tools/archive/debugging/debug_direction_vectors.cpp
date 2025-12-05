/**
 * @file debug_direction_vectors.cpp
 * @brief Debug tool to compare direction vectors and frames between legacy and modern
 */

#include <iostream>
#include <iomanip>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 980 997\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);

    std::cout << "============================================================\n";
    std::cout << "Direction Vector Debug Tool\n";
    std::cout << "============================================================\n";
    std::cout << "PDB file: " << pdb_file << "\n";
    std::cout << "Pair: (" << idx1 << ", " << idx2 << ")\n\n";

    // Parse PDB
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);

    // Find residues by legacy index
    const Residue* res1 = nullptr;
    const Residue* res2 = nullptr;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx == idx1) {
                    res1 = &residue;
                }
                if (legacy_idx == idx2) {
                    res2 = &residue;
                }
            }
        }
    }

    if (!res1 || !res2) {
        std::cerr << "Error: Could not find residues " << idx1 << " and/or " << idx2 << "\n";
        return 1;
    }

    std::cout << "Residue 1 (legacy_idx=" << idx1 << "):\n";
    std::cout << "  Name: " << res1->name() << "\n";
    std::cout << "  Chain: " << res1->chain_id() << "\n";
    std::cout << "  Seq: " << res1->seq_num() << "\n";

    std::cout << "\nResidue 2 (legacy_idx=" << idx2 << "):\n";
    std::cout << "  Name: " << res2->name() << "\n";
    std::cout << "  Chain: " << res2->chain_id() << "\n";
    std::cout << "  Seq: " << res2->seq_num() << "\n";

    // Calculate frames
    BaseFrameCalculator calculator("data/templates");
    bool is_rna = false;
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == " O2'") {
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
    // Note: BaseFrameCalculator detects RNA automatically from templates

    // Calculate frames for both residues
    Structure mutable_structure = structure; // Need mutable for calculate_frame
    Residue* mutable_res1 = nullptr;
    Residue* mutable_res2 = nullptr;

    for (auto& chain : mutable_structure.chains()) {
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
        std::cerr << "Error: Could not find mutable residues\n";
        return 1;
    }

    std::cout << "\n============================================================\n";
    std::cout << "STEP 1: Calculate frames\n";
    std::cout << "============================================================\n";

    FrameCalculationResult frame1_result = calculator.calculate_frame(*mutable_res1);
    FrameCalculationResult frame2_result = calculator.calculate_frame(*mutable_res2);

    if (!frame1_result.is_valid || !frame2_result.is_valid) {
        std::cerr << "Error: Frame calculation failed\n";
        return 1;
    }

    // Store frames on residues (required for validation)
    mutable_res1->set_reference_frame(frame1_result.frame);
    mutable_res2->set_reference_frame(frame2_result.frame);

    ReferenceFrame frame1 = frame1_result.frame;
    ReferenceFrame frame2 = frame2_result.frame;

    std::cout << "Frame 1 origin: [" << frame1.origin().x() << ", " << frame1.origin().y() << ", "
              << frame1.origin().z() << "]\n";
    std::cout << "Frame 1 rotation matrix:\n";
    Matrix3D rot1 = frame1.rotation();
    for (int i = 0; i < 3; i++) {
        Vector3D col = rot1.column(i);
        std::cout << "  [" << col.x() << ", " << col.y() << ", " << col.z() << "]\n";
    }

    std::cout << "\nFrame 2 origin: [" << frame2.origin().x() << ", " << frame2.origin().y() << ", "
              << frame2.origin().z() << "]\n";
    std::cout << "Frame 2 rotation matrix:\n";
    Matrix3D rot2 = frame2.rotation();
    for (int i = 0; i < 3; i++) {
        Vector3D col = rot2.column(i);
        std::cout << "  [" << col.x() << ", " << col.y() << ", " << col.z() << "]\n";
    }

    std::cout << "\n============================================================\n";
    std::cout << "STEP 2: Calculate direction vectors\n";
    std::cout << "============================================================\n";

    // Calculate direction vectors (i, j) order
    double dir_x_ij, dir_y_ij, dir_z_ij;
    dir_x_ij = frame1.x_axis().dot(frame2.x_axis());
    dir_y_ij = frame1.y_axis().dot(frame2.y_axis());
    dir_z_ij = frame1.z_axis().dot(frame2.z_axis());

    std::cout << "Direction vectors (res1, res2) order:\n";
    std::cout << "  dir_x = " << std::fixed << std::setprecision(6) << dir_x_ij << "\n";
    std::cout << "  dir_y = " << dir_y_ij << "\n";
    std::cout << "  dir_z = " << dir_z_ij << "\n";

    // Calculate direction vectors (j, i) order
    double dir_x_ji, dir_y_ji, dir_z_ji;
    dir_x_ji = frame2.x_axis().dot(frame1.x_axis());
    dir_y_ji = frame2.y_axis().dot(frame1.y_axis());
    dir_z_ji = frame2.z_axis().dot(frame1.z_axis());

    std::cout << "\nDirection vectors (res2, res1) order:\n";
    std::cout << "  dir_x = " << dir_x_ji << "\n";
    std::cout << "  dir_y = " << dir_y_ji << "\n";
    std::cout << "  dir_z = " << dir_z_ji << "\n";

    std::cout << "\nNote: Direction vectors should be symmetric (same for both orders)\n";

    std::cout << "\n============================================================\n";
    std::cout << "STEP 3: Run validation\n";
    std::cout << "============================================================\n";

    BasePairValidator validator;
    // Use mutable_res1 and mutable_res2 which have frames stored
    ValidationResult result_ij = validator.validate(*mutable_res1, *mutable_res2);
    ValidationResult result_ji = validator.validate(*mutable_res2, *mutable_res1);

    std::cout << "Validation (res1, res2):\n";
    std::cout << "  dir_x = " << result_ij.dir_x << "\n";
    std::cout << "  dir_y = " << result_ij.dir_y << "\n";
    std::cout << "  dir_z = " << result_ij.dir_z << "\n";
    std::cout << "  is_valid = " << (result_ij.is_valid ? "YES" : "NO") << "\n";

    std::cout << "\nValidation (res2, res1):\n";
    std::cout << "  dir_x = " << result_ji.dir_x << "\n";
    std::cout << "  dir_y = " << result_ji.dir_y << "\n";
    std::cout << "  dir_z = " << result_ji.dir_z << "\n";
    std::cout << "  is_valid = " << (result_ji.is_valid ? "YES" : "NO") << "\n";

    std::cout << "\n============================================================\n";
    std::cout << "STEP 4: Compare with legacy JSON\n";
    std::cout << "============================================================\n";
    std::cout << "Check data/json_legacy/pair_validation/<PDB_ID>.json for legacy values\n";

    return 0;
}
