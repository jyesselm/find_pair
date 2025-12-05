/**
 * @file debug_dorg_discrepancy.cpp
 * @brief Debug tool to investigate dorg discrepancy between frame_calc.json and validation
 *
 * Issue: For pair (495, 498) involving PSU (seq=516) and C (seq=519):
 * - Frame origins in frame_calc.json are ~6.46 Å apart
 * - Validation record shows dorg = 17.25 Å
 *
 * This tool traces through the entire flow to find where the discrepancy occurs.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna;
using json = nlohmann::json;

void print_separator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << title << "\n";
    std::cout << std::string(60, '=') << "\n";
}

void print_vector(const std::string& name, const geometry::Vector3D& v) {
    std::cout << "  " << name << ": [" << std::fixed << std::setprecision(4) << v.x() << ", "
              << v.y() << ", " << v.z() << "]\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> [legacy_idx1 legacy_idx2]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 495 498\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int target_idx1 = (argc > 2) ? std::stoi(argv[2]) : 495;
    int target_idx2 = (argc > 3) ? std::stoi(argv[3]) : 498;

    std::cout << "Debug dorg discrepancy for pair (" << target_idx1 << ", " << target_idx2 << ")\n";
    std::cout << "PDB file: " << pdb_file << "\n";

    // Step 1: Parse PDB
    print_separator("STEP 1: Parse PDB and build residue mapping");

    io::PdbParser parser;
    core::Structure structure = parser.parse_file(pdb_file);

    std::cout << "Total atoms: " << structure.num_atoms() << "\n";
    std::cout << "Total chains: " << structure.chains().size() << "\n";

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

    std::cout << "Residues mapped: " << residue_by_legacy_idx.size() << "\n";
    std::cout << "Max legacy idx: " << max_legacy_idx << "\n";

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
    std::cout << "  Seq: " << res1->seq_num() << "\n";
    std::cout << "  Chain: " << res1->chain_id() << "\n";
    std::cout << "  Num atoms: " << res1->num_atoms() << "\n";
    std::cout << "  Has frame before calc: " << (res1->reference_frame().has_value() ? "YES" : "NO")
              << "\n";

    std::cout << "Residue 2 (legacy_idx=" << target_idx2 << "):\n";
    std::cout << "  Name: " << res2->name() << "\n";
    std::cout << "  Seq: " << res2->seq_num() << "\n";
    std::cout << "  Chain: " << res2->chain_id() << "\n";
    std::cout << "  Num atoms: " << res2->num_atoms() << "\n";
    std::cout << "  Has frame before calc: " << (res2->reference_frame().has_value() ? "YES" : "NO")
              << "\n";

    // Step 3: Calculate frames
    print_separator("STEP 3: Calculate frames for target residues");

    algorithms::BaseFrameCalculator calculator;

    auto frame_result1 = calculator.calculate_frame(*res1);
    std::cout << "Frame calculation for residue 1:\n";
    std::cout << "  Is valid: " << (frame_result1.is_valid ? "YES" : "NO") << "\n";
    std::cout << "  RMS fit: " << frame_result1.rms_fit << "\n";
    std::cout << "  Template: " << frame_result1.template_file << "\n";
    if (frame_result1.is_valid) {
        print_vector("Translation (frame origin)", frame_result1.translation);
    }

    auto frame_result2 = calculator.calculate_frame(*res2);
    std::cout << "Frame calculation for residue 2:\n";
    std::cout << "  Is valid: " << (frame_result2.is_valid ? "YES" : "NO") << "\n";
    std::cout << "  RMS fit: " << frame_result2.rms_fit << "\n";
    std::cout << "  Template: " << frame_result2.template_file << "\n";
    if (frame_result2.is_valid) {
        print_vector("Translation (frame origin)", frame_result2.translation);
    }

    // Step 4: Check stored frames on residue objects
    print_separator("STEP 4: Check frames stored on residue objects");

    std::cout << "Residue 1 reference_frame:\n";
    if (res1->reference_frame().has_value()) {
        auto frame1 = res1->reference_frame().value();
        print_vector("Origin", frame1.origin());
        std::cout << "  Rotation matrix:\n";
        for (int i = 0; i < 3; ++i) {
            std::cout << "    [" << frame1.rotation().at(i, 0) << ", " << frame1.rotation().at(i, 1)
                      << ", " << frame1.rotation().at(i, 2) << "]\n";
        }
    } else {
        std::cout << "  NO FRAME STORED!\n";
    }

    std::cout << "Residue 2 reference_frame:\n";
    if (res2->reference_frame().has_value()) {
        auto frame2 = res2->reference_frame().value();
        print_vector("Origin", frame2.origin());
        std::cout << "  Rotation matrix:\n";
        for (int i = 0; i < 3; ++i) {
            std::cout << "    [" << frame2.rotation().at(i, 0) << ", " << frame2.rotation().at(i, 1)
                      << ", " << frame2.rotation().at(i, 2) << "]\n";
        }
    } else {
        std::cout << "  NO FRAME STORED!\n";
    }

    // Step 5: Calculate dorg from stored frames
    print_separator("STEP 5: Calculate dorg from stored frames");

    if (res1->reference_frame().has_value() && res2->reference_frame().has_value()) {
        auto frame1 = res1->reference_frame().value();
        auto frame2 = res2->reference_frame().value();

        geometry::Vector3D diff = frame1.origin() - frame2.origin();
        double dorg = diff.length();

        std::cout << "Origin 1: [" << frame1.origin().x() << ", " << frame1.origin().y() << ", "
                  << frame1.origin().z() << "]\n";
        std::cout << "Origin 2: [" << frame2.origin().x() << ", " << frame2.origin().y() << ", "
                  << frame2.origin().z() << "]\n";
        std::cout << "Difference: [" << diff.x() << ", " << diff.y() << ", " << diff.z() << "]\n";
        std::cout << "dorg (calculated from stored frames): " << dorg << " Å\n";
    }

    // Step 6: Run validation and compare
    print_separator("STEP 6: Run validation and compare");

    algorithms::BasePairValidator validator;
    algorithms::ValidationResult result = validator.validate(*res1, *res2);

    std::cout << "Validation result:\n";
    std::cout << "  dorg (from validation): " << result.dorg << " Å\n";
    std::cout << "  d_v: " << result.d_v << "\n";
    std::cout << "  plane_angle: " << result.plane_angle << "\n";
    std::cout << "  dNN: " << result.dNN << "\n";
    std::cout << "  quality_score: " << result.quality_score << "\n";
    std::cout << "  is_valid: " << (result.is_valid ? "YES" : "NO") << "\n";
    std::cout << "\n  Validation checks:\n";
    std::cout << "    distance_check: " << (result.distance_check ? "PASS" : "FAIL") << "\n";
    std::cout << "    d_v_check: " << (result.d_v_check ? "PASS" : "FAIL") << "\n";
    std::cout << "    plane_angle_check: " << (result.plane_angle_check ? "PASS" : "FAIL") << "\n";
    std::cout << "    dNN_check: " << (result.dNN_check ? "PASS" : "FAIL") << "\n";

    // Step 7: Compare with frame_result translations
    print_separator("STEP 7: Compare calculations");

    if (frame_result1.is_valid && frame_result2.is_valid) {
        geometry::Vector3D diff_from_results =
            frame_result1.translation - frame_result2.translation;
        double dorg_from_results = diff_from_results.length();

        std::cout << "dorg calculated from frame_result translations: " << dorg_from_results
                  << " Å\n";
        std::cout << "dorg from validation: " << result.dorg << " Å\n";
        std::cout << "Difference: " << std::abs(dorg_from_results - result.dorg) << " Å\n";

        if (std::abs(dorg_from_results - result.dorg) > 0.01) {
            std::cout << "\n*** DISCREPANCY DETECTED! ***\n";
            std::cout << "The frame_result translations don't match what validation is using.\n";
            std::cout << "This suggests the frames stored on residues are different from "
                         "frame_results.\n";
        }
    }

    // Step 8: Check if there's an indexing issue
    print_separator("STEP 8: Check for indexing issues");

    // List nearby residues to check if there's an off-by-one error
    std::cout << "Residues near target indices:\n";
    for (int idx = target_idx1 - 2; idx <= target_idx2 + 2; ++idx) {
        auto it = residue_by_legacy_idx.find(idx);
        if (it != residue_by_legacy_idx.end()) {
            core::Residue* res = it->second;
            std::cout << "  idx=" << idx << ": " << res->name() << " seq=" << res->seq_num()
                      << " chain=" << res->chain_id();
            if (res->reference_frame().has_value()) {
                auto origin = res->reference_frame().value().origin();
                std::cout << " origin=[" << origin.x() << "," << origin.y() << "," << origin.z()
                          << "]";
            }
            std::cout << "\n";
        }
    }

    print_separator("INVESTIGATION COMPLETE");

    return 0;
}
