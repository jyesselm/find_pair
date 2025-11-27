/**
 * @file compare_quality_score_components.cpp
 * @brief Compare all components of quality score calculation between legacy and modern
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <nlohmann/json.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;
using json = nlohmann::json;

void print_quality_breakdown(const ValidationResult& result, int idx1, int idx2) {
    std::cout << "\n============================================================\n";
    std::cout << "QUALITY SCORE BREAKDOWN: Pair (" << idx1 << ", " << idx2 << ")\n";
    std::cout << "============================================================\n\n";
    
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "BASE QUALITY SCORE COMPONENTS:\n";
    std::cout << "  dorg:        " << result.dorg << " Å\n";
    std::cout << "  d_v:          " << result.d_v << " Å\n";
    std::cout << "  plane_angle:  " << result.plane_angle << " degrees\n";
    std::cout << "  Base formula: dorg + 2.0 * d_v + plane_angle / 20.0\n";
    double base_quality = result.dorg + 2.0 * result.d_v + result.plane_angle / 20.0;
    std::cout << "  Base quality: " << base_quality << "\n";
    std::cout << "  Recorded quality_score: " << result.quality_score << "\n";
    
    std::cout << "\nHYDROGEN BONDS:\n";
    std::cout << "  num_base_hb: " << result.num_base_hb << "\n";
    std::cout << "  num_o2_hb:   " << result.num_o2_hb << "\n";
    std::cout << "  Total H-bonds: " << result.hbonds.size() << "\n";
    if (!result.hbonds.empty()) {
        std::cout << "  H-bond details:\n";
        int num_good_hb = 0; // Count H-bonds with type='-' and distance in [2.5, 3.5]
        for (size_t i = 0; i < result.hbonds.size() && i < 10; ++i) {
            const auto& hb = result.hbonds[i];
            bool is_good = (hb.type == '-' && hb.distance >= 2.5 && hb.distance <= 3.5);
            if (is_good) {
                num_good_hb++;
            }
            std::cout << "    " << (i+1) << ". " << hb.donor_atom << " -> " << hb.acceptor_atom 
                      << " (distance: " << hb.distance << " Å, type: '" << hb.type << "'";
            if (is_good) {
                std::cout << " [GOOD - counts for adjust_pairQuality]";
            }
            std::cout << ")\n";
        }
        if (result.hbonds.size() > 10) {
            std::cout << "    ... and " << (result.hbonds.size() - 10) << " more\n";
        }
        std::cout << "  Good H-bonds (type='-' and 2.5-3.5 Å): " << num_good_hb << "\n";
        double adjust_pairQuality = (num_good_hb >= 2) ? -3.0 : -static_cast<double>(num_good_hb);
        std::cout << "  adjust_pairQuality: " << adjust_pairQuality << "\n";
        std::cout << "  Adjusted quality (base + adjust): " << (result.quality_score + adjust_pairQuality) << "\n";
    }
    
    std::cout << "\nVALIDATION CHECKS:\n";
    std::cout << "  distance_check:    " << (result.distance_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  d_v_check:          " << (result.d_v_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  plane_angle_check:  " << (result.plane_angle_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  dNN_check:          " << (result.dNN_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  overlap_check:      " << (result.overlap_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  hbond_check:        " << (result.hbond_check ? "PASS" : "FAIL") << "\n";
    std::cout << "  is_valid:           " << (result.is_valid ? "YES" : "NO") << "\n";
    
    std::cout << "\nBASE PAIR TYPE:\n";
    std::cout << "  bp_type: " << static_cast<int>(result.bp_type) << "\n";
    
    std::cout << "\nDIRECTION VECTORS:\n";
    std::cout << "  dir_x: " << result.dir_x << "\n";
    std::cout << "  dir_y: " << result.dir_y << "\n";
    std::cout << "  dir_z: " << result.dir_z << "\n";
    
    std::cout << "\nGEOMETRIC PARAMETERS:\n";
    std::cout << "  dNN: " << result.dNN << " Å\n";
    std::cout << "  overlap_area: " << result.overlap_area << " Å²\n";
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <residue1_idx> <residue2_idx>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 968 1024\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int idx1 = std::stoi(argv[2]);
    int idx2 = std::stoi(argv[3]);

    std::cout << "============================================================\n";
    std::cout << "Quality Score Component Analysis Tool\n";
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

    std::cout << "Residue 1 (legacy_idx=" << idx1 << "): " << mutable_res1->name() 
              << " Chain " << mutable_res1->chain_id() << " Seq " << mutable_res1->seq_num() << "\n";
    std::cout << "Residue 2 (legacy_idx=" << idx2 << "): " << mutable_res2->name() 
              << " Chain " << mutable_res2->chain_id() << " Seq " << mutable_res2->seq_num() << "\n";

    // Calculate and store frames
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
            if (is_rna) break;
        }
        if (is_rna) break;
    }

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

    // Print breakdown
    print_quality_breakdown(result, idx1, idx2);

    // Now calculate adjusted quality score (as done in find_best_partner)
    std::cout << "\n============================================================\n";
    std::cout << "ADJUSTED QUALITY SCORE CALCULATION\n";
    std::cout << "============================================================\n\n";
    
    // Calculate bp_type_id (simplified - would need full BasePairFinder context)
    std::cout << "Note: Full adjusted quality score calculation requires BasePairFinder context\n";
    std::cout << "This includes:\n";
    std::cout << "  1. adjust_pairQuality() based on H-bonds\n";
    std::cout << "  2. bp_type_id calculation (requires step parameters)\n";
    std::cout << "  3. bp_type_id == 2 adjustment (-2.0)\n";
    
    return 0;
}

