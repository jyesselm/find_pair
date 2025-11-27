/**
 * @file trace_pair_selection.cpp
 * @brief Trace the pair selection process for specific residues
 */

#include <iostream>
#include <iomanip>
#include <x3dna/core/structure.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/io/json_writer.hpp>

using namespace x3dna::core;
using namespace x3dna::io;
using namespace x3dna::algorithms;
using namespace x3dna::geometry;

void trace_residue_selection(int legacy_idx,
                             const std::map<int, const Residue*>& residue_by_legacy_idx,
                             int max_legacy_idx) {
    std::cout << "\n============================================================\n";
    std::cout << "TRACING: Residue " << legacy_idx << " selection\n";
    std::cout << "============================================================\n\n";
    
    // Simulate find_best_partner for this residue
    std::vector<bool> matched_indices(max_legacy_idx + 1, false);
    std::map<std::pair<int, int>, ValidationResult> phase1_validation_results;
    
    auto it = residue_by_legacy_idx.find(legacy_idx);
    if (it == residue_by_legacy_idx.end() || !it->second) {
        std::cout << "Residue " << legacy_idx << " not found\n";
        return;
    }
    
    const Residue* res1 = it->second;
    BasePairValidator validator;
    
    std::cout << "Residue " << legacy_idx << ": " << res1->name() 
              << " Chain " << res1->chain_id() << " Seq " << res1->seq_num() << "\n\n";
    
    std::cout << "Checking all potential partners:\n";
    std::cout << std::fixed << std::setprecision(6);
    
    double best_score = std::numeric_limits<double>::max();
    int best_partner = -1;
    
    // Check all other residues
    for (int legacy_idx2 = 1; legacy_idx2 <= max_legacy_idx; ++legacy_idx2) {
        if (legacy_idx2 == legacy_idx || matched_indices[legacy_idx2]) {
            continue;
        }
        
        auto it2 = residue_by_legacy_idx.find(legacy_idx2);
        if (it2 == residue_by_legacy_idx.end() || !it2->second) {
            continue;
        }
        
        const Residue* res2 = it2->second;
        
        // Validate pair
        ValidationResult result;
        if (legacy_idx < legacy_idx2) {
            result = validator.validate(*res1, *res2);
        } else {
            result = validator.validate(*res2, *res1);
        }
        
        if (!result.is_valid) {
            continue;
        }
        
        // Calculate adjusted quality score (simplified - without bp_type_id adjustment)
        double quality_adjustment = 0.0;
        int num_good_hb = 0;
        for (const auto& hb : result.hbonds) {
            if (hb.type == '-' && hb.distance >= 2.5 && hb.distance <= 3.5) {
                num_good_hb++;
            }
        }
        if (num_good_hb >= 2) {
            quality_adjustment = -3.0;
        } else {
            quality_adjustment = -static_cast<double>(num_good_hb);
        }
        
        double adjusted_quality = result.quality_score + quality_adjustment;
        
        std::cout << "  Partner " << legacy_idx2 << " (" << res2->name() << "): "
                  << "base=" << result.quality_score << ", adjust=" << quality_adjustment
                  << ", adjusted=" << adjusted_quality;
        
        if (adjusted_quality < best_score) {
            best_score = adjusted_quality;
            best_partner = legacy_idx2;
            std::cout << " [NEW BEST]";
        }
        std::cout << "\n";
    }
    
    std::cout << "\nBEST PARTNER: " << best_partner;
    if (best_partner > 0) {
        std::cout << " (adjusted quality: " << best_score << ")";
    }
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file>\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb\n";
        return 1;
    }

    std::string pdb_file = argv[1];

    std::cout << "============================================================\n";
    std::cout << "Pair Selection Tracing Tool\n";
    std::cout << "============================================================\n";
    std::cout << "PDB file: " << pdb_file << "\n\n";

    // Parse PDB
    PdbParser parser;
    parser.set_include_hetatm(true);
    parser.set_include_waters(true);
    Structure structure = parser.parse_file(pdb_file);

    // Build residue map and calculate frames
    std::map<int, const Residue*> residue_by_legacy_idx;
    int max_legacy_idx = 0;
    
    BaseFrameCalculator calculator("data/templates");
    
    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    // Calculate and store frame
                    FrameCalculationResult frame_result = calculator.calculate_frame(residue);
                    if (frame_result.is_valid) {
                        residue.set_reference_frame(frame_result.frame);
                    }
                    residue_by_legacy_idx[legacy_idx] = &residue;
                    if (legacy_idx > max_legacy_idx) {
                        max_legacy_idx = legacy_idx;
                    }
                }
            }
        }
    }

    // Trace selection for key residues
    std::vector<int> residues_to_trace = {968, 980, 998, 1024, 1188};
    
    for (int legacy_idx : residues_to_trace) {
        if (residue_by_legacy_idx.find(legacy_idx) != residue_by_legacy_idx.end()) {
            trace_residue_selection(legacy_idx, residue_by_legacy_idx, max_legacy_idx);
        }
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "MUTUAL BEST MATCH ANALYSIS\n";
    std::cout << "============================================================\n\n";
    std::cout << "For a pair to be selected, both residues must select each other\n";
    std::cout << "as their best partner (mutual best match).\n";

    return 0;
}

