/**
 * @file investigate_missing_pairs.cpp
 * @brief Investigate why specific pairs are not found in validation
 *
 * This tool checks why pairs that exist in legacy validation are missing
 * from modern validation. It checks:
 * 1. Residue recognition (is_nucleotide)
 * 2. Frame availability
 * 3. Residue indices
 * 4. Early rejection reasons
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/base_pair_finder.hpp>
#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/io/residue_index_fixer.hpp>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>

using json = nlohmann::json;

struct PairInfo {
    int idx1;
    int idx2;
};

void check_residue(const x3dna::core::Residue& residue, int legacy_idx, const std::string& label) {
    std::cout << "  " << label << " (legacy_idx=" << legacy_idx << "):\n";
    std::cout << "    ResName: " << residue.name() << "\n";
    std::cout << "    ChainID: " << residue.chain_id() << "\n";
    std::cout << "    ResSeq: " << residue.seq_num() << "\n";
    char ins_code = residue.insertion();
    std::cout << "    Insertion: " << (ins_code != ' ' ? std::string(1, ins_code) : "none") << "\n";
    std::cout << "    Num atoms: " << residue.atoms().size() << "\n";

    // Check if nucleotide using BasePairFinder's is_nucleotide function (static)
    bool is_nuc = x3dna::algorithms::BasePairFinder::is_nucleotide(residue);
    std::cout << "    Is nucleotide: " << (is_nuc ? "YES" : "NO") << "\n";
    std::cout << "    ResidueType: " << static_cast<int>(residue.residue_type()) << "\n";

    // Check frame
    if (residue.reference_frame().has_value()) {
        std::cout << "    Frame: AVAILABLE\n";
        auto frame = residue.reference_frame().value();
        std::cout << "      Origin: (" << frame.origin().x() << ", " << frame.origin().y() << ", "
                  << frame.origin().z() << ")\n";
    } else {
        std::cout << "    Frame: MISSING\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <pdb_file> <legacy_idx1> <legacy_idx2> [legacy_json_file]\n";
        std::cerr << "Example: " << argv[0] << " data/pdb/6CAQ.pdb 495 498\n";
        std::cerr << "         " << argv[0]
                  << " data/pdb/6CAQ.pdb 495 498 data/json_legacy/base_frame_calc/6CAQ.json\n";
        return 1;
    }

    std::string pdb_file = argv[1];
    int legacy_idx1 = std::stoi(argv[2]);
    int legacy_idx2 = std::stoi(argv[3]);
    std::string legacy_json_file = (argc > 4) ? argv[4] : "";

    std::cout << "Investigating pair (" << legacy_idx1 << ", " << legacy_idx2 << ") in " << pdb_file
              << "\n";
    std::cout << "=" << 60 << "\n\n";

    // Parse PDB
    x3dna::io::PdbParser parser;
    auto structure = parser.parse_file(pdb_file);

    // Fix indices if legacy JSON provided
    if (!legacy_json_file.empty() && std::filesystem::exists(legacy_json_file)) {
        std::cout << "Fixing residue indices from: " << legacy_json_file << "\n";
        x3dna::io::fix_residue_indices_from_json(structure, legacy_json_file);
    }

    // Build residue map
    std::map<int, const x3dna::core::Residue*> residue_by_legacy_idx;
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residue_by_legacy_idx[legacy_idx] = &residue;
                }
            }
        }
    }

    // Find residues
    auto it1 = residue_by_legacy_idx.find(legacy_idx1);
    auto it2 = residue_by_legacy_idx.find(legacy_idx2);

    if (it1 == residue_by_legacy_idx.end()) {
        std::cout << "ERROR: Residue " << legacy_idx1 << " not found!\n";
        std::cout << "Available indices: ";
        int count = 0;
        for (const auto& [idx, res] : residue_by_legacy_idx) {
            if (count++ < 10) {
                std::cout << idx << " ";
            }
        }
        std::cout << "...\n";
        return 1;
    }

    if (it2 == residue_by_legacy_idx.end()) {
        std::cout << "ERROR: Residue " << legacy_idx2 << " not found!\n";
        return 1;
    }

    const x3dna::core::Residue* res1 = it1->second;
    const x3dna::core::Residue* res2 = it2->second;

    std::cout << "1. Residue Information:\n";
    std::cout << "-" << 60 << "\n";
    check_residue(*res1, legacy_idx1, "Residue 1");
    std::cout << "\n";
    check_residue(*res2, legacy_idx2, "Residue 2");
    std::cout << "\n";

    // Check if both are nucleotides using BasePairFinder's is_nucleotide function (static)
    bool res1_is_nuc = x3dna::algorithms::BasePairFinder::is_nucleotide(*res1);
    bool res2_is_nuc = x3dna::algorithms::BasePairFinder::is_nucleotide(*res2);

    std::cout << "2. Validation Checks:\n";
    std::cout << "-" << 60 << "\n";

    if (!res1_is_nuc) {
        std::cout << "❌ Residue 1 is NOT recognized as nucleotide - pair will be skipped\n";
    }
    if (!res2_is_nuc) {
        std::cout << "❌ Residue 2 is NOT recognized as nucleotide - pair will be skipped\n";
    }
    if (!res1->reference_frame().has_value()) {
        std::cout << "❌ Residue 1 has NO frame - pair will be skipped\n";
    }
    if (!res2->reference_frame().has_value()) {
        std::cout << "❌ Residue 2 has NO frame - pair will be skipped\n";
    }

    if (res1_is_nuc && res2_is_nuc && res1->reference_frame().has_value() &&
        res2->reference_frame().has_value()) {
        std::cout << "✅ Both residues are nucleotides with frames - attempting validation\n\n";

        // Try validation
        x3dna::algorithms::BasePairValidator validator;
        auto result = validator.validate(*res1, *res2);

        std::cout << "3. Validation Result:\n";
        std::cout << "-" << 60 << "\n";
        std::cout << "  is_valid: " << (result.is_valid ? "YES" : "NO") << "\n";
        std::cout << "  dorg: " << result.dorg << "\n";
        std::cout << "  d_v: " << result.d_v << "\n";
        std::cout << "  dNN: " << result.dNN << "\n";
        std::cout << "  plane_angle: " << result.plane_angle << "\n";
        std::cout << "  overlap_area: " << result.overlap_area << "\n";
        std::cout << "  num_base_hb: " << result.num_base_hb << "\n";
        std::cout << "  quality_score: " << result.quality_score << "\n";

        if (!result.is_valid) {
            std::cout << "\n  ❌ Pair FAILED validation\n";
            std::cout << "  This explains why it's not in validation records\n";
        } else {
            std::cout << "\n  ✅ Pair PASSED validation\n";
            std::cout
                << "  ⚠️  But it's missing from validation records - check Phase 1 iteration\n";
        }
    } else {
        std::cout << "\n❌ Pair cannot be validated due to missing requirements\n";
    }

    return 0;
}
