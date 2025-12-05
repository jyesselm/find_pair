/**
 * @file debug_protocol_6v9q.cpp
 * @brief Debug tool to investigate why FindPairProtocol finds 0 pairs for 6V9Q
 */

#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/core/structure.hpp>
#include <iostream>
#include <filesystem>

using namespace x3dna;

int main(int argc, char* argv[]) {
    std::string pdb_file = "data/pdb/6V9Q.pdb";
    if (argc > 1) {
        pdb_file = argv[1];
    }

    std::cout << "=== Debugging FindPairProtocol for " << pdb_file << " ===" << std::endl;

    // Parse PDB
    io::PdbParser parser;
    core::Structure structure;
    try {
        structure = parser.parse_file(pdb_file);
        std::cout << "✓ PDB parsed successfully" << std::endl;
        std::cout << "  Chains: " << structure.num_chains() << std::endl;
        std::cout << "  Residues: " << structure.num_residues() << std::endl;
        std::cout << "  Atoms: " << structure.num_atoms() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "✗ Failed to parse PDB: " << e.what() << std::endl;
        return 1;
    }

    // Check legacy residue indices
    size_t residues_with_legacy_idx = 0;
    size_t residues_with_frames = 0;
    size_t nucleotide_residues = 0;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (!residue.atoms().empty()) {
                int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                if (legacy_idx > 0) {
                    residues_with_legacy_idx++;
                }
            }
            if (residue.reference_frame().has_value()) {
                residues_with_frames++;
            }
            // Check if nucleotide (includes modified nucleotides via is_nucleotide check)
            // For modified nucleotides, we check for nitrogen atoms
            bool is_nuc = residue.is_nucleotide();
            if (!is_nuc) {
                // Check for modified nucleotides by looking for nitrogen atoms
                for (const auto& atom : residue.atoms()) {
                    std::string name = atom.name();
                    if (name == " N1" || name == " N9" || name == " N3") {
                        is_nuc = true;
                        break;
                    }
                }
            }
            if (is_nuc) {
                nucleotide_residues++;
            }
        }
    }

    std::cout << "\n=== Structure Analysis ===" << std::endl;
    std::cout << "Residues with legacy_idx > 0: " << residues_with_legacy_idx << std::endl;
    std::cout << "Residues with frames: " << residues_with_frames << std::endl;
    std::cout << "Nucleotide residues: " << nucleotide_residues << std::endl;

    // Create protocol
    protocols::FindPairProtocol protocol;
    auto& config = config::ConfigManager::instance();
    config.set_defaults(); // Ensure defaults are set
    protocol.set_config_manager(config);

    // Print config parameters
    std::cout << "\n=== Config Parameters ===" << std::endl;
    const auto& thresholds = config.thresholds();
    std::cout << "min_dorg: " << thresholds.min_dorg << std::endl;
    std::cout << "max_dorg: " << thresholds.max_dorg << std::endl;
    std::cout << "min_dv: " << thresholds.min_dv << std::endl;
    std::cout << "max_dv: " << thresholds.max_dv << std::endl;
    std::cout << "min_dNN: " << thresholds.min_dNN << std::endl;
    std::cout << "max_dNN: " << thresholds.max_dNN << std::endl;
    std::cout << "min_plane_angle: " << thresholds.min_plane_angle << std::endl;
    std::cout << "max_plane_angle: " << thresholds.max_plane_angle << std::endl;
    std::cout << "min_base_hb: " << thresholds.min_base_hb << std::endl;
    std::cout << "hb_lower: " << thresholds.hb_lower << std::endl;
    std::cout << "hb_dist1: " << thresholds.hb_dist1 << std::endl;
    std::cout << "overlap_threshold: " << thresholds.overlap_threshold << std::endl;

    // Execute protocol
    std::cout << "\n=== Executing Protocol ===" << std::endl;
    try {
        protocol.execute(structure);
        std::cout << "✓ Protocol executed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "✗ Protocol execution failed: " << e.what() << std::endl;
        return 1;
    }

    // Check results
    const auto& base_pairs = protocol.base_pairs();
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Base pairs found: " << base_pairs.size() << std::endl;

    if (base_pairs.empty()) {
        std::cout << "\n⚠ No base pairs found! Investigating..." << std::endl;

        // Check frames after protocol execution
        size_t residues_with_frames_after = 0;
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                if (residue.reference_frame().has_value()) {
                    residues_with_frames_after++;
                }
            }
        }
        std::cout << "Residues with frames after protocol: " << residues_with_frames_after
                  << std::endl;

        // Check if any residues are nucleotides with frames and legacy_idx
        size_t nuc_with_frames = 0;
        std::vector<int> legacy_indices_with_frames;
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                bool is_nuc = residue.is_nucleotide();
                if (!is_nuc) {
                    // Check for modified nucleotides
                    for (const auto& atom : residue.atoms()) {
                        std::string name = atom.name();
                        if (name == " N1" || name == " N9" || name == " N3") {
                            is_nuc = true;
                            break;
                        }
                    }
                }
                if (is_nuc && residue.reference_frame().has_value()) {
                    nuc_with_frames++;
                    if (!residue.atoms().empty()) {
                        int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                        if (legacy_idx > 0) {
                            legacy_indices_with_frames.push_back(legacy_idx);
                        }
                    }
                }
            }
        }
        std::cout << "Nucleotide residues with frames: " << nuc_with_frames << std::endl;
        std::cout << "Legacy indices with frames (first 20): ";
        std::sort(legacy_indices_with_frames.begin(), legacy_indices_with_frames.end());
        for (size_t i = 0; i < std::min(legacy_indices_with_frames.size(), size_t(20)); ++i) {
            std::cout << legacy_indices_with_frames[i] << " ";
        }
        std::cout << std::endl;
        if (!legacy_indices_with_frames.empty()) {
            std::cout << "Max legacy index with frame: " << legacy_indices_with_frames.back()
                      << std::endl;
        }

        // Find max legacy index overall
        int max_legacy_idx = 0;
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                if (!residue.atoms().empty()) {
                    int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                    if (legacy_idx > max_legacy_idx) {
                        max_legacy_idx = legacy_idx;
                    }
                }
            }
        }
        std::cout << "Max legacy index in structure: " << max_legacy_idx << std::endl;

        // Check legacy expected pairs
        std::cout << "\nExpected legacy pairs (from data/json_legacy/base_pair/6V9Q.json):"
                  << std::endl;
        std::cout << "  Residue indices: 42, 44, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60"
                  << std::endl;
        std::cout << "  Checking if these have frames and legacy_idx..." << std::endl;
        for (int expected_idx : {42, 44, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60}) {
            bool found = false;
            for (const auto& chain : structure.chains()) {
                for (const auto& residue : chain.residues()) {
                    if (!residue.atoms().empty()) {
                        int legacy_idx = residue.atoms()[0].legacy_residue_idx();
                        if (legacy_idx == expected_idx) {
                            found = true;
                            bool has_frame = residue.reference_frame().has_value();
                            bool is_nuc = residue.is_nucleotide();
                            std::cout << "  Residue " << expected_idx << ": found=" << found
                                      << ", has_frame=" << has_frame << ", is_nucleotide=" << is_nuc
                                      << std::endl;
                            break;
                        }
                    }
                }
                if (found)
                    break;
            }
            if (!found) {
                std::cout << "  Residue " << expected_idx << ": NOT FOUND" << std::endl;
            }
        }
    } else {
        std::cout << "\n✓ Found " << base_pairs.size() << " base pairs:" << std::endl;
        for (size_t i = 0; i < base_pairs.size(); ++i) {
            const auto& pair = base_pairs[i];
            // BasePair uses 0-based indices that are legacy_idx - 1
            // So to get legacy indices, we add 1
            int legacy_idx1 = static_cast<int>(pair.residue_idx1()) + 1;
            int legacy_idx2 = static_cast<int>(pair.residue_idx2()) + 1;
            std::cout << "  " << (i + 1) << ": legacy " << legacy_idx1 << " <-> " << legacy_idx2
                      << " (bp_type=" << pair.bp_type() << ")" << std::endl;
        }
    }

    return 0;
}
