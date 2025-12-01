/**
 * @file analyze_protocol.cpp
 * @brief AnalyzeProtocol implementation
 */

#include <x3dna/protocols/analyze_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/structure_legacy_order.hpp>
#include <iostream>
#include <algorithm>
#include <map>

namespace x3dna {
namespace protocols {

AnalyzeProtocol::AnalyzeProtocol(const std::filesystem::path& template_path)
    : frame_calculator_(template_path)
{
}

void AnalyzeProtocol::execute(const std::filesystem::path& input_file) {
    // Step 1: Parse .inp file
    input_data_ = io::InputFileParser::parse(input_file);

    // Step 2: Load PDB structure
    core::Structure structure = load_structure(input_data_.pdb_file);

    // Step 3: Store base pairs from input file
    base_pairs_ = input_data_.base_pairs;

    // Step 3.5: Convert atom indices to residue indices if needed
    // Legacy input files may contain atom indices instead of residue indices
    convert_atom_indices_to_residue_indices(structure);

    // Step 4: Execute protocol on structure
    execute(structure);
}

void AnalyzeProtocol::execute(core::Structure& structure) {
    // Check legacy mode from config if not explicitly set
    if (!legacy_mode_ && config_) {
        legacy_mode_ = config_->legacy_mode();
    }

    // Step 1: Recalculate frames for all residues
    recalculate_frames(structure);

    // Step 2: Calculate step and helical parameters
    calculate_parameters(structure);
}

void AnalyzeProtocol::recalculate_frames(core::Structure& structure) {
    // Check if frames already exist on residues (from find_pair phase)
    // Only recalculate if frames are missing
    bool need_recalculate = false;
    size_t frames_found = 0;
    size_t frames_missing = 0;
    
    // First pass: check which frames exist
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (residue.reference_frame().has_value()) {
                frames_found++;
            } else {
                frames_missing++;
                need_recalculate = true;
            }
        }
    }
    
    // Store original frames for verification if they exist
    std::map<int, core::ReferenceFrame> original_frames;
    if (frames_found > 0) {
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                if (residue.reference_frame().has_value()) {
                    int legacy_idx = structure.get_legacy_idx_for_residue(&residue);
                    if (legacy_idx > 0) {
                        original_frames[legacy_idx] = residue.reference_frame().value();
                    }
                }
            }
        }
    }
    
    if (need_recalculate) {
        // Some frames are missing, recalculate all frames
        // Detect RNA by checking for O2' atoms
        bool is_rna = false;
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " O2'") {
                        is_rna = true;
                        break;
                    }
                }
                if (is_rna) break;
            }
            if (is_rna) break;
        }

        frame_calculator_.set_is_rna(is_rna);
        frame_calculator_.calculate_all_frames(structure);
        
        if (frames_found > 0) {
            std::cout << "Recalculated " << frames_missing << " missing frames "
                      << "(reused " << frames_found << " existing frames from find_pair)\n";
        } else {
            std::cout << "Calculated frames for all " << frames_missing << " residues\n";
        }
    } else {
        // All frames already exist - reuse frames from find_pair phase
        std::cout << "Reusing " << frames_found << " frames from find_pair phase\n";
    }

    // Verify frames match (if we had original frames)
    if (!original_frames.empty()) {
        size_t frames_verified = 0;
        size_t frames_differ = 0;
        const double tolerance = 1e-6;
        
        for (const auto& chain : structure.chains()) {
            for (const auto& residue : chain.residues()) {
                if (residue.reference_frame().has_value()) {
                    int legacy_idx = structure.get_legacy_idx_for_residue(&residue);
                    if (legacy_idx > 0) {
                        auto it = original_frames.find(legacy_idx);
                        if (it != original_frames.end()) {
                            // Compare frames
                            const auto& original = it->second;
                            // Store current frame to avoid temporary
                            auto current_frame_opt = residue.reference_frame();
                            if (!current_frame_opt.has_value()) {
                                continue;
                            }
                            const auto& current = current_frame_opt.value();
                            
                            // Compare origins
                            auto org_diff = current.origin() - original.origin();
                            double org_dist = std::sqrt(org_diff.x() * org_diff.x() + 
                                                       org_diff.y() * org_diff.y() + 
                                                       org_diff.z() * org_diff.z());
                            
                            // Compare rotation matrices (element-wise)
                            const auto& rot_orig = original.rotation();
                            const auto& rot_curr = current.rotation();
                            double max_rot_diff = 0.0;
                            for (int i = 0; i < 3; ++i) {
                                for (int j = 0; j < 3; ++j) {
                                    double diff = std::abs(rot_curr.at(i, j) - rot_orig.at(i, j));
                                    if (diff > max_rot_diff) {
                                        max_rot_diff = diff;
                                    }
                                }
                            }
                            
                            if (org_dist < tolerance && max_rot_diff < tolerance) {
                                frames_verified++;
                            } else {
                                frames_differ++;
                                if (frames_differ <= 5) {  // Report first 5 differences
                                    std::cerr << "Warning: Frame mismatch for residue " << legacy_idx << "\n";
                                    std::cerr << "  Origin diff: " << org_dist << " (tolerance: " << tolerance << ")\n";
                                    std::cerr << "  Rotation max diff: " << max_rot_diff << " (tolerance: " << tolerance << ")\n";
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if (frames_verified > 0) {
            std::cout << "Verified " << frames_verified << " frames match find_pair phase";
            if (frames_differ > 0) {
                std::cout << " (" << frames_differ << " differ)";
            }
            std::cout << "\n";
        }
    }

    // Update base pairs with frames from residues
    // For each base pair, get the frames from the corresponding residues
    // The base pair indices from .inp file are 0-based but represent legacy residue indices
    // So we convert back to 1-based to use get_residue_by_legacy_idx()
    for (auto& pair : base_pairs_) {
        // Convert from 0-based to 1-based legacy index
        int legacy_idx1 = static_cast<int>(pair.residue_idx1() + 1);
        int legacy_idx2 = static_cast<int>(pair.residue_idx2() + 1);

        const auto* res1 = core::get_residue_by_legacy_idx(structure, legacy_idx1);
        const auto* res2 = core::get_residue_by_legacy_idx(structure, legacy_idx2);

        if (res1 && res1->reference_frame().has_value()) {
            pair.set_frame1(res1->reference_frame().value());
        }

        if (res2 && res2->reference_frame().has_value()) {
            pair.set_frame2(res2->reference_frame().value());
        }
    }
}

void AnalyzeProtocol::calculate_parameters(core::Structure& /* structure */) {
    // Clear previous results
    step_parameters_.clear();
    helical_parameters_.clear();

    if (base_pairs_.size() < 2) {
        // Need at least 2 pairs to calculate step parameters
        return;
    }

    // Determine range of pairs to process (based on -S option)
    // step_start_ and step_size_ are 1-based (matching legacy)
    // Legacy: for (j = istart; j <= nbpm1; j += istep) where j is 1-based base pair index
    size_t start_idx = (step_start_ > 0) ? (step_start_ - 1) : 0;  // Convert to 0-based
    if (start_idx >= base_pairs_.size()) {
        return;  // Start index out of range
    }

    // Calculate parameters for consecutive pairs
    // Apply step_size_ by skipping pairs
    // Use 1-based base pair indices for JSON recording (matching legacy)
    for (size_t i = start_idx; i + 1 < base_pairs_.size(); i += step_size_) {
        const auto& pair1 = base_pairs_[i];
        const auto& pair2 = base_pairs_[i + 1];

        // Verify both pairs have frames
        if (!pair1.frame1().has_value() || !pair1.frame2().has_value() ||
            !pair2.frame1().has_value() || !pair2.frame2().has_value()) {
            continue;  // Skip pairs without frames
        }

        // Legacy get_parameters() uses refs_i_j(j, j+1, orien[i], org[i], r1, o1, r2, o2)
        // which extracts frames from the same strand (orien[i] for duplex i)
        // For duplex 1: uses strand 1 frames (pair_num[1][j])
        // For duplex 2: uses strand 2 frames (pair_num[2][j])
        // Modern uses frame1() from each pair, which should be strand 1 (first residue)
        // This matches legacy duplex 1 behavior
        
        // Use frame1() from each pair (matching legacy strand 1 frames)
        // Note: Legacy get_parameters() does NOT apply frame reversals
        // (Frame reversals are only in analyze.c Rotmat path, not get_parameters path)
        core::ReferenceFrame frame1 = pair1.frame1().value();
        core::ReferenceFrame frame2 = pair2.frame1().value();

        // Calculate step parameters using frames (no reversals - matching get_parameters)
        auto step_params = param_calculator_.calculate_step_parameters(frame1, frame2);
        step_parameters_.push_back(step_params);

        // Record to JSON if writer provided
        // Use 1-based base pair indices (matching legacy: j, j+1 where j is 1-based)
        if (json_writer_) {
            size_t bp_idx1 = i + 1;  // Convert 0-based vector index to 1-based base pair index
            size_t bp_idx2 = i + 2;  // Next pair index (1-based)
            json_writer_->record_bpstep_params(bp_idx1, bp_idx2, step_params);
        }

        // Calculate helical parameters
        auto helical_params = param_calculator_.calculate_helical_parameters(pair1, pair2);
        helical_parameters_.push_back(helical_params);

        // Record to JSON if writer provided
        // Use 1-based base pair indices (matching legacy: j, j+1 where j is 1-based)
        if (json_writer_) {
            size_t bp_idx1 = i + 1;  // Convert 0-based vector index to 1-based base pair index
            size_t bp_idx2 = i + 2;  // Next pair index (1-based)
            json_writer_->record_helical_params(bp_idx1, bp_idx2, helical_params);
        }
    }

    // Handle circular structure (last pair with first pair)
    if (circular_structure_ && base_pairs_.size() >= 2) {
        const auto& pair1 = base_pairs_.back();
        const auto& pair2 = base_pairs_.front();

        if (pair1.frame1().has_value() && pair1.frame2().has_value() &&
            pair2.frame1().has_value() && pair2.frame2().has_value()) {
            auto step_params = param_calculator_.calculate_step_parameters(pair1, pair2);
            step_parameters_.push_back(step_params);

            // Record to JSON if writer provided
            if (json_writer_) {
                size_t bp_idx1 = base_pairs_.size();  // Last pair (1-based)
                size_t bp_idx2 = 1;  // First pair (1-based)
                json_writer_->record_bpstep_params(bp_idx1, bp_idx2, step_params);
            }

            auto helical_params = param_calculator_.calculate_helical_parameters(pair1, pair2);
            helical_parameters_.push_back(helical_params);

            // Record to JSON if writer provided
            if (json_writer_) {
                size_t bp_idx1 = base_pairs_.size();  // Last pair (1-based)
                size_t bp_idx2 = 1;  // First pair (1-based)
                json_writer_->record_helical_params(bp_idx1, bp_idx2, helical_params);
            }
        }
    }
}

core::Structure AnalyzeProtocol::load_structure(const std::filesystem::path& pdb_file) {
    // Configure parser based on input data
    // Note: input_data_.flags may contain information about HETATM inclusion
    // For now, we'll use default settings matching legacy behavior

    // Parse PDB file
    return pdb_parser_.parse_file(pdb_file);
}

void AnalyzeProtocol::convert_atom_indices_to_residue_indices(const core::Structure& structure) {
    // Get number of residues in legacy order
    auto residues_legacy = core::get_residues_in_legacy_order(structure);
    size_t num_residues = residues_legacy.size();

    // Build map from atom index to residue index
    // Legacy atom indices are 1-based, stored in atom.legacy_atom_idx()
    std::map<int, int> atom_idx_to_residue_idx;
    
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            // Get legacy residue index for this residue
            int legacy_residue_idx = structure.get_legacy_idx_for_residue(&residue);
            if (legacy_residue_idx == 0) {
                continue;  // Residue not found in legacy order
            }
            
            // Map all atoms in this residue to the residue index
            for (const auto& atom : residue.atoms()) {
                int legacy_atom_idx = atom.legacy_atom_idx();
                if (legacy_atom_idx > 0) {
                    atom_idx_to_residue_idx[legacy_atom_idx] = legacy_residue_idx;
                }
            }
        }
    }

    // Debug: Check if we have atom indices mapped
    if (atom_idx_to_residue_idx.empty() && !base_pairs_.empty()) {
        std::cerr << "Warning: No legacy atom indices found in structure. "
                  << "Atom index conversion may not work correctly.\n";
    }

    // Convert base pair indices if they appear to be atom indices
    // Detection: if index > num_residues, it's likely an atom index
    int converted_count = 0;
    for (auto& pair : base_pairs_) {
        size_t idx1 = pair.residue_idx1();
        size_t idx2 = pair.residue_idx2();
        
        // Convert from 0-based to 1-based for comparison
        int idx1_1based = static_cast<int>(idx1 + 1);
        int idx2_1based = static_cast<int>(idx2 + 1);
        
        // Check if indices are likely atom indices (larger than residue count)
        bool idx1_is_atom = (idx1_1based > static_cast<int>(num_residues));
        bool idx2_is_atom = (idx2_1based > static_cast<int>(num_residues));
        
        if (idx1_is_atom) {
            auto it = atom_idx_to_residue_idx.find(idx1_1based);
            if (it != atom_idx_to_residue_idx.end() && it->second > 0) {
                // Convert to 0-based residue index
                pair.set_residue_idx1(static_cast<size_t>(it->second - 1));
                converted_count++;
            }
        }
        
        if (idx2_is_atom) {
            auto it = atom_idx_to_residue_idx.find(idx2_1based);
            if (it != atom_idx_to_residue_idx.end() && it->second > 0) {
                // Convert to 0-based residue index
                pair.set_residue_idx2(static_cast<size_t>(it->second - 1));
                converted_count++;
            }
        }
    }
    
    if (converted_count > 0) {
        std::cout << "Converted " << converted_count << " atom indices to residue indices\n";
    }
}

} // namespace protocols
} // namespace x3dna

