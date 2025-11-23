/**
 * @file base_frame_calculator.cpp
 * @brief Implementation of base frame calculator
 */

#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace x3dna {
namespace algorithms {

BaseFrameCalculator::BaseFrameCalculator(const std::filesystem::path& template_path)
    : templates_(template_path) {
}

FrameCalculationResult BaseFrameCalculator::calculate_frame(core::Residue& residue) {
    FrameCalculationResult result = calculate_frame_impl(residue);
    
    // Store frame in residue if valid
    if (result.is_valid) {
        residue.set_reference_frame(result.frame);
    }
    
    return result;
}

FrameCalculationResult BaseFrameCalculator::calculate_frame_const(const core::Residue& residue) const {
    return calculate_frame_impl(residue);
}

FrameCalculationResult BaseFrameCalculator::calculate_frame_impl(const core::Residue& residue) const {
    FrameCalculationResult result;
    result.is_valid = false;
    
    // Get residue type
    core::ResidueType residue_type = residue.residue_type();
    
    // DEBUG: Log residue info (can be enabled with -DDEBUG_FRAME_CALC)
    #ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: Calculating frame for residue: " << residue.name()
              << " " << residue.chain_id() << ":" << residue.seq_num()
              << " (type=" << static_cast<int>(residue_type) 
              << ", one_letter=" << residue.one_letter_code() << ")\n";
    #endif
    
    // Check if valid nucleotide
    // Legacy includes modified nucleotides that have ring atoms but aren't in NT_LIST
    // Check if residue has ring atoms as a fallback (similar to legacy residue_ident)
    bool has_ring_atoms = false;
    int ring_atom_count = 0;
    bool has_purine_atoms = false;
    
    if (residue_type == core::ResidueType::UNKNOWN ||
        residue_type == core::ResidueType::AMINO_ACID) {
        // Check for ring atoms (C4, N3, C2, N1, C6, C5 are common to all)
        static const std::vector<std::string> common_ring_atoms = {
            " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "
        };
        static const std::vector<std::string> purine_ring_atoms = {
            " N7 ", " C8 ", " N9 "
        };
        
        for (const auto& atom_name : common_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    ring_atom_count++;
                    break;
                }
            }
        }
        
        // Check for purine-specific atoms
        for (const auto& atom_name : purine_ring_atoms) {
            for (const auto& atom : residue.atoms()) {
                if (atom.name() == atom_name) {
                    has_purine_atoms = true;
                    break;
                }
            }
        }
        
        has_ring_atoms = (ring_atom_count >= 3);
        
        #ifdef DEBUG_FRAME_CALC
        if (has_ring_atoms) {
            std::cerr << "DEBUG: Residue has " << ring_atom_count 
                      << " ring atoms - treating as nucleotide\n";
        }
        #endif
        
        // If it has ring atoms, determine type from atoms
        if (has_ring_atoms) {
            if (has_purine_atoms) {
                // Determine purine type (A vs G) by checking for characteristic atoms
                bool has_o6 = false, has_n6 = false, has_n2 = false;
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " O6 ") has_o6 = true;
                    if (atom.name() == " N6 ") has_n6 = true;
                    if (atom.name() == " N2 ") has_n2 = true;
                }
                // G has O6, or N2 without N6; A has N6 without O6
                residue_type = (has_o6 || (!has_n6 && has_n2)) ? 
                               core::ResidueType::GUANINE : 
                               core::ResidueType::ADENINE;
            } else {
                // Determine pyrimidine type (C vs T vs U) by checking for characteristic atoms
                bool has_n4 = false, has_c5m = false;
                for (const auto& atom : residue.atoms()) {
                    if (atom.name() == " N4 ") has_n4 = true;
                    if (atom.name() == " C5M" || atom.name() == " C7 ") has_c5m = true;
                }
                // C has N4, T has C5M/C7, otherwise U (U has O4 but we default to U if no N4 or C5M)
                if (has_n4) {
                    residue_type = core::ResidueType::CYTOSINE;
                } else if (has_c5m) {
                    residue_type = core::ResidueType::THYMINE;
                } else {
                    residue_type = core::ResidueType::URACIL;  // Default for pyrimidines (includes 5MU, PSU)
                }
            }
        } else {
            #ifdef DEBUG_FRAME_CALC
            std::cerr << "DEBUG: Skipping - not a nucleotide (type=" 
                      << static_cast<int>(residue_type) << ", no ring atoms)\n";
            #endif
            return result; // Cannot calculate frame for non-nucleotides
        }
    }
    
    // Load standard template
    core::Structure standard_template;
    try {
        standard_template = templates_.load_template(residue_type);
        result.template_file = templates_.get_template_path(residue_type);
        
        #ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Template loaded: " << result.template_file << "\n";
        std::cerr << "DEBUG: Template has " << standard_template.num_atoms() << " atoms\n";
        #endif
    } catch (const std::exception& e) {
        // Template not found or failed to load
        #ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Template loading failed: " << e.what() << "\n";
        std::cerr << "DEBUG: Template path: " << templates_.template_path() << "\n";
        std::cerr << "DEBUG: Residue type: " << static_cast<int>(residue_type) << "\n";
        #endif
        return result;
    }
    
    // Match ring atoms
    MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template, is_rna_, legacy_mode_);
    
    #ifdef DEBUG_FRAME_CALC
    std::cerr << "DEBUG: Matched " << matched.num_matched << " atoms\n";
    std::cerr << "DEBUG: Matched atom names: ";
    for (const auto& name : matched.atom_names) {
        std::cerr << name << " ";
    }
    std::cerr << "\n";
    #endif
    
    if (!matched.is_valid()) {
        #ifdef DEBUG_FRAME_CALC
        std::cerr << "DEBUG: Not enough matched atoms (need >= 3, got " 
                  << matched.num_matched << ")\n";
        std::cerr << "DEBUG: Residue has " << residue.num_atoms() << " total atoms\n";
        #endif
        // Not enough matched atoms (need at least 3)
        return result;
    }
    
    result.num_matched = matched.num_matched;
    result.matched_atoms = matched.atom_names;
    
    // Extract coordinates from matched atoms
    std::vector<geometry::Vector3D> standard_coords;
    std::vector<geometry::Vector3D> experimental_coords;
    
    for (size_t i = 0; i < matched.num_matched; ++i) {
        standard_coords.push_back(matched.standard[i].position());
        experimental_coords.push_back(matched.experimental[i].position());
    }
    
    // Perform least-squares fitting
    geometry::LeastSquaresFitter fitter;
    geometry::LeastSquaresFitter::FitResult fit_result = fitter.fit(standard_coords, experimental_coords);
    
    // Extract rotation and translation
    result.rotation_matrix = fit_result.rotation;
    result.translation = fit_result.translation;
    result.rms_fit = fit_result.rms;
    
    // Create reference frame from rotation and translation
    result.frame = core::ReferenceFrame(result.rotation_matrix, result.translation);
    result.is_valid = true;
    
    return result;
}

void BaseFrameCalculator::calculate_all_frames(core::Structure& structure) {
    // Iterate through chains and residues, calculating frames
    for (auto& chain : structure.chains()) {
        for (auto& residue : chain.residues()) {
            // Skip non-nucleotides
            if (residue.residue_type() == core::ResidueType::UNKNOWN ||
                residue.residue_type() == core::ResidueType::AMINO_ACID) {
                continue;
            }
            
            // Calculate frame (modifies residue to store frame)
            calculate_frame(residue);
        }
    }
}

void BaseFrameCalculator::set_template_path(const std::filesystem::path& template_path) {
    templates_.set_template_path(template_path);
}

} // namespace algorithms
} // namespace x3dna

