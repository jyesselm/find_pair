/**
 * @file base_frame_calculator.cpp
 * @brief Implementation of base frame calculator
 */

#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <stdexcept>
#include <cmath>

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
    
    // Check if valid nucleotide
    if (residue_type == core::ResidueType::UNKNOWN ||
        residue_type == core::ResidueType::AMINO_ACID) {
        return result; // Cannot calculate frame for non-nucleotides
    }
    
    // Load standard template
    core::Structure standard_template;
    try {
        standard_template = templates_.load_template(residue_type);
        result.template_file = templates_.get_template_path(residue_type);
    } catch (const std::exception&) {
        // Template not found or failed to load
        return result;
    }
    
    // Match ring atoms
    MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template, is_rna_);
    
    if (!matched.is_valid()) {
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

