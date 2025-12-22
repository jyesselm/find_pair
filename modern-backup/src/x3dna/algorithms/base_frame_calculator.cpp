/**
 * @file base_frame_calculator.cpp
 * @brief Implementation of base frame calculator
 */

#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/validation_constants.hpp>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace algorithms {

namespace {
// Use constants from centralized constants.hpp
using constants::nucleotide_lists::EXCLUDED_MOLECULES;
using constants::nucleotide_lists::NT_LIST;
using constants::ring_data::COMMON_RING_ATOMS;
using constants::ring_data::PURINE_RING_ATOMS;
using constants::ring_data::RING_ATOM_NAMES;
using constants::ring_data::STANDARD_RING_GEOMETRY;

/**
 * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
 */
RmsdCheckResult check_nt_type_by_rmsd(const core::Residue& residue) {
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int ring_atom_count = 0;
    bool has_c1_prime = false;
    int purine_atom_count = 0;

    // Try to match ALL ring atoms
    for (size_t i = 0; i < RING_ATOM_NAMES.size(); ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                const auto& pos = atom.position();
                experimental_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));
                standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                             STANDARD_RING_GEOMETRY[i][2]));
                if (i == 1 || i == 3 || i == 6 || i == 8) {
                    ring_atom_count++;
                }
                if (i >= 6) {
                    purine_atom_count++;
                }
                break;
            }
        }
    }

    // Check for C1' or C1R atom
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " C1'" || atom.name() == " C1R") {
            has_c1_prime = true;
            break;
        }
    }

    if (ring_atom_count == 0 && !has_c1_prime) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }
    if (experimental_coords.size() < 3) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }

    // Collect matched atom names
    std::vector<std::string> matched_names;
    for (size_t i = 0; i < RING_ATOM_NAMES.size(); ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                matched_names.push_back(std::string(atom_name));
                break;
            }
        }
    }

    // Perform least-squares fitting
    geometry::LeastSquaresFitter fitter;
    try {
        auto fit_result = fitter.fit(standard_coords, experimental_coords);
        return {fit_result.rms, purine_atom_count > 0, matched_names, experimental_coords, standard_coords};
    } catch (const std::exception&) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }
}

// Known non-nucleotide molecules to exclude
bool is_excluded_molecule(const std::string& res_name) {
    for (const auto& excluded : EXCLUDED_MOLECULES) {
        if (res_name == excluded) {
            return true;
        }
    }
    return false;
}

// Standard nucleotide list
bool is_in_nt_list(const std::string& res_name) {
    std::string upper = res_name;
    for (char& c : upper) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    for (const auto& nt : NT_LIST) {
        if (upper == nt) {
            return true;
        }
    }
    return false;
}

// Count ring atoms in a residue
std::tuple<int, bool> count_ring_atoms(const core::Residue& residue) {
    int count = 0;
    bool has_purine = false;

    for (const auto& atom_name : COMMON_RING_ATOMS) {
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                count++;
                break;
            }
        }
    }
    for (const auto& atom_name : PURINE_RING_ATOMS) {
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                count++;
                has_purine = true;
                break;
            }
        }
    }
    return {count, has_purine};
}

// Detect purine atoms
bool detect_purine_atoms(const core::Residue& residue) {
    bool has_n7 = false, has_c8 = false, has_n9 = false;
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " N7 ")
            has_n7 = true;
        if (atom.name() == " C8 ")
            has_c8 = true;
        if (atom.name() == " N9 ")
            has_n9 = true;
    }
    return has_n7 || has_c8 || has_n9;
}

// Determine purine type (A vs G)
core::ResidueType determine_purine_type(const core::Residue& residue) {
    bool has_o6 = false, has_n6 = false, has_n2 = false;
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " O6 ")
            has_o6 = true;
        if (atom.name() == " N6 ")
            has_n6 = true;
        if (atom.name() == " N2 ")
            has_n2 = true;
    }
    return (has_o6 || (!has_n6 && has_n2)) ? core::ResidueType::GUANINE : core::ResidueType::ADENINE;
}

// Determine pyrimidine type
core::ResidueType determine_pyrimidine_type(const core::Residue& residue, char one_letter) {
    bool has_n4 = false, has_c5m = false;
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " N4 ")
            has_n4 = true;
        if (atom.name() == " C5M" || atom.name() == " C7 ")
            has_c5m = true;
    }

    // Check for pseudouridine
    auto c1p = residue.find_atom(" C1'");
    auto n1 = residue.find_atom(" N1 ");
    auto c5 = residue.find_atom(" C5 ");
    if (c1p && n1 && c5) {
        double dist_n1 = (c1p->position() - n1->position()).length();
        double dist_c5 = (c1p->position() - c5->position()).length();
        if (dist_c5 <= 2.0 && dist_n1 > 2.0) {
            return core::ResidueType::PSEUDOURIDINE;
        }
    }

    if (has_n4)
        return core::ResidueType::CYTOSINE;
    if (has_c5m && one_letter != 'u')
        return core::ResidueType::THYMINE;
    return core::ResidueType::URACIL;
}

// Forward declaration
std::optional<double> try_pyrimidine_rmsd(const core::Residue& residue);

/**
 * @brief Result of RMSD validation with pyrimidine fallback
 */
struct RmsdValidationResult {
    bool is_valid = false;
    std::optional<double> rmsd;
    bool used_pyrimidine_fallback = false;
    bool should_clear_purine_flag = false;
    RmsdCheckResult rmsd_check;
};

/**
 * @brief Perform RMSD validation with pyrimidine fallback
 *
 * @param residue Residue to validate
 * @param res_name Trimmed residue name
 * @return RmsdValidationResult with validation status
 */
RmsdValidationResult validate_with_rmsd(const core::Residue& residue, const std::string& res_name) {
    constexpr double RMSD_THRESHOLD = validation_constants::NT_RMSD_CUTOFF;
    RmsdValidationResult result;

    result.rmsd_check = check_nt_type_by_rmsd(residue);
    result.rmsd = result.rmsd_check.rmsd;

    // Check if RMSD is within threshold
    if (!result.rmsd.has_value() || *result.rmsd > RMSD_THRESHOLD) {
        // Try pyrimidine fallback if purine atoms were found
        if (result.rmsd_check.found_purine_atoms) {
            auto pyr_rmsd = try_pyrimidine_rmsd(residue);
            if (!pyr_rmsd.has_value() || *pyr_rmsd > RMSD_THRESHOLD) {
                return result; // is_valid remains false
            }
            // Pyrimidine fallback succeeded
            if (!core::ModifiedNucleotideRegistry::contains(res_name)) {
                result.should_clear_purine_flag = true;
            }
            result.used_pyrimidine_fallback = true;
            result.rmsd = pyr_rmsd;
            result.is_valid = true;
        }
        return result; // is_valid remains false if no purine atoms
    }

    result.is_valid = true;
    return result;
}

/**
 * @brief Apply RMSD fallback atoms to matched atoms structure
 *
 * When template matching fails but RMSD check found atoms, use those atoms.
 *
 * @param residue Residue to extract atoms from
 * @param rmsd_check RMSD check result with matched atoms
 * @param matched MatchedAtoms structure to populate
 * @return true if fallback was applied successfully
 */
bool apply_rmsd_fallback_atoms(const core::Residue& residue, const RmsdCheckResult& rmsd_check, MatchedAtoms& matched) {
    if (rmsd_check.matched_atom_names.empty() || rmsd_check.matched_atom_names.size() < 3) {
        return false;
    }

    matched.num_matched = rmsd_check.matched_atom_names.size();
    matched.atom_names = rmsd_check.matched_atom_names;
    matched.experimental.clear();
    matched.standard.clear();

    for (size_t i = 0; i < rmsd_check.matched_atom_names.size(); ++i) {
        const std::string& atom_name = rmsd_check.matched_atom_names[i];
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                matched.experimental.push_back(atom);
                break;
            }
        }
        const auto& std_coord = rmsd_check.matched_standard_coords[i];
        core::Atom std_atom(atom_name, geometry::Vector3D(std_coord.x(), std_coord.y(), std_coord.z()), "", ' ', 0,
                            'A');
        matched.standard.push_back(std_atom);
    }

    return matched.experimental.size() >= 3 && matched.standard.size() >= 3;
}

// Try pyrimidine-only RMSD check
std::optional<double> try_pyrimidine_rmsd(const core::Residue& residue) {
    std::vector<geometry::Vector3D> exp_coords;
    std::vector<geometry::Vector3D> std_coords;

    for (size_t i = 0; i < 6; ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                const auto& pos = atom.position();
                exp_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));
                std_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                        STANDARD_RING_GEOMETRY[i][2]));
                break;
            }
        }
    }

    if (exp_coords.size() < 3)
        return std::nullopt;

    geometry::LeastSquaresFitter fitter;
    try {
        auto result = fitter.fit(std_coords, exp_coords);
        return result.rms;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

} // namespace

BaseFrameCalculator::BaseFrameCalculator(const std::filesystem::path& template_path) : templates_(template_path) {}

FrameCalculationResult BaseFrameCalculator::calculate_frame(core::Residue& residue) {
    FrameCalculationResult result = calculate_frame_impl(residue);
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

    // Get residue info
    core::ResidueType residue_type = residue.residue_type();
    std::string res_name = residue.name();
    while (!res_name.empty() && res_name[0] == ' ')
        res_name.erase(0, 1);
    while (!res_name.empty() && res_name.back() == ' ')
        res_name.pop_back();

    char one_letter = residue.one_letter_code();
    if (one_letter == ' ' || is_excluded_molecule(res_name)) {
        return result;
    }

    bool needs_rmsd_check = !is_in_nt_list(res_name);
    bool has_ring_atoms = false;
    bool has_purine_atoms = false;
    bool used_pyrimidine_fallback = false;

    // Check for ring atoms if needed
    bool should_check_rings = (residue_type == core::ResidueType::UNKNOWN ||
                               residue_type == core::ResidueType::AMINO_ACID ||
                               residue_type == core::ResidueType::NONCANONICAL_RNA || needs_rmsd_check);

    if (should_check_rings) {
        auto [ring_count, has_purine] = count_ring_atoms(residue);
        has_ring_atoms = (ring_count >= 3);
        has_purine_atoms = has_purine || detect_purine_atoms(residue);
    } else {
        has_ring_atoms = true;
    }

    // RMSD validation with pyrimidine fallback
    RmsdValidationResult rmsd_validation;
    if (has_ring_atoms) {
        rmsd_validation = validate_with_rmsd(residue, res_name);
        if (!rmsd_validation.is_valid) {
            return result;
        }
        if (rmsd_validation.should_clear_purine_flag) {
            has_purine_atoms = false;
        }
        used_pyrimidine_fallback = rmsd_validation.used_pyrimidine_fallback;
    } else {
        return result;
    }

    // Determine residue type
    bool is_registry_nucleotide = core::ModifiedNucleotideRegistry::contains(res_name);

    if (!is_registry_nucleotide &&
        (residue_type == core::ResidueType::UNKNOWN || residue_type == core::ResidueType::AMINO_ACID ||
         residue_type == core::ResidueType::NONCANONICAL_RNA || needs_rmsd_check)) {
        if (has_ring_atoms) {
            bool is_purine_by_letter = (one_letter == 'A' || one_letter == 'G' || one_letter == 'a' ||
                                        one_letter == 'g');
            if (has_purine_atoms || is_purine_by_letter) {
                residue_type = determine_purine_type(residue);
            } else {
                residue_type = determine_pyrimidine_type(residue, one_letter);
            }
        } else {
            return result;
        }
    }

    // Load template
    bool is_modified = std::islower(static_cast<unsigned char>(one_letter));
    core::Structure standard_template;
    try {
        standard_template = templates_.load_template(residue_type, is_modified);
        result.template_file = templates_.get_template_path(residue_type, is_modified);
    } catch (const std::exception&) {
        return result;
    }

    // Match ring atoms
    core::ResidueType matching_type = residue_type;
    if (used_pyrimidine_fallback &&
        (residue_type == core::ResidueType::ADENINE || residue_type == core::ResidueType::GUANINE)) {
        matching_type = core::ResidueType::URACIL;
    }
    MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template,
                                                  std::optional<core::ResidueType>(matching_type));

    // Fallback to RMSD check atoms if template matching failed
    if (!matched.is_valid() && rmsd_validation.rmsd.has_value()) {
        if (!apply_rmsd_fallback_atoms(residue, rmsd_validation.rmsd_check, matched)) {
            return result;
        }
    } else if (!matched.is_valid()) {
        return result;
    }

    result.num_matched = matched.num_matched;
    result.matched_atoms = matched.atom_names;

    // Extract coordinates
    std::vector<geometry::Vector3D> standard_coords;
    std::vector<geometry::Vector3D> experimental_coords;
    for (size_t i = 0; i < matched.num_matched; ++i) {
        standard_coords.push_back(matched.standard[i].position());
        experimental_coords.push_back(matched.experimental[i].position());
    }
    result.matched_standard_coords = standard_coords;
    result.matched_experimental_coords = experimental_coords;

    // Perform least-squares fitting
    geometry::LeastSquaresFitter fitter;
    auto fit_result = fitter.fit(standard_coords, experimental_coords);

    result.rotation_matrix = fit_result.rotation;
    result.translation = fit_result.translation;
    result.rms_fit = fit_result.rms;
    // Pass num_matched to ReferenceFrame so validator can determine purine/pyrimidine fit
    result.frame = core::ReferenceFrame(result.rotation_matrix, result.translation, result.num_matched);
    result.is_valid = true;

    return result;
}

void BaseFrameCalculator::calculate_all_frames(core::Structure& structure) {
    std::vector<core::Residue*> residues;
    for (const auto* ptr : structure.residues_in_legacy_order()) {
        residues.push_back(const_cast<core::Residue*>(ptr));
    }

    for (auto* residue : residues) {
        if (residue->residue_type() == core::ResidueType::AMINO_ACID) {
            continue;
        }
        // Frame is stored in residue as side effect; result intentionally discarded
        (void)calculate_frame(*residue);
    }
}

void BaseFrameCalculator::set_template_path(const std::filesystem::path& template_path) {
    templates_.set_template_path(template_path);
}

bool BaseFrameCalculator::detect_rna(const core::Structure& structure) {
    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (residue.find_atom(" O2'").has_value() || residue.find_atom(" O2*").has_value()) {
                return true;
            }
        }
    }
    return false;
}

} // namespace algorithms
} // namespace x3dna
