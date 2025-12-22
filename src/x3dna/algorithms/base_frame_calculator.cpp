/**
 * @file base_frame_calculator.cpp
 * @brief Implementation of base frame calculator
 */

#include <x3dna/algorithms/base_frame_calculator.hpp>
#include <x3dna/algorithms/validation_constants.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <x3dna/core/typing/type_registry.hpp>
#include <x3dna/core/typing/atom_type.hpp>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace algorithms {

using core::AtomType;
using core::RING_ATOM_TYPES;
using core::NUM_RING_ATOM_TYPES;

namespace {
// Standard nucleotide ring geometry (from legacy xyz_ring array)
// Matches RA_LIST order: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
constexpr std::array<std::array<double, 3>, 9> STANDARD_RING_GEOMETRY = {{
    {{-1.265, 3.177, 0.000}}, // C4
    {{-2.342, 2.364, 0.001}}, // N3
    {{-1.999, 1.087, 0.000}}, // C2
    {{-0.700, 0.641, 0.000}}, // N1
    {{0.424, 1.460, 0.000}},  // C6
    {{0.071, 2.833, 0.000}},  // C5
    {{0.870, 3.969, 0.000}},  // N7 (purine)
    {{0.023, 4.962, 0.000}},  // C8 (purine)
    {{-1.289, 4.551, 0.000}}, // N9 (purine)
}};

// Ring atom names are now accessed via RING_ATOM_TYPES from atom_type.hpp
// Use core::typing::to_string(AtomType) to get the string representation

/**
 * @brief Check nucleotide type by RMSD (matches legacy check_nt_type_by_rmsd)
 */
RmsdCheckResult check_nt_type_by_rmsd(const core::Residue& residue) {
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int nN = 0;
    bool has_c1_prime = false;
    int purine_atom_count = 0;

    // Try to match ALL ring atoms using AtomType for O(1) comparison
    for (size_t i = 0; i < NUM_RING_ATOM_TYPES; ++i) {
        AtomType target_type = RING_ATOM_TYPES[i];
        const core::Atom* atom = residue.find_atom_by_type(target_type);
        if (atom != nullptr) {
            const auto& pos = atom->position();
            experimental_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));
            standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                         STANDARD_RING_GEOMETRY[i][2]));
            // N3, N1, N7, N9 are nitrogen atoms (indices 1, 3, 6, 8)
            if (i == 1 || i == 3 || i == 6 || i == 8) {
                nN++;
            }
            if (i >= 6) {
                purine_atom_count++;
            }
        }
    }

    // Check for C1' atom using AtomType (O(1) comparison)
    // Also check for C1R via string (not in AtomType enum - rare alternative name)
    if (residue.has_atom_type(AtomType::C1_PRIME)) {
        has_c1_prime = true;
    } else {
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == "C1R") {
                has_c1_prime = true;
                break;
            }
        }
    }

    if (nN == 0 && !has_c1_prime) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }
    if (experimental_coords.size() < 3) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }

    // Collect matched atom names using AtomType
    std::vector<std::string> matched_names;
    for (size_t i = 0; i < NUM_RING_ATOM_TYPES; ++i) {
        AtomType target_type = RING_ATOM_TYPES[i];
        if (residue.has_atom_type(target_type)) {
            matched_names.push_back(std::string(core::typing::to_string(target_type)));
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
    static const std::vector<std::string> excluded = {"MES", "HEPES", "TRIS", "EDO", "GOL", "SO4",
                                                      "PO4", "ACT",   "FMT",  "EFZ", "LYA"};
    return std::find(excluded.begin(), excluded.end(), res_name) != excluded.end();
}

// Standard nucleotide list
bool is_in_nt_list(const std::string& res_name) {
    static const std::vector<std::string> NT_LIST = {"A",   "C",   "G",   "T",   "U",   "PSU", "P5P", "PU", "I",  "DI",
                                                     "ADP", "GDP", "CDP", "UDP", "TDP", "DA",  "DC",  "DG", "DT", "DU"};
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

// Count ring atoms in a residue using AtomType
std::tuple<int, bool> count_ring_atoms(const core::Residue& residue) {
    // Common ring atoms (indices 0-5 in RING_ATOM_TYPES)
    static constexpr AtomType common_ring[] = {
        AtomType::C4, AtomType::N3, AtomType::C2, AtomType::N1, AtomType::C6, AtomType::C5
    };
    // Purine-only ring atoms (indices 6-8 in RING_ATOM_TYPES)
    static constexpr AtomType purine_ring[] = {AtomType::N7, AtomType::C8, AtomType::N9};

    int count = 0;
    bool has_purine = false;

    for (auto type : common_ring) {
        if (residue.has_atom_type(type)) {
            count++;
        }
    }
    for (auto type : purine_ring) {
        if (residue.has_atom_type(type)) {
            count++;
            has_purine = true;
        }
    }
    return {count, has_purine};
}

// Check if one-letter code represents a purine base (A, G, I)
[[nodiscard]] bool is_purine_letter(char letter) {
    const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(letter)));
    return upper == 'A' || upper == 'G' || upper == 'I';
}

// Detect purine atoms using AtomType (O(1) comparison)
bool detect_purine_atoms(const core::Residue& residue) {
    return residue.has_atom_type(AtomType::N7) ||
           residue.has_atom_type(AtomType::C8) ||
           residue.has_atom_type(AtomType::N9);
}

// Determine purine type (A vs G) using AtomType
core::typing::BaseType determine_purine_type(const core::Residue& residue) {
    bool has_o6 = residue.has_atom_type(AtomType::O6);
    bool has_n6 = residue.has_atom_type(AtomType::N6);
    bool has_n2 = residue.has_atom_type(AtomType::N2);

    if (has_o6 || (!has_n6 && has_n2)) {
        return core::typing::BaseType::GUANINE;
    }
    return core::typing::BaseType::ADENINE;
}

// Determine pyrimidine type using AtomType
core::typing::BaseType determine_pyrimidine_type(const core::Residue& residue, char one_letter) {
    bool has_n4 = residue.has_atom_type(AtomType::N4);
    bool has_c5m = residue.has_atom_type(AtomType::C5M) || residue.has_atom_type(AtomType::C7);

    // Check for pseudouridine
    const core::Atom* c1p = residue.find_atom_by_type(AtomType::C1_PRIME);
    const core::Atom* n1 = residue.find_atom_by_type(AtomType::N1);
    const core::Atom* c5 = residue.find_atom_by_type(AtomType::C5);
    if (c1p && n1 && c5) {
        double dist_n1 = (c1p->position() - n1->position()).length();
        double dist_c5 = (c1p->position() - c5->position()).length();
        if (dist_c5 <= 2.0 && dist_n1 > 2.0) {
            return core::typing::BaseType::PSEUDOURIDINE;
        }
    }

    if (has_n4)
        return core::typing::BaseType::CYTOSINE;
    if (has_c5m && one_letter != 'u')
        return core::typing::BaseType::THYMINE;
    return core::typing::BaseType::URACIL;
}

// Try pyrimidine-only RMSD check using AtomType
std::optional<double> try_pyrimidine_rmsd(const core::Residue& residue) {
    std::vector<geometry::Vector3D> exp_coords;
    std::vector<geometry::Vector3D> std_coords;

    // Only check first 6 ring atoms (pyrimidine ring, excluding purine-only N7, C8, N9)
    for (size_t i = 0; i < 6; ++i) {
        AtomType target_type = RING_ATOM_TYPES[i];
        const core::Atom* atom = residue.find_atom_by_type(target_type);
        if (atom != nullptr) {
            const auto& pos = atom->position();
            exp_coords.push_back(geometry::Vector3D(pos.x(), pos.y(), pos.z()));
            std_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                    STANDARD_RING_GEOMETRY[i][2]));
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

// Check if residue needs special type detection
bool needs_type_detection(const core::Residue& residue) {
    auto mol_type = residue.molecule_type();
    auto base_type = residue.base_type();
    // Unknown molecule type or unknown base type for a nucleotide
    return mol_type == core::typing::MoleculeType::UNKNOWN ||
           mol_type == core::typing::MoleculeType::PROTEIN ||
           (mol_type == core::typing::MoleculeType::NUCLEIC_ACID && base_type == core::typing::BaseType::UNKNOWN);
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
    core::typing::BaseType base_type = residue.base_type();
    std::string res_name = residue.name();
    while (!res_name.empty() && res_name[0] == ' ')
        res_name.erase(0, 1);
    while (!res_name.empty() && res_name.back() == ' ')
        res_name.pop_back();

    char one_letter = core::one_letter_code(residue);
    if (one_letter == ' ' || is_excluded_molecule(res_name)) {
        return result;
    }

    bool needs_rmsd_check = !is_in_nt_list(res_name);
    bool has_ring_atoms = false;
    bool has_purine_atoms = false;
    bool used_pyrimidine_fallback = false;

    // Check for ring atoms if needed
    bool should_check_rings = needs_type_detection(residue) || needs_rmsd_check;

    if (should_check_rings) {
        auto [ring_count, has_purine] = count_ring_atoms(residue);
        has_ring_atoms = (ring_count >= 3);
        has_purine_atoms = has_purine || detect_purine_atoms(residue);
    } else {
        has_ring_atoms = true;
    }

    // RMSD check
    RmsdCheckResult rmsd_check;
    std::optional<double> rmsd_result;
    bool found_purine_atoms = false;
    constexpr double RMSD_THRESHOLD = validation_constants::NT_RMSD_CUTOFF;

    if (has_ring_atoms) {
        rmsd_check = check_nt_type_by_rmsd(residue);
        rmsd_result = rmsd_check.rmsd;
        found_purine_atoms = rmsd_check.found_purine_atoms;

        if (!rmsd_result.has_value() || *rmsd_result > RMSD_THRESHOLD) {
            if (found_purine_atoms) {
                auto pyr_rmsd = try_pyrimidine_rmsd(residue);
                if (!pyr_rmsd.has_value() || *pyr_rmsd > RMSD_THRESHOLD) {
                    return result;
                }
                if (!core::TypeRegistry::instance().is_nucleotide(res_name)) {
                    has_purine_atoms = false;
                }
                used_pyrimidine_fallback = true;
                rmsd_result = pyr_rmsd;
            } else {
                return result;
            }
        }
    }

    if (!has_ring_atoms) {
        return result;
    }

    // Determine base type if unknown
    bool is_registry_nucleotide = core::TypeRegistry::instance().is_nucleotide(res_name);

    if (!is_registry_nucleotide && (needs_type_detection(residue) || needs_rmsd_check)) {
        if (has_ring_atoms) {
            if (has_purine_atoms || is_purine_letter(one_letter)) {
                base_type = determine_purine_type(residue);
            } else {
                base_type = determine_pyrimidine_type(residue, one_letter);
            }
        } else {
            return result;
        }
    }

    // Load template
    bool is_modified = std::islower(static_cast<unsigned char>(one_letter));
    core::Structure standard_template;
    try {
        standard_template = templates_.load_template(base_type, is_modified);
        result.template_file = templates_.get_template_path(base_type, is_modified);
    } catch (const std::exception&) {
        return result;
    }

    // Match ring atoms
    core::typing::BaseType matching_type = base_type;
    if (used_pyrimidine_fallback &&
        (base_type == core::typing::BaseType::ADENINE || base_type == core::typing::BaseType::GUANINE)) {
        matching_type = core::typing::BaseType::URACIL;
    }
    MatchedAtoms matched = RingAtomMatcher::match(residue, standard_template,
                                                  std::optional<core::typing::BaseType>(matching_type));

    // Fallback to RMSD check atoms if template matching failed
    if (!matched.is_valid() && has_ring_atoms && rmsd_result.has_value() && !rmsd_check.matched_atom_names.empty() &&
        rmsd_check.matched_atom_names.size() >= 3) {
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
            core::Atom std_atom = core::Atom::create(atom_name, geometry::Vector3D(std_coord.x(), std_coord.y(), std_coord.z()))
                                  .build();
            matched.standard.push_back(std_atom);
        }

        if (matched.experimental.size() < 3 || matched.standard.size() < 3) {
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
    result.frame = core::ReferenceFrame(result.rotation_matrix, result.translation);
    result.is_valid = true;

    return result;
}

void BaseFrameCalculator::calculate_all_frames(core::Structure& structure) {
    std::vector<core::Residue*> residues;
    for (const auto* ptr : structure.residues_in_legacy_order()) {
        residues.push_back(const_cast<core::Residue*>(ptr));
    }

    for (auto* residue : residues) {
        if (residue->is_protein()) {
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
            // Use AtomType for O2' check (fast O(1) comparison)
            // Also check O2* via string (alternative name not in AtomType enum)
            if (residue.has_atom_type(AtomType::O2_PRIME) || residue.find_atom("O2*").has_value()) {
                return true;
            }
        }
    }
    return false;
}

} // namespace algorithms
} // namespace x3dna
