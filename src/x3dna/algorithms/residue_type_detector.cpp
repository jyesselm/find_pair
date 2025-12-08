/**
 * @file residue_type_detector.cpp
 * @brief Implementation of residue type detector
 */

#include <x3dna/algorithms/residue_type_detector.hpp>
#include <x3dna/geometry/least_squares_fitter.hpp>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace algorithms {

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

// Legacy RA_LIST order for ring atoms
constexpr std::array<const char*, 9> RING_ATOM_NAMES = {" C4 ", " N3 ", " C2 ", " N1 ", " C6 ",
                                                        " C5 ", " N7 ", " C8 ", " N9 "};

} // namespace

RmsdCheckResult ResidueTypeDetector::check_by_rmsd(const core::Residue& residue) {
    std::vector<geometry::Vector3D> experimental_coords;
    std::vector<geometry::Vector3D> standard_coords;
    int purine_atom_count = 0;
    bool has_c1_prime = false;

    // Match ring atoms (matches legacy logic in find_pair)
    // Legacy RA_LIST: " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
    for (size_t i = 0; i < RING_ATOM_NAMES.size(); ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];

        // Find matching atom in residue
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                experimental_coords.push_back(atom.position());
                standard_coords.push_back(geometry::Vector3D(STANDARD_RING_GEOMETRY[i][0], STANDARD_RING_GEOMETRY[i][1],
                                                             STANDARD_RING_GEOMETRY[i][2]));

                // Count purine atoms (N7, C8, N9)
                if (atom_name == std::string(" N7 ") || atom_name == std::string(" C8 ") ||
                    atom_name == std::string(" N9 ")) {
                    purine_atom_count++;
                }

                break;
            }
        }
    }

    // Check for C1' or C1R (sugar atom)
    // Some nucleotides like NMN use C1R instead of C1'
    for (const auto& atom : residue.atoms()) {
        if (atom.name() == " C1'" || atom.name() == " C1R") {
            has_c1_prime = true;
            break;
        }
    }

    // Legacy requires: (!nN && !C1_prime) -> return DUMMY
    int nN = static_cast<int>(experimental_coords.size());
    if (nN == 0 && !has_c1_prime) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }

    // Need at least 3 atoms for RMSD calculation
    if (experimental_coords.size() < 3) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }

    // Collect matched atom names in the order they were found
    std::vector<std::string> matched_names;
    for (size_t i = 0; i < RING_ATOM_NAMES.size(); ++i) {
        const char* atom_name = RING_ATOM_NAMES[i];
        // Check if this atom was found in residue
        for (const auto& atom : residue.atoms()) {
            if (atom.name() == atom_name) {
                matched_names.push_back(std::string(atom_name));
                break;
            }
        }
    }

    // Perform least-squares fitting (matches legacy ls_fitting)
    geometry::LeastSquaresFitter fitter;
    try {
        auto fit_result = fitter.fit(standard_coords, experimental_coords);
        return {fit_result.rms, purine_atom_count > 0, matched_names, experimental_coords, standard_coords};
    } catch (const std::exception&) {
        return {std::nullopt, purine_atom_count > 0, {}, {}, {}};
    }
}

bool ResidueTypeDetector::is_in_nt_list(const std::string& res_name) {
    static const std::vector<std::string> NT_LIST = {"A",   "C",   "G",   "T",   "U",   "PSU", "P5P", "PU", "I",  "DI",
                                                     "ADP", "GDP", "CDP", "UDP", "TDP", "DA",  "DC",  "DG", "DT", "DU"};

    std::string res_upper = res_name;
    for (char& c : res_upper)
        c = std::toupper(c);

    return std::find(NT_LIST.begin(), NT_LIST.end(), res_upper) != NT_LIST.end();
}

bool ResidueTypeDetector::is_purine(core::ResidueType type) {
    return type == core::ResidueType::ADENINE || type == core::ResidueType::GUANINE ||
           type == core::ResidueType::INOSINE;
}

TypeDetectionResult ResidueTypeDetector::detect_type(const core::Residue& residue) {
    TypeDetectionResult result;
    result.detected_type = residue.residue_type();
    result.used_fallback = false;
    result.detection_method = "standard";

    // If already classified and not UNKNOWN, use it
    if (result.detected_type != core::ResidueType::UNKNOWN) {
        return result;
    }

    // Get residue name
    std::string res_name = residue.name();
    while (!res_name.empty() && res_name[0] == ' ')
        res_name.erase(0, 1);
    while (!res_name.empty() && res_name.back() == ' ')
        res_name.pop_back();

    // Check if in NT_LIST - if so, no RMSD check needed
    if (is_in_nt_list(res_name)) {
        result.detection_method = "nt_list";
        return result;
    }

    // Not in NT_LIST - try RMSD check
    auto rmsd_result = check_by_rmsd(residue);
    result.rmsd = rmsd_result.rmsd;

    if (rmsd_result.rmsd.has_value()) {
        result.used_fallback = true;
        result.detection_method = "rmsd";

        // Use purine atom presence to determine type
        if (rmsd_result.found_purine_atoms) {
            result.detected_type = core::ResidueType::ADENINE; // Default purine
        } else {
            result.detected_type = core::ResidueType::CYTOSINE; // Default pyrimidine
        }
    }

    return result;
}

} // namespace algorithms
} // namespace x3dna
