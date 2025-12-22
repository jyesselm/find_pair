/**
 * @file base_pair_validator.cpp
 * @brief Implementation of base pair validation (matches legacy check_pair)
 */

#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/validation/overlap_calculator.hpp>
#include <x3dna/core/atom.hpp>
#include <x3dna/core/typing.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <x3dna/algorithms/ring_atom_matcher.hpp>
#include <x3dna/algorithms/hydrogen_bond.hpp>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cctype>
#include <vector>
#include <set>
#include <map>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>

namespace x3dna {
namespace algorithms {

using namespace x3dna::core;
using namespace x3dna::geometry;

ValidationResult BasePairValidator::validate(const Residue& res1, const Residue& res2) const {
    ValidationResult result;

    // Skip if same residue
    if (&res1 == &res2) {
        return result;
    }

    // Both residues must have reference frames
    auto frame1_opt = res1.reference_frame();
    auto frame2_opt = res2.reference_frame();
    if (!frame1_opt.has_value() || !frame2_opt.has_value()) {
        return result;
    }

    const ReferenceFrame& frame1 = frame1_opt.value();
    const ReferenceFrame& frame2 = frame2_opt.value();

    // Calculate average z-axis and origin
    Vector3D oave, zave;
    get_bp_zoave(frame1, frame2, oave, zave);

    // Calculate dorg (distance between origins)
    // Legacy: ddxyz(org[i], org[j], dorg) means dorg = org[i] - org[j]
    // Note: Legacy uses 1-based indexing, so org[i] is first residue, org[j] is second
    // Modern: frame1 is first (i), frame2 is second (j)
    Vector3D dorg = frame1.origin() - frame2.origin();
    result.dorg = dorg.length();

    // Calculate direction vectors (dot products of frame axes)
    calculate_direction_vectors(frame1, frame2, result.dir_x, result.dir_y, result.dir_z);

    // Calculate d_v (vertical distance)
    result.d_v = std::abs(dorg.dot(zave));

    // Calculate plane angle (angle between z-axes, 0-90 degrees)
    result.plane_angle = z1_z2_angle_in_0_to_90(frame1.z_axis(), frame2.z_axis());

    // Calculate dNN (distance between N1/N9 atoms)
    auto n1_n9_1 = find_n1_n9_position(res1);
    auto n1_n9_2 = find_n1_n9_position(res2);
    if (n1_n9_1.has_value() && n1_n9_2.has_value()) {
        Vector3D dNN_vec = n1_n9_1.value() - n1_n9_2.value();
        result.dNN = dNN_vec.length();
    } else {
        result.dNN = validation_constants::DNN_FALLBACK; // Large value if N1/N9 not found
    }

    // Calculate quality score (matches rtn_val[5])
    result.quality_score = result.dorg + validation_constants::D_V_WEIGHT * result.d_v +
                           result.plane_angle / validation_constants::PLANE_ANGLE_DIVISOR;

    // Perform validation checks
    result.distance_check = in_range(result.dorg, params_.min_dorg, params_.max_dorg);
    result.d_v_check = in_range(result.d_v, params_.min_dv, params_.max_dv);
    result.plane_angle_check = in_range(result.plane_angle, params_.min_plane_angle, params_.max_plane_angle);
    result.dNN_check = in_range(result.dNN, params_.min_dNN, params_.max_dNN);

    // Check overlap area
    result.overlap_area = calculate_overlap_area(res1, res2, oave, zave);
    result.overlap_check = (result.overlap_area < params_.overlap_threshold);

    // If all distance/angle checks pass and overlap is acceptable, check hydrogen bonds
    bool cdns = result.distance_check && result.d_v_check && result.plane_angle_check && result.dNN_check &&
                result.overlap_check;


    if (cdns) {
        // Count H-bonds simply (BEFORE validation) - matches legacy check_pair behavior
        // This is the key fix: legacy counts H-bonds before validation for pair validation
        using namespace x3dna::algorithms::hydrogen_bond;
        HBondDetector detector(HBondDetectionParams::legacy_compatible());
        detector.count_potential_hbonds(res1, res2, result.num_base_hb, result.num_o2_hb);

        // Check H-bond requirement (matches legacy lines 4616-4617)
        if (params_.min_base_hb > 0) {
            result.hbond_check = (result.num_base_hb >= params_.min_base_hb);
        } else {
            result.hbond_check = (result.num_o2_hb > 0 || result.num_base_hb > 0);
        }

        // Pair is valid if all checks pass
        result.is_valid = result.hbond_check;

        // Find validated H-bonds (AFTER validation) - used for adjust_pairQuality
        // This matches legacy hb_numlist behavior which uses validated H-bonds
        if (result.is_valid) {
            result.hbonds = find_hydrogen_bonds(res1, res2);
        }

        // Determine base pair type (simplified - would need calculate_more_bppars)
        if (result.is_valid) {
            // For now, set to UNKNOWN - will be determined by calculate_more_bppars equivalent
            result.bp_type = BasePairType::UNKNOWN;
        }
    }

    return result;
}

void BasePairValidator::calculate_direction_vectors(const ReferenceFrame& frame1, const ReferenceFrame& frame2,
                                                    double& dir_x, double& dir_y, double& dir_z) {

    // Dot products of corresponding axes
    dir_x = frame1.x_axis().dot(frame2.x_axis());
    dir_y = frame1.y_axis().dot(frame2.y_axis());
    dir_z = frame1.z_axis().dot(frame2.z_axis());
}

void BasePairValidator::get_bp_zoave(const ReferenceFrame& frame1, const ReferenceFrame& frame2, Vector3D& oave,
                                     Vector3D& zave) {

    // Average origin (matches legacy avexyz)
    oave = (frame1.origin() + frame2.origin()) * 0.5;

    // Average z-axis (matches legacy get_bp_zoave exactly)
    // Legacy: d = dot(z1, z2)
    //         if d > 0: zave = z2 + z1 (sumxyz)
    //         else: zave = z2 - z1 (ddxyz)
    //         then normalize
    Vector3D z1 = frame1.z_axis();
    Vector3D z2 = frame2.z_axis();
    double d = z1.dot(z2);

    if (d > 0.0) {
        // z-axes point in same direction: add them
        zave = z2 + z1;
    } else {
        // z-axes point in opposite directions: subtract them
        zave = z2 - z1;
    }

    // Normalize (matches legacy vec_norm)
    double len = zave.length();
    if (len > 1e-10) {
        zave = zave / len;
    } else {
        // Fallback if vectors are opposite and cancel out
        zave = z1;
    }
}

double BasePairValidator::z1_z2_angle_in_0_to_90(const Vector3D& z1, const Vector3D& z2) {

    double dot_product = z1.dot(z2);
    // Clamp to [-1, 1] to avoid numerical issues
    dot_product = std::max(-1.0, std::min(1.0, dot_product));

    double angle_rad = std::acos(dot_product);
    double angle_deg = angle_rad * 180.0 / M_PI;

    // Return angle in 0-90 degree range
    if (angle_deg > 90.0) {
        return 180.0 - angle_deg;
    }
    return angle_deg;
}

std::optional<Vector3D> BasePairValidator::find_n1_n9_position(const Residue& residue) {
    // Match legacy glyco_N logic exactly (org/src/cmn_fncs.c lines 4680-4730)
    // Legacy uses bseq (base type letter) to determine purine/pyrimidine:
    //   - isR=1 for purines (A, G, I and their lowercase modified forms)
    //   - isR=0 for pyrimidines (C, T, U, P and their lowercase modified forms)
    // Then glyco_N looks for N9 (purine) or N1 (pyrimidine)
    // If not found, legacy falls back to finding any atom with '9' or '1' in name

    // CRITICAL: Use one_letter_code (which matches legacy bseq) to determine purine/pyrimidine
    // Do NOT use atom presence (e.g., C8) because some modified pyrimidines like 70U
    // have C8 as part of their modification but should still use N1
    char one_letter = core::one_letter_code(residue);

    // Purines: A, G, I (and their lowercase modified forms a, g, i)
    // Pyrimidines: C, T, U, P (and their lowercase modified forms c, t, u, p)
    bool is_purine = typing::is_purine_letter(one_letter);

    if (is_purine) {
        // Purine: find N9 using AtomType (O(1) comparison)
        const core::Atom* n9 = residue.find_atom_by_type(core::AtomType::N9);
        if (n9 != nullptr) {
            return n9->position();
        }
        // N9 not found - legacy fallback: find atom with '9' in name
        // This handles modified nucleotides like 8B4 that have C9 but not N9
        for (const auto& atom : residue.atoms()) {
            const std::string& name = atom.name();
            if (name.find('9') != std::string::npos) {
                return atom.position();
            }
        }
    } else {
        // Pyrimidine: find N1
        // Special case: for P/p bases, legacy uses C5
        char base = core::one_letter_code(residue);
        if (base == 'P' || base == 'p') {
            const core::Atom* c5 = residue.find_atom_by_type(core::AtomType::C5);
            if (c5 != nullptr) {
                return c5->position();
            }
        }
        // Otherwise, use N1 using AtomType (O(1) comparison)
        const core::Atom* n1 = residue.find_atom_by_type(core::AtomType::N1);
        if (n1 != nullptr) {
            return n1->position();
        }
        // N1 not found - legacy fallback: find atom with '1' in name
        for (const auto& atom : residue.atoms()) {
            const std::string& name = atom.name();
            if (name.find('1') != std::string::npos) {
                return atom.position();
            }
        }
    }

    return std::nullopt;
}

double BasePairValidator::calculate_overlap_area(const Residue& res1, const Residue& res2, const Vector3D& oave,
                                                 const Vector3D& zave) const {
    // Delegate to OverlapCalculator with cached ring data
    return validation::OverlapCalculator::calculate(res1, res2, oave, zave, ring_data_cache_);
}

std::vector<core::hydrogen_bond> BasePairValidator::find_hydrogen_bonds(const Residue& res1,
                                                                        const Residue& res2) const {
    // Use new HBondDetector with legacy-compatible parameters
    using namespace x3dna::algorithms::hydrogen_bond;

    HBondDetector detector(HBondDetectionParams::legacy_compatible());

    // Detect ALL H-bonds (not just base-base) to match baseline behavior
    auto result = detector.detect_all_hbonds_detailed(res1, res2,
        core::typing::MoleculeType::NUCLEIC_ACID, core::typing::MoleculeType::NUCLEIC_ACID);

    // Helper to pad atom name to 4 characters (matches legacy " O2 " format)
    auto pad_atom_name = [](const std::string& name) -> std::string {
        if (name.size() >= 4) return name;
        std::string padded = " " + name;
        while (padded.size() < 4) padded += " ";
        return padded;
    };

    // Return all classified bonds (including INVALID) to match baseline
    std::vector<core::hydrogen_bond> all_hbonds;
    for (const auto& hbond : result.all_classified_bonds) {
        core::hydrogen_bond h;
        h.donor_atom = pad_atom_name(hbond.donor_atom_name);
        h.acceptor_atom = pad_atom_name(hbond.acceptor_atom_name);
        h.distance = hbond.distance;
        h.type = hbond.legacy_type_char();
        all_hbonds.push_back(h);
    }

    return all_hbonds;
}

char BasePairValidator::donor_acceptor(char base1, char base2, const std::string& atom1, const std::string& atom2) {
    // Matches legacy donor_acceptor function
    // CB_LIST = "ACGITU" (A=0, C=1, G=2, I=3, T=4, U=5)
    static const char* CB_LIST = "ACGITU";
    static const char* da_types[] = {"AD", "AX", "XD", "XX", "DA", "DX", "XA"};
    static const char* bb_da[] = {"O1P", "O2P", "O5'", "O4'", "O3'", "O2'"};
    static const char bb_da_type[] = {'A', 'A', 'A', 'A', 'A', 'X'};

    // Base-specific atom patterns: [base_index][atom_index] = "atom_name_type"
    // base_da[base_index][atom_index] format: "N9" with type '?' (either)
    static const char* base_da[6][6] = {// A (Adenine, index 0)
                                        {"N9", "N7", "N6", "N1", "N3", nullptr},
                                        // C (Cytosine, index 1)
                                        {"N1", "O2", "N3", "N4", nullptr, nullptr},
                                        // G (Guanine, index 2)
                                        {"N9", "N7", "O6", "N1", "N2", "N3"},
                                        // I (Inosine, index 3)
                                        {"N9", "N7", "O6", "N1", "N3", nullptr},
                                        // T (Thymine, index 4)
                                        {"N1", "O2", "N3", "O4", nullptr, nullptr},
                                        // U (Uracil, index 5)
                                        {"N1", "O2", "N3", "O4", nullptr, nullptr}};

    // Base-specific atom types: [base_index][atom_index] = type
    static const char base_da_type[6][6] = {// A
                                            {'?', 'A', 'D', 'A', 'A', '\0'},
                                            // C
                                            {'?', 'A', 'A', 'D', '\0', '\0'},
                                            // G
                                            {'?', 'A', 'A', 'D', 'D', 'A'},
                                            // I
                                            {'?', 'A', 'A', 'D', 'A', '\0'},
                                            // T
                                            {'?', 'A', 'D', 'A', '\0', '\0'},
                                            // U
                                            {'?', 'A', 'D', 'A', '\0', '\0'}};

    char hbatom_type = '*';
    char ia = '\0', ja = '\0';

    // Find base indices in CB_LIST
    int inum = -1, jnum = -1;
    const char* pchar = strchr(CB_LIST, std::toupper(static_cast<unsigned char>(base1)));
    if (pchar)
        inum = pchar - CB_LIST;
    pchar = strchr(CB_LIST, std::toupper(static_cast<unsigned char>(base2)));
    if (pchar)
        jnum = pchar - CB_LIST;

    if (inum < 0 || jnum < 0) {
        return hbatom_type; // Invalid base
    }

    // Check backbone atoms first (names are now trimmed)
    for (int i = 0; i < 6; ++i) {
        if (atom1 == bb_da[i]) {
            ia = bb_da_type[i];
        }
        if (atom2 == bb_da[i]) {
            ja = bb_da_type[i];
        }
    }

    // Check base-specific atoms if not found in backbone
    if (!ia && inum >= 0 && inum < 6) {
        for (int i = 0; i < 6; ++i) {
            if (base_da[inum][i] == nullptr)
                break;
            if (atom1 == base_da[inum][i]) {
                ia = base_da_type[inum][i];
                break;
            }
        }
    }

    if (!ja && jnum >= 0 && jnum < 6) {
        for (int i = 0; i < 6; ++i) {
            if (base_da[jnum][i] == nullptr)
                break;
            if (atom2 == base_da[jnum][i]) {
                ja = base_da_type[jnum][i];
                break;
            }
        }
    }

    // Check if donor-acceptor combination is valid
    if (ia && ja) {
        char da[3] = {ia, ja, '\0'};
        // Check against valid da_types: "AD", "AX", "XD", "XX", "DA", "DX", "XA"
        for (int i = 0; i < 7; ++i) {
            if (std::strcmp(da, da_types[i]) == 0) {
                hbatom_type = '-'; // Valid standard H-bond
                break;
            }
        }
    }

    return hbatom_type;
}

bool BasePairValidator::pattern_match(const std::string& str, const std::string& pattern) {
    // Matches legacy str_pmatch: pattern matching where '.' matches any character
    if (str.length() != pattern.length()) {
        return false;
    }
    for (size_t i = 0; i < str.length() && i < pattern.length(); i++) {
        if (pattern[i] != '.' && pattern[i] != str[i]) {
            return false;
        }
    }
    return true;
}

} // namespace algorithms
} // namespace x3dna
