/**
 * @file quality_score_calculator.cpp
 * @brief Implementation of QualityScoreCalculator
 */

#include <x3dna/algorithms/quality_score_calculator.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/geometry/matrix3d.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <cmath>

namespace x3dna {
namespace algorithms {

// Watson-Crick pair list (matches legacy WC_LIST)
const std::vector<std::string> QualityScoreCalculator::WC_LIST = {
    "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
};

double QualityScoreCalculator::calculate_selection_score(
    const ValidationResult& result,
    const core::Residue& res1,
    const core::Residue& res2
) const {
    // Start with raw quality score
    double adjusted_score = result.quality_score;

    // Apply H-bond quality adjustment
    adjusted_score += adjust_pair_quality(result.hbonds);

    // Apply bp_type_id == 2 bonus
    int bp_type_id = calculate_bp_type_id(res1, res2, result);
    if (bp_type_id == 2) {
        adjusted_score -= 2.0;
    }

    return adjusted_score;
}

double QualityScoreCalculator::adjust_pair_quality(
    const std::vector<core::hydrogen_bond>& hbonds
) const {
    // Match legacy adjust_pairQuality logic
    // Count good hydrogen bonds (distance in [2.5, 3.5] Angstroms)
    // Legacy: dval = num_list[k][3] / MFACTOR, then checks if in [2.5, 3.5]
    // Modern: hbonds already have distance calculated
    //
    // Legacy skips h-bonds with type '*' (num_list[k][0] == 1)
    // In modern, we have three types:
    //   '-' = standard (good) h-bond
    //   '*' = non-standard h-bond (SKIP)
    //   ' ' = initially unvalidated but can still be counted
    // Legacy ONLY skips '*' types, all others (including ' ') are counted
    int num_good_hb = 0;
    for (const auto& hbond : hbonds) {
        // Legacy flow: hb_info string excludes type ' ' h-bonds (see get_hbond_ij)
        // Then adjust_pairQuality skips type '*' via num_list[k][0]
        // Net result: only type '-' h-bonds are counted for quality adjustment
        if (hbond.type != '-') {
            continue;
        }
        // Check if distance is in good range [2.5, 3.5]
        // CRITICAL: Legacy uses %4.2f format in hb_info string, which rounds to 2 decimals
        // Then hb_numlist parses this string, so 2.4995 becomes 2.50
        // To match legacy, round distance to 2 decimal places before range check
        double rounded_dist = std::round(hbond.distance * 100.0) / 100.0;
        if (rounded_dist >= 2.5 && rounded_dist <= 3.5) {
            num_good_hb++;
        }
    }

    // Legacy: if (num_good_hb >= 2) return -3.0; else return -num_good_hb;
    if (num_good_hb >= 2) {
        return -3.0;
    } else {
        return -static_cast<double>(num_good_hb);
    }
}

int QualityScoreCalculator::calculate_bp_type_id(
    const core::Residue& res1,
    const core::Residue& res2,
    const ValidationResult& result
) const {
    // Match legacy check_wc_wobble_pair logic
    // Legacy: *bpid = -1 initially in calculate_more_bppars
    // Then: if (dir_x > 0.0 && dir_y < 0.0 && dir_z < 0.0) {
    //           check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6]);
    //           if (*bpid == 2) rtn_val[5] -= 2.0;
    //       }
    // check_wc_wobble_pair sets:
    //   *bpid = 1 if fabs(shear) in [1.8, 2.8]
    //   *bpid = 2 if fabs(shear) <= 1.8 && base pair matches WC_LIST

    // Start with -1 (legacy initial value)
    int bp_type_id = -1;

    // Only set to 0 for invalid pairs
    if (!result.is_valid) {
        return 0;
    }

    // Check direction vector condition (matches legacy)
    if (result.dir_x > 0.0 && result.dir_y < 0.0 && result.dir_z < 0.0) {
        // Get reference frames from residues
        if (!res1.reference_frame().has_value() || !res2.reference_frame().has_value()) {
            return bp_type_id; // Keep -1 if frames not available
        }

        // Calculate step parameters for this base pair
        // Legacy: bpstep_par(r2, org[j], r1, org[i], ...)
        // CRITICAL: Legacy reverses y and z columns of r2 when dir_z <= 0
        // See legacy code: r2[l][k] = (k == 1 || dir_z > 0) ? orien[j][koffset + l] :
        // -orien[j][koffset + l];
        core::ReferenceFrame frame1 = res1.reference_frame().value();
        core::ReferenceFrame frame2 = res2.reference_frame().value();

        // Apply frame reversal if dir_z <= 0 (matches legacy logic)
        if (result.dir_z <= 0.0) {
            // Reverse y and z columns (columns 1 and 2 in 0-based indexing) of frame2
            geometry::Matrix3D rot2 = frame2.rotation();
            geometry::Vector3D y_col = rot2.column(1);
            geometry::Vector3D z_col = rot2.column(2);
            rot2.set_column(1, -y_col);
            rot2.set_column(2, -z_col);
            frame2 = core::ReferenceFrame(rot2, frame2.origin());
        }

        // Use frame2 first, frame1 second (matches legacy order)
        core::BasePairStepParameters params = param_calculator_.calculate_step_parameters(frame2, frame1);

        // CRITICAL: Legacy has a bug - it passes wrong parameters to check_wc_wobble_pair!
        // Legacy calls: check_wc_wobble_pair(bpid, bpi, pars[1], pars[2], pars[6])
        // Where pars[1]=Shift, pars[2]=Slide, pars[6]=Twist
        // But function expects: (shear, stretch, opening)
        // So legacy incorrectly uses:
        //   pars[1] (Shift) as shear ❌
        //   pars[2] (Slide) as stretch ❌
        //   pars[6] (Twist) as opening ✅
        // To match legacy output exactly, we must replicate this bug:
        double shear = params.shift;   // BUG: Should be params.slide
        double stretch = params.slide; // BUG: Should be params.rise
        double opening = params.twist; // Correct

        // Get base pair type string (e.g., "AT", "GC")
        // CRITICAL: Use ResidueType to get base letter (more reliable than one_letter_code)
        // Legacy uses: sprintf(bpi, "%c%c", toupper((int) bseq[i]), toupper((int) bseq[j]));
        char base1 = get_base_letter(res1.residue_type());
        char base2 = get_base_letter(res2.residue_type());
        std::string bp_type = std::string(1, base1) + std::string(1, base2);

        // Check stretch and opening thresholds (matches legacy: fabs(stretch) > 2.0 ||
        // fabs(opening) > 60) CRITICAL: Legacy uses fabs(opening) > 60 (not >=), and opening is in
        // degrees
        if (std::abs(stretch) > 2.0 || std::abs(opening) > 60.0) {
            return bp_type_id; // Keep -1
        }

        // Check for wobble pair (fabs(shear) in [1.8, 2.8])
        // CRITICAL: Legacy checks this first, then WC check can overwrite if both conditions met
        if (std::abs(shear) >= 1.8 && std::abs(shear) <= 2.8) {
            bp_type_id = 1; // Wobble
        }

        // Check for Watson-Crick pair (fabs(shear) <= 1.8 AND in WC_LIST)
        // CRITICAL: This can overwrite wobble assignment if shear <= 1.8
        if (std::abs(shear) <= 1.8) {
            for (const auto& wc : WC_LIST) {
                if (bp_type == wc) {
                    bp_type_id = 2; // Watson-Crick
                    break;
                }
            }
            // If not in WC_LIST, keep previous assignment (wobble if set, otherwise -1)
        }
    }

    return bp_type_id;
}

char QualityScoreCalculator::get_base_letter(core::ResidueType type) {
    // Convert ResidueType to one-letter code (matches legacy bseq character)
    switch (type) {
        case core::ResidueType::ADENINE:
            return 'A';
        case core::ResidueType::CYTOSINE:
            return 'C';
        case core::ResidueType::GUANINE:
            return 'G';
        case core::ResidueType::THYMINE:
            return 'T';
        case core::ResidueType::URACIL:
            return 'U';
        case core::ResidueType::INOSINE:
            return 'I';
        case core::ResidueType::PSEUDOURIDINE:
            return 'P';
        default:
            // For modified nucleotides or unknown, try one_letter_code as fallback
            return '?';
    }
}

}  // namespace algorithms
}  // namespace x3dna

