/**
 * @file parameters.hpp
 * @brief Centralized parameter management - single source of truth
 *
 * Parameters can be:
 * 1. Used with compile-time defaults (always available)
 * 2. Loaded from resources/config/parameters.json at runtime
 *
 * The JSON file is also readable by Python for consistent values across languages.
 */

#pragma once

#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <optional>

namespace x3dna {
namespace config {

/**
 * @struct ValidationDistanceParams
 * @brief Distance thresholds for base pair validation
 */
struct ValidationDistanceParams {
    double min_dorg = 0.0;
    double max_dorg = 15.0;
    double min_dv = 0.0;
    double max_dv = 2.5;
    double min_dNN = 4.5;
    double max_dNN = 1e18;
};

/**
 * @struct ValidationAngleParams
 * @brief Angle thresholds for base pair validation
 */
struct ValidationAngleParams {
    double min_plane_angle = 0.0;
    double max_plane_angle = 65.0;
};

/**
 * @struct HydrogenBondDetectionParams
 * @brief H-bond detection parameters
 */
struct HydrogenBondDetectionParams {
    double hb_lower = 1.8;
    double hb_dist1 = 4.0;
    double hb_dist2 = 0.0;
    std::string hb_atoms = ".O.N";
};

/**
 * @struct HydrogenBondThresholds
 * @brief H-bond validation thresholds
 */
struct HydrogenBondThresholds {
    double good_min = 2.5;
    double good_max = 3.5;
    double filter_max = 3.6;
    double nonstandard_min = 2.6;
    double nonstandard_max = 3.2;
    double default_dist2 = 4.5;
};

/**
 * @struct QualityScoreParams
 * @brief Quality score calculation parameters
 */
struct QualityScoreParams {
    double d_v_weight = 2.0;
    double plane_angle_divisor = 20.0;
    double wc_bonus = 2.0;
};

/**
 * @struct NucleotideParams
 * @brief Nucleotide identification parameters
 */
struct NucleotideParams {
    double rmsd_cutoff = 0.2618;
    double dnn_fallback = 1e10;
    double bond_distance = 2.0;
    double min_atom_distance = 0.1;
};

/**
 * @struct HelixParams
 * @brief Helix organization parameters
 */
struct HelixParams {
    double helix_break = 7.8;
    double end_stack_xang = 125.0;
    double std_curved = 0.6;
};

/**
 * @struct MiscParams
 * @brief Miscellaneous parameters
 */
struct MiscParams {
    std::string alt_list = "A1";
    double o3p_dist = 4.5;
    double xbig = 1e18;
    double gamut = 5e8;
    double overlap_threshold = 0.01;
};

/**
 * @class Parameters
 * @brief Singleton for accessing all algorithm parameters
 *
 * Usage:
 *   auto& params = Parameters::instance();
 *   double threshold = params.validation_distance().max_dorg;
 *   double hb_min = params.hbond_thresholds().good_min;
 */
class Parameters {
public:
    [[nodiscard]] static Parameters& instance() {
        static Parameters instance;
        return instance;
    }

    // Load from JSON file (optional - defaults are always available)
    bool load_from_file(const std::filesystem::path& path) {
        try {
            std::ifstream f(path);
            if (!f.is_open()) return false;

            nlohmann::json j;
            f >> j;
            load_from_json(j);
            loaded_from_file_ = true;
            return true;
        } catch (...) {
            return false;
        }
    }

    void load_from_json(const nlohmann::json& j) {
        // Validation distance
        if (j.contains("validation") && j["validation"].contains("distance")) {
            const auto& d = j["validation"]["distance"];
            validation_distance_.min_dorg = d.value("min_dorg", validation_distance_.min_dorg);
            validation_distance_.max_dorg = d.value("max_dorg", validation_distance_.max_dorg);
            validation_distance_.min_dv = d.value("min_dv", validation_distance_.min_dv);
            validation_distance_.max_dv = d.value("max_dv", validation_distance_.max_dv);
            validation_distance_.min_dNN = d.value("min_dNN", validation_distance_.min_dNN);
            validation_distance_.max_dNN = d.value("max_dNN", validation_distance_.max_dNN);
        }

        // Validation angle
        if (j.contains("validation") && j["validation"].contains("angle")) {
            const auto& a = j["validation"]["angle"];
            validation_angle_.min_plane_angle = a.value("min_plane_angle", validation_angle_.min_plane_angle);
            validation_angle_.max_plane_angle = a.value("max_plane_angle", validation_angle_.max_plane_angle);
        }

        // Overlap threshold
        if (j.contains("validation")) {
            misc_.overlap_threshold = j["validation"].value("overlap_threshold", misc_.overlap_threshold);
        }

        // H-bond detection
        if (j.contains("hydrogen_bond") && j["hydrogen_bond"].contains("detection")) {
            const auto& d = j["hydrogen_bond"]["detection"];
            hbond_detection_.hb_lower = d.value("hb_lower", hbond_detection_.hb_lower);
            hbond_detection_.hb_dist1 = d.value("hb_dist1", hbond_detection_.hb_dist1);
            hbond_detection_.hb_dist2 = d.value("hb_dist2", hbond_detection_.hb_dist2);
            hbond_detection_.hb_atoms = d.value("hb_atoms", hbond_detection_.hb_atoms);
        }

        // H-bond thresholds
        if (j.contains("hydrogen_bond") && j["hydrogen_bond"].contains("thresholds")) {
            const auto& t = j["hydrogen_bond"]["thresholds"];
            hbond_thresholds_.good_min = t.value("good_min", hbond_thresholds_.good_min);
            hbond_thresholds_.good_max = t.value("good_max", hbond_thresholds_.good_max);
            hbond_thresholds_.filter_max = t.value("filter_max", hbond_thresholds_.filter_max);
            hbond_thresholds_.nonstandard_min = t.value("nonstandard_min", hbond_thresholds_.nonstandard_min);
            hbond_thresholds_.nonstandard_max = t.value("nonstandard_max", hbond_thresholds_.nonstandard_max);
            hbond_thresholds_.default_dist2 = t.value("default_dist2", hbond_thresholds_.default_dist2);
        }

        // Quality score
        if (j.contains("quality_score")) {
            const auto& q = j["quality_score"];
            quality_score_.d_v_weight = q.value("d_v_weight", quality_score_.d_v_weight);
            quality_score_.plane_angle_divisor = q.value("plane_angle_divisor", quality_score_.plane_angle_divisor);
            quality_score_.wc_bonus = q.value("wc_bonus", quality_score_.wc_bonus);
        }

        // Nucleotide
        if (j.contains("nucleotide")) {
            const auto& n = j["nucleotide"];
            nucleotide_.rmsd_cutoff = n.value("rmsd_cutoff", nucleotide_.rmsd_cutoff);
            nucleotide_.dnn_fallback = n.value("dnn_fallback", nucleotide_.dnn_fallback);
            nucleotide_.bond_distance = n.value("bond_distance", nucleotide_.bond_distance);
            nucleotide_.min_atom_distance = n.value("min_atom_distance", nucleotide_.min_atom_distance);
        }

        // Helix
        if (j.contains("helix")) {
            const auto& h = j["helix"];
            helix_.helix_break = h.value("helix_break", helix_.helix_break);
            helix_.end_stack_xang = h.value("end_stack_xang", helix_.end_stack_xang);
            helix_.std_curved = h.value("std_curved", helix_.std_curved);
        }

        // Misc
        if (j.contains("misc")) {
            const auto& m = j["misc"];
            misc_.alt_list = m.value("alt_list", misc_.alt_list);
            misc_.o3p_dist = m.value("o3p_dist", misc_.o3p_dist);
            misc_.xbig = m.value("xbig", misc_.xbig);
            misc_.gamut = m.value("gamut", misc_.gamut);
        }
    }

    // Accessors
    [[nodiscard]] const ValidationDistanceParams& validation_distance() const { return validation_distance_; }
    [[nodiscard]] const ValidationAngleParams& validation_angle() const { return validation_angle_; }
    [[nodiscard]] const HydrogenBondDetectionParams& hbond_detection() const { return hbond_detection_; }
    [[nodiscard]] const HydrogenBondThresholds& hbond_thresholds() const { return hbond_thresholds_; }
    [[nodiscard]] const QualityScoreParams& quality_score() const { return quality_score_; }
    [[nodiscard]] const NucleotideParams& nucleotide() const { return nucleotide_; }
    [[nodiscard]] const HelixParams& helix() const { return helix_; }
    [[nodiscard]] const MiscParams& misc() const { return misc_; }

    [[nodiscard]] bool loaded_from_file() const { return loaded_from_file_; }

    // Delete copy/move
    Parameters(const Parameters&) = delete;
    Parameters& operator=(const Parameters&) = delete;

private:
    Parameters() = default;

    ValidationDistanceParams validation_distance_;
    ValidationAngleParams validation_angle_;
    HydrogenBondDetectionParams hbond_detection_;
    HydrogenBondThresholds hbond_thresholds_;
    QualityScoreParams quality_score_;
    NucleotideParams nucleotide_;
    HelixParams helix_;
    MiscParams misc_;
    bool loaded_from_file_ = false;
};

// Convenience namespace for compile-time defaults (backwards compatible)
namespace defaults {
    // Validation
    constexpr double MAX_DORG = 15.0;
    constexpr double MAX_DV = 2.5;
    constexpr double MAX_PLANE_ANGLE = 65.0;
    constexpr double MIN_DNN = 4.5;
    constexpr double OVERLAP_THRESHOLD = 0.01;

    // H-bond detection
    constexpr double HB_LOWER = 1.8;
    constexpr double HB_DIST1 = 4.0;
    constexpr double HB_DEFAULT_DIST2 = 4.5;  // Default for find_hydrogen_bonds

    // H-bond validation thresholds
    constexpr double HB_GOOD_MIN = 2.5;
    constexpr double HB_GOOD_MAX = 3.5;
    constexpr double HB_FILTER_MAX = 3.6;
    constexpr double HB_NONSTANDARD_MIN = 2.6;
    constexpr double HB_NONSTANDARD_MAX = 3.2;
    constexpr int HB_LINKAGE_CONFLICT = 18;

    // Quality score
    constexpr double D_V_WEIGHT = 2.0;
    constexpr double PLANE_ANGLE_DIVISOR = 20.0;
    constexpr double WC_QUALITY_BONUS = 2.0;

    // Nucleotide
    constexpr double NT_RMSD_CUTOFF = 0.2618;
    constexpr double DNN_FALLBACK = 1e10;

    // Helix
    constexpr double HELIX_BREAK = 7.8;
    constexpr double END_STACK_XANG = 125.0;

    // Misc
    constexpr double XBIG = 1e18;
    constexpr double GAMUT = 5e8;
}

} // namespace config
} // namespace x3dna
