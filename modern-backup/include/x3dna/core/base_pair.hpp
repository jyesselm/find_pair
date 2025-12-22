/**
 * @file base_pair.hpp
 * @brief BasePair class representing a base pair between two residues
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <optional>  // Still needed for basepair_idx_
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/hydrogen_bond.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

namespace x3dna {
namespace core {

using geometry::Vector3D;

/**
 * @enum BasePairType
 * @brief Type of base pair
 */
enum class BasePairType { WATSON_CRICK, WOBBLE, HOOGSTEEN, UNKNOWN };

// Type alias for backward compatibility with code using the old struct name
using hydrogen_bond = HydrogenBond;

/**
 * @class BasePair
 * @brief Represents a base pair between two nucleotide residues
 *
 * Reference frames are always present - a BasePair cannot exist without frames.
 * Use the constructor with frames:
 *   BasePair(idx1, idx2, frame1, frame2)
 *
 * For JSON deserialization, identity frames are used as defaults.
 */
class BasePair {
private:
    // Private member variables (declared first to help IntelliSense)
    size_t residue_idx1_ = 0;                   // Index of first residue
    size_t residue_idx2_ = 0;                   // Index of second residue
    BasePairType type_ = BasePairType::UNKNOWN; // Base pair type
    std::string bp_type_;                       // Base pair type string (e.g., "CG", "AT")
    ReferenceFrame frame1_;                     // Reference frame for first residue (identity if not set)
    ReferenceFrame frame2_;                     // Reference frame for second residue (identity if not set)
    std::vector<HydrogenBond> hbonds_;          // Hydrogen bonds
    std::optional<size_t> basepair_idx_;        // Optional index for tracking (assigned when recording)
    bool finding_order_swapped_ = false; // True if indices were swapped during normalization (finding order was j,i not
                                         // i,j)

    /**
     * @brief Set base pair type from string and update enum
     *
     * Note: This is a simplified classification based on the bp_type string.
     * The original code (check_wc_wobble_pair) also uses geometric parameters
     * (shear, stretch, opening) for classification. This function matches
     * the WC_LIST from the original code: "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
     *
     * Wobble pairs (GT, TG, GU, UG) are not in WC_LIST and would be classified
     * as wobble based on geometry in the original code.
     */
    void update_type_from_bp_type() {
        // Legacy uppercases base types before comparing (cmn_fncs.c:4529)
        // So "Gc" should match "GC", "Ug" should match "UG"
        std::string upper_bp;
        for (char c : bp_type_) {
            upper_bp += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }

        // Watson-Crick pairs from WC_LIST (excluding "XX" placeholder)
        if (upper_bp == "AT" || upper_bp == "TA" || upper_bp == "AU" || upper_bp == "UA" || upper_bp == "GC" ||
            upper_bp == "CG" || upper_bp == "IC" || upper_bp == "CI") {
            type_ = BasePairType::WATSON_CRICK;
        }
        // Wobble pairs (not in WC_LIST, but commonly occur)
        // Note: Original code classifies these based on geometry, not just string
        else if (upper_bp == "GT" || upper_bp == "TG" || upper_bp == "GU" || upper_bp == "UG") {
            type_ = BasePairType::WOBBLE;
        } else {
            type_ = BasePairType::UNKNOWN;
        }
    }

public:
    /**
     * @brief Constructor with residue indices and required frames
     * @param idx1 Index of first residue
     * @param idx2 Index of second residue
     * @param frame1 Reference frame for first residue
     * @param frame2 Reference frame for second residue
     * @param type Base pair type (default: UNKNOWN)
     */
    BasePair(size_t idx1, size_t idx2, const ReferenceFrame& frame1, const ReferenceFrame& frame2,
             BasePairType type = BasePairType::UNKNOWN)
        : residue_idx1_(idx1), residue_idx2_(idx2), type_(type), frame1_(frame1), frame2_(frame2) {}

private:
    // Private default constructor for JSON deserialization only
    BasePair() = default;
    friend class nlohmann::basic_json<>;  // Allow nlohmann::json to use default constructor

public:

    // Getters
    [[nodiscard]] size_t residue_idx1() const {
        return residue_idx1_;
    }
    [[nodiscard]] size_t residue_idx2() const {
        return residue_idx2_;
    }
    [[nodiscard]] BasePairType type() const {
        return type_;
    }
    [[nodiscard]] const std::string& bp_type() const {
        return bp_type_;
    }
    [[nodiscard]] const std::vector<hydrogen_bond>& hydrogen_bonds() const {
        return hbonds_;
    }
    [[nodiscard]] std::optional<size_t> basepair_idx() const {
        return basepair_idx_;
    }

    /**
     * @brief Check if the original finding order was swapped during normalization
     * @return True if pair was found in (j,i) order but stored as (i,j) where i < j
     */
    [[nodiscard]] bool finding_order_swapped() const {
        return finding_order_swapped_;
    }

    /**
     * @brief Get reference frame for first residue
     */
    [[nodiscard]] const ReferenceFrame& frame1() const {
        return frame1_;
    }

    /**
     * @brief Get reference frame for second residue
     */
    [[nodiscard]] const ReferenceFrame& frame2() const {
        return frame2_;
    }

    /**
     * @brief Get the reference frame to use for step parameter calculation
     * @param strand_swapped Whether the five2three algorithm swapped this pair's strand direction
     * @return The appropriate frame for step calculation
     *
     * This encapsulates the frame selection logic that matches legacy behavior:
     * - Legacy stores pairs in finding order (searching_residue, best_partner)
     * - Modern normalizes to (smaller_index, larger_index) and tracks finding_order_swapped
     * - Legacy's five2three may swap frames based on helix direction (strand_swapped)
     * - The correct frame is determined by XOR of these two flags
     */
    [[nodiscard]] const ReferenceFrame& get_step_frame(bool strand_swapped) const {
        bool use_larger_index_frame = (finding_order_swapped_ != strand_swapped);
        return use_larger_index_frame ? frame2_ : frame1_;
    }

    // Modification methods
    // BasePair is constructed incrementally during pair finding workflow.
    // Data comes from multiple sources (validation, frame calculation, H-bond detection)
    // so setters are required for the incremental construction pattern.

    void set_residue_idx1(size_t idx) {
        residue_idx1_ = idx;
    }
    void set_residue_idx2(size_t idx) {
        residue_idx2_ = idx;
    }
    void set_type(BasePairType type) {
        type_ = type;
    }
    void set_bp_type(const std::string& bp_type) {
        bp_type_ = bp_type;
        update_type_from_bp_type();
    }
    void set_basepair_idx(size_t idx) {
        basepair_idx_ = idx;
    }
    void set_finding_order_swapped(bool swapped) {
        finding_order_swapped_ = swapped;
    }
    void set_frame1(const ReferenceFrame& frame) {
        frame1_ = frame;
    }
    void set_frame2(const ReferenceFrame& frame) {
        frame2_ = frame;
    }
    void add_hydrogen_bond(const hydrogen_bond& hbond) {
        hbonds_.push_back(hbond);
    }
    void set_hydrogen_bonds(const std::vector<hydrogen_bond>& hbonds) {
        hbonds_ = hbonds;
    }

    /**
     * @brief Calculate distance between origins of the two reference frames
     * @return Distance in Angstroms
     */
    [[nodiscard]] double origin_distance() const {
        return frame1_.origin().distance_to(frame2_.origin());
    }

    /**
     * @brief Calculate plane angle between the two base planes
     * @return Angle in radians
     */
    [[nodiscard]] double plane_angle() const {
        Vector3D z1 = frame1_.z_axis();
        Vector3D z2 = frame2_.z_axis();
        double dot = std::clamp(z1.dot(z2), -1.0, 1.0);
        return std::acos(dot);
    }

    /**
     * @brief Calculate N-N distance (distance between N1/N9 atoms)
     * Note: This is a placeholder - actual implementation requires atom access
     * @return Distance in Angstroms
     */
    [[nodiscard]] double n_n_distance() const {
        // Will be implemented when we have residue access
        return 0.0;
    }

    /**
     * @brief Get direction vector (z-axis dot product)
     * @return Dot product of z-axes (negative for valid base pairs)
     */
    [[nodiscard]] double direction_dot_product() const {
        return frame1_.direction_dot_product(frame2_);
    }

    /**
     * @brief Convert to legacy JSON format (base_pair record)
     */
    [[nodiscard]] nlohmann::json to_json_legacy() const {
        nlohmann::json j;
        j["type"] = "base_pair";
        j["base_i"] = static_cast<long>(residue_idx1_);
        j["base_j"] = static_cast<long>(residue_idx2_);
        j["bp_type"] = bp_type_;

        // Frame 1
        j["orien_i"] = frame1_.rotation().to_json_legacy();
        j["org_i"] = frame1_.origin().to_json();

        // For frame2 (orien_j), legacy code applies a sign flip when dir_z <= 0
        // Legacy: r2[l][k] = (k == 1 || dir_z > 0) ? orien[j][...] : -orien[j][...]
        // This negates columns 2 and 3 (y and z axes) when dir_z <= 0
        double dir_z = frame1_.z_axis().dot(frame2_.z_axis());

        if (dir_z <= 0.0) {
            // Apply legacy sign flip: negate y and z columns
            geometry::Matrix3D rot2 = frame2_.rotation();
            geometry::Vector3D y_col = rot2.column(1);
            geometry::Vector3D z_col = rot2.column(2);
            rot2.set_column(1, -y_col);
            rot2.set_column(2, -z_col);
            j["orien_j"] = rot2.to_json_legacy();
        } else {
            j["orien_j"] = frame2_.rotation().to_json_legacy();
        }
        j["org_j"] = frame2_.origin().to_json();

        // Direction vector (dot products of corresponding frame axes)
        // NOTE: Legacy has a bug - stores [dir_y, dir_z, 0.0] instead of [dir_x, dir_y, dir_z]
        double dir_y = frame1_.y_axis().dot(frame2_.y_axis());
        j["dir_xyz"] = nlohmann::json::array({dir_y, dir_z, 0.0});

        // Base pair index (if set)
        if (basepair_idx_.has_value()) {
            j["basepair_idx"] = static_cast<long>(basepair_idx_.value());
        }

        return j;
    }

    /**
     * @brief Create BasePair from legacy JSON format
     */
    [[nodiscard]] static BasePair from_json_legacy(const nlohmann::json& j) {
        size_t idx1 = j.value("base_i", 0);
        size_t idx2 = j.value("base_j", 0);
        std::string bp_type_str = j.value("bp_type", "");

        // Parse reference frames (use identity if not present)
        ReferenceFrame frame1, frame2;
        if (j.contains("orien_i") && j.contains("org_i")) {
            nlohmann::json frame1_json;
            frame1_json["orien"] = j["orien_i"];
            frame1_json["org"] = j["org_i"];
            frame1 = ReferenceFrame::from_json_legacy(frame1_json);
        }
        if (j.contains("orien_j") && j.contains("org_j")) {
            nlohmann::json frame2_json;
            frame2_json["orien"] = j["orien_j"];
            frame2_json["org"] = j["org_j"];
            frame2 = ReferenceFrame::from_json_legacy(frame2_json);
        }

        BasePair bp(idx1, idx2, frame1, frame2);
        bp.set_bp_type(bp_type_str);

        // Parse base pair index (if present)
        if (j.contains("basepair_idx")) {
            bp.set_basepair_idx(j.value("basepair_idx", 0));
        }

        // Parse hydrogen bonds
        if (j.contains("hbonds") && j["hbonds"].is_array()) {
            for (const auto& hb_json : j["hbonds"]) {
                hydrogen_bond hbond;
                hbond.donor_atom = hb_json.value("donor_atom", "");
                hbond.acceptor_atom = hb_json.value("acceptor_atom", "");
                hbond.distance = hb_json.value("distance", 0.0);
                std::string type_str = hb_json.value("type", " ");
                hbond.type = type_str.empty() ? ' ' : type_str[0];
                if (hb_json.contains("hbond_idx")) {
                    hbond.hbond_idx = hb_json.value("hbond_idx", 0);
                }
                bp.add_hydrogen_bond(hbond);
            }
        }

        return bp;
    }

    /**
     * @brief Convert to modern JSON format
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;
        j["residue_idx1"] = residue_idx1_;
        j["residue_idx2"] = residue_idx2_;
        j["bp_type"] = bp_type_;
        if (basepair_idx_.has_value()) {
            j["basepair_idx"] = basepair_idx_.value();
        }
        j["frame1"] = frame1_.to_json();
        j["frame2"] = frame2_.to_json();
        j["hydrogen_bonds"] = nlohmann::json::array();
        for (size_t i = 0; i < hbonds_.size(); ++i) {
            const auto& hbond = hbonds_[i];
            nlohmann::json hb_json;
            hb_json["donor_atom"] = hbond.donor_atom;
            hb_json["acceptor_atom"] = hbond.acceptor_atom;
            hb_json["distance"] = hbond.distance;
            hb_json["type"] = std::string(1, hbond.type);
            if (hbond.hbond_idx.has_value()) {
                hb_json["hbond_idx"] = hbond.hbond_idx.value();
            } else {
                hb_json["hbond_idx"] = i;
            }
            j["hydrogen_bonds"].push_back(hb_json);
        }
        return j;
    }

    /**
     * @brief Create BasePair from modern JSON format
     */
    [[nodiscard]] static BasePair from_json(const nlohmann::json& j) {
        size_t idx1 = j.value("residue_idx1", 0);
        size_t idx2 = j.value("residue_idx2", 0);
        std::string bp_type_str = j.value("bp_type", "");

        // Parse frames (use identity if not present)
        ReferenceFrame frame1, frame2;
        if (j.contains("frame1")) {
            frame1 = ReferenceFrame::from_json(j["frame1"]);
        }
        if (j.contains("frame2")) {
            frame2 = ReferenceFrame::from_json(j["frame2"]);
        }

        BasePair bp(idx1, idx2, frame1, frame2);
        bp.set_bp_type(bp_type_str);

        if (j.contains("basepair_idx")) {
            bp.set_basepair_idx(j.value("basepair_idx", 0));
        }

        if (j.contains("hydrogen_bonds") && j["hydrogen_bonds"].is_array()) {
            for (const auto& hb_json : j["hydrogen_bonds"]) {
                hydrogen_bond hbond;
                hbond.donor_atom = hb_json.value("donor_atom", "");
                hbond.acceptor_atom = hb_json.value("acceptor_atom", "");
                hbond.distance = hb_json.value("distance", 0.0);
                std::string type_str = hb_json.value("type", " ");
                hbond.type = type_str.empty() ? ' ' : type_str[0];
                if (hb_json.contains("hbond_idx")) {
                    hbond.hbond_idx = hb_json.value("hbond_idx", 0);
                }
                bp.add_hydrogen_bond(hbond);
            }
        }

        return bp;
    }
};

} // namespace core
} // namespace x3dna
