/**
 * @file base_pair.hpp
 * @brief BasePair class representing a base pair between two residues
 */

#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <optional>
#include <nlohmann/json.hpp>
#include <x3dna/core/reference_frame.hpp>
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

/**
 * @struct hydrogen_bond
 * @brief Represents a hydrogen bond in a base pair
 */
struct hydrogen_bond {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance;
    char type;                       // '-' for standard, ' ' for non-standard
    std::optional<size_t> hbond_idx; // Optional index for tracking (assigned when recording)
};

/**
 * @class BasePair
 * @brief Represents a base pair between two nucleotide residues
 */
class BasePair {
private:
    // Private member variables (declared first to help IntelliSense)
    size_t residue_idx1_ = 0;                   // Index of first residue
    size_t residue_idx2_ = 0;                   // Index of second residue
    BasePairType type_ = BasePairType::UNKNOWN; // Base pair type
    std::string bp_type_;                       // Base pair type string (e.g., "CG", "AT")
    std::optional<ReferenceFrame> frame1_;      // Reference frame for first residue
    std::optional<ReferenceFrame> frame2_;      // Reference frame for second residue
    std::vector<hydrogen_bond> hbonds_;         // Hydrogen bonds
    std::optional<size_t> basepair_idx_;        // Optional index for tracking (assigned when recording)
    bool finding_order_swapped_ = false;        // True if indices were swapped during normalization (finding order was j,i not i,j)

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
        if (upper_bp == "AT" || upper_bp == "TA" || upper_bp == "AU" || upper_bp == "UA" ||
            upper_bp == "GC" || upper_bp == "CG" || upper_bp == "IC" || upper_bp == "CI") {
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
     * @brief Default constructor
     */
    BasePair() = default;

    /**
     * @brief Constructor with residue indices and type
     * @param idx1 Index of first residue
     * @param idx2 Index of second residue
     * @param type Base pair type
     */
    BasePair(size_t idx1, size_t idx2, BasePairType type) : residue_idx1_(idx1), residue_idx2_(idx2), type_(type) {}

    // Getters
    size_t residue_idx1() const {
        return residue_idx1_;
    }
    size_t residue_idx2() const {
        return residue_idx2_;
    }
    BasePairType type() const {
        return type_;
    }
    const std::string& bp_type() const {
        return bp_type_;
    }
    const std::vector<hydrogen_bond>& hydrogen_bonds() const {
        return hbonds_;
    }
    std::optional<size_t> basepair_idx() const {
        return basepair_idx_;
    }

    /**
     * @brief Check if the original finding order was swapped during normalization
     * @return True if pair was found in (j,i) order but stored as (i,j) where i < j
     */
    bool finding_order_swapped() const {
        return finding_order_swapped_;
    }

    /**
     * @brief Get reference frame for first residue
     */
    std::optional<ReferenceFrame> frame1() const {
        return frame1_;
    }

    /**
     * @brief Get reference frame for second residue
     */
    std::optional<ReferenceFrame> frame2() const {
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
    std::optional<ReferenceFrame> get_step_frame(bool strand_swapped) const {
        bool use_larger_index_frame = (finding_order_swapped_ != strand_swapped);
        return use_larger_index_frame ? frame2_ : frame1_;
    }

    // Setters
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

    /**
     * @brief Set finding order swapped flag
     * @param swapped True if pair was found in (j,i) order but stored as (i,j) where i < j
     */
    void set_finding_order_swapped(bool swapped) {
        finding_order_swapped_ = swapped;
    }

    /**
     * @brief Set reference frame for first residue
     */
    void set_frame1(const ReferenceFrame& frame) {
        frame1_ = frame;
    }

    /**
     * @brief Set reference frame for second residue
     */
    void set_frame2(const ReferenceFrame& frame) {
        frame2_ = frame;
    }

    /**
     * @brief Add a hydrogen bond
     */
    void add_hydrogen_bond(const hydrogen_bond& hbond) {
        hbonds_.push_back(hbond);
    }

    /**
     * @brief Set all hydrogen bonds at once
     */
    void set_hydrogen_bonds(const std::vector<hydrogen_bond>& hbonds) {
        hbonds_ = hbonds;
    }

    /**
     * @brief Calculate distance between origins of the two reference frames
     * @return Distance in Angstroms
     */
    double origin_distance() const {
        if (!frame1_.has_value() || !frame2_.has_value()) {
            return 0.0;
        }
        return frame1_->origin().distance_to(frame2_->origin());
    }

    /**
     * @brief Calculate plane angle between the two base planes
     * @return Angle in radians
     */
    double plane_angle() const {
        if (!frame1_.has_value() || !frame2_.has_value()) {
            return 0.0;
        }
        Vector3D z1 = frame1_->z_axis();
        Vector3D z2 = frame2_->z_axis();
        double dot = z1.dot(z2);
        // Clamp to [-1, 1] for acos
        if (dot > 1.0)
            dot = 1.0;
        if (dot < -1.0)
            dot = -1.0;
        return std::acos(dot);
    }

    /**
     * @brief Calculate N-N distance (distance between N1/N9 atoms)
     * Note: This is a placeholder - actual implementation requires atom access
     * @return Distance in Angstroms
     */
    double n_n_distance() const {
        // Will be implemented when we have residue access
        return 0.0;
    }

    /**
     * @brief Get direction vector (z-axis dot product)
     * @return Dot product of z-axes (negative for valid base pairs)
     */
    double direction_dot_product() const {
        if (!frame1_.has_value() || !frame2_.has_value()) {
            return 0.0;
        }
        return frame1_->direction_dot_product(frame2_.value());
    }

    /**
     * @brief Convert to legacy JSON format (base_pair record)
     */
    nlohmann::json to_json_legacy() const {
        nlohmann::json j;
        j["type"] = "base_pair";
        j["base_i"] = static_cast<long>(residue_idx1_);
        j["base_j"] = static_cast<long>(residue_idx2_);
        j["bp_type"] = bp_type_;

        if (frame1_.has_value()) {
            j["orien_i"] = frame1_->rotation().to_json_legacy();
            j["org_i"] = frame1_->origin().to_json();
        }

        // For frame2 (orien_j), legacy code applies a sign flip when dir_z <= 0
        // Legacy: r2[l][k] = (k == 1 || dir_z > 0) ? orien[j][...] : -orien[j][...]
        // This negates columns 2 and 3 (y and z axes) when dir_z <= 0
        if (frame2_.has_value()) {
            double dir_z = 0.0;
            if (frame1_.has_value()) {
                dir_z = frame1_->z_axis().dot(frame2_->z_axis());
            }

            if (dir_z <= 0.0 && frame1_.has_value()) {
                // Apply legacy sign flip: negate y and z columns
                geometry::Matrix3D rot2 = frame2_->rotation();
                geometry::Vector3D y_col = rot2.column(1);
                geometry::Vector3D z_col = rot2.column(2);
                rot2.set_column(1, -y_col);
                rot2.set_column(2, -z_col);
                j["orien_j"] = rot2.to_json_legacy();
            } else {
                j["orien_j"] = frame2_->rotation().to_json_legacy();
            }
            j["org_j"] = frame2_->origin().to_json();
        }

        // Direction vector (dot products of corresponding frame axes)
        // Legacy: dir_x = dot(&orien[i][0], &orien[j][0])
        //         dir_y = dot(&orien[i][3], &orien[j][3])
        //         dir_z = dot(&orien[i][6], &orien[j][6])
        //
        // NOTE: Legacy has a bug in json_writer_record_base_pair:
        //   It declares: double dir_xyz_arr[4] = {dir_x, dir_y, dir_z};
        //   But uses: dir_xyz[1], dir_xyz[2], dir_xyz[3] (1-based indexing)
        //   So it actually stores: [dir_y, dir_z, 0.0] (skipping dir_x!)
        //
        // To match legacy exactly, we need to replicate this bug:
        if (frame1_.has_value() && frame2_.has_value()) {
            // Calculate all three direction components (even though we only use y and z)
            double dir_y = frame1_->y_axis().dot(frame2_->y_axis());
            double dir_z = frame1_->z_axis().dot(frame2_->z_axis());
            // Match legacy bug: store [dir_y, dir_z, 0.0] instead of [dir_x, dir_y, dir_z]
            j["dir_xyz"] = nlohmann::json::array({dir_y, dir_z, 0.0});
        }

        // Base pair index (if set)
        if (basepair_idx_.has_value()) {
            j["basepair_idx"] = static_cast<long>(basepair_idx_.value());
        }

        // Hydrogen bonds
        // NOTE: Legacy does NOT store hbonds in base_pair records - they are in separate hbond_list
        // records For exact legacy match, we should not include them in base_pair (Legacy stores
        // hbonds separately in hbond_list records) Commented out to match legacy exactly: if
        // (!hbonds_.empty()) {
        //     j["num_hbonds"] = static_cast<long>(hbonds_.size());
        //     j["hbonds"] = nlohmann::json::array();
        //     ...
        // }

        return j;
    }

    /**
     * @brief Create BasePair from legacy JSON format
     */
    static BasePair from_json_legacy(const nlohmann::json& j) {
        size_t idx1 = j.value("base_i", 0);
        size_t idx2 = j.value("base_j", 0);
        std::string bp_type = j.value("bp_type", "");

        BasePair bp(idx1, idx2, BasePairType::UNKNOWN);
        bp.set_bp_type(bp_type);

        // Parse reference frames
        if (j.contains("orien_i") && j.contains("org_i")) {
            nlohmann::json frame1_json;
            frame1_json["orien"] = j["orien_i"];
            frame1_json["org"] = j["org_i"];
            bp.set_frame1(ReferenceFrame::from_json_legacy(frame1_json));
        }

        if (j.contains("orien_j") && j.contains("org_j")) {
            nlohmann::json frame2_json;
            frame2_json["orien"] = j["orien_j"];
            frame2_json["org"] = j["org_j"];
            bp.set_frame2(ReferenceFrame::from_json_legacy(frame2_json));
        }

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
    nlohmann::json to_json() const {
        nlohmann::json j;
        j["residue_idx1"] = residue_idx1_;
        j["residue_idx2"] = residue_idx2_;
        j["bp_type"] = bp_type_;
        if (basepair_idx_.has_value()) {
            j["basepair_idx"] = basepair_idx_.value();
        }
        if (frame1_.has_value()) {
            j["frame1"] = frame1_->to_json();
        }
        if (frame2_.has_value()) {
            j["frame2"] = frame2_->to_json();
        }
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
    static BasePair from_json(const nlohmann::json& j) {
        size_t idx1 = j.value("residue_idx1", 0);
        size_t idx2 = j.value("residue_idx2", 0);
        std::string bp_type = j.value("bp_type", "");

        BasePair bp(idx1, idx2, BasePairType::UNKNOWN);
        bp.set_bp_type(bp_type);

        if (j.contains("frame1")) {
            bp.set_frame1(ReferenceFrame::from_json(j["frame1"]));
        }
        if (j.contains("frame2")) {
            bp.set_frame2(ReferenceFrame::from_json(j["frame2"]));
        }

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
