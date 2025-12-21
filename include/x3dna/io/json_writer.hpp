/**
 * @file json_writer.hpp
 * @brief Modern C++ JSON writer for X3DNA calculation records
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>
#include <map>
#include <array>
#include <set>
#include <nlohmann/json.hpp>
#include <x3dna/core/structure.hpp>
#include <x3dna/core/base_pair.hpp>
#include <x3dna/core/reference_frame.hpp>
#include <x3dna/core/parameters.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <x3dna/algorithms/helix_organizer.hpp>
#include <x3dna/geometry/vector3d.hpp>
#include <x3dna/geometry/matrix3d.hpp>

namespace x3dna {
namespace io {

/**
 * @class JsonWriter
 * @brief Writes calculation records to JSON format (legacy and modern)
 *
 * This class provides a modern C++ interface for writing X3DNA calculation
 * records to JSON format. It matches the legacy JSON format exactly while
 * providing a clean, type-safe API.
 */
class JsonWriter {
public:
    /**
     * @brief Constructor
     * @param pdb_file Path to the PDB file being processed
     */
    explicit JsonWriter(const std::filesystem::path& pdb_file);

    /**
     * @brief Destructor - automatically finalizes if not already done
     */
    ~JsonWriter();

    /**
     * @brief Get the JSON object (can be called before finalize)
     * @return Reference to the JSON object
     */
    nlohmann::json& json() {
        return json_;
    }
    const nlohmann::json& json() const {
        return json_;
    }

    /**
     * @brief Get base pairs in the order they were recorded
     * @return Vector of base pairs in recording order (for step calculation)
     */
    const std::vector<core::BasePair>& ordered_base_pairs() const {
        return ordered_base_pairs_;
    }

    /**
     * @brief Get the JSON as a string
     * @param pretty_print Whether to format with indentation
     * @return JSON string
     */
    std::string to_string(bool pretty_print = true) const;

    /**
     * @brief Write JSON to file
     * @param output_path Path to output JSON file or directory
     * @param pretty_print Whether to format with indentation
     */
    void write_to_file(const std::filesystem::path& output_path, bool pretty_print = true) const;

    /**
     * @brief Write split files to record-type-specific directories
     * @param output_dir Base output directory (e.g., data/json)
     * @param pretty_print Whether to format with indentation
     */
    void write_split_files(const std::filesystem::path& output_dir, bool pretty_print = true) const;

    // Record writing methods - matching legacy format

    /**
     * @brief Record PDB atom data
     * @param structure Structure containing atoms
     */
    void record_pdb_atoms(const core::Structure& structure);

    /**
     * @brief Record residue indices (seidx) - maps residues to atom ranges
     * @param structure Structure containing residues and atoms
     */
    void record_residue_indices(const core::Structure& structure);

    /**
     * @brief Record base frame calculation
     * @param residue_idx Residue index (0-based)
     * @param base_type Base type character (A, C, G, T, U)
     * @param standard_template Path to standard template file
     * @param rms_fit RMS fit value
     * @param matched_atoms Vector of matched atom names
     * @param residue_name Residue name (e.g., "  G")
     * @param chain_id Chain identifier
     * @param residue_seq Residue sequence number
     * @param insertion Insertion code (default "")
     */
    void record_base_frame_calc(size_t residue_idx, char base_type, const std::filesystem::path& standard_template,
                                double rms_fit, const std::vector<std::string>& matched_atoms,
                                const std::string& residue_name = "", const std::string& chain_id = "",
                                int residue_seq = 0, const std::string& insertion = "");

    /**
     * @brief Record least-squares fitting details
     * @param residue_idx Residue index (0-based)
     * @param num_points Number of points used in fitting
     * @param rms_fit RMS fit value
     * @param rotation_matrix 3x3 rotation matrix
     * @param translation Translation vector
     * @param residue_name Residue name (e.g., "  G")
     * @param chain_id Chain identifier
     * @param residue_seq Residue sequence number
     * @param insertion Insertion code (default "")
     */
    void record_ls_fitting(size_t residue_idx, size_t num_points, double rms_fit,
                           const geometry::Matrix3D& rotation_matrix, const geometry::Vector3D& translation,
                           const std::string& residue_name = "", const std::string& chain_id = "", int residue_seq = 0,
                           const std::string& insertion = "");

    /**
     * @brief Record frame calculation with matched coordinates
     * @param residue_idx Residue index (0-based)
     * @param base_type Base type character
     * @param template_file Path to template file
     * @param rms_fit RMS fit value
     * @param matched_std_xyz Standard template coordinates (num_matched x 3)
     * @param matched_exp_xyz Experimental coordinates (num_matched x 3)
     * @param residue_name Residue name (e.g., "  G")
     * @param chain_id Chain identifier
     * @param residue_seq Residue sequence number
     * @param insertion Insertion code (default "")
     */
    void record_frame_calc(size_t residue_idx, char base_type, const std::filesystem::path& template_file,
                           double rms_fit, const std::vector<geometry::Vector3D>& matched_std_xyz,
                           const std::vector<geometry::Vector3D>& matched_exp_xyz, const std::string& residue_name = "",
                           const std::string& chain_id = "", int residue_seq = 0, const std::string& insertion = "");

    /**
     * @brief Record base pair information
     * @param pair Base pair object
     */
    void record_base_pair(const core::BasePair& pair);

    /**
     * @brief Record base pair step parameters
     * @param bp_idx1 First base pair index (0-based)
     * @param bp_idx2 Second base pair index (0-based)
     * @param params Base pair step parameters
     * @param pair1 Optional first base pair (to extract res_ids)
     * @param pair2 Optional second base pair (to extract res_ids)
     */
    void record_bpstep_params(size_t bp_idx1, size_t bp_idx2, const core::BasePairStepParameters& params,
                              const core::BasePair* pair1 = nullptr, const core::BasePair* pair2 = nullptr);

    /**
     * @brief Record helical parameters
     * @param bp_idx1 First base pair index (0-based)
     * @param bp_idx2 Second base pair index (0-based)
     * @param params Helical parameters
     * @param pair1 Optional first base pair (to extract res_ids)
     * @param pair2 Optional second base pair (to extract res_ids)
     */
    void record_helical_params(size_t bp_idx1, size_t bp_idx2, const core::HelicalParameters& params,
                               const core::BasePair* pair1 = nullptr, const core::BasePair* pair2 = nullptr);

    /**
     * @brief Record all reference frames
     * @param structure Structure containing residues with frames
     */
    void record_all_ref_frames(const core::Structure& structure);

    /**
     * @brief Record removed atom during parsing
     * @param pdb_line Original PDB line
     * @param reason Reason for removal
     * @param atom_serial Atom serial number
     * @param atom_name Atom name
     * @param residue_name Residue name
     * @param chain_id Chain ID
     * @param residue_seq Residue sequence number
     * @param xyz Coordinates (can be nullptr)
     * @param model_num Model number
     */
    void record_removed_atom(const std::string& pdb_line, const std::string& reason, int atom_serial,
                             const std::string& atom_name, const std::string& residue_name,
                             const std::string& chain_id, int residue_seq, const geometry::Vector3D* xyz, int model_num);

    /**
     * @brief Record summary of removed atoms
     * @param num_removed Total number of removed atoms
     */
    void record_removed_atoms_summary(size_t num_removed);

    /**
     * @brief Record pair validation results (matches legacy format)
     * @param base_i First base index (1-based for legacy format)
     * @param base_j Second base index (1-based for legacy format)
     * @param is_valid Whether pair is valid
     * @param bp_type_id Base pair type ID
     * @param dir_x X component of direction vector
     * @param dir_y Y component of direction vector
     * @param dir_z Z component of direction vector
     * @param rtn_val Return values array (5 values: [0]=dorg, [1]=d_v, [2]=plane_angle, [3]=dNN,
     * [4]=quality_score)
     * @param params Validation parameters (for thresholds)
     * @param res_id_i Unique residue identifier for first base (optional)
     * @param res_id_j Unique residue identifier for second base (optional)
     */
    void record_pair_validation(size_t base_i, size_t base_j, bool is_valid, int bp_type_id, double dir_x, double dir_y,
                                double dir_z, const std::array<double, 5>& rtn_val,
                                const algorithms::ValidationParameters& params,
                                const std::string& res_id_i = "", const std::string& res_id_j = "");

    /**
     * @brief Record distance checks (matches legacy format)
     * @param base_i First base index (1-based for legacy format)
     * @param base_j Second base index (1-based for legacy format)
     * @param dorg Distance between origins
     * @param dNN Distance between N1/N9 atoms
     * @param plane_angle Angle between z-axes
     * @param d_v Vertical distance
     * @param overlap_area Overlap area
     * @param res_id_i Unique residue identifier for first base (optional)
     * @param res_id_j Unique residue identifier for second base (optional)
     */
    void record_distance_checks(size_t base_i, size_t base_j, double dorg, double dNN, double plane_angle, double d_v,
                                double overlap_area, const std::string& res_id_i = "", const std::string& res_id_j = "");

    /**
     * @brief Record hydrogen bond list
     * @param base_i First base index (0-based)
     * @param base_j Second base index (0-based)
     * @param hbonds Vector of hydrogen bond information
     * @param res_id_i Unique residue identifier for first base (optional)
     * @param res_id_j Unique residue identifier for second base (optional)
     */
    void record_hbond_list(size_t base_i, size_t base_j, const std::vector<core::hydrogen_bond>& hbonds,
                           const std::string& res_id_i = "", const std::string& res_id_j = "");

    /**
     * @brief Record the original base pair selection from find_bestpair
     * @param selected_pairs Vector of (base_i, base_j) pairs (1-based legacy indices)
     *
     * This records the pairs actually selected by find_bestpair (mutual best matches),
     * which is the original base pair identification before any reordering.
     */
    void record_find_bestpair_selection(const std::vector<std::pair<size_t, size_t>>& selected_pairs);

    /**
     * @brief Record best partner candidates for debugging (all candidates considered)
     * @param res_i Residue index (legacy 1-based)
     * @param candidates Vector of candidate information: (res_j, is_eligible, score, bp_type_id)
     * @param best_j Best partner residue index (legacy 1-based)
     * @param best_score Best score
     */
    void record_best_partner_candidates(int res_i, const std::vector<std::tuple<int, bool, double, int>>& candidates,
                                        int best_j, double best_score);

    /**
     * @brief Record mutual best decision
     * @param res_i First residue index (legacy 1-based)
     * @param res_j Second residue index (legacy 1-based)
     * @param best_j_for_i Best partner for res_i
     * @param best_i_for_j Best partner for res_j
     * @param is_mutual Whether they are mutual best
     * @param was_selected Whether pair was selected
     */
    void record_mutual_best_decision(int res_i, int res_j, int best_j_for_i, int best_i_for_j, bool is_mutual,
                                     bool was_selected);

    /**
     * @brief Record iteration state after each find_bestpair iteration
     * @param iteration_num Iteration number (1-based)
     * @param num_matched Number of matched residues
     * @param num_total Total number of residues
     * @param matched_indices Vector indicating which residues are matched (size = num_total + 1,
     * 1-based)
     * @param pairs Vector of pairs found so far (each pair is (res_i, res_j) in legacy 1-based
     * indices)
     */
    void record_iteration_state(int iteration_num, int num_matched, int num_total,
                                const std::vector<bool>& matched_indices,
                                const std::vector<std::pair<int, int>>& pairs);

    /**
     * @brief Record helix organization decisions from five2three algorithm
     * @param helix_num Helix number (1-based)
     * @param helix Helix segment information
     * @param pair_order Ordered pair indices
     * @param pairs Base pairs
     * @param strand_swapped Swap status for each pair
     */
    void record_helix_organization(size_t helix_num, const algorithms::HelixSegment& helix,
                                   const std::vector<size_t>& pair_order, const std::vector<core::BasePair>& pairs,
                                   const std::vector<bool>& strand_swapped);

    /**
     * @brief Record pair neighbor context (equivalent to legacy bp_order)
     * @param pairs Base pairs
     * @param context Vector of PairContextInfo structs
     */
    void record_bp_context(const std::vector<core::BasePair>& pairs,
                           const std::vector<algorithms::PairContextInfo>& context);

private:
    std::filesystem::path pdb_file_;
    std::string pdb_name_;
    nlohmann::json json_;

    // Split file storage: map from calculation type to array of records
    std::map<std::string, nlohmann::json> split_records_;

    // Index counters for tracking basepairs and hbonds
    // Legacy uses 1-based bp_idx, so start at 1
    size_t basepair_idx_counter_ = 1;
    size_t hbond_idx_counter_ = 0;

    // Track recorded base pairs to avoid duplicates (normalized by (min, max))
    std::set<std::pair<size_t, size_t>> recorded_base_pairs_;

    // Store base pairs in order for step calculation (matches legacy base_pair order)
    std::vector<core::BasePair> ordered_base_pairs_;

    /**
     * @brief Initialize JSON structure
     */
    void initialize_json();

    /**
     * @brief Add a calculation record to the JSON (also stores in split_records_)
     * @param record Record to add
     */
    void add_calculation_record(const nlohmann::json& record);

    /**
     * @brief Escape string for JSON
     * @param str String to escape
     * @return Escaped string
     */
    static std::string escape_string(const std::string& str);

    /**
     * @brief Format double for JSON (6 decimal places, or null if empty)
     * @param value Double value
     * @return JSON value (number or null)
     */
    static nlohmann::json format_double(double value);
};

} // namespace io
} // namespace x3dna
