/**
 * @file json_writer.cpp
 * @brief Implementation of modern C++ JSON writer
 */

#include <x3dna/io/json_writer.hpp>
#include <x3dna/algorithms/base_pair_validator.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <iostream>

namespace x3dna {
namespace io {

// Constants
constexpr double EMPTY_CRITERION = 1e-10;

JsonWriter::JsonWriter(const std::filesystem::path& pdb_file) : pdb_file_(pdb_file) {
    initialize_json();
}

JsonWriter::~JsonWriter() {
    // Destructor - no explicit finalization needed
    // JSON is available via json() method at any time
}

void JsonWriter::initialize_json() {
    json_ = nlohmann::json::object();
    json_["pdb_file"] = pdb_file_.string();

    // Extract PDB name from path (stem without extension)
    pdb_name_ = pdb_file_.stem().string();
    json_["pdb_name"] = pdb_name_;

    // Use split files format (like org code)
    json_["calculations"] = nlohmann::json::array();
    nlohmann::json split_note;
    split_note["_note"] = "Calculations are split into separate files: " + pdb_name_ + "_*.json";
    split_note["_split_files"] = true;
    json_["calculations"].push_back(split_note);

    // Store metadata
    json_["metadata"] = nlohmann::json::object();
    json_["metadata"]["version"] = "Modern X3DNA C++ Implementation";
}

std::string JsonWriter::to_string(bool pretty_print) const {
    if (pretty_print) {
        return json_.dump(2); // 2 spaces indentation
    }
    return json_.dump();
}

void JsonWriter::write_to_file(const std::filesystem::path& output_path, bool pretty_print) const {
    // If output_path is a directory, write split files directly
    if (std::filesystem::is_directory(output_path) || output_path.extension().empty()) {
        write_split_files(output_path, pretty_print);
        return;
    }

    // If output_path is a file, treat parent as directory and write split files
    // (no main file - only split files in record-type directories)
    write_split_files(output_path.parent_path(), pretty_print);
}

void JsonWriter::add_calculation_record(const nlohmann::json& record) {
    json_["calculations"].push_back(record);

    // Also store in split_records_ for split file output
    if (record.contains("type") && record["type"].is_string()) {
        std::string calc_type = record["type"];
        // Do NOT remove type field from record for split file (comparison script expects it)
        // nlohmann::json split_record = record;
        // split_record.erase("type"); // Removed - keep type field for comparison

        if (split_records_.find(calc_type) == split_records_.end()) {
            split_records_[calc_type] = nlohmann::json::array();
        }
        split_records_[calc_type].push_back(record); // Push original record with "type"
    }
}

void JsonWriter::write_split_files(const std::filesystem::path& output_dir,
                                   bool pretty_print) const {
    if (split_records_.empty()) {
        return;
    }

    // Map record types to directory names
    std::map<std::string, std::string> type_to_dir = {
        {"pdb_atoms", "pdb_atoms"},
        {"base_frame_calc", "base_frame_calc"},
        {"frame_calc", "frame_calc"},
        {"ls_fitting", "ls_fitting"}, // ls_fitting has its own directory (matches legacy)
        {"base_pair", "base_pair"},
        {"pair_validation", "pair_validation"},
        {"distance_checks", "distance_checks"},
        {"hbond_list", "hbond_list"},
        {"find_bestpair_selection", "find_bestpair_selection"},
        {"bpstep_params", "bpstep_params"},
        {"helical_params", "helical_params"},
        {"best_partner_candidates", "best_partner_candidates"},
        {"mutual_best_decision", "mutual_best_decisions"},
        {"iteration_states", "iteration_states"},
    };

    for (const auto& [calc_type, records] : split_records_) {
        // Get directory name for this record type
        std::string dir_name = type_to_dir.count(calc_type) ? type_to_dir[calc_type] : calc_type;

        // Create record-type-specific directory
        std::filesystem::path record_dir = output_dir / dir_name;
        std::filesystem::create_directories(record_dir);

        // Write file: <PDB_ID>.json in the record-type directory
        std::filesystem::path split_file = record_dir / (pdb_name_ + ".json");
        std::ofstream file(split_file);
        if (!file.is_open()) {
            continue;
        }

        if (pretty_print) {
            file << records.dump(2);
        } else {
            file << records.dump();
        }
    }
}

void JsonWriter::load_pdb_lines() const {
    if (pdb_lines_loaded_) {
        return;
    }

    std::ifstream file(pdb_file_);
    if (!file.is_open()) {
        pdb_lines_loaded_ = true; // Mark as loaded even if failed
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        pdb_lines_.push_back(line);
    }

    pdb_lines_loaded_ = true;
}

std::string JsonWriter::get_pdb_line(size_t line_number) const {
    load_pdb_lines();

    // line_number is 1-based
    if (line_number > 0 && line_number <= pdb_lines_.size()) {
        return pdb_lines_[line_number - 1];
    }
    return "";
}

void JsonWriter::record_pdb_atoms(const core::Structure& structure) {
    nlohmann::json record;
    record["type"] = "pdb_atoms";

    // Count total atoms
    size_t total_atoms = structure.num_atoms();
    record["num_atoms"] = total_atoms;

    // Build atoms array
    nlohmann::json atoms_array = nlohmann::json::array();
    size_t sequential_idx = 1; // Fallback for atoms without legacy mapping

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                nlohmann::json atom_json;

                // Use legacy_atom_idx for atom_idx to match legacy exactly (1-based)
                int legacy_atom_idx = atom.legacy_atom_idx();

                if (legacy_atom_idx > 0) {
                    atom_json["atom_idx"] = legacy_atom_idx;
                } else {
                    // Fallback for atoms without legacy mapping
                    atom_json["atom_idx"] = sequential_idx++;
                }

                // Core fields (match legacy exactly)
                atom_json["atom_name"] = atom.name();
                atom_json["residue_name"] = atom.residue_name();
                atom_json["chain_id"] = std::string(1, atom.chain_id());
                atom_json["residue_seq"] = atom.residue_seq();

                // Coordinates
                const auto& pos = atom.position();
                atom_json["xyz"] = nlohmann::json::array({pos.x(), pos.y(), pos.z()});

                // Record type (A for ATOM, H for HETATM)
                atom_json["record_type"] = std::string(1, atom.record_type());

                // Optional metadata (only if present)
                if (atom.line_number() > 0) {
                    atom_json["line_number"] = atom.line_number();

                    // Add PDB line for debugging (like legacy code)
                    std::string pdb_line = get_pdb_line(atom.line_number());
                    if (!pdb_line.empty()) {
                        atom_json["pdb_line"] = pdb_line;
                    }
                }

                atoms_array.push_back(atom_json);
            }
        }
    }

    record["atoms"] = atoms_array;

    // add_calculation_record() will handle adding to split_records_
    // No need to manually add here (was causing duplicates)
    add_calculation_record(record);
}

void JsonWriter::record_residue_indices(const core::Structure& structure) {
    // Get residues in legacy order (PDB file order, grouped by ResName+ChainID+ResSeq+insertion)
    auto residues = structure.residues_in_legacy_order();

    if (residues.empty()) {
        return;
    }

    nlohmann::json record;
    record["type"] = "residue_indices";
    record["num_residue"] = residues.size();

    nlohmann::json seidx_array = nlohmann::json::array();

    // Build seidx array with 1-based residue indices
    for (size_t i = 0; i < residues.size(); i++) {
        const core::Residue* residue = residues[i];

        // Get atom range for this residue (start_atom, end_atom)
        auto [start_atom, end_atom] = residue->atom_range();

        // Skip residues with no valid atom indices
        if (start_atom == 0 && end_atom == 0) {
            continue;
        }

        nlohmann::json entry;
        entry["residue_idx"] = static_cast<int>(i + 1); // 1-based residue index
        entry["start_atom"] = start_atom;
        entry["end_atom"] = end_atom;
        seidx_array.push_back(entry);
    }

    record["seidx"] = seidx_array;
    record["num_residue"] = seidx_array.size(); // Update count after filtering

    add_calculation_record(record);
}

void JsonWriter::record_base_frame_calc(size_t residue_idx, char base_type,
                                        const std::filesystem::path& standard_template,
                                        double rms_fit,
                                        const std::vector<std::string>& matched_atoms,
                                        const std::string& residue_name, char chain_id,
                                        int residue_seq, char insertion) {
    nlohmann::json record;
    record["type"] = "base_frame_calc";
    record["residue_idx"] = residue_idx;

    // Note: residue_idx is already 1-based (legacy format) when passed from frame_json_recorder
    record["legacy_residue_idx"] = static_cast<int>(residue_idx);

    record["base_type"] = std::string(1, base_type);

    // Add residue identification information (matching legacy format)
    if (!residue_name.empty()) {
        record["residue_name"] = residue_name;
    }
    record["chain_id"] = std::string(1, chain_id);
    record["residue_seq"] = residue_seq;
    if (insertion != ' ') {
        record["insertion"] = std::string(1, insertion);
    }

    record["standard_template"] = standard_template.string();
    record["rms_fit"] = format_double(rms_fit);
    record["num_matched_atoms"] = matched_atoms.size();

    // Convert matched atoms to JSON array
    nlohmann::json atoms_array = nlohmann::json::array();
    for (const auto& atom_name : matched_atoms) {
        atoms_array.push_back(atom_name);
    }
    record["matched_atoms"] = atoms_array;

    add_calculation_record(record);
}

void JsonWriter::record_ls_fitting(size_t residue_idx, size_t num_points, double rms_fit,
                                   const geometry::Matrix3D& rotation_matrix,
                                   const geometry::Vector3D& translation,
                                   const std::string& residue_name, char chain_id, int residue_seq,
                                   char insertion) {
    nlohmann::json record;
    record["type"] = "ls_fitting";
    record["residue_idx"] = residue_idx;

    // Note: residue_idx is already 1-based (legacy format) when passed from frame_json_recorder
    record["legacy_residue_idx"] = static_cast<int>(residue_idx);

    // Add residue identification information (matching legacy format)
    if (!residue_name.empty()) {
        record["residue_name"] = residue_name;
    }
    record["chain_id"] = std::string(1, chain_id);
    record["residue_seq"] = residue_seq;
    if (insertion != ' ') {
        record["insertion"] = std::string(1, insertion);
    }

    record["num_points"] = num_points;
    record["rms_fit"] = format_double(rms_fit);

    // Rotation matrix as 3x3 nested array
    nlohmann::json rot_array = nlohmann::json::array();
    for (int i = 0; i < 3; ++i) {
        nlohmann::json row = nlohmann::json::array();
        for (int j = 0; j < 3; ++j) {
            row.push_back(format_double(rotation_matrix.at(i, j)));
        }
        rot_array.push_back(row);
    }
    record["rotation_matrix"] = rot_array;

    // Translation as 3-element array
    record["translation"] =
        nlohmann::json::array({format_double(translation.x()), format_double(translation.y()),
                               format_double(translation.z())});

    add_calculation_record(record);
}

void JsonWriter::record_frame_calc(size_t residue_idx, char base_type,
                                   const std::filesystem::path& template_file, double rms_fit,
                                   const std::vector<geometry::Vector3D>& matched_std_xyz,
                                   const std::vector<geometry::Vector3D>& matched_exp_xyz,
                                   const std::string& residue_name, char chain_id, int residue_seq,
                                   char insertion) {
    if (matched_std_xyz.size() != matched_exp_xyz.size()) {
        throw std::invalid_argument("Matched coordinate arrays must have same size");
    }

    nlohmann::json record;
    record["type"] = "frame_calc";
    record["residue_idx"] = residue_idx;

    // Note: residue_idx is already 1-based (legacy format) when passed from frame_json_recorder
    record["legacy_residue_idx"] = static_cast<int>(residue_idx);
    record["base_type"] = std::string(1, base_type);

    // Add residue identification information (matching legacy format)
    if (!residue_name.empty()) {
        record["residue_name"] = residue_name;
    }
    record["chain_id"] = std::string(1, chain_id);
    record["residue_seq"] = residue_seq;
    if (insertion != ' ') {
        record["insertion"] = std::string(1, insertion);
    }

    record["template_file"] = template_file.string();
    record["rms_fit"] = format_double(rms_fit);
    record["num_matched_atoms"] = matched_std_xyz.size();

    // Matched coordinates array
    nlohmann::json coords_array = nlohmann::json::array();
    for (size_t i = 0; i < matched_std_xyz.size(); ++i) {
        nlohmann::json coord_json;
        coord_json["atom_idx"] = static_cast<int>(i + 1); // 1-based

        // Standard coordinates
        const auto& std_xyz = matched_std_xyz[i];
        coord_json["std_xyz"] = nlohmann::json::array(
            {format_double(std_xyz.x()), format_double(std_xyz.y()), format_double(std_xyz.z())});

        // Experimental coordinates
        const auto& exp_xyz = matched_exp_xyz[i];
        coord_json["exp_xyz"] = nlohmann::json::array(
            {format_double(exp_xyz.x()), format_double(exp_xyz.y()), format_double(exp_xyz.z())});

        coords_array.push_back(coord_json);
    }
    record["matched_coordinates"] = coords_array;

    add_calculation_record(record);
}

void JsonWriter::record_base_pair(const core::BasePair& pair) {
    // Convert residue indices to legacy format (1-based, counting all residues)
    // BasePair stores 0-based indices, but legacy uses 1-based
    size_t base_i = pair.residue_idx1() + 1; // Convert to 1-based
    size_t base_j = pair.residue_idx2() + 1; // Convert to 1-based

    // Normalize pair key to avoid duplicates (always use (min, max) order)
    std::pair<size_t, size_t> pair_key =
        (base_i < base_j) ? std::make_pair(base_i, base_j) : std::make_pair(base_j, base_i);

    // Check if this pair has already been recorded
    if (recorded_base_pairs_.find(pair_key) != recorded_base_pairs_.end()) {
        // Skip duplicate - already recorded this pair
        return;
    }

    // Mark as recorded
    recorded_base_pairs_.insert(pair_key);

    // Create a mutable copy to assign index
    core::BasePair pair_with_idx = pair;
    pair_with_idx.set_basepair_idx(basepair_idx_counter_++);

    // Assign hbond indices
    std::vector<core::hydrogen_bond> hbonds = pair_with_idx.hydrogen_bonds();
    for (size_t i = 0; i < hbonds.size(); ++i) {
        hbonds[i].hbond_idx = hbond_idx_counter_++;
    }
    pair_with_idx.set_hydrogen_bonds(hbonds);

    // Use BasePair's to_json_legacy() to ensure exact format match with legacy JSON
    // This ensures base_i, base_j, orien_i, orien_j, org_i, org_j, dir_xyz format
    nlohmann::json record = pair_with_idx.to_json_legacy();

    // Note: to_json_legacy() already sets "type" = "base_pair", "base_i", "base_j", "bp_type",
    // "orien_i", "orien_j", "org_i", "org_j", "dir_xyz", "basepair_idx", and "hbonds" if present

    // Set the indices (already converted above)
    record["base_i"] = static_cast<long>(base_i);
    record["base_j"] = static_cast<long>(base_j);

    add_calculation_record(record);
}

void JsonWriter::record_bpstep_params(size_t bp_idx1, size_t bp_idx2,
                                      const core::BasePairStepParameters& params) {
    nlohmann::json record;
    record["type"] = "bpstep_params";
    record["bp_idx1"] = bp_idx1;
    record["bp_idx2"] = bp_idx2;

    // 6 parameters: Shift, Slide, Rise, Tilt, Roll, Twist
    record["shift"] = format_double(params.shift);
    record["slide"] = format_double(params.slide);
    record["rise"] = format_double(params.rise);
    record["tilt"] = format_double(params.tilt);
    record["roll"] = format_double(params.roll);
    record["twist"] = format_double(params.twist);

    // Midstep frame if available
    if (params.midstep_frame.has_value()) {
        record["midstep_frame"] = params.midstep_frame->to_json_legacy();
    }

    add_calculation_record(record);
}

void JsonWriter::record_helical_params(size_t bp_idx1, size_t bp_idx2,
                                       const core::HelicalParameters& params) {
    nlohmann::json record;
    record["type"] = "helical_params";
    record["bp_idx1"] = bp_idx1;
    record["bp_idx2"] = bp_idx2;

    // 6 parameters: x_displacement, y_displacement, rise, inclination, tip, twist
    record["x_displacement"] = format_double(params.x_displacement);
    record["y_displacement"] = format_double(params.y_displacement);
    record["rise"] = format_double(params.rise);
    record["inclination"] = format_double(params.inclination);
    record["tip"] = format_double(params.tip);
    record["twist"] = format_double(params.twist);

    // Helical midstep frame if available
    if (params.midstep_frame.has_value()) {
        record["midstep_frame"] = params.midstep_frame->to_json_legacy();
    }

    add_calculation_record(record);
}

void JsonWriter::record_all_ref_frames(const core::Structure& structure) {
    nlohmann::json record;
    record["type"] = "all_ref_frames";

    nlohmann::json frames_array = nlohmann::json::array();
    size_t residue_idx = 0;

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            if (residue.reference_frame().has_value()) {
                nlohmann::json frame_json;
                frame_json["residue_idx"] = residue_idx;
                frame_json["frame"] = residue.reference_frame()->to_json_legacy();
                frames_array.push_back(frame_json);
            }
            residue_idx++;
        }
    }

    record["frames"] = frames_array;
    record["num_frames"] = frames_array.size();

    add_calculation_record(record);
}

void JsonWriter::record_removed_atom(const std::string& pdb_line, const std::string& reason,
                                     int atom_serial, const std::string& atom_name,
                                     const std::string& residue_name, char chain_id,
                                     int residue_seq, const geometry::Vector3D* xyz,
                                     int model_num) {
    nlohmann::json record;
    record["type"] = "removed_atom";

    if (!reason.empty()) {
        record["reason"] = reason;
    }
    if (!pdb_line.empty()) {
        record["pdb_line"] = pdb_line;
    }
    if (atom_serial > 0) {
        record["atom_serial"] = atom_serial;
    }
    if (!atom_name.empty()) {
        record["atom_name"] = atom_name;
    }
    if (!residue_name.empty()) {
        record["residue_name"] = residue_name;
    }
    if (chain_id != ' ') {
        record["chain_id"] = std::string(1, chain_id);
    }
    if (residue_seq > 0) {
        record["residue_seq"] = residue_seq;
    }
    if (xyz != nullptr) {
        record["xyz"] = nlohmann::json::array(
            {format_double(xyz->x()), format_double(xyz->y()), format_double(xyz->z())});
    }
    record["model_num"] = model_num;

    add_calculation_record(record);
}

void JsonWriter::record_removed_atoms_summary(size_t num_removed) {
    nlohmann::json record;
    record["type"] = "removed_atoms_summary";
    record["num_removed"] = num_removed;

    add_calculation_record(record);
}

void JsonWriter::record_pair_validation(size_t base_i, size_t base_j, bool is_valid, int bp_type_id,
                                        double dir_x, double dir_y, double dir_z,
                                        const std::array<double, 5>& rtn_val,
                                        const algorithms::ValidationParameters& params) {
    // NOTE: We receive 0-based indices, but need to output 1-based for legacy compatibility
    // Legacy pair_validation records use 1-based indices (e.g., base_i=1 to 20 for 20 residues)
    if ((base_i + 1 == 93 && base_j + 1 == 130) || (base_i + 1 == 130 && base_j + 1 == 93)) {
        std::cerr << "[DEBUG JSON writer] pair_validation(" << base_i+1 << "," << base_j+1 
                  << ") bp_type_id=" << bp_type_id << " qs=" << rtn_val[4] << std::endl;
    }
    nlohmann::json record;
    record["type"] = "pair_validation";
    record["base_i"] = static_cast<long>(base_i + 1); // Convert to 1-based for legacy
    record["base_j"] = static_cast<long>(base_j + 1); // Convert to 1-based for legacy
    record["is_valid"] = static_cast<long>(is_valid ? 1 : 0); // Legacy uses long
    record["bp_type_id"] = static_cast<long>(bp_type_id);

    // Direction vectors (nested object, matches legacy format)
    nlohmann::json dir_vectors;
    dir_vectors["dir_x"] = format_double(dir_x);
    dir_vectors["dir_y"] = format_double(dir_y);
    dir_vectors["dir_z"] = format_double(dir_z);
    record["direction_vectors"] = dir_vectors;

    // Calculated values (nested object, matches legacy format)
    // rtn_val[0] = dorg, [1] = d_v, [2] = plane_angle, [3] = dNN, [4] = quality_score
    nlohmann::json calc_values;
    calc_values["dorg"] = format_double(rtn_val[0]);
    calc_values["d_v"] = format_double(rtn_val[1]);
    calc_values["plane_angle"] = format_double(rtn_val[2]);
    calc_values["dNN"] = format_double(rtn_val[3]);
    calc_values["quality_score"] = format_double(rtn_val[4]);
    record["calculated_values"] = calc_values;

    // Validation checks (nested object, matches legacy format)
    nlohmann::json validation_checks;
    validation_checks["distance_check"] =
        (rtn_val[0] >= params.min_dorg && rtn_val[0] <= params.max_dorg);
    validation_checks["d_v_check"] = (rtn_val[1] >= params.min_dv && rtn_val[1] <= params.max_dv);
    validation_checks["plane_angle_check"] =
        (rtn_val[2] >= params.min_plane_angle && rtn_val[2] <= params.max_plane_angle);
    validation_checks["dNN_check"] = (rtn_val[3] >= params.min_dNN && rtn_val[3] <= params.max_dNN);
    record["validation_checks"] = validation_checks;

    // Thresholds (nested object, matches legacy format)
    nlohmann::json thresholds;
    thresholds["min_dorg"] = format_double(params.min_dorg);
    thresholds["max_dorg"] = format_double(params.max_dorg);
    thresholds["min_dv"] = format_double(params.min_dv);
    thresholds["max_dv"] = format_double(params.max_dv);
    thresholds["min_plane_angle"] = format_double(params.min_plane_angle);
    thresholds["max_plane_angle"] = format_double(params.max_plane_angle);
    thresholds["min_dNN"] = format_double(params.min_dNN);
    thresholds["max_dNN"] = format_double(params.max_dNN);
    record["thresholds"] = thresholds;

    add_calculation_record(record);
}

void JsonWriter::record_distance_checks(size_t base_i, size_t base_j, double dorg, double dNN,
                                        double plane_angle, double d_v, double overlap_area) {
    // NOTE: We receive 0-based indices, but need to output 1-based for legacy compatibility
    // Also, legacy outputs BOTH (i,j) and (j,i) - so we need to output both pairs

    // Output first pair (i, j) with 1-based indices
    nlohmann::json record;
    record["type"] = "distance_checks";
    record["base_i"] = static_cast<long>(base_i + 1); // Convert to 1-based
    record["base_j"] = static_cast<long>(base_j + 1); // Convert to 1-based

    // Values (nested object, matches legacy format)
    nlohmann::json values;
    values["dorg"] = format_double(dorg);
    values["dNN"] = format_double(dNN);
    values["plane_angle"] = format_double(plane_angle);
    values["d_v"] = format_double(d_v);
    // Legacy outputs 0.0 instead of null for zero overlap
    if (overlap_area == 0.0 || std::isnan(overlap_area)) {
        values["overlap_area"] = 0.0;
    } else {
        values["overlap_area"] = format_double(overlap_area);
    }
    record["values"] = values;

    add_calculation_record(record);

    // Legacy also outputs the reversed pair (j, i)
    nlohmann::json record_reversed;
    record_reversed["type"] = "distance_checks";
    record_reversed["base_i"] = static_cast<long>(base_j + 1); // Swap: j becomes i
    record_reversed["base_j"] = static_cast<long>(base_i + 1); // Swap: i becomes j
    record_reversed["values"] = values;                        // Same values

    add_calculation_record(record_reversed);
}

void JsonWriter::record_hbond_list(size_t base_i, size_t base_j,
                                   const std::vector<core::hydrogen_bond>& hbonds) {
    // NOTE: We receive 0-based indices, but need to output 1-based for legacy compatibility
    nlohmann::json record;
    record["type"] = "hbond_list";
    record["base_i"] = static_cast<long>(base_i + 1); // Convert to 1-based for legacy
    record["base_j"] = static_cast<long>(base_j + 1); // Convert to 1-based for legacy
    record["num_hbonds"] = hbonds.size();

    nlohmann::json hbonds_array = nlohmann::json::array();
    // Legacy records ALL H-bonds including invalid ones (type=' ')
    // Legacy uses 1-based per-pair indexing: for (i = 1; i <= num_hbonds; i++) hbond_idx = i
    for (size_t i = 0; i < hbonds.size(); ++i) {
        const auto& hbond = hbonds[i];
        nlohmann::json hbond_json;
        hbond_json["donor_atom"] = hbond.donor_atom;
        hbond_json["acceptor_atom"] = hbond.acceptor_atom;
        hbond_json["distance"] = format_double(hbond.distance);
        hbond_json["type"] = std::string(1, hbond.type);
        // Use existing hbond_idx if set (for legacy comparison), otherwise use 1-based per-pair
        // indexing Legacy: records ALL H-bonds (including type=' ') with sequential hbond_idx
        // starting at 1
        if (hbond.hbond_idx.has_value()) {
            hbond_json["hbond_idx"] = static_cast<long>(hbond.hbond_idx.value());
        } else {
            hbond_json["hbond_idx"] =
                static_cast<long>(i + 1); // 1-based per-pair indexing (matches legacy)
        }
        hbonds_array.push_back(hbond_json);
    }
    record["hbonds"] = hbonds_array;

    add_calculation_record(record);
}

std::string JsonWriter::escape_string(const std::string& str) {
    std::ostringstream escaped;
    for (char c : str) {
        if (c == '"' || c == '\\') {
            escaped << '\\' << c;
        } else if (c == '\n') {
            escaped << "\\n";
        } else {
            escaped << c;
        }
    }
    return escaped.str();
}

nlohmann::json JsonWriter::format_double(double value) {
    // Check for NaN or Inf first
    if (std::isnan(value) || std::isinf(value)) {
        return nlohmann::json(nullptr);
    }
    // Return null if value is effectively zero (matches legacy behavior)
    if (std::abs(value) < EMPTY_CRITERION) {
        return nlohmann::json(nullptr);
    }
    // Format with 6 decimal places (matches legacy format)
    return nlohmann::json(std::round(value * 1000000.0) / 1000000.0);
}

void JsonWriter::record_find_bestpair_selection(
    const std::vector<std::pair<size_t, size_t>>& selected_pairs) {
    nlohmann::json record;
    record["type"] = "find_bestpair_selection";
    record["num_bp"] = selected_pairs.size();

    nlohmann::json pairs_array = nlohmann::json::array();
    for (const auto& pair : selected_pairs) {
        pairs_array.push_back(nlohmann::json::array({pair.first, pair.second}));
    }
    record["pairs"] = pairs_array;

    add_calculation_record(record);
}

void JsonWriter::record_best_partner_candidates(
    int res_i, const std::vector<std::tuple<int, bool, double, int>>& candidates, int best_j,
    double best_score) {
    nlohmann::json record;
    record["type"] = "best_partner_candidates";
    record["res_i"] = res_i;
    record["num_candidates"] = candidates.size();
    record["best_partner"] = best_j;
    record["best_score"] = best_score;

    nlohmann::json candidates_array = nlohmann::json::array();
    for (const auto& cand : candidates) {
        int res_j = std::get<0>(cand);
        bool is_eligible = std::get<1>(cand);
        double score = std::get<2>(cand);
        int bp_type_id = std::get<3>(cand);

        nlohmann::json cand_json;
        cand_json["res_j"] = res_j;
        cand_json["is_eligible"] = is_eligible ? 1 : 0;
        cand_json["score"] = score;
        cand_json["bp_type_id"] = bp_type_id;
        cand_json["is_best"] = (res_j == best_j) ? 1 : 0;

        candidates_array.push_back(cand_json);
    }
    record["candidates"] = candidates_array;

    add_calculation_record(record);
}

void JsonWriter::record_mutual_best_decision(int res_i, int res_j, int best_j_for_i,
                                             int best_i_for_j, bool is_mutual, bool was_selected) {
    nlohmann::json record;
    record["type"] = "mutual_best_decision";
    record["res1"] = res_i;
    record["res2"] = res_j;
    record["best_partner_for_res1"] = best_j_for_i;
    record["best_partner_for_res2"] = best_i_for_j;
    record["is_mutual"] = is_mutual ? 1 : 0;
    record["was_selected"] = was_selected ? 1 : 0;

    add_calculation_record(record);
}

void JsonWriter::record_iteration_state(int iteration_num, int num_matched, int num_total,
                                        const std::vector<bool>& matched_indices,
                                        const std::vector<std::pair<int, int>>& pairs) {
    nlohmann::json record;
    record["type"] = "iteration_states";
    record["iteration_num"] = iteration_num;
    record["num_matched"] = num_matched;
    record["num_total"] = num_total;

    // Collect pairs found in this iteration
    nlohmann::json pairs_array = nlohmann::json::array();
    for (const auto& pair : pairs) {
        pairs_array.push_back(nlohmann::json::array({pair.first, pair.second}));
    }
    record["pairs_found_in_iteration"] = pairs_array;

    // Collect matched residues
    nlohmann::json matched_array = nlohmann::json::array();
    for (size_t i = 1; i < matched_indices.size() && i <= static_cast<size_t>(num_total); ++i) {
        if (matched_indices[i]) {
            matched_array.push_back(static_cast<int>(i));
        }
    }
    record["matched_residues"] = matched_array;

    add_calculation_record(record);
}

} // namespace io
} // namespace x3dna
