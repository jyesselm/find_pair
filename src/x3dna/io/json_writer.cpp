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

JsonWriter::JsonWriter(const std::filesystem::path& pdb_file,
                       const std::filesystem::path& legacy_json_file) 
    : pdb_file_(pdb_file), legacy_mappings_loaded_(false) {
    initialize_json();
    
    // Load legacy mappings if legacy JSON file is provided
    if (!legacy_json_file.empty() && std::filesystem::exists(legacy_json_file)) {
        load_legacy_mappings(legacy_json_file);
    }
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
    // Write main metadata file
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + output_path.string());
    }

    if (pretty_print) {
        file << json_.dump(2);
    } else {
        file << json_.dump();
    }
    
    // Write split files for each calculation type
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

void JsonWriter::write_split_files(const std::filesystem::path& output_dir, bool pretty_print) const {
    if (split_records_.empty()) {
        std::cerr << "[JSON_WRITER] Warning: No split records to write (split_records_ is empty)\n";
        return;
    }
    
    std::cerr << "[JSON_WRITER] Writing " << split_records_.size() << " split files to " << output_dir << "\n";
    
    for (const auto& [calc_type, records] : split_records_) {
        std::filesystem::path split_file = output_dir / (pdb_name_ + "_" + calc_type + ".json");
        std::ofstream file(split_file);
        if (!file.is_open()) {
            std::cerr << "[JSON_WRITER] Warning: Could not open split file for writing: " << split_file << "\n";
            continue;
        }
        
        if (pretty_print) {
            file << records.dump(2);
        } else {
            file << records.dump();
        }
        
        std::cerr << "[JSON_WRITER] Wrote split file: " << split_file << " (" << records.size() << " entries)\n";
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

void JsonWriter::load_legacy_mappings(const std::filesystem::path& legacy_json_file) {
    if (legacy_mappings_loaded_) {
        return;
    }
    
    try {
        std::ifstream file(legacy_json_file);
        if (!file.is_open()) {
            std::cerr << "[JSON_WRITER] Warning: Could not open legacy JSON file: " << legacy_json_file << "\n";
            legacy_mappings_loaded_ = true;
            return;
        }
        
        nlohmann::json legacy_data;
        file >> legacy_data;
        
        // Extract atoms from legacy JSON (handle both split files and grouped format)
        nlohmann::json legacy_atoms;
        nlohmann::json calculations = legacy_data.value("calculations", nlohmann::json::object());
        
        if (calculations.is_object()) {
            // Grouped format: calculations is a dict
            auto pdb_atoms_group = calculations.value("pdb_atoms", nlohmann::json::array());
            if (pdb_atoms_group.is_array() && !pdb_atoms_group.empty()) {
                legacy_atoms = pdb_atoms_group[0].value("atoms", nlohmann::json::array());
            }
        } else if (calculations.is_array()) {
            // Array format: find pdb_atoms entry
            for (const auto& calc : calculations) {
                if (calc.value("type", std::string("")) == "pdb_atoms") {
                    legacy_atoms = calc.value("atoms", nlohmann::json::array());
                    break;
                }
            }
        }
        
        // Try split file if not found in main file
        if (legacy_atoms.empty()) {
            std::filesystem::path split_file = legacy_json_file.parent_path() / 
                (legacy_json_file.stem().string() + "_pdb_atoms.json");
            if (std::filesystem::exists(split_file)) {
                std::ifstream split_stream(split_file);
                if (split_stream.is_open()) {
                    nlohmann::json split_data;
                    split_stream >> split_data;
                    if (split_data.is_array() && !split_data.empty()) {
                        legacy_atoms = split_data[0].value("atoms", nlohmann::json::array());
                    }
                }
            }
        }
        
        // Build atom index mapping
        for (const auto& atom : legacy_atoms) {
            if (!atom.is_object()) continue;
            
            std::string chain_id_str = atom.value("chain_id", std::string(""));
            char chain_id = chain_id_str.empty() ? ' ' : chain_id_str[0];
            int residue_seq = atom.value("residue_seq", 0);
            std::string insertion_str = atom.value("insertion", std::string(" "));
            char insertion = insertion_str.empty() ? ' ' : insertion_str[0];
            std::string atom_name = atom.value("atom_name", std::string(""));
            int atom_idx = atom.value("atom_idx", 0);
            
            if (atom_idx > 0 && !atom_name.empty()) {
                legacy_atom_idx_map_[std::make_tuple(chain_id, residue_seq, insertion, atom_name)] = atom_idx;
            }
        }
        
        // Build residue index mapping from atoms
        // Legacy residue indices are assigned sequentially as residues are encountered
        // We infer residue boundaries by detecting changes in (chain_id, residue_seq, insertion)
        int residue_idx = 1; // 1-based
        char last_chain_id = '\0';
        int last_residue_seq = 0;
        char last_insertion = '\0';
        
        for (const auto& atom : legacy_atoms) {
            if (!atom.is_object()) continue;
            
            std::string chain_id_str = atom.value("chain_id", std::string(""));
            char chain_id = chain_id_str.empty() ? ' ' : chain_id_str[0];
            int residue_seq = atom.value("residue_seq", 0);
            std::string insertion_str = atom.value("insertion", std::string(" "));
            char insertion = insertion_str.empty() ? ' ' : insertion_str[0];
            
            // Check if this is a new residue (different from previous)
            if (residue_seq != last_residue_seq || chain_id != last_chain_id || insertion != last_insertion) {
                // New residue - add to mapping
                auto key = std::make_tuple(chain_id, residue_seq, insertion);
                if (legacy_residue_idx_map_.find(key) == legacy_residue_idx_map_.end()) {
                    legacy_residue_idx_map_[key] = residue_idx;
                    residue_idx++;
                }
                
                last_chain_id = chain_id;
                last_residue_seq = residue_seq;
                last_insertion = insertion;
            }
        }
        
        std::cerr << "[JSON_WRITER] Loaded " << legacy_atom_idx_map_.size() 
                  << " legacy atom indices and " << legacy_residue_idx_map_.size() 
                  << " legacy residue indices\n";
        
        legacy_mappings_loaded_ = true;
    } catch (const std::exception& e) {
        std::cerr << "[JSON_WRITER] Error loading legacy mappings: " << e.what() << "\n";
        legacy_mappings_loaded_ = true;
    }
}

int JsonWriter::get_legacy_atom_idx(char chain_id, int residue_seq, char insertion, const std::string& atom_name) const {
    auto key = std::make_tuple(chain_id, residue_seq, insertion, atom_name);
    auto it = legacy_atom_idx_map_.find(key);
    if (it != legacy_atom_idx_map_.end()) {
        return it->second;
    }
    return 0;
}

int JsonWriter::get_legacy_residue_idx(char chain_id, int residue_seq, char insertion) const {
    auto key = std::make_tuple(chain_id, residue_seq, insertion);
    auto it = legacy_residue_idx_map_.find(key);
    if (it != legacy_residue_idx_map_.end()) {
        return it->second;
    }
    return 0;
}

void JsonWriter::set_legacy_indices_on_structure(core::Structure& structure) {
    // Set legacy indices on all atoms in the structure
    // This makes legacy indices available to all algorithms that use the structure
    structure.set_legacy_indices(legacy_atom_idx_map_, legacy_residue_idx_map_);
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
                
                // Use the atom's stored legacy_atom_idx (set when structure was created)
                // This ensures modern atom_idx matches legacy atom_idx for direct comparison
                int legacy_atom_idx = atom.legacy_atom_idx();
                
                if (legacy_atom_idx > 0) {
                    atom_json["atom_idx"] = legacy_atom_idx;
                    atom_json["legacy_atom_idx"] = legacy_atom_idx;
                } else {
                    // For atoms not in legacy, use sequential indexing
                    atom_json["atom_idx"] = sequential_idx++;
                    // No legacy_atom_idx for atoms not in legacy
                }
                atom_json["atom_name"] = atom.name();
                atom_json["residue_name"] = atom.residue_name();
                atom_json["chain_id"] = std::string(1, atom.chain_id());
                atom_json["residue_seq"] = atom.residue_seq();

                // Coordinates
                const auto& pos = atom.position();
                atom_json["xyz"] = nlohmann::json::array({pos.x(), pos.y(), pos.z()});

                // Record type (A for ATOM, H for HETATM)
                atom_json["record_type"] = std::string(1, atom.record_type());

                // Additional metadata
                if (atom.alt_loc() != ' ' && atom.alt_loc() != '\0') {
                    atom_json["alt_loc"] = std::string(1, atom.alt_loc());
                }
                if (atom.insertion() != ' ' && atom.insertion() != '\0') {
                    atom_json["insertion"] = std::string(1, atom.insertion());
                }
                if (atom.line_number() > 0) {
                    atom_json["line_number"] = atom.line_number();
                    
                    // Add PDB line for debugging (like org code)
                    std::string pdb_line = get_pdb_line(atom.line_number());
                    if (!pdb_line.empty()) {
                        atom_json["pdb_line"] = pdb_line;
                    }
                }
                if (atom.atom_serial() > 0) {
                    atom_json["atom_serial"] = atom.atom_serial();
                }
                
                // Use the atom's stored legacy_residue_idx (set when structure was created)
                int legacy_residue_idx = atom.legacy_residue_idx();
                if (legacy_residue_idx > 0) {
                    atom_json["legacy_residue_idx"] = legacy_residue_idx;
                }

                atoms_array.push_back(atom_json);
            }
        }
    }

    record["atoms"] = atoms_array;
    
    // Store in split_records_ for split file output
    nlohmann::json split_record = record;
    split_record.erase("type");
    split_records_["pdb_atoms"] = nlohmann::json::array();
    split_records_["pdb_atoms"].push_back(split_record);
    
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
    
    // Add legacy residue index if available
    int legacy_residue_idx = get_legacy_residue_idx(chain_id, residue_seq, insertion);
    if (legacy_residue_idx > 0) {
        record["legacy_residue_idx"] = legacy_residue_idx;
    }
    
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
    
    // Add legacy residue index if available
    int legacy_residue_idx = get_legacy_residue_idx(chain_id, residue_seq, insertion);
    if (legacy_residue_idx > 0) {
        record["legacy_residue_idx"] = legacy_residue_idx;
    }

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
    
    // Add legacy residue index if available
    int legacy_residue_idx = get_legacy_residue_idx(chain_id, residue_seq, insertion);
    if (legacy_residue_idx > 0) {
        record["legacy_residue_idx"] = legacy_residue_idx;
    }
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
    std::pair<size_t, size_t> pair_key = (base_i < base_j) ? 
        std::make_pair(base_i, base_j) : std::make_pair(base_j, base_i);
    
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
    nlohmann::json record;
    record["type"] = "pair_validation";
    record["base_i"] = static_cast<long>(base_i);  // Legacy uses long
    record["base_j"] = static_cast<long>(base_j);
    record["is_valid"] = static_cast<long>(is_valid ? 1 : 0);  // Legacy uses long
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
    validation_checks["distance_check"] = (rtn_val[0] >= params.min_dorg && rtn_val[0] <= params.max_dorg);
    validation_checks["d_v_check"] = (rtn_val[1] >= params.min_dv && rtn_val[1] <= params.max_dv);
    validation_checks["plane_angle_check"] = (rtn_val[2] >= params.min_plane_angle && rtn_val[2] <= params.max_plane_angle);
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

void JsonWriter::record_distance_checks(size_t base_i, size_t base_j,
                                        double dorg, double dNN, double plane_angle, double d_v,
                                        double overlap_area) {
    nlohmann::json record;
    record["type"] = "distance_checks";
    record["base_i"] = static_cast<long>(base_i);  // Legacy uses long
    record["base_j"] = static_cast<long>(base_j);

    // Values (nested object, matches legacy format)
    nlohmann::json values;
    values["dorg"] = format_double(dorg);
    values["dNN"] = format_double(dNN);
    values["plane_angle"] = format_double(plane_angle);
    values["d_v"] = format_double(d_v);
    values["overlap_area"] = format_double(overlap_area);
    record["values"] = values;

    add_calculation_record(record);
}

void JsonWriter::record_hbond_list(size_t base_i, size_t base_j,
                                   const std::vector<core::hydrogen_bond>& hbonds) {
    nlohmann::json record;
    record["type"] = "hbond_list";
    record["base_i"] = base_i;
    record["base_j"] = base_j;
    record["num_hbonds"] = hbonds.size();

    nlohmann::json hbonds_array = nlohmann::json::array();
    for (size_t i = 0; i < hbonds.size(); ++i) {
        const auto& hbond = hbonds[i];
        nlohmann::json hbond_json;
        hbond_json["donor_atom"] = hbond.donor_atom;
        hbond_json["acceptor_atom"] = hbond.acceptor_atom;
        hbond_json["distance"] = format_double(hbond.distance);
        hbond_json["type"] = std::string(1, hbond.type);
        // Assign hbond_idx if not already set, otherwise use existing or position
        if (hbond.hbond_idx.has_value()) {
            hbond_json["hbond_idx"] = static_cast<long>(hbond.hbond_idx.value());
        } else {
            hbond_json["hbond_idx"] = static_cast<long>(hbond_idx_counter_++);
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

void JsonWriter::record_find_bestpair_selection(const std::vector<std::pair<size_t, size_t>>& selected_pairs) {
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

} // namespace io
} // namespace x3dna
