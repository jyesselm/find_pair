/**
 * @file json_writer.cpp
 * @brief Implementation of modern C++ JSON writer
 */

#include <x3dna/io/json_writer.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

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

    json_["calculations"] = nlohmann::json::array();
}

std::string JsonWriter::to_string(bool pretty_print) const {
    if (pretty_print) {
        return json_.dump(2); // 2 spaces indentation
    }
    return json_.dump();
}

void JsonWriter::write_to_file(const std::filesystem::path& output_path, bool pretty_print) const {
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + output_path.string());
    }

    if (pretty_print) {
        file << json_.dump(2);
    } else {
        file << json_.dump();
    }
}

void JsonWriter::add_calculation_record(const nlohmann::json& record) {
    json_["calculations"].push_back(record);
}

void JsonWriter::record_pdb_atoms(const core::Structure& structure) {
    nlohmann::json record;
    record["type"] = "pdb_atoms";

    // Count total atoms
    size_t total_atoms = structure.num_atoms();
    record["num_atoms"] = total_atoms;

    // Build atoms array
    nlohmann::json atoms_array = nlohmann::json::array();
    size_t atom_idx = 1; // 1-based indexing for legacy format

    for (const auto& chain : structure.chains()) {
        for (const auto& residue : chain.residues()) {
            for (const auto& atom : residue.atoms()) {
                nlohmann::json atom_json;
                atom_json["atom_idx"] = atom_idx++;
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
                }
                if (atom.atom_serial() > 0) {
                    atom_json["atom_serial"] = atom.atom_serial();
                }

                atoms_array.push_back(atom_json);
            }
        }
    }

    record["atoms"] = atoms_array;
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
    nlohmann::json record;
    record["type"] = "base_pair";
    record["residue_idx1"] = pair.residue_idx1();
    record["residue_idx2"] = pair.residue_idx2();
    record["bp_type"] = pair.bp_type();

    // Base pair type enum
    std::string type_str;
    switch (pair.type()) {
        case core::BasePairType::WATSON_CRICK:
            type_str = "WATSON_CRICK";
            break;
        case core::BasePairType::WOBBLE:
            type_str = "WOBBLE";
            break;
        case core::BasePairType::HOOGSTEEN:
            type_str = "HOOGSTEEN";
            break;
        default:
            type_str = "UNKNOWN";
            break;
    }
    record["pair_type"] = type_str;

    // Reference frames if available
    if (pair.frame1().has_value()) {
        record["frame1"] = pair.frame1()->to_json_legacy();
    }
    if (pair.frame2().has_value()) {
        record["frame2"] = pair.frame2()->to_json_legacy();
    }

    // Hydrogen bonds
    const auto& hbonds = pair.hydrogen_bonds();
    if (!hbonds.empty()) {
        nlohmann::json hbonds_array = nlohmann::json::array();
        for (const auto& hbond : hbonds) {
            nlohmann::json hbond_json;
            hbond_json["donor_atom"] = hbond.donor_atom;
            hbond_json["acceptor_atom"] = hbond.acceptor_atom;
            hbond_json["distance"] = format_double(hbond.distance);
            hbond_json["type"] = std::string(1, hbond.type);
            hbonds_array.push_back(hbond_json);
        }
        record["hydrogen_bonds"] = hbonds_array;
    }

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
                                        const std::array<double, 5>& rtn_val) {
    nlohmann::json record;
    record["type"] = "pair_validation";
    record["base_i"] = base_i;
    record["base_j"] = base_j;
    record["is_valid"] = is_valid;
    record["bp_type_id"] = bp_type_id;
    record["dir_x"] = format_double(dir_x);
    record["dir_y"] = format_double(dir_y);
    record["dir_z"] = format_double(dir_z);

    // Return values array
    nlohmann::json rtn_array = nlohmann::json::array();
    for (double val : rtn_val) {
        rtn_array.push_back(format_double(val));
    }
    record["rtn_val"] = rtn_array;

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
    for (const auto& hbond : hbonds) {
        nlohmann::json hbond_json;
        hbond_json["donor_atom"] = hbond.donor_atom;
        hbond_json["acceptor_atom"] = hbond.acceptor_atom;
        hbond_json["distance"] = format_double(hbond.distance);
        hbond_json["type"] = std::string(1, hbond.type);
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
    // Return null if value is effectively zero (matches legacy behavior)
    if (std::abs(value) < EMPTY_CRITERION) {
        return nlohmann::json(nullptr);
    }
    // Format with 6 decimal places (matches legacy format)
    return nlohmann::json(std::round(value * 1000000.0) / 1000000.0);
}

} // namespace io
} // namespace x3dna
