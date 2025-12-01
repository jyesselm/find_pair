/**
 * @file input_file_writer.cpp
 * @brief InputFileWriter implementation
 */

#include <x3dna/io/input_file_writer.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace x3dna {
namespace io {

void InputFileWriter::write(const std::filesystem::path& output_path,
                             const std::filesystem::path& pdb_file,
                             const std::vector<core::BasePair>& base_pairs,
                             int duplex_number,
                             int flags) {
    std::string output_file_name = default_output_filename(pdb_file);
    write(output_path, pdb_file, output_file_name, base_pairs, duplex_number, flags);
}

void InputFileWriter::write(const std::filesystem::path& output_path,
                             const std::filesystem::path& pdb_file,
                             const std::string& output_file_name,
                             const std::vector<core::BasePair>& base_pairs,
                             int duplex_number,
                             int flags) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path.string());
    }

    // Line 1: PDB file path
    // Write the path as a string - use the original path string as provided
    // This preserves relative paths and matches legacy behavior
    out << pdb_file.string() << "\n";

    // Line 2: Output file name
    out << output_file_name << "\n";

    // Line 3: Duplex number
    out << "    " << duplex_number << "         # duplex\n";

    // Line 4: Number of base pairs
    out << std::setw(5) << base_pairs.size() << "         # number of base-pairs\n";

    // Line 5: Flags
    out << "    " << flags << " " << std::setw(5) << 0
        << "    # explicit bp numbering/hetero atoms\n";

    // Base pair lines
    // Format: bp_num res1 res2 flag # comment
    // Note: Convert from 0-based to 1-based for residue indices
    for (size_t i = 0; i < base_pairs.size(); ++i) {
        const auto& pair = base_pairs[i];
        size_t bp_num = i + 1;  // 1-based base pair number
        size_t res1 = pair.residue_idx1() + 1;  // Convert to 1-based
        size_t res2 = pair.residue_idx2() + 1;  // Convert to 1-based
        int flag = 0;  // Default flag value

        // Get base pair type string if available
        std::string bp_type = pair.bp_type();
        std::string comment = "";
        if (!bp_type.empty()) {
            comment = " # " + bp_type;
        }

        out << std::setw(5) << bp_num << " " << std::setw(5) << res1 << " "
            << std::setw(5) << res2 << " " << std::setw(5) << flag << comment << "\n";
    }

    out.close();
}

void InputFileWriter::write_ref_frames(const std::filesystem::path& output_path,
                                        const std::vector<core::BasePair>& base_pairs,
                                        const core::Structure& structure) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path.string());
    }

    // Set precision for floating point output
    out << std::fixed << std::setprecision(4);

    // Line 1: Number of base pairs
    out << std::setw(5) << base_pairs.size() << " base-pairs\n";

    // Get residues for descriptions
    auto residues = structure.all_residues();
    
    // Create parameter calculator to compute proper midstep frames
    algorithms::ParameterCalculator calc;

    // For each base pair, output the mid-frame using bpstep_par algorithm
    for (size_t i = 0; i < base_pairs.size(); ++i) {
        const auto& bp = base_pairs[i];
        size_t bp_num = i + 1;  // 1-based

        // Get bp_type (e.g., "GC" -> "G-C")
        std::string bp_type = bp.bp_type();
        std::string formatted_bp_type;
        if (bp_type.size() >= 2) {
            formatted_bp_type = std::string(1, bp_type[0]) + "-" + std::string(1, bp_type[1]);
        } else {
            formatted_bp_type = bp_type;
        }

        // Get residue descriptions
        std::string res1_desc = "unknown";
        std::string res2_desc = "unknown";
        if (bp.residue_idx1() < residues.size()) {
            res1_desc = format_residue_description(*residues[bp.residue_idx1()]);
        }
        if (bp.residue_idx2() < residues.size()) {
            res2_desc = format_residue_description(*residues[bp.residue_idx2()]);
        }

        // Line: ...     N bp_type   # res1_desc - res2_desc
        out << "..." << std::setw(6) << bp_num << " " << formatted_bp_type
            << "   # " << res1_desc << " - " << res2_desc << "\n";

        // Get the reference frames
        auto frame1 = bp.frame1();
        auto frame2 = bp.frame2();

        if (frame1.has_value() && frame2.has_value()) {
            // Calculate middle frame using cehs_average/bpstep_par algorithm
            // 
            // Note: The y/z axis orientation depends on residue ordering within
            // the pair. Legacy uses variable ordering based on find_bestpair scan
            // order, while modern consistently uses smaller-residue-first.
            // This causes y/z axes to differ by a 180° rotation around x-axis
            // for ~35% of pairs. The origin and x-axis always match exactly.
            auto mid_frame = calc.calculate_pair_frame(frame1.value(), frame2.value());
            auto mid_org = mid_frame.origin();
            auto mid_x = mid_frame.x_axis();
            auto mid_y = mid_frame.y_axis();
            auto mid_z = mid_frame.z_axis();

            // Output origin
            out << std::setw(10) << mid_org.x() << std::setw(10) << mid_org.y()
                << std::setw(10) << mid_org.z() << "  # origin\n";

            // Output x-axis
            out << std::setw(10) << mid_x.x() << std::setw(10) << mid_x.y()
                << std::setw(10) << mid_x.z() << "  # x-axis\n";

            // Output y-axis
            out << std::setw(10) << mid_y.x() << std::setw(10) << mid_y.y()
                << std::setw(10) << mid_y.z() << "  # y-axis\n";

            // Output z-axis
            out << std::setw(10) << mid_z.x() << std::setw(10) << mid_z.y()
                << std::setw(10) << mid_z.z() << "  # z-axis\n";
        } else {
            // No frames available - output identity frame
            out << std::setw(10) << 0.0 << std::setw(10) << 0.0
                << std::setw(10) << 0.0 << "  # origin\n";
            out << std::setw(10) << 1.0 << std::setw(10) << 0.0
                << std::setw(10) << 0.0 << "  # x-axis\n";
            out << std::setw(10) << 0.0 << std::setw(10) << 1.0
                << std::setw(10) << 0.0 << "  # y-axis\n";
            out << std::setw(10) << 0.0 << std::setw(10) << 0.0
                << std::setw(10) << 1.0 << "  # z-axis\n";
        }
    }

    out.close();
}

std::string InputFileWriter::format_residue_description(const core::Residue& residue) {
    // Format: chain:seq:[name]letter
    // Example: T:...4_:[..G]G
    std::ostringstream oss;
    oss << residue.chain_id() << ":";

    // Format sequence number (right-aligned in 4 chars, with dots for padding)
    std::ostringstream seq_ss;
    seq_ss << residue.seq_num();
    std::string seq_str = seq_ss.str();
    if (residue.insertion() != ' ' && residue.insertion() != '\0') {
        seq_str += residue.insertion();
    } else {
        seq_str += "_";
    }

    // Pad with dots to make it 4 chars + insertion
    while (seq_str.size() < 5) {
        seq_str = "." + seq_str;
    }
    oss << seq_str << ":";

    // Format residue name (right-padded with dots in brackets)
    std::string name = residue.name();
    while (name.size() < 3) {
        name = "." + name;
    }
    oss << "[" << name << "]";

    // Add one-letter code
    oss << residue.one_letter_code();

    return oss.str();
}

std::string InputFileWriter::default_output_filename(const std::filesystem::path& pdb_file) {
    // Default output file name: remove extension and add .outp
    std::string stem = pdb_file.stem().string();
    return stem + ".outp";
}

std::map<std::pair<int,int>, int> InputFileWriter::parse_legacy_inp_ordering(
    const std::filesystem::path& inp_file) {
    
    std::map<std::pair<int,int>, int> ordering;
    
    std::ifstream in(inp_file);
    if (!in.is_open()) {
        return ordering;  // Return empty map on error
    }
    
    std::string line;
    int line_num = 0;
    
    // Skip first 5 header lines
    while (std::getline(in, line) && line_num < 5) {
        ++line_num;
    }
    
    // Parse base pair lines: "  res1  res2  flag #..."
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        int res1, res2;
        
        // Try to parse as .inp format (without bp number prefix first)
        // Legacy format: "  947  1013   0 #    1 | ..."
        if (iss >> res1 >> res2) {
            // Create canonical key (min, max)
            int min_res = std::min(res1, res2);
            int max_res = std::max(res1, res2);
            std::pair<int,int> key(min_res, max_res);
            
            // Store which residue was listed first
            ordering[key] = res1;
        }
    }
    
    return ordering;
}

void InputFileWriter::write_ref_frames(const std::filesystem::path& output_path,
                                        const std::vector<core::BasePair>& base_pairs,
                                        const core::Structure& structure,
                                        const std::map<std::pair<int,int>, int>& legacy_pair_ordering) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path.string());
    }

    // Set precision for floating point output
    out << std::fixed << std::setprecision(4);

    // Line 1: Number of base pairs
    out << std::setw(5) << base_pairs.size() << " base-pairs\n";

    // Get residues for descriptions
    auto residues = structure.all_residues();
    
    // Create parameter calculator to compute proper midstep frames
    algorithms::ParameterCalculator calc;

    // For each base pair, output the mid-frame using bpstep_par algorithm
    for (size_t i = 0; i < base_pairs.size(); ++i) {
        const auto& bp = base_pairs[i];
        size_t bp_num = i + 1;  // 1-based

        // Get bp_type (e.g., "GC" -> "G-C")
        std::string bp_type = bp.bp_type();
        std::string formatted_bp_type;
        if (bp_type.size() >= 2) {
            formatted_bp_type = std::string(1, bp_type[0]) + "-" + std::string(1, bp_type[1]);
        } else {
            formatted_bp_type = bp_type;
        }

        // Get residue descriptions
        std::string res1_desc = "unknown";
        std::string res2_desc = "unknown";
        if (bp.residue_idx1() < residues.size()) {
            res1_desc = format_residue_description(*residues[bp.residue_idx1()]);
        }
        if (bp.residue_idx2() < residues.size()) {
            res2_desc = format_residue_description(*residues[bp.residue_idx2()]);
        }

        // Line: ...     N bp_type   # res1_desc - res2_desc
        out << "..." << std::setw(6) << bp_num << " " << formatted_bp_type
            << "   # " << res1_desc << " - " << res2_desc << "\n";

        // Get the reference frames
        auto frame1 = bp.frame1();
        auto frame2 = bp.frame2();

        if (frame1.has_value() && frame2.has_value()) {
            // Check legacy ordering for this pair
            // Convert to 1-based indices for lookup (legacy uses 1-based)
            int res1_1based = bp.residue_idx1() + 1;
            int res2_1based = bp.residue_idx2() + 1;
            int min_res = std::min(res1_1based, res2_1based);
            int max_res = std::max(res1_1based, res2_1based);
            std::pair<int,int> key(min_res, max_res);
            
            // Determine frame order based on legacy ordering
            core::ReferenceFrame f1 = frame1.value();
            core::ReferenceFrame f2 = frame2.value();
            
            // If legacy had larger residue first, swap frames for calculation
            auto it = legacy_pair_ordering.find(key);
            if (it != legacy_pair_ordering.end() && it->second == max_res) {
                // Legacy had larger-first: swap frame order
                std::swap(f1, f2);
            }
            
            // Calculate middle frame using cehs_average/bpstep_par algorithm
            auto mid_frame = calc.calculate_pair_frame(f1, f2);
            auto mid_org = mid_frame.origin();
            auto mid_x = mid_frame.x_axis();
            auto mid_y = mid_frame.y_axis();
            auto mid_z = mid_frame.z_axis();

            // Output origin
            out << std::setw(10) << mid_org.x() << std::setw(10) << mid_org.y()
                << std::setw(10) << mid_org.z() << "  # origin\n";

            // Output x-axis
            out << std::setw(10) << mid_x.x() << std::setw(10) << mid_x.y()
                << std::setw(10) << mid_x.z() << "  # x-axis\n";

            // Output y-axis
            out << std::setw(10) << mid_y.x() << std::setw(10) << mid_y.y()
                << std::setw(10) << mid_y.z() << "  # y-axis\n";

            // Output z-axis
            out << std::setw(10) << mid_z.x() << std::setw(10) << mid_z.y()
                << std::setw(10) << mid_z.z() << "  # z-axis\n";
        } else {
            // No frames available - output identity frame
            out << std::setw(10) << 0.0 << std::setw(10) << 0.0
                << std::setw(10) << 0.0 << "  # origin\n";
            out << std::setw(10) << 1.0 << std::setw(10) << 0.0
                << std::setw(10) << 0.0 << "  # x-axis\n";
            out << std::setw(10) << 0.0 << std::setw(10) << 1.0
                << std::setw(10) << 0.0 << "  # y-axis\n";
            out << std::setw(10) << 0.0 << std::setw(10) << 0.0
                << std::setw(10) << 1.0 << "  # z-axis\n";
        }
    }

    out.close();
}

} // namespace io
} // namespace x3dna

