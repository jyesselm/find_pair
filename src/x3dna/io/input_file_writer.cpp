/**
 * @file input_file_writer.cpp
 * @brief InputFileWriter implementation
 */

#include <x3dna/io/input_file_writer.hpp>
#include <x3dna/core/nucleotide_utils.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace x3dna {
namespace io {

void InputFileWriter::write(const std::filesystem::path& output_path, const std::filesystem::path& pdb_file,
                            const std::vector<core::BasePair>& base_pairs, int duplex_number, int flags) {
    std::string output_file_name = default_output_filename(pdb_file);
    write(output_path, pdb_file, output_file_name, base_pairs, duplex_number, flags);
}

void InputFileWriter::write(const std::filesystem::path& output_path, const std::filesystem::path& pdb_file,
                            const std::string& output_file_name, const std::vector<core::BasePair>& base_pairs,
                            int duplex_number, int flags) {
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
    out << "    " << flags << " " << std::setw(5) << 0 << "    # explicit bp numbering/hetero atoms\n";

    // Base pair lines
    // Format: bp_num res1 res2 flag # comment
    // Note: Convert from 0-based to 1-based for residue indices
    for (size_t i = 0; i < base_pairs.size(); ++i) {
        const auto& pair = base_pairs[i];
        size_t bp_num = i + 1;                 // 1-based base pair number
        size_t res1 = pair.residue_idx1() + 1; // Convert to 1-based
        size_t res2 = pair.residue_idx2() + 1; // Convert to 1-based
        int flag = 0;                          // Default flag value

        // Get base pair type string if available
        std::string bp_type = pair.bp_type();
        std::string comment = "";
        if (!bp_type.empty()) {
            comment = " # " + bp_type;
        }

        out << std::setw(5) << bp_num << " " << std::setw(5) << res1 << " " << std::setw(5) << res2 << " "
            << std::setw(5) << flag << comment << "\n";
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
        size_t bp_num = i + 1; // 1-based

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
        out << "..." << std::setw(6) << bp_num << " " << formatted_bp_type << "   # " << res1_desc << " - " << res2_desc
            << "\n";

        // Get the reference frames
        const auto& frame1 = bp.frame1();
        const auto& frame2 = bp.frame2();

        // Calculate middle frame using cehs_average/bpstep_par algorithm
        //
        // Note: Without legacy ordering information, we cannot determine
        // which residue is strand 1 vs strand 2. Legacy uses strand 2 first,
        // strand 1 second (refs_right_left). For exact legacy matching,
        // use the version with legacy_pair_ordering parameter.
        //
        // Default assumption: residue_idx1 = strand 1, residue_idx2 = strand 2
        // Legacy uses strand 2 first, strand 1 second (refs_right_left)
        // So we need: frame2 (strand 2) first, frame1 (strand 1) second
        // This matches legacy refs_right_left() behavior
        auto mid_frame = calc.calculate_pair_frame(frame2, frame1);
        auto mid_org = mid_frame.origin();
        auto mid_x = mid_frame.x_axis();
        auto mid_y = mid_frame.y_axis();
        auto mid_z = mid_frame.z_axis();

        // Output origin
        out << std::setw(10) << mid_org.x() << std::setw(10) << mid_org.y() << std::setw(10) << mid_org.z()
            << "  # origin\n";

        // Output x-axis
        out << std::setw(10) << mid_x.x() << std::setw(10) << mid_x.y() << std::setw(10) << mid_x.z()
            << "  # x-axis\n";

        // Output y-axis
        out << std::setw(10) << mid_y.x() << std::setw(10) << mid_y.y() << std::setw(10) << mid_y.z()
            << "  # y-axis\n";

        // Output z-axis
        out << std::setw(10) << mid_z.x() << std::setw(10) << mid_z.y() << std::setw(10) << mid_z.z()
            << "  # z-axis\n";
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
    if (!residue.insertion().empty()) {
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
    oss << core::one_letter_code(residue);

    return oss.str();
}

std::string InputFileWriter::default_output_filename(const std::filesystem::path& pdb_file) {
    // Default output file name: remove extension and add .outp
    std::string stem = pdb_file.stem().string();
    return stem + ".outp";
}

std::map<std::pair<int, int>, int> InputFileWriter::parse_legacy_inp_ordering(const std::filesystem::path& inp_file) {

    std::map<std::pair<int, int>, int> ordering;

    std::ifstream in(inp_file);
    if (!in.is_open()) {
        return ordering; // Return empty map on error
    }

    std::string line;
    int line_num = 0;

    // Skip first 5 header lines
    while (std::getline(in, line) && line_num < 5) {
        ++line_num;
    }

    // Parse base pair lines
    // Legacy format: "  res1  res2  flag #    1 | ..." (no bp_num prefix)
    // Modern format: "  bp_num  res1  res2  flag # type" (has bp_num prefix)
    // In both: res1 = strand 1, res2 = strand 2
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        int first, second, third;

        // Try to read first three integers
        if (iss >> first >> second >> third) {
            int res1, res2;

            // Check if first number is small (likely a bp_num 1-1000)
            // If second and third are also small, it's probably legacy format
            // If first is small but second/third are larger (residue indices), it's modern format
            // Heuristic: if first < 100 and second > first, likely modern format with bp_num
            if (first < 100 && second > first && second > 10) {
                // Modern format: bp_num res1 res2
                res1 = second;
                res2 = third;
            } else {
                // Legacy format: res1 res2 flag
                res1 = first;
                res2 = second;
            }

            // Create canonical key (min, max)
            int min_res = std::min(res1, res2);
            int max_res = std::max(res1, res2);
            std::pair<int, int> key(min_res, max_res);

            // Store which residue was listed first (strand 1)
            // This tells us the ordering in legacy .inp file
            ordering[key] = res1;
        }
    }

    return ordering;
}

void InputFileWriter::write_ref_frames(const std::filesystem::path& output_path,
                                       const std::vector<core::BasePair>& base_pairs, const core::Structure& structure,
                                       const std::map<std::pair<int, int>, int>& legacy_pair_ordering) {
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
        size_t bp_num = i + 1; // 1-based

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
        out << "..." << std::setw(6) << bp_num << " " << formatted_bp_type << "   # " << res1_desc << " - " << res2_desc
            << "\n";

        // Get the reference frames
        const auto& frame1 = bp.frame1();
        const auto& frame2 = bp.frame2();

        // Check legacy ordering for this pair
        // Convert to 1-based indices for lookup (legacy uses 1-based)
        int res1_1based = bp.residue_idx1() + 1;
        int res2_1based = bp.residue_idx2() + 1;
        int min_res = std::min(res1_1based, res2_1based);
        int max_res = std::max(res1_1based, res2_1based);
        std::pair<int, int> key(min_res, max_res);

        // Determine frame order based on legacy strand assignment
        // Legacy: refs_right_left uses strand 2 first, strand 1 second
        // In legacy .inp: res1 = strand 1, res2 = strand 2
        // So we need: strand 2 frame first, strand 1 frame second
        core::ReferenceFrame f_strand2; // Will be passed first to calculate_pair_frame
        core::ReferenceFrame f_strand1; // Will be passed second to calculate_pair_frame

        auto it = legacy_pair_ordering.find(key);
        if (it != legacy_pair_ordering.end()) {
            // We have legacy ordering information
            // it->second is the residue that was listed first in legacy .inp (strand 1)
            int legacy_strand1_res = it->second;

            // Determine which modern residue corresponds to strand 1
            if (res1_1based == legacy_strand1_res) {
                // res1 is strand 1, res2 is strand 2
                f_strand1 = frame1;
                f_strand2 = frame2;
            } else if (res2_1based == legacy_strand1_res) {
                // res2 is strand 1, res1 is strand 2
                f_strand1 = frame2;
                f_strand2 = frame1;
            } else {
                // Can't determine strand assignment, use default order
                // (This shouldn't happen if legacy ordering is correct)
                f_strand2 = frame2;
                f_strand1 = frame1;
            }
        } else {
            // No legacy ordering available - use default assumption
            // In modern code: residue_idx1 is typically strand 1, residue_idx2 is strand 2
            // But legacy uses strand 2 first, strand 1 second
            // So we use frame2 (strand 2) first, frame1 (strand 1) second
            f_strand1 = frame1;
            f_strand2 = frame2;
        }

        // Calculate middle frame using cehs_average/bpstep_par algorithm
        // Legacy uses strand 2 first, strand 1 second (refs_right_left)
        // Always use this order to match legacy behavior
        auto mid_frame = calc.calculate_pair_frame(f_strand2, f_strand1);
        auto mid_org = mid_frame.origin();
        auto mid_x = mid_frame.x_axis();
        auto mid_y = mid_frame.y_axis();
        auto mid_z = mid_frame.z_axis();

        // Output origin
        out << std::setw(10) << mid_org.x() << std::setw(10) << mid_org.y() << std::setw(10) << mid_org.z()
            << "  # origin\n";

        // Output x-axis
        out << std::setw(10) << mid_x.x() << std::setw(10) << mid_x.y() << std::setw(10) << mid_x.z()
            << "  # x-axis\n";

        // Output y-axis
        out << std::setw(10) << mid_y.x() << std::setw(10) << mid_y.y() << std::setw(10) << mid_y.z()
            << "  # y-axis\n";

        // Output z-axis
        out << std::setw(10) << mid_z.x() << std::setw(10) << mid_z.y() << std::setw(10) << mid_z.z()
            << "  # z-axis\n";
    }

    out.close();
}

void InputFileWriter::write_step_params(const std::filesystem::path& output_path,
                                        const std::vector<core::BasePairStepParameters>& step_params,
                                        const std::vector<core::BasePair>& base_pairs,
                                        const core::Structure& structure) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path.string());
    }

    // Helper to get one-letter code for a residue using legacy index
    // Base pair indices are 0-based but correspond to legacy 1-based indices
    auto get_base_code = [&structure](size_t res_idx) -> char {
        // Convert 0-based index to 1-based legacy index
        int legacy_idx = static_cast<int>(res_idx + 1);
        auto residue = structure.get_residue_by_legacy_idx(legacy_idx);
        if (residue) {
            return core::one_letter_code(*residue);
        }
        return '-';
    };

    // Header line
    out << "#        Shift    Slide     Rise     Tilt     Roll    Twist\n";

    // Output each step
    for (size_t i = 0; i < step_params.size() && i + 1 < base_pairs.size(); ++i) {
        const auto& params = step_params[i];
        const auto& bp1 = base_pairs[i];
        const auto& bp2 = base_pairs[i + 1];

        // Get base types from residues
        char bp1_base1 = get_base_code(bp1.residue_idx1());
        char bp1_base2 = get_base_code(bp1.residue_idx2());
        char bp2_base1 = get_base_code(bp2.residue_idx1());
        char bp2_base2 = get_base_code(bp2.residue_idx2());

        // Format step label: e.g., "UA/UA" for step between U-A and U-A pairs
        std::string step_label;
        step_label = std::string(1, bp1_base1) + std::string(1, bp2_base1) + "/" + std::string(1, bp1_base2) +
                     std::string(1, bp2_base2);

        out << std::fixed << std::setprecision(2);
        out << std::setw(5) << step_label << " ";
        out << std::setw(8) << params.shift << " ";
        out << std::setw(8) << params.slide << " ";
        out << std::setw(8) << params.rise << " ";
        out << std::setw(8) << params.tilt << " ";
        out << std::setw(8) << params.roll << " ";
        out << std::setw(8) << params.twist << "\n";
    }

    out.close();
}

void InputFileWriter::write_helical_params(const std::filesystem::path& output_path,
                                           const std::vector<core::HelicalParameters>& helical_params,
                                           const std::vector<core::BasePair>& base_pairs,
                                           const core::Structure& structure) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path.string());
    }

    // Helper to get one-letter code for a residue using legacy index
    // Base pair indices are 0-based but correspond to legacy 1-based indices
    auto get_base_code = [&structure](size_t res_idx) -> char {
        // Convert 0-based index to 1-based legacy index
        int legacy_idx = static_cast<int>(res_idx + 1);
        auto residue = structure.get_residue_by_legacy_idx(legacy_idx);
        if (residue) {
            return core::one_letter_code(*residue);
        }
        return '-';
    };

    // Header line
    out << "#        X-disp   Y-disp   h-Rise    Incl.     Tip   h-Twist\n";

    // Output each step
    for (size_t i = 0; i < helical_params.size() && i + 1 < base_pairs.size(); ++i) {
        const auto& params = helical_params[i];
        const auto& bp1 = base_pairs[i];
        const auto& bp2 = base_pairs[i + 1];

        // Get base types from residues
        char bp1_base1 = get_base_code(bp1.residue_idx1());
        char bp1_base2 = get_base_code(bp1.residue_idx2());
        char bp2_base1 = get_base_code(bp2.residue_idx1());
        char bp2_base2 = get_base_code(bp2.residue_idx2());

        // Format step label: e.g., "UA/UA" for step between U-A and U-A pairs
        std::string step_label;
        step_label = std::string(1, bp1_base1) + std::string(1, bp2_base1) + "/" + std::string(1, bp1_base2) +
                     std::string(1, bp2_base2);

        out << std::fixed << std::setprecision(2);
        out << std::setw(5) << step_label << " ";
        out << std::setw(8) << params.x_displacement << " ";
        out << std::setw(8) << params.y_displacement << " ";
        out << std::setw(8) << params.rise << " ";
        out << std::setw(8) << params.inclination << " ";
        out << std::setw(8) << params.tip << " ";
        out << std::setw(8) << params.twist << "\n";
    }

    out.close();
}

} // namespace io
} // namespace x3dna
