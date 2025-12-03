/**
 * @file command_line_parser.hpp
 * @brief Command-line argument parser for x3dna applications
 */

#pragma once

#include <string>
#include <filesystem>
#include <vector>

namespace x3dna {
namespace apps {

/**
 * @struct FindPairOptions
 * @brief Command-line options for find_pair application
 */
struct FindPairOptions {
    std::filesystem::path pdb_file;
    std::filesystem::path output_file;
    bool single_strand = false;                // -S or -1 flag
    bool find_all_pairs = false;               // -P flag
    bool divide_helices = false;               // -D flag
    bool curves = false;                       // -C flag
    bool curves_plus = false;                  // -c+ flag
    bool hetatm = false;                       // -T flag
    bool detailed = false;                     // -Z flag
    bool waters = false;                       // -W flag
    bool hjb = false;                          // -hjb flag
    std::string map_file;                      // -m flag
    bool legacy_mode = false;                  // --legacy-mode flag
    bool fix_indices_from_legacy_json = false; // --fix-indices flag
    std::string legacy_json_file = "";         // --fix-indices=FILE (auto-detected if empty)
    std::string legacy_inp_file = "";          // --legacy-inp=FILE for pair ordering

    /**
     * @brief Check if any option is set
     */
    bool has_options() const {
        return single_strand || find_all_pairs || divide_helices || curves || curves_plus ||
               hetatm || detailed || waters || hjb || !map_file.empty() || legacy_mode ||
               fix_indices_from_legacy_json;
    }
};

/**
 * @struct AnalyzeOptions
 * @brief Command-line options for analyze application
 */
struct AnalyzeOptions {
    std::filesystem::path input_file; // .inp file
    bool calculate_torsions = false;  // -t flag
    std::string torsion_file;         // -t=filename
    bool bz = true;                   // -bz flag (default: true)
    bool ring = false;                // -ri flag
    bool simple_pars = true;          // -si flag (default: true)
    bool abi = false;                 // -abi flag
    bool circular = false;            // -circ flag
    bool icnt = false;                // -C flag
    bool waters = false;              // -W flag
    size_t step_start = 1;            // -S=step,start
    size_t step_size = 1;             // -S=step,start
    bool legacy_mode = false;         // --legacy-mode flag
};

/**
 * @class CommandLineParser
 * @brief Parses command-line arguments for x3dna applications
 */
class CommandLineParser {
public:
    /**
     * @brief Parse find_pair command-line arguments
     * @param argc Argument count
     * @param argv Argument vector
     * @return FindPairOptions structure
     * @throws std::runtime_error on parsing error
     */
    static FindPairOptions parse_find_pair(int argc, char* argv[]);

    /**
     * @brief Parse analyze command-line arguments
     * @param argc Argument count
     * @param argv Argument vector
     * @return AnalyzeOptions structure
     * @throws std::runtime_error on parsing error
     */
    static AnalyzeOptions parse_analyze(int argc, char* argv[]);

    /**
     * @brief Print find_pair usage message
     */
    static void print_find_pair_usage(const char* program_name);

    /**
     * @brief Print analyze usage message
     */
    static void print_analyze_usage(const char* program_name);

private:
    /**
     * @brief Check if argument is a global option (like --help, --version)
     */
    static bool is_global_option(const std::string& arg);

    /**
     * @brief Check if argument is legacy mode flag
     */
    static bool is_legacy_mode(const std::string& arg);

    /**
     * @brief Extract value from option (e.g., "-m=value" -> "value")
     */
    static std::string extract_option_value(const std::string& arg);

    /**
     * @brief Check if option has value
     */
    static bool option_has_value(const std::string& arg);
};

} // namespace apps
} // namespace x3dna
