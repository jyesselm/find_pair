/**
 * @file command_line_parser.cpp
 * @brief CommandLineParser implementation
 */

#include <x3dna/apps/command_line_parser.hpp>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cstring>

namespace x3dna {
namespace apps {

FindPairOptions CommandLineParser::parse_find_pair(int argc, char* argv[]) {
    FindPairOptions options;

    if (argc < 2) {
        print_find_pair_usage(argv[0]);
        throw std::runtime_error("Insufficient arguments");
    }

    int arg_idx = 1;
    std::vector<std::string> positional_args;

    // Parse arguments
    while (arg_idx < argc) {
        std::string arg = argv[arg_idx];

        // Check for legacy mode
        if (is_legacy_mode(arg)) {
            options.legacy_mode = true;
            arg_idx++;
            continue;
        }

        // Check for global options
        if (is_global_option(arg)) {
            arg_idx++;
            continue;
        }

        // Check if it's an option (starts with -)
        if (arg[0] == '-') {
            // Handle multi-character flags like -SDC
            for (size_t i = 1; i < arg.length(); ++i) {
                char flag = arg[i];
                switch (flag) {
                    case 'S':
                    case '1':
                        options.single_strand = true;
                        break;
                    case 'P':
                        options.find_all_pairs = true;
                        break;
                    case 'D':
                        options.divide_helices = true;
                        break;
                    case 'C':
                        options.curves = true;
                        break;
                    case 'T':
                        options.hetatm = true;
                        break;
                    case 'Z':
                        options.detailed = true;
                        break;
                    case 'W':
                        options.waters = true;
                        break;
                    default:
                        // Unknown flag, skip
                        break;
                }
            }

            // Handle special options
            if (arg == "-c+" || arg == "--c+") {
                options.curves_plus = true;
            } else if (arg == "-hjb" || arg == "--hjb") {
                options.hjb = true;
            } else if (arg.find("-m") == 0) {
                // Map file option: -m or -m=filename
                if (option_has_value(arg)) {
                    options.map_file = extract_option_value(arg);
                } else {
                    options.map_file = "Gaussian";
                }
            }

            arg_idx++;
        } else {
            // Positional argument (file path)
            positional_args.push_back(arg);
            arg_idx++;
        }
    }

    // Parse positional arguments
    if (positional_args.empty()) {
        print_find_pair_usage(argv[0]);
        throw std::runtime_error("PDB file not specified");
    }

    options.pdb_file = positional_args[0];

    // Optional output file
    if (positional_args.size() > 1) {
        options.output_file = positional_args[1];
    } else {
        // Default output file: remove extension and add .inp
        std::string stem = options.pdb_file.stem().string();
        options.output_file = options.pdb_file.parent_path() / (stem + ".inp");
    }

    return options;
}

AnalyzeOptions CommandLineParser::parse_analyze(int argc, char* argv[]) {
    AnalyzeOptions options;

    if (argc < 2) {
        print_analyze_usage(argv[0]);
        throw std::runtime_error("Insufficient arguments");
    }

    int arg_idx = 1;
    std::vector<std::string> positional_args;

    // Parse arguments
    while (arg_idx < argc) {
        std::string arg = argv[arg_idx];

        // Check for legacy mode
        if (is_legacy_mode(arg)) {
            options.legacy_mode = true;
            arg_idx++;
            continue;
        }

        // Check for global options
        if (is_global_option(arg)) {
            arg_idx++;
            continue;
        }

        // Check if it's an option
        if (arg[0] == '-') {
            if (arg == "-bz" || arg == "--bz") {
                options.bz = true;
            } else if (arg == "-no-bz" || arg == "--no-bz") {
                options.bz = false;
            } else if (arg == "-ri" || arg == "--ri") {
                options.ring = true;
            } else if (arg == "-si" || arg == "--si") {
                options.simple_pars = true;
            } else if (arg == "-no-si" || arg == "--no-si") {
                options.simple_pars = false;
            } else if (arg == "-abi" || arg == "--abi") {
                options.abi = true;
            } else if (arg == "-circ" || arg == "--circ") {
                options.circular = true;
            } else if (arg.find("-t") == 0) {
                // Torsion option: -t or -t=filename
                options.calculate_torsions = true;
                if (option_has_value(arg)) {
                    options.torsion_file = extract_option_value(arg);
                }
            } else if (arg.find("-S=") == 0) {
                // Step option: -S=step,start
                std::string value = extract_option_value(arg);
                size_t comma_pos = value.find(',');
                if (comma_pos != std::string::npos) {
                    options.step_size = std::stoul(value.substr(0, comma_pos));
                    options.step_start = std::stoul(value.substr(comma_pos + 1));
                } else {
                    options.step_size = std::stoul(value);
                    options.step_start = 1;
                }
            } else {
                // Handle single-character flags
                for (size_t i = 1; i < arg.length(); ++i) {
                    char flag = arg[i];
                    switch (flag) {
                        case 'C':
                            options.icnt = true;
                            break;
                        case 'W':
                            options.waters = true;
                            break;
                        default:
                            break;
                    }
                }
            }

            arg_idx++;
        } else {
            // Positional argument (input file)
            positional_args.push_back(arg);
            arg_idx++;
        }
    }

    // Parse positional arguments
    if (positional_args.empty()) {
        print_analyze_usage(argv[0]);
        throw std::runtime_error("Input file (.inp) not specified");
    }

    options.input_file = positional_args[0];

    return options;
}

void CommandLineParser::print_find_pair_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " [options] <pdb_file> [output_file]\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -S, -1          Single strand mode\n";
    std::cerr << "  -P               Find all pairs\n";
    std::cerr << "  -D               Divide helices\n";
    std::cerr << "  -C               Curves output\n";
    std::cerr << "  -c+              Curves+ output\n";
    std::cerr << "  -T               Include HETATM records\n";
    std::cerr << "  -Z               Detailed output\n";
    std::cerr << "  -W               Include waters\n";
    std::cerr << "  -hjb             HJB option\n";
    std::cerr << "  -m[=filename]    Map file (default: Gaussian)\n";
    std::cerr << "  --legacy-mode    Enable legacy compatibility mode\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << program_name << " 1H4S.pdb\n";
    std::cerr << "  " << program_name << " --legacy-mode 1H4S.pdb output.inp\n";
}

void CommandLineParser::print_analyze_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " [options] <input_file.inp>\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -t[=filename]    Calculate torsions\n";
    std::cerr << "  -bz, --bz        BZ option (default: on)\n";
    std::cerr << "  -ri, --ri        Ring option\n";
    std::cerr << "  -si, --si        Simple parameters (default: on)\n";
    std::cerr << "  -abi, --abi      ABI option\n";
    std::cerr << "  -circ, --circ    Circular structure\n";
    std::cerr << "  -C               ICNT option\n";
    std::cerr << "  -W               Include waters\n";
    std::cerr << "  -S=step,start    Step parameters\n";
    std::cerr << "  --legacy-mode    Enable legacy compatibility mode\n";
    std::cerr << "\nExample:\n";
    std::cerr << "  " << program_name << " input.inp\n";
    std::cerr << "  " << program_name << " --legacy-mode -S=1,1 input.inp\n";
}

bool CommandLineParser::is_global_option(const std::string& arg) {
    return arg == "--help" || arg == "-h" || arg == "--version" || arg == "-v";
}

bool CommandLineParser::is_legacy_mode(const std::string& arg) {
    return arg == "--legacy-mode" || arg == "--legacy";
}

std::string CommandLineParser::extract_option_value(const std::string& arg) {
    size_t eq_pos = arg.find('=');
    if (eq_pos != std::string::npos) {
        return arg.substr(eq_pos + 1);
    }
    return "";
}

bool CommandLineParser::option_has_value(const std::string& arg) {
    return arg.find('=') != std::string::npos;
}

} // namespace apps
} // namespace x3dna

