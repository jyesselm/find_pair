/**
 * @file analyze_app.cpp
 * @brief Main executable for analyze application
 */

#include <x3dna/apps/command_line_parser.hpp>
#include <x3dna/protocols/analyze_protocol.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/io/input_file_parser.hpp>
#include <x3dna/io/json_writer.hpp>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
#include <exception>

int main(int argc, char* argv[]) {
    try {
        // Parse command-line arguments
        auto options = x3dna::apps::CommandLineParser::parse_analyze(argc, argv);

        // Set legacy mode in config manager if requested
        auto& config = x3dna::config::ConfigManager::instance();
        if (options.legacy_mode) {
            config.set_legacy_mode(true);
        }

        // Parse input file to get PDB file path for JSON writer
        auto input_data = x3dna::io::InputFileParser::parse(options.input_file);
        
        // Create JSON writer (optional - only if we want JSON output)
        // For now, we'll create it to enable step parameter recording
        x3dna::io::JsonWriter json_writer(input_data.pdb_file);
        
        // Create protocol
        x3dna::protocols::AnalyzeProtocol protocol;
        protocol.set_config_manager(config);
        protocol.set_calculate_torsions(options.calculate_torsions);
        protocol.set_simple_parameters(options.simple_pars);
        protocol.set_circular_structure(options.circular);
        protocol.set_step_start(options.step_start);
        protocol.set_step_size(options.step_size);
        protocol.set_legacy_mode(options.legacy_mode);
        protocol.set_json_writer(&json_writer);

        // Execute protocol
        std::cout << "Analyzing input file: " << options.input_file << "\n";
        protocol.execute(options.input_file);

        // Get results
        const auto& step_params = protocol.step_parameters();
        const auto& helical_params = protocol.helical_parameters();
        
        std::cout << "Calculated " << step_params.size() << " step parameters\n";
        std::cout << "Calculated " << helical_params.size() << " helical parameters\n";

        // Print detailed parameter comparison
        if (!step_params.empty()) {
            std::cout << "\n=== Step Parameters ===\n";
            std::cout << "#    Shift    Slide    Rise     Tilt     Roll     Twist\n";
            for (size_t i = 0; i < step_params.size(); ++i) {
                const auto& params = step_params[i];
                std::cout << std::fixed << std::setprecision(2);
                std::cout << std::setw(3) << (i + 1) << "  ";
                std::cout << std::setw(7) << params.shift << "  ";
                std::cout << std::setw(7) << params.slide << "  ";
                std::cout << std::setw(7) << params.rise << "  ";
                std::cout << std::setw(7) << params.tilt << "  ";
                std::cout << std::setw(7) << params.roll << "  ";
                std::cout << std::setw(7) << params.twist << "\n";
            }
        }

        if (!helical_params.empty()) {
            std::cout << "\n=== Helical Parameters ===\n";
            std::cout << "#  X-disp   Y-disp   h-Rise   Incl.    Tip      h-Twist\n";
            for (size_t i = 0; i < helical_params.size(); ++i) {
                const auto& params = helical_params[i];
                std::cout << std::fixed << std::setprecision(2);
                std::cout << std::setw(3) << (i + 1) << "  ";
                std::cout << std::setw(7) << params.x_displacement << "  ";
                std::cout << std::setw(7) << params.y_displacement << "  ";
                std::cout << std::setw(7) << params.rise << "  ";
                std::cout << std::setw(7) << params.inclination << "  ";
                std::cout << std::setw(7) << params.tip << "  ";
                std::cout << std::setw(7) << params.twist << "\n";
            }
        }

        // Write JSON files to data/json directory
        std::filesystem::path json_output_dir = "data/json";
        std::filesystem::create_directories(json_output_dir);
        json_writer.write_split_files(json_output_dir, true);
        
        // TODO: Write output files (.par files, .outp file)
        std::cout << "\nDone!\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred\n";
        return 1;
    }
}

