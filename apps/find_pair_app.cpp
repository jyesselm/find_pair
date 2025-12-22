/**
 * @file find_pair_app.cpp
 * @brief Main executable for find_pair application
 */

#include <x3dna/apps/command_line_parser.hpp>
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/protocols/analyze_protocol.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/input_file_writer.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/config/config_manager.hpp>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <chrono>
#include <iomanip>

namespace {
    // Simple timer for profiling
    class Timer {
    public:
        void start() { start_ = std::chrono::high_resolution_clock::now(); }
        double elapsed_ms() const {
            auto end = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double, std::milli>(end - start_).count();
        }
    private:
        std::chrono::high_resolution_clock::time_point start_;
    };

    bool g_show_timing = false;

    void print_timing(const char* label, double ms) {
        if (g_show_timing) {
            std::cout << "[TIMING] " << std::setw(30) << std::left << label
                      << std::fixed << std::setprecision(1) << ms << " ms\n";
        }
    }
}

int main(int argc, char* argv[]) {
    try {
        Timer total_timer, step_timer;
        total_timer.start();

        // Parse command-line arguments
        auto options = x3dna::apps::CommandLineParser::parse_find_pair(argc, argv);

        // Check for timing flag (--timing or -t) and --no-json flag
        bool skip_json = false;
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--timing" || std::string(argv[i]) == "-t") {
                g_show_timing = true;
            }
            if (std::string(argv[i]) == "--no-json") {
                skip_json = true;
            }
        }

        // Set legacy mode in config manager if requested
        auto& config = x3dna::config::ConfigManager::instance();
        if (options.legacy_mode) {
            config.set_legacy_mode(true);
        }

        // Parse PDB file
        x3dna::io::PdbParser parser;
        if (options.hetatm) {
            parser.set_include_hetatm(true);
        }
        if (options.waters) {
            parser.set_include_waters(true);
        }

        step_timer.start();
        std::cout << "Parsing PDB file: " << options.pdb_file << "\n";
        auto structure = parser.parse_file(options.pdb_file);
        print_timing("PDB parsing", step_timer.elapsed_ms());

        // Create JSON writer for step-by-step debugging (if enabled)
        std::unique_ptr<x3dna::io::JsonWriter> json_writer;
        if (!skip_json) {
            json_writer = std::make_unique<x3dna::io::JsonWriter>(options.pdb_file);
        }

        // Create protocol
        x3dna::protocols::FindPairProtocol protocol;
        protocol.set_config_manager(config);
        protocol.set_single_strand_mode(options.single_strand);
        protocol.set_find_all_pairs(options.find_all_pairs);
        protocol.set_divide_helices(options.divide_helices);
        protocol.set_legacy_mode(options.legacy_mode);
        if (json_writer) {
            protocol.set_json_writer(json_writer.get());
        }
        // Execute protocol
        step_timer.start();
        std::cout << "Finding base pairs...\n";
        protocol.execute(structure);
        print_timing("Find pairs (total)", step_timer.elapsed_ms());

        // Get results
        const auto& base_pairs = protocol.base_pairs();
        std::cout << "Found " << base_pairs.size() << " base pairs\n";

        // Write JSON output (if JsonWriter was enabled)
        if (json_writer) {
            step_timer.start();
            std::filesystem::path json_output_dir = "data/json";
            json_writer->write_to_file(json_output_dir, true);
            std::cout << "JSON debug output written to " << json_output_dir << "\n";
            print_timing("JSON writing", step_timer.elapsed_ms());
        }

        // Write output file (.inp format)
        if (!base_pairs.empty()) {
            // Use the original PDB file path (as provided on command line)
            // Convert to string to preserve relative paths
            std::string pdb_file_str = options.pdb_file.string();
            x3dna::io::InputFileWriter::write(options.output_file, std::filesystem::path(pdb_file_str), base_pairs,
                                              2, // duplex_number
                                              1  // flags (explicit bp numbering)
            );
            std::cout << "Output file written: " << options.output_file << "\n";

            // Write ref_frames_modern.dat
            if (!options.legacy_inp_file.empty()) {
                // Use legacy pair ordering for exact frame matching
                auto legacy_ordering = x3dna::io::InputFileWriter::parse_legacy_inp_ordering(options.legacy_inp_file);
                if (!legacy_ordering.empty()) {
                    x3dna::io::InputFileWriter::write_ref_frames("ref_frames_modern.dat", base_pairs, structure,
                                                                 legacy_ordering);
                    std::cout << "Reference frames written: ref_frames_modern.dat "
                              << "(using legacy ordering from " << options.legacy_inp_file << ")\n";
                } else {
                    std::cerr << "[WARNING] Could not parse legacy inp file: " << options.legacy_inp_file << "\n";
                    x3dna::io::InputFileWriter::write_ref_frames("ref_frames_modern.dat", base_pairs, structure);
                    std::cout << "Reference frames written: ref_frames_modern.dat\n";
                }
            } else {
                x3dna::io::InputFileWriter::write_ref_frames("ref_frames_modern.dat", base_pairs, structure);
                std::cout << "Reference frames written: ref_frames_modern.dat\n";
            }

            // Calculate and write step/helical parameters
            if (base_pairs.size() >= 2) {
                step_timer.start();
                std::cout << "Calculating step and helical parameters...\n";

                // Create AnalyzeProtocol to calculate parameters
                x3dna::protocols::AnalyzeProtocol analyze_protocol;
                analyze_protocol.set_config_manager(config);
                analyze_protocol.set_legacy_mode(options.legacy_mode);

                // Execute on the .inp file we just wrote
                analyze_protocol.execute(options.output_file);
                print_timing("Analyze protocol", step_timer.elapsed_ms());

                // Get results
                const auto& step_params = analyze_protocol.step_parameters();
                const auto& helical_params = analyze_protocol.helical_parameters();
                const auto& analyze_base_pairs = analyze_protocol.base_pairs();

                std::cout << "Calculated " << step_params.size() << " step parameters\n";
                std::cout << "Calculated " << helical_params.size() << " helical parameters\n";

                // Write .par files
                if (!step_params.empty()) {
                    x3dna::io::InputFileWriter::write_step_params("bp_step.par", step_params, analyze_base_pairs,
                                                                  structure);
                    std::cout << "Step parameters written: bp_step.par\n";
                }

                if (!helical_params.empty()) {
                    x3dna::io::InputFileWriter::write_helical_params("bp_helical.par", helical_params,
                                                                     analyze_base_pairs, structure);
                    std::cout << "Helical parameters written: bp_helical.par\n";
                }
            }
        } else {
            std::cout << "No base pairs found - no output file written\n";
        }

        std::cout << "Done!\n";
        print_timing("TOTAL TIME", total_timer.elapsed_ms());

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred\n";
        return 1;
    }
}
