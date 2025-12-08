/**
 * @file find_pair_app.cpp
 * @brief Main executable for find_pair application
 */

#include <x3dna/apps/command_line_parser.hpp>
#include <x3dna/protocols/find_pair_protocol.hpp>
#include <x3dna/io/pdb_parser.hpp>
#include <x3dna/io/input_file_writer.hpp>
#include <x3dna/io/json_writer.hpp>
#include <x3dna/config/config_manager.hpp>
#include <iostream>
#include <stdexcept>
#include <exception>

int main(int argc, char* argv[]) {
    try {
        // Parse command-line arguments
        auto options = x3dna::apps::CommandLineParser::parse_find_pair(argc, argv);

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

        std::cout << "Parsing PDB file: " << options.pdb_file << "\n";
        auto structure = parser.parse_file(options.pdb_file);

        // Create JSON writer for step-by-step debugging (if enabled)
        std::unique_ptr<x3dna::io::JsonWriter> json_writer;
        // Always enable JSON output for debugging
        json_writer = std::make_unique<x3dna::io::JsonWriter>(options.pdb_file);

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
        std::cout << "Finding base pairs...\n";
        protocol.execute(structure);

        // Get results
        const auto& base_pairs = protocol.base_pairs();
        std::cout << "Found " << base_pairs.size() << " base pairs\n";

        // Write JSON output (if JsonWriter was enabled)
        if (json_writer) {
            std::filesystem::path json_output_dir = "data/json";
            json_writer->write_to_file(json_output_dir, true);
            std::cout << "JSON debug output written to " << json_output_dir << "\n";
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
        } else {
            std::cout << "No base pairs found - no output file written\n";
        }

        std::cout << "Done!\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred\n";
        return 1;
    }
}
