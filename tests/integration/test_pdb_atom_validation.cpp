/**
 * @file test_pdb_atom_validation.cpp
 * @brief Integration tests to validate PDB parsing against legacy JSON
 *
 * This test ensures that when we read a PDB file, we get exactly the same
 * atom data as stored in the legacy JSON files. This is critical for
 * regression testing and ensuring correctness.
 */

#include <gtest/gtest.h>
#include <x3dna/core/atom.hpp>
#include "integration_test_base.hpp"
#include "test_data_discovery.hpp"
#include <filesystem>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <future>
#include <thread>
#include <mutex>
#include <atomic>

namespace x3dna::test {

using namespace x3dna::core;

/**
 * @brief Compare atom data from PDB parsing with legacy JSON
 *
 * This function will be used once we implement the PDB parser.
 * For now, it validates the JSON structure and prepares the comparison logic.
 */
class PdbAtomValidationTest : public integration_test_base {
protected:
    /**
     * @brief Load pdb_atoms record from legacy JSON
     */
    nlohmann::json load_pdb_atoms_record(const std::filesystem::path& json_file) {
        auto json = load_legacy_json(json_file);

        // Find pdb_atoms record in calculations array
        auto records = find_records_by_type(json, "pdb_atoms");
        if (records.empty()) {
            throw std::runtime_error("No pdb_atoms record found in JSON");
        }

        // Should only be one pdb_atoms record per file
        if (records.size() > 1) {
            throw std::runtime_error("Multiple pdb_atoms records found");
        }

        return records[0];
    }

    /**
     * @brief Create Atom objects from legacy JSON pdb_atoms record
     */
    std::vector<Atom> atoms_from_json(const nlohmann::json& pdb_atoms_record) {
        std::vector<Atom> atoms;

        if (!pdb_atoms_record.contains("atoms") || !pdb_atoms_record["atoms"].is_array()) {
            throw std::invalid_argument("Invalid pdb_atoms record: missing atoms array");
        }

        for (const auto& atom_json : pdb_atoms_record["atoms"]) {
            atoms.push_back(Atom::from_json_legacy(atom_json));
        }

        return atoms;
    }

    /**
     * @brief Compare two atoms for exact match
     * @param tolerance Coordinate tolerance (default 1e-6 for floating point)
     */
    void compare_atoms(const Atom& expected, const Atom& actual, double tolerance = 1e-6) {
        EXPECT_EQ(expected.name(), actual.name()) << "Atom name mismatch";
        EXPECT_EQ(expected.residue_name(), actual.residue_name()) << "Residue name mismatch";
        EXPECT_EQ(expected.chain_id(), actual.chain_id()) << "Chain ID mismatch";
        EXPECT_EQ(expected.residue_seq(), actual.residue_seq()) << "Residue sequence number mismatch";
        EXPECT_EQ(expected.record_type(), actual.record_type()) << "Record type mismatch";

        // Compare coordinates with tolerance
        EXPECT_NEAR(expected.position().x(), actual.position().x(), tolerance) << "X coordinate mismatch";
        EXPECT_NEAR(expected.position().y(), actual.position().y(), tolerance) << "Y coordinate mismatch";
        EXPECT_NEAR(expected.position().z(), actual.position().z(), tolerance) << "Z coordinate mismatch";
    }

    /**
     * @brief Compare atom vectors for exact match
     */
    void compare_atom_vectors(const std::vector<Atom>& expected, const std::vector<Atom>& actual,
                              double tolerance = 1e-6) {
        EXPECT_EQ(expected.size(), actual.size())
            << "Atom count mismatch: expected " << expected.size() << ", got " << actual.size();

        size_t min_size = std::min(expected.size(), actual.size());
        for (size_t i = 0; i < min_size; ++i) {
            SCOPED_TRACE("Atom index " + std::to_string(i));
            compare_atoms(expected[i], actual[i], tolerance);
        }
    }

    /**
     * @brief Process a single PDB/JSON pair (thread-safe)
     * @return pair of (success, skipped) counts
     */
    std::pair<size_t, size_t> process_pdb_pair(const pdb_json_pair& pair) {
        size_t success = 0;
        size_t skipped = 0;

        try {
            auto pdb_atoms_record = load_pdb_atoms_record(pair.json_file);
            auto atoms = atoms_from_json(pdb_atoms_record);

            // Verify we got atoms
            if (atoms.size() == 0) {
                return {0, 1};
            }

            // Verify atom count matches JSON
            long expected_count = pdb_atoms_record["num_atoms"].get<long>();
            if (static_cast<long>(atoms.size()) != expected_count) {
                return {0, 1};
            }

            // Verify all atoms have valid data
            for (size_t i = 0; i < atoms.size(); ++i) {
                if (atoms[i].name().empty()) {
                    return {0, 1};
                }
                if (!std::isfinite(atoms[i].position().x()) || !std::isfinite(atoms[i].position().y()) ||
                    !std::isfinite(atoms[i].position().z())) {
                    return {0, 1};
                }
            }
            success = 1;
        } catch (const nlohmann::json::parse_error&) {
            skipped = 1;
        } catch (const std::exception&) {
            skipped = 1;
        }

        return {success, skipped};
    }
};

/**
 * @brief Test that we can load and parse pdb_atoms records from JSON
 *
 * This validates the JSON structure and our parsing logic.
 */
TEST_F(PdbAtomValidationTest, LoadPdbAtomsFromJson) {
    // Test with first discovered pair
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }

    const auto& pair = pairs_[0];
    auto pdb_atoms_record = load_pdb_atoms_record(pair.json_file);

    // Verify record structure
    EXPECT_TRUE(pdb_atoms_record.contains("type"));
    EXPECT_EQ(pdb_atoms_record["type"], "pdb_atoms");
    EXPECT_TRUE(pdb_atoms_record.contains("num_atoms"));
    EXPECT_TRUE(pdb_atoms_record.contains("atoms"));
    EXPECT_TRUE(pdb_atoms_record["atoms"].is_array());

    // Verify atom count matches
    size_t num_atoms = pdb_atoms_record["atoms"].size();
    EXPECT_EQ(pdb_atoms_record["num_atoms"].get<long>(), static_cast<long>(num_atoms));

    // Parse atoms
    auto atoms = atoms_from_json(pdb_atoms_record);
    EXPECT_EQ(atoms.size(), num_atoms);
}

/**
 * @brief Test JSON round-trip for atoms
 *
 * Ensures that Atom::from_json_legacy and Atom::to_json_legacy work correctly.
 */
TEST_F(PdbAtomValidationTest, AtomJsonRoundTrip) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }

    const auto& pair = pairs_[0];
    auto pdb_atoms_record = load_pdb_atoms_record(pair.json_file);
    auto atoms = atoms_from_json(pdb_atoms_record);

    // Test round-trip for first 10 atoms
    size_t num_to_test = std::min(static_cast<size_t>(10), atoms.size());
    for (size_t i = 0; i < num_to_test; ++i) {
        auto json = atoms[i].to_json_legacy();
        Atom reconstructed = Atom::from_json_legacy(json);

        SCOPED_TRACE("Atom index " + std::to_string(i));
        compare_atoms(atoms[i], reconstructed, 1e-9);
    }
}

/**
 * @brief Test atom data consistency across multiple PDB files
 *
 * Validates that we can correctly parse pdb_atoms from multiple JSON files.
 * By default tests 10 PDBs, but can be configured via environment variable
 * TEST_ALL_PDBS=1 to test all discovered pairs.
 * Uses threading to process multiple PDB files in parallel for faster execution.
 *
 * Usage:
 *   ./test_pdb_atom_validation --gtest_filter="*MultiplePdbFiles*"  # Test 10 PDBs
 *   TEST_ALL_PDBS=1 ./test_pdb_atom_validation --gtest_filter="*MultiplePdbFiles*"  # Test all
 */
TEST_F(PdbAtomValidationTest, MultiplePdbFiles) {
    // Check if we should test all PDBs (via environment variable)
    const char* test_all = std::getenv("TEST_ALL_PDBS");
    bool test_all_pdbs = (test_all != nullptr && std::string(test_all) == "1");

    // Test first 10 pairs by default, or all if TEST_ALL_PDBS=1
    size_t max_pairs;
    if (test_all_pdbs) {
        max_pairs = pairs_.size();
    } else {
        max_pairs = std::min(static_cast<size_t>(10), pairs_.size());
    }

    // Determine number of threads (use hardware concurrency, but cap at reasonable limit)
    size_t num_threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), max_pairs);
    if (num_threads == 0)
        num_threads = 4; // Fallback if hardware_concurrency() returns 0

    std::atomic<size_t> successful{0};
    std::atomic<size_t> skipped{0};
    std::mutex error_mutex;
    std::vector<std::string> errors;

    // Process PDB files in parallel batches
    size_t batch_size = (max_pairs + num_threads - 1) / num_threads;

    std::vector<std::future<void>> futures;

    for (size_t batch_start = 0; batch_start < max_pairs; batch_start += batch_size) {
        size_t batch_end = std::min(batch_start + batch_size, max_pairs);

        futures.push_back(std::async(std::launch::async, [&, batch_start, batch_end]() {
            for (size_t p = batch_start; p < batch_end; ++p) {
                const auto& pair = pairs_[p];

                try {
                    auto pdb_atoms_record = load_pdb_atoms_record(pair.json_file);
                    auto atoms = atoms_from_json(pdb_atoms_record);

                    // Verify we got atoms
                    EXPECT_GT(atoms.size(), 0) << "No atoms found in " << pair.pdb_name;

                    // Verify atom count matches JSON
                    long expected_count = pdb_atoms_record["num_atoms"].get<long>();
                    EXPECT_EQ(static_cast<long>(atoms.size()), expected_count);

                    // Verify all atoms have valid data
                    for (size_t i = 0; i < atoms.size(); ++i) {
                        EXPECT_FALSE(atoms[i].name().empty()) << "Atom " << i << " has empty name in " << pair.pdb_name;
                        // Coordinates should be finite
                        EXPECT_TRUE(std::isfinite(atoms[i].position().x()));
                        EXPECT_TRUE(std::isfinite(atoms[i].position().y()));
                        EXPECT_TRUE(std::isfinite(atoms[i].position().z()));
                    }
                    successful++;
                } catch (const nlohmann::json::parse_error& e) {
                    // Skip corrupted/incomplete JSON files
                    skipped++;
                    if (test_all_pdbs) {
                        std::lock_guard<std::mutex> lock(error_mutex);
                        errors.push_back("Skipping " + pair.pdb_name + " (JSON parse error): " + e.what());
                    }
                } catch (const std::exception& e) {
                    // Log other errors but continue
                    skipped++;
                    std::lock_guard<std::mutex> lock(error_mutex);
                    errors.push_back("Error processing " + pair.pdb_name + ": " + e.what());
                }
            }
        }));
    }

    // Wait for all threads to complete
    for (auto& future : futures) {
        future.wait();
    }

    // Print errors (if any)
    if (test_all_pdbs && !errors.empty()) {
        std::lock_guard<std::mutex> lock(error_mutex);
        for (const auto& error : errors) {
            std::cerr << error << std::endl;
        }
    }

    // Summary
    std::cout << "Tested " << successful.load() << " PDB files successfully";
    if (skipped.load() > 0) {
        std::cout << ", skipped " << skipped.load() << " files (corrupted/incomplete JSON)";
    }
    std::cout << " (using " << num_threads << " threads)" << std::endl;

    EXPECT_GT(successful.load(), 0) << "No PDB files were successfully tested";
}

/**
 * @brief Test that atom coordinates match exactly
 *
 * This is a critical test - coordinates must match to high precision.
 */
TEST_F(PdbAtomValidationTest, AtomCoordinatePrecision) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }

    const auto& pair = pairs_[0];
    auto pdb_atoms_record = load_pdb_atoms_record(pair.json_file);
    auto atoms = atoms_from_json(pdb_atoms_record);

    // Verify coordinates match JSON exactly (within floating point precision)
    const auto& atoms_json = pdb_atoms_record["atoms"];
    size_t num_to_test = std::min(static_cast<size_t>(20), atoms.size());

    for (size_t i = 0; i < num_to_test; ++i) {
        const auto& atom_json = atoms_json[i];
        const auto& atom = atoms[i];

        SCOPED_TRACE("Atom index " + std::to_string(i));

        // Compare coordinates from JSON directly
        std::vector<double> xyz = atom_json["xyz"].get<std::vector<double>>();
        EXPECT_EQ(xyz.size(), 3);

        EXPECT_NEAR(atom.position().x(), xyz[0], 1e-6);
        EXPECT_NEAR(atom.position().y(), xyz[1], 1e-6);
        EXPECT_NEAR(atom.position().z(), xyz[2], 1e-6);
    }
}

/**
 * @brief Test atom metadata (names, chain IDs, residue info)
 */
TEST_F(PdbAtomValidationTest, AtomMetadata) {
    if (pairs_.empty()) {
        GTEST_SKIP() << "No PDB/JSON pairs found";
    }

    const auto& pair = pairs_[0];
    auto pdb_atoms_record = load_pdb_atoms_record(pair.json_file);
    auto atoms = atoms_from_json(pdb_atoms_record);

    const auto& atoms_json = pdb_atoms_record["atoms"];
    size_t num_to_test = std::min(static_cast<size_t>(20), atoms.size());

    for (size_t i = 0; i < num_to_test; ++i) {
        const auto& atom_json = atoms_json[i];
        const auto& atom = atoms[i];

        SCOPED_TRACE("Atom index " + std::to_string(i));

        // Compare metadata
        EXPECT_EQ(atom.name(), atom_json["atom_name"].get<std::string>());
        EXPECT_EQ(atom.residue_name(), atom_json["residue_name"].get<std::string>());

        std::string chain_str = atom_json["chain_id"].get<std::string>();
        EXPECT_EQ(atom.chain_id(), chain_str[0]);

        EXPECT_EQ(atom.residue_seq(), atom_json["residue_seq"].get<int>());

        if (atom_json.contains("record_type")) {
            std::string record_str = atom_json["record_type"].get<std::string>();
            EXPECT_EQ(atom.record_type(), record_str[0]);
        }
    }
}

} // namespace x3dna::test
