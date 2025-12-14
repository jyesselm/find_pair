/**
 * @file pair_validation_debugger.hpp
 * @brief Debug infrastructure for comparing pair validation between legacy and modern code
 *
 * Enable debugging by setting environment variable:
 *   X3DNA_DEBUG_PAIRS=1         - Enable all pair debugging
 *   X3DNA_DEBUG_PAIRS=1EHZ      - Debug specific PDB
 *   X3DNA_DEBUG_PAIRS=1EHZ:1,72 - Debug specific pair in specific PDB
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <optional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <nlohmann/json.hpp>

namespace x3dna {
namespace debug {

/**
 * @brief Stores validation details for a single pair
 */
struct PairValidationDetails {
    int base_i = 0;
    int base_j = 0;

    // Geometric checks
    double dorg = 0.0;
    double d_v = 0.0;
    double plane_angle = 0.0;
    double dNN = 0.0;
    double overlap_area = 0.0;

    // Direction vectors
    double dir_x = 0.0;
    double dir_y = 0.0;
    double dir_z = 0.0;

    // Validation results
    bool distance_check = false;
    bool d_v_check = false;
    bool plane_angle_check = false;
    bool dNN_check = false;
    bool overlap_check = false;
    bool hbond_check = false;
    bool is_valid = false;

    // Quality scoring
    double quality_score = 0.0;
    double hbond_adjustment = 0.0;
    double adjusted_quality = 0.0;
    int bp_type_id = -1;

    // H-bond details
    int num_base_hbonds = 0;
    int num_o2_hbonds = 0;
    int num_good_hbonds = 0;

    // Convert to JSON for comparison
    nlohmann::json to_json() const;

    // Load from legacy JSON file
    static std::optional<PairValidationDetails> from_legacy_json(
        const nlohmann::json& legacy_data, int base_i, int base_j);
};

/**
 * @brief Comparison result between legacy and modern validation
 */
struct ValidationComparison {
    PairValidationDetails legacy;
    PairValidationDetails modern;

    std::vector<std::string> differences;
    bool matches = true;

    // Generate detailed comparison report
    std::string generate_report() const;
};

/**
 * @brief Debug infrastructure for pair validation
 *
 * Enabled via X3DNA_DEBUG_PAIRS environment variable.
 * Provides:
 * 1. Side-by-side comparison of validation results
 * 2. Field-level diff highlighting
 * 3. JSON export for analysis
 */
class PairValidationDebugger {
public:
    /**
     * @brief Get singleton instance
     */
    static PairValidationDebugger& instance();

    /**
     * @brief Check if debugging is enabled
     */
    bool is_enabled() const { return enabled_; }

    /**
     * @brief Check if debugging is enabled for specific PDB
     */
    bool is_enabled_for_pdb(const std::string& pdb_id) const;

    /**
     * @brief Check if debugging is enabled for specific pair
     */
    bool should_debug_pair(const std::string& pdb_id, int base_i, int base_j) const;

    /**
     * @brief Set current PDB being processed
     */
    void set_current_pdb(const std::string& pdb_id);

    /**
     * @brief Record modern validation result
     */
    void record_modern_validation(const PairValidationDetails& details);

    /**
     * @brief Load legacy validation results from JSON directory
     */
    bool load_legacy_results(const std::string& json_dir);

    /**
     * @brief Compare specific pair between legacy and modern
     */
    ValidationComparison compare_pair(int base_i, int base_j) const;

    /**
     * @brief Compare all recorded pairs
     */
    std::vector<ValidationComparison> compare_all_pairs() const;

    /**
     * @brief Print comparison report to stderr
     */
    void print_comparison_report() const;

    /**
     * @brief Export comparison to JSON file
     */
    void export_comparison_json(const std::string& output_path) const;

    /**
     * @brief Log a debug message (only if enabled)
     */
    void log(const std::string& message) const;

    /**
     * @brief Log a debug message for specific pair
     */
    void log_pair(int base_i, int base_j, const std::string& message) const;

private:
    PairValidationDebugger();
    ~PairValidationDebugger() = default;

    // Non-copyable
    PairValidationDebugger(const PairValidationDebugger&) = delete;
    PairValidationDebugger& operator=(const PairValidationDebugger&) = delete;

    // Parse environment variable
    void parse_env_config();

    // Check if pair key matches filter
    bool matches_pair_filter(int base_i, int base_j) const;

    bool enabled_ = false;
    std::string filter_pdb_;
    std::vector<std::pair<int, int>> filter_pairs_;

    std::string current_pdb_;
    std::map<std::pair<int, int>, PairValidationDetails> modern_results_;
    std::map<std::pair<int, int>, PairValidationDetails> legacy_results_;

    mutable std::ofstream log_file_;
};

/**
 * @brief RAII helper for scoped debugging
 */
class ScopedPairDebug {
public:
    ScopedPairDebug(const std::string& pdb_id, const std::string& json_dir);
    ~ScopedPairDebug();

    bool is_enabled() const;

private:
    std::string pdb_id_;
    bool was_enabled_;
};

// Convenience macros for conditional debugging
#define X3DNA_DEBUG_PAIR_IF_ENABLED(base_i, base_j, msg) \
    do { \
        auto& dbg = x3dna::debug::PairValidationDebugger::instance(); \
        if (dbg.should_debug_pair(dbg.current_pdb(), base_i, base_j)) { \
            dbg.log_pair(base_i, base_j, msg); \
        } \
    } while(0)

#define X3DNA_DEBUG_LOG_IF_ENABLED(msg) \
    do { \
        auto& dbg = x3dna::debug::PairValidationDebugger::instance(); \
        if (dbg.is_enabled()) { \
            dbg.log(msg); \
        } \
    } while(0)

} // namespace debug
} // namespace x3dna
