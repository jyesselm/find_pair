/**
 * @file pair_validation_debugger.cpp
 * @brief Implementation of pair validation debugging infrastructure
 */

#include <x3dna/debug/pair_validation_debugger.hpp>
#include <x3dna/config/config_manager.hpp>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <set>

namespace x3dna {
namespace debug {

// PairValidationDetails implementation

nlohmann::json PairValidationDetails::to_json() const {
    nlohmann::json j;
    j["base_i"] = base_i;
    j["base_j"] = base_j;

    j["geometry"]["dorg"] = dorg;
    j["geometry"]["d_v"] = d_v;
    j["geometry"]["plane_angle"] = plane_angle;
    j["geometry"]["dNN"] = dNN;
    j["geometry"]["overlap_area"] = overlap_area;

    j["direction"]["dir_x"] = dir_x;
    j["direction"]["dir_y"] = dir_y;
    j["direction"]["dir_z"] = dir_z;

    j["checks"]["distance_check"] = distance_check;
    j["checks"]["d_v_check"] = d_v_check;
    j["checks"]["plane_angle_check"] = plane_angle_check;
    j["checks"]["dNN_check"] = dNN_check;
    j["checks"]["overlap_check"] = overlap_check;
    j["checks"]["hbond_check"] = hbond_check;
    j["checks"]["is_valid"] = is_valid;

    j["quality"]["quality_score"] = quality_score;
    j["quality"]["hbond_adjustment"] = hbond_adjustment;
    j["quality"]["adjusted_quality"] = adjusted_quality;
    j["quality"]["bp_type_id"] = bp_type_id;

    j["hbonds"]["num_base_hbonds"] = num_base_hbonds;
    j["hbonds"]["num_o2_hbonds"] = num_o2_hbonds;
    j["hbonds"]["num_good_hbonds"] = num_good_hbonds;

    return j;
}

std::optional<PairValidationDetails> PairValidationDetails::from_legacy_json(
    const nlohmann::json& legacy_data, int base_i, int base_j) {

    PairValidationDetails details;
    details.base_i = base_i;
    details.base_j = base_j;

    // Search through legacy data for this pair
    for (const auto& record : legacy_data) {
        int rec_i = record.value("base_i", 0);
        int rec_j = record.value("base_j", 0);

        // Normalize pair key (min, max)
        int norm_i = std::min(rec_i, rec_j);
        int norm_j = std::max(rec_i, rec_j);
        int search_i = std::min(base_i, base_j);
        int search_j = std::max(base_i, base_j);

        if (norm_i == search_i && norm_j == search_j) {
            // Found the pair - extract fields
            details.is_valid = record.value("is_valid", 0) == 1;
            details.bp_type_id = record.value("bp_type_id", -1);

            // Direction vectors
            if (record.contains("direction_vectors")) {
                auto& dir = record["direction_vectors"];
                details.dir_x = dir.value("dir_x", 0.0);
                details.dir_y = dir.value("dir_y", 0.0);
                details.dir_z = dir.value("dir_z", 0.0);
            }

            // Calculated values (from distance_checks usually)
            if (record.contains("calculated_values")) {
                auto& calc = record["calculated_values"];
                details.dorg = calc.value("dorg", 0.0);
                details.d_v = calc.value("d_v", 0.0);
                details.plane_angle = calc.value("plane_angle", 0.0);
                details.dNN = calc.value("dNN", 0.0);
                details.quality_score = calc.value("quality_score", 0.0);
            }

            // Validation checks
            if (record.contains("validation_checks")) {
                auto& checks = record["validation_checks"];
                details.distance_check = checks.value("distance_check", false);
                details.d_v_check = checks.value("d_v_check", false);
                details.plane_angle_check = checks.value("plane_angle_check", false);
                details.dNN_check = checks.value("dNN_check", false);
            }

            return details;
        }
    }

    return std::nullopt;
}

// ValidationComparison implementation

std::string ValidationComparison::generate_report() const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);

    ss << "=== Pair (" << legacy.base_i << ", " << legacy.base_j << ") Comparison ===\n\n";

    auto compare_field = [&](const std::string& name, double leg, double mod, double tol = 1e-6) {
        double diff = std::abs(leg - mod);
        bool match = diff <= tol;
        ss << "  " << std::setw(20) << std::left << name << ": "
           << std::setw(12) << leg << " vs " << std::setw(12) << mod;
        if (!match) {
            ss << " [DIFF: " << diff << "]";
        } else {
            ss << " [OK]";
        }
        ss << "\n";
    };

    auto compare_bool = [&](const std::string& name, bool leg, bool mod) {
        bool match = (leg == mod);
        ss << "  " << std::setw(20) << std::left << name << ": "
           << (leg ? "true" : "false") << " vs " << (mod ? "true" : "false");
        if (!match) {
            ss << " [MISMATCH]";
        } else {
            ss << " [OK]";
        }
        ss << "\n";
    };

    auto compare_int = [&](const std::string& name, int leg, int mod) {
        bool match = (leg == mod);
        ss << "  " << std::setw(20) << std::left << name << ": "
           << leg << " vs " << mod;
        if (!match) {
            ss << " [MISMATCH]";
        } else {
            ss << " [OK]";
        }
        ss << "\n";
    };

    ss << "--- Geometry ---\n";
    compare_field("dorg", legacy.dorg, modern.dorg);
    compare_field("d_v", legacy.d_v, modern.d_v);
    compare_field("plane_angle", legacy.plane_angle, modern.plane_angle);
    compare_field("dNN", legacy.dNN, modern.dNN);
    compare_field("overlap_area", legacy.overlap_area, modern.overlap_area);

    ss << "\n--- Direction Vectors ---\n";
    compare_field("dir_x", legacy.dir_x, modern.dir_x);
    compare_field("dir_y", legacy.dir_y, modern.dir_y);
    compare_field("dir_z", legacy.dir_z, modern.dir_z);

    ss << "\n--- Validation Checks ---\n";
    compare_bool("distance_check", legacy.distance_check, modern.distance_check);
    compare_bool("d_v_check", legacy.d_v_check, modern.d_v_check);
    compare_bool("plane_angle_check", legacy.plane_angle_check, modern.plane_angle_check);
    compare_bool("dNN_check", legacy.dNN_check, modern.dNN_check);
    compare_bool("overlap_check", legacy.overlap_check, modern.overlap_check);
    compare_bool("hbond_check", legacy.hbond_check, modern.hbond_check);
    compare_bool("is_valid", legacy.is_valid, modern.is_valid);

    ss << "\n--- Quality Scoring ---\n";
    compare_field("quality_score", legacy.quality_score, modern.quality_score);
    compare_field("hbond_adjustment", legacy.hbond_adjustment, modern.hbond_adjustment);
    compare_field("adjusted_quality", legacy.adjusted_quality, modern.adjusted_quality);
    compare_int("bp_type_id", legacy.bp_type_id, modern.bp_type_id);

    ss << "\n--- H-bond Counts ---\n";
    compare_int("num_base_hbonds", legacy.num_base_hbonds, modern.num_base_hbonds);
    compare_int("num_o2_hbonds", legacy.num_o2_hbonds, modern.num_o2_hbonds);
    compare_int("num_good_hbonds", legacy.num_good_hbonds, modern.num_good_hbonds);

    ss << "\n--- Summary ---\n";
    if (matches) {
        ss << "  RESULT: MATCH\n";
    } else {
        ss << "  RESULT: MISMATCH\n";
        ss << "  Differences:\n";
        for (const auto& diff : differences) {
            ss << "    - " << diff << "\n";
        }
    }

    return ss.str();
}

// PairValidationDebugger implementation

PairValidationDebugger::PairValidationDebugger() {
    parse_env_config();
}

PairValidationDebugger& PairValidationDebugger::instance() {
    static PairValidationDebugger instance;
    return instance;
}

void PairValidationDebugger::parse_env_config() {
    // Use ConfigManager to get debug settings (which reads from env vars)
    auto& cfg = config::ConfigManager::instance();
    cfg.init_debug_from_environment();
    const auto& debug_cfg = cfg.debug_config();

    if (!debug_cfg.debug_pairs) {
        enabled_ = false;
        return;
    }

    enabled_ = true;
    std::string config = debug_cfg.debug_pairs_filter;

    // If no filter, debug all pairs
    if (config.empty()) {
        filter_pdb_.clear();
        filter_pairs_.clear();
        return;
    }

    // Parse format: "PDB_ID" or "PDB_ID:i,j" or "PDB_ID:i,j;i2,j2"
    // Check for colon separator (PDB:pairs format)
    size_t colon_pos = config.find(':');
    if (colon_pos != std::string::npos) {
        filter_pdb_ = config.substr(0, colon_pos);
        std::string pairs_str = config.substr(colon_pos + 1);

        // Parse pairs (format: "i,j" or "i,j;i2,j2")
        std::istringstream iss(pairs_str);
        std::string pair_token;
        while (std::getline(iss, pair_token, ';')) {
            size_t comma_pos = pair_token.find(',');
            if (comma_pos != std::string::npos) {
                try {
                    int i = std::stoi(pair_token.substr(0, comma_pos));
                    int j = std::stoi(pair_token.substr(comma_pos + 1));
                    filter_pairs_.push_back({std::min(i, j), std::max(i, j)});
                } catch (...) {
                    // Invalid pair format, skip
                }
            }
        }
    } else {
        // Just PDB ID
        filter_pdb_ = config;
    }

    std::cerr << "[X3DNA_DEBUG] Pair validation debugging enabled";
    if (!filter_pdb_.empty()) {
        std::cerr << " for PDB: " << filter_pdb_;
    }
    if (!filter_pairs_.empty()) {
        std::cerr << " pairs: ";
        for (const auto& p : filter_pairs_) {
            std::cerr << "(" << p.first << "," << p.second << ") ";
        }
    }
    std::cerr << "\n";
}

bool PairValidationDebugger::is_enabled_for_pdb(const std::string& pdb_id) const {
    if (!enabled_) return false;
    if (filter_pdb_.empty()) return true;

    // Case-insensitive comparison
    std::string pdb_upper = pdb_id;
    std::string filter_upper = filter_pdb_;
    std::transform(pdb_upper.begin(), pdb_upper.end(), pdb_upper.begin(), ::toupper);
    std::transform(filter_upper.begin(), filter_upper.end(), filter_upper.begin(), ::toupper);

    return pdb_upper == filter_upper;
}

bool PairValidationDebugger::should_debug_pair(const std::string& pdb_id, int base_i, int base_j) const {
    if (!is_enabled_for_pdb(pdb_id)) return false;
    if (filter_pairs_.empty()) return true;

    return matches_pair_filter(base_i, base_j);
}

bool PairValidationDebugger::matches_pair_filter(int base_i, int base_j) const {
    int norm_i = std::min(base_i, base_j);
    int norm_j = std::max(base_i, base_j);

    for (const auto& p : filter_pairs_) {
        if (p.first == norm_i && p.second == norm_j) {
            return true;
        }
    }
    return false;
}

void PairValidationDebugger::set_current_pdb(const std::string& pdb_id) {
    current_pdb_ = pdb_id;
    modern_results_.clear();
    legacy_results_.clear();

    if (is_enabled_for_pdb(pdb_id)) {
        std::cerr << "[X3DNA_DEBUG] Processing PDB: " << pdb_id << "\n";
    }
}

void PairValidationDebugger::record_modern_validation(const PairValidationDetails& details) {
    if (!is_enabled_for_pdb(current_pdb_)) return;

    int norm_i = std::min(details.base_i, details.base_j);
    int norm_j = std::max(details.base_i, details.base_j);
    modern_results_[{norm_i, norm_j}] = details;
}

bool PairValidationDebugger::load_legacy_results(const std::string& json_dir) {
    if (!is_enabled_for_pdb(current_pdb_)) return false;

    // Try to load pair_validation JSON
    std::string path = json_dir + "/pair_validation/" + current_pdb_ + ".json";
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "[X3DNA_DEBUG] Could not open legacy JSON: " << path << "\n";
        return false;
    }

    try {
        nlohmann::json legacy_data = nlohmann::json::parse(file);
        for (const auto& record : legacy_data) {
            int base_i = record.value("base_i", 0);
            int base_j = record.value("base_j", 0);

            auto details = PairValidationDetails::from_legacy_json(legacy_data, base_i, base_j);
            if (details) {
                int norm_i = std::min(base_i, base_j);
                int norm_j = std::max(base_i, base_j);
                legacy_results_[{norm_i, norm_j}] = *details;
            }
        }
        std::cerr << "[X3DNA_DEBUG] Loaded " << legacy_results_.size()
                  << " legacy pair validation records\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[X3DNA_DEBUG] Error parsing legacy JSON: " << e.what() << "\n";
        return false;
    }
}

ValidationComparison PairValidationDebugger::compare_pair(int base_i, int base_j) const {
    ValidationComparison result;

    int norm_i = std::min(base_i, base_j);
    int norm_j = std::max(base_i, base_j);

    auto leg_it = legacy_results_.find({norm_i, norm_j});
    auto mod_it = modern_results_.find({norm_i, norm_j});

    if (leg_it == legacy_results_.end()) {
        result.differences.push_back("Pair not found in legacy results");
        result.matches = false;
        return result;
    }
    if (mod_it == modern_results_.end()) {
        result.differences.push_back("Pair not found in modern results");
        result.matches = false;
        return result;
    }

    result.legacy = leg_it->second;
    result.modern = mod_it->second;

    // Compare all fields
    constexpr double tol = 1e-6;

    auto check_float = [&](const std::string& name, double leg, double mod) {
        if (std::abs(leg - mod) > tol) {
            result.differences.push_back(name + ": " + std::to_string(leg) + " vs " + std::to_string(mod));
            result.matches = false;
        }
    };

    auto check_bool = [&](const std::string& name, bool leg, bool mod) {
        if (leg != mod) {
            result.differences.push_back(name + ": " + (leg ? "true" : "false") + " vs " + (mod ? "true" : "false"));
            result.matches = false;
        }
    };

    auto check_int = [&](const std::string& name, int leg, int mod) {
        if (leg != mod) {
            result.differences.push_back(name + ": " + std::to_string(leg) + " vs " + std::to_string(mod));
            result.matches = false;
        }
    };

    // Geometry
    check_float("dorg", result.legacy.dorg, result.modern.dorg);
    check_float("d_v", result.legacy.d_v, result.modern.d_v);
    check_float("plane_angle", result.legacy.plane_angle, result.modern.plane_angle);
    check_float("dNN", result.legacy.dNN, result.modern.dNN);
    check_float("overlap_area", result.legacy.overlap_area, result.modern.overlap_area);

    // Direction vectors
    check_float("dir_x", result.legacy.dir_x, result.modern.dir_x);
    check_float("dir_y", result.legacy.dir_y, result.modern.dir_y);
    check_float("dir_z", result.legacy.dir_z, result.modern.dir_z);

    // Checks
    check_bool("is_valid", result.legacy.is_valid, result.modern.is_valid);
    check_int("bp_type_id", result.legacy.bp_type_id, result.modern.bp_type_id);

    // Quality
    check_float("quality_score", result.legacy.quality_score, result.modern.quality_score);
    check_float("adjusted_quality", result.legacy.adjusted_quality, result.modern.adjusted_quality);

    return result;
}

std::vector<ValidationComparison> PairValidationDebugger::compare_all_pairs() const {
    std::vector<ValidationComparison> comparisons;

    // Collect all unique pair keys
    std::set<std::pair<int, int>> all_keys;
    for (const auto& [key, _] : legacy_results_) {
        all_keys.insert(key);
    }
    for (const auto& [key, _] : modern_results_) {
        all_keys.insert(key);
    }

    for (const auto& key : all_keys) {
        if (filter_pairs_.empty() || matches_pair_filter(key.first, key.second)) {
            comparisons.push_back(compare_pair(key.first, key.second));
        }
    }

    return comparisons;
}

void PairValidationDebugger::print_comparison_report() const {
    if (!is_enabled_for_pdb(current_pdb_)) return;

    auto comparisons = compare_all_pairs();

    int total = comparisons.size();
    int matches = 0;
    int mismatches = 0;

    for (const auto& comp : comparisons) {
        if (comp.matches) {
            ++matches;
        } else {
            ++mismatches;
            std::cerr << comp.generate_report() << "\n";
        }
    }

    std::cerr << "\n=== SUMMARY ===\n";
    std::cerr << "Total pairs compared: " << total << "\n";
    std::cerr << "Matches: " << matches << "\n";
    std::cerr << "Mismatches: " << mismatches << "\n";
}

void PairValidationDebugger::export_comparison_json(const std::string& output_path) const {
    nlohmann::json output;
    output["pdb_id"] = current_pdb_;

    auto comparisons = compare_all_pairs();
    for (const auto& comp : comparisons) {
        nlohmann::json pair_json;
        pair_json["base_i"] = comp.legacy.base_i;
        pair_json["base_j"] = comp.legacy.base_j;
        pair_json["matches"] = comp.matches;
        pair_json["differences"] = comp.differences;
        pair_json["legacy"] = comp.legacy.to_json();
        pair_json["modern"] = comp.modern.to_json();
        output["pairs"].push_back(pair_json);
    }

    std::ofstream file(output_path);
    file << output.dump(2);
}

void PairValidationDebugger::log(const std::string& message) const {
    if (!enabled_) return;
    std::cerr << "[X3DNA_DEBUG] " << message << "\n";
}

void PairValidationDebugger::log_pair(int base_i, int base_j, const std::string& message) const {
    if (!enabled_) return;
    std::cerr << "[X3DNA_DEBUG] (" << base_i << "," << base_j << ") " << message << "\n";
}

// ScopedPairDebug implementation

ScopedPairDebug::ScopedPairDebug(const std::string& pdb_id, const std::string& json_dir)
    : pdb_id_(pdb_id) {
    auto& dbg = PairValidationDebugger::instance();
    was_enabled_ = dbg.is_enabled();
    if (was_enabled_) {
        dbg.set_current_pdb(pdb_id);
        dbg.load_legacy_results(json_dir);
    }
}

ScopedPairDebug::~ScopedPairDebug() {
    if (was_enabled_) {
        auto& dbg = PairValidationDebugger::instance();
        dbg.print_comparison_report();
    }
}

bool ScopedPairDebug::is_enabled() const {
    return was_enabled_;
}

} // namespace debug
} // namespace x3dna
