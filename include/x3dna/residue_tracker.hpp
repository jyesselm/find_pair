#pragma once

#include <vector>
#include <string>
#include <optional>

namespace x3dna {

/**
 * @brief Record for tracking a single residue through the parsing and filtering process
 *
 * Tracks:
 * - Order residue was read from PDB file
 * - Whether it was filtered out (and why)
 * - Final modern index (0-based)
 * - Legacy index from JSON (1-based)
 * - PDB properties for matching
 */
struct ResidueRecord {
    // Tracking indices
    int read_index;            // Order read from PDB (0-based)
    int legacy_index;          // Index from legacy JSON (1-based, -1 if not set)
    int modern_index;          // Final index after filtering (0-based, -1 if filtered)
    bool filtered;             // Was this residue filtered out?
    std::string filter_reason; // Why was it filtered (empty if not filtered)

    // PDB properties for matching
    std::string chain_id;
    int residue_seq;
    std::string insertion;
    std::string residue_name;

    ResidueRecord(int read_idx, const std::string& chain, int seq, const std::string& ins, const std::string& name)
        : read_index(read_idx), legacy_index(-1), modern_index(-1), filtered(false), filter_reason(""), chain_id(chain),
          residue_seq(seq), insertion(ins), residue_name(name) {}
};

/**
 * @brief Tracks residues as they are read and filtered
 *
 * This is THE CRITICAL FIX for index matching.
 *
 * Purpose:
 * - Track every residue in the order it's read from PDB
 * - Track which residues get filtered out (and why)
 * - Map modern indices to legacy indices
 * - Validate 100% match before allowing any comparisons
 *
 * Usage:
 * 1. As each residue is read from PDB: add_residue()
 * 2. When a residue is filtered: mark_filtered()
 * 3. After final residue list is created: assign_modern_index()
 * 4. Load legacy indices from JSON: load_legacy_indices()
 * 5. Validate match: validate()
 * 6. Export for debugging: export_mapping()
 */
class ResidueTracker {
public:
    ResidueTracker() = default;

    /**
     * @brief Add residue as it's read from PDB (in order)
     *
     * Call this for EVERY residue as it's read, before any filtering.
     * Order matters - this creates the read_index.
     */
    void add_residue(const std::string& chain_id, int residue_seq, const std::string& insertion,
                     const std::string& residue_name);

    /**
     * @brief Mark residue as filtered (won't appear in final list)
     *
     * @param read_index The read_index of the residue to mark as filtered
     * @param reason Why it was filtered (e.g., "no ring atoms", "failed frame calc")
     */
    void mark_filtered(int read_index, const std::string& reason);

    /**
     * @brief Assign final modern index to a residue
     *
     * Call this after the final residue list is created, for each residue
     * that made it through filtering.
     *
     * @param read_index The read_index of the residue
     * @param modern_index The final 0-based index in the modern code
     */
    void assign_modern_index(int read_index, int modern_index);

    /**
     * @brief Load legacy indices from JSON file
     *
     * Reads legacy base_frame_calc JSON and matches residues by PDB properties
     * (chain_id, residue_seq, insertion).
     *
     * @param legacy_json_path Path to legacy base_frame_calc JSON file
     * @return true if loaded successfully, false otherwise
     */
    bool load_legacy_indices(const std::string& legacy_json_path);

    /**
     * @brief Validation result structure
     */
    struct ValidationResult {
        bool success;
        int num_residues_read;
        int num_legacy;
        int num_modern;
        int num_filtered;
        int num_matched;
        int num_unmatched;
        std::vector<std::string> errors;

        std::string to_string() const;
    };

    /**
     * @brief Validate that our indices match legacy exactly
     *
     * This is THE critical check. Must return success=true before
     * any comparisons are allowed.
     *
     * Checks:
     * - num_modern == num_legacy
     * - All non-filtered residues have both modern and legacy indices
     * - No duplicate indices
     *
     * @return ValidationResult with detailed information
     */
    ValidationResult validate() const;

    /**
     * @brief Export mapping to JSON for debugging
     *
     * Creates a JSON file with all residue tracking information.
     * Useful for investigating validation failures.
     *
     * @param output_path Path to write JSON file
     */
    void export_mapping(const std::string& output_path) const;

    /**
     * @brief Get legacy index for a modern index
     *
     * @param modern_index 0-based modern index
     * @return Legacy index (1-based) or nullopt if not found
     */
    std::optional<int> get_legacy_index(int modern_index) const;

    /**
     * @brief Get modern index for a legacy index
     *
     * @param legacy_index 1-based legacy index
     * @return Modern index (0-based) or nullopt if not found
     */
    std::optional<int> get_modern_index(int legacy_index) const;

    /**
     * @brief Get all residues (for debugging)
     */
    const std::vector<ResidueRecord>& get_residues() const {
        return residues_;
    }

    /**
     * @brief Get number of residues read
     */
    size_t size() const {
        return residues_.size();
    }

    /**
     * @brief Clear all tracking data
     */
    void clear() {
        residues_.clear();
    }

private:
    std::vector<ResidueRecord> residues_;

    /**
     * @brief Helper to find residue by PDB properties
     *
     * @return read_index if found, -1 otherwise
     */
    int find_by_pdb_props(const std::string& chain, int seq, const std::string& ins) const;
};

} // namespace x3dna
