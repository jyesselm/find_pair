/**
 * @file structure.hpp
 * @brief Structure class representing a complete PDB structure
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <x3dna/core/chain.hpp>
#include <x3dna/core/residue.hpp>
#include <x3dna/core/structure_legacy_order.hpp>

namespace x3dna {
namespace core {

/**
 * @class Structure
 * @brief Represents a complete PDB structure with chains, residues, and atoms
 */
class Structure {
public:
    /**
     * @brief Default constructor
     */
    Structure() = default;

    /**
     * @brief Constructor with PDB ID
     * @param pdb_id PDB identifier (e.g., "157D", "100D")
     */
    explicit Structure(const std::string& pdb_id) : pdb_id_(pdb_id) {}

    // Getters
    [[nodiscard]] const std::string& pdb_id() const {
        return pdb_id_;
    }
    [[nodiscard]] const std::vector<Chain>& chains() const {
        return chains_;
    }
    std::vector<Chain>& chains() {
        return chains_;
    }
    [[nodiscard]] size_t num_chains() const {
        return chains_.size();
    }

    // Container-like interface for chains
    [[nodiscard]] auto begin() const {
        return chains_.begin();
    }
    [[nodiscard]] auto end() const {
        return chains_.end();
    }
    [[nodiscard]] auto begin() {
        return chains_.begin();
    }
    [[nodiscard]] auto end() {
        return chains_.end();
    }
    [[nodiscard]] size_t size() const {
        return chains_.size();
    }
    [[nodiscard]] bool empty() const {
        return chains_.empty();
    }
    [[nodiscard]] const Chain& operator[](size_t idx) const {
        return chains_[idx];
    }
    [[nodiscard]] Chain& operator[](size_t idx) {
        return chains_[idx];
    }

    /**
     * @brief Get total number of residues in all chains
     */
    [[nodiscard]] size_t num_residues() const {
        size_t count = 0;
        for (const auto& chain : chains_) {
            count += chain.num_residues();
        }
        return count;
    }

    /**
     * @brief Get total number of atoms in all chains
     */
    [[nodiscard]] size_t num_atoms() const {
        size_t count = 0;
        for (const auto& chain : chains_) {
            count += chain.num_atoms();
        }
        return count;
    }

    // Setters
    void set_pdb_id(const std::string& pdb_id) {
        pdb_id_ = pdb_id;
    }

    /**
     * @brief Add a chain to this structure
     */
    void add_chain(const Chain& chain) {
        chains_.push_back(chain);
    }

    /**
     * @brief Set legacy indices on all atoms in this structure
     * @param atom_idx_map Map from (chain_id, residue_seq, insertion, atom_name) -> legacy_atom_idx
     * @param residue_idx_map Map from (chain_id, residue_seq, insertion) -> legacy_residue_idx
     */
    void set_legacy_indices(const std::map<std::tuple<std::string, int, std::string, std::string>, int>& atom_idx_map,
                            const std::map<std::tuple<std::string, int, std::string>, int>& residue_idx_map) {
        for (auto& chain : chains_) {
            for (auto& residue : chain.residues()) {
                const std::string& chain_id = residue.chain_id();
                int residue_seq = residue.seq_num();
                const std::string& insertion = residue.insertion();

                // Get legacy residue index for this residue
                auto residue_key = std::make_tuple(chain_id, residue_seq, insertion);
                auto residue_it = residue_idx_map.find(residue_key);
                int legacy_residue_idx = (residue_it != residue_idx_map.end()) ? residue_it->second : 0;

                // Set legacy_residue_idx on the residue itself
                if (legacy_residue_idx > 0) {
                    residue.set_legacy_residue_idx(legacy_residue_idx);
                }

                // Set legacy atom indices
                for (auto& atom : residue.atoms()) {
                    auto atom_key = std::make_tuple(chain_id, residue_seq, insertion, atom.name());
                    auto atom_it = atom_idx_map.find(atom_key);
                    if (atom_it != atom_idx_map.end()) {
                        atom.set_legacy_atom_idx(atom_it->second);
                    }
                }
            }
        }
    }

    /**
     * @brief Get all residues from all chains
     * @return Vector of residue pointers (non-owning)
     */
    [[nodiscard]] std::vector<const Residue*> all_residues() const {
        std::vector<const Residue*> residues;
        for (const auto& chain : chains_) {
            for (const auto& residue : chain.residues()) {
                residues.push_back(&residue);
            }
        }
        return residues;
    }

    /**
     * @brief Get all nucleotide residues from all chains
     * @return Vector of nucleotide residue pointers (non-owning)
     */
    [[nodiscard]] std::vector<const Residue*> nucleotides() const {
        std::vector<const Residue*> nts;
        for (const auto& chain : chains_) {
            for (const auto& residue : chain.residues()) {
                if (residue.is_nucleotide()) {
                    nts.push_back(&residue);
                }
            }
        }
        return nts;
    }

    /**
     * @brief Get all residues in legacy order (PDB file order)
     *
     * Returns residues in the same order as legacy's residue_idx() function:
     * - Processes atoms in PDB file order (by line_number)
     * - Groups by (ResName, ChainID, ResSeq, insertion)
     * - Returns unique residues in order they first appear
     *
     * This matches the legacy residue indexing used throughout the codebase.
     *
     * @return Vector of residue pointers in legacy order (non-owning)
     */
    [[nodiscard]] std::vector<const Residue*> residues_in_legacy_order() const;

    /**
     * @brief Get residue by legacy index (1-based)
     *
     * Finds the residue that would be at the given legacy index when counting
     * in legacy order (PDB file order).
     *
     * @param legacy_idx Legacy residue index (1-based)
     * @return Pointer to residue, or nullptr if not found
     */
    [[nodiscard]] const Residue* get_residue_by_legacy_idx(int legacy_idx) const;

    /**
     * @brief Get legacy index for a residue
     *
     * Returns the legacy index (1-based) for the given residue when counting
     * in legacy order (PDB file order).
     *
     * @param residue The residue to find the index for
     * @return Legacy residue index (1-based), or 0 if not found
     */
    [[nodiscard]] int get_legacy_idx_for_residue(const Residue* residue) const;

    /**
     * @brief Find chain by ID
     * @param chain_id Chain identifier
     * @return Optional chain if found
     */
    [[nodiscard]] std::optional<Chain> find_chain(const std::string& chain_id) const {
        for (const auto& chain : chains_) {
            if (chain.chain_id() == chain_id) {
                return chain;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Set record type for a residue
     * @param chain_id Chain identifier
     * @param seq_num Sequence number
     * @param insertion Insertion code
     * @param record_type 'A' for ATOM, 'H' for HETATM
     */
    void set_residue_record_type(const std::string& chain_id, int seq_num, const std::string& insertion,
                                 char record_type) {
        residue_record_types_[std::make_tuple(chain_id, seq_num, insertion)] = record_type;
    }

    /**
     * @brief Get record type for a residue
     * @param chain_id Chain identifier
     * @param seq_num Sequence number
     * @param insertion Insertion code
     * @return 'A' for ATOM, 'H' for HETATM (default 'A' if not found)
     */
    [[nodiscard]] char get_residue_record_type(const std::string& chain_id, int seq_num,
                                               const std::string& insertion) const {
        auto key = std::make_tuple(chain_id, seq_num, insertion);
        auto it = residue_record_types_.find(key);
        return (it != residue_record_types_.end()) ? it->second : 'A';
    }

private:
    std::string pdb_id_;        // PDB identifier
    std::vector<Chain> chains_; // Chains in this structure
    // Note: Base pairs are NOT part of Structure - they are derived/calculated data
    // and should be managed separately (e.g., in algorithms/protocols)

    // Map from (chain_id, seq_num, insertion) to record_type ('A' or 'H')
    std::map<std::tuple<std::string, int, std::string>, char> residue_record_types_;

    // Structure resolution in Angstroms (0.0 = unknown/not applicable)
    double resolution_ = 0.0;

public:
    // Resolution accessors
    /**
     * @brief Get structure resolution
     * @return Resolution in Angstroms (0.0 if unknown)
     */
    [[nodiscard]] double resolution() const { return resolution_; }

    /**
     * @brief Set structure resolution
     * @param res Resolution in Angstroms
     */
    void set_resolution(double res) { resolution_ = res; }

    /**
     * @brief Check if resolution is known
     * @return true if resolution > 0
     */
    [[nodiscard]] bool has_resolution() const { return resolution_ > 0.0; }
};

} // namespace core
} // namespace x3dna
