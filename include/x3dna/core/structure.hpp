/**
 * @file structure.hpp
 * @brief Structure class representing a complete PDB structure
 */

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <tuple>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
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
    const std::string& pdb_id() const {
        return pdb_id_;
    }
    const std::vector<Chain>& chains() const {
        return chains_;
    }
    std::vector<Chain>& chains() {
        return chains_;
    }
    size_t num_chains() const {
        return chains_.size();
    }

    /**
     * @brief Get total number of residues in all chains
     */
    size_t num_residues() const {
        size_t count = 0;
        for (const auto& chain : chains_) {
            count += chain.num_residues();
        }
        return count;
    }

    /**
     * @brief Get total number of atoms in all chains
     */
    size_t num_atoms() const {
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
    void set_legacy_indices(
        const std::map<std::tuple<char, int, char, std::string>, int>& atom_idx_map,
        const std::map<std::tuple<char, int, char>, int>& residue_idx_map) {
        for (auto& chain : chains_) {
            for (auto& residue : chain.residues()) {
                char chain_id = residue.chain_id();
                int residue_seq = residue.seq_num();
                char insertion = residue.insertion();
                
                // Get legacy residue index for this residue
                auto residue_key = std::make_tuple(chain_id, residue_seq, insertion);
                auto residue_it = residue_idx_map.find(residue_key);
                int legacy_residue_idx = (residue_it != residue_idx_map.end()) ? residue_it->second : 0;
                
                // Set legacy indices on all atoms in this residue
                for (auto& atom : residue.atoms()) {
                    auto atom_key = std::make_tuple(chain_id, residue_seq, insertion, atom.name());
                    auto atom_it = atom_idx_map.find(atom_key);
                    if (atom_it != atom_idx_map.end()) {
                        atom.set_legacy_atom_idx(atom_it->second);
                    }
                    if (legacy_residue_idx > 0) {
                        atom.set_legacy_residue_idx(legacy_residue_idx);
                    }
                }
            }
        }
    }

    /**
     * @brief Get all residues from all chains
     * @return Vector of residue pointers (non-owning)
     */
    std::vector<const Residue*> all_residues() const {
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
    std::vector<const Residue*> nucleotides() const {
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
    std::vector<const Residue*> residues_in_legacy_order() const;

    /**
     * @brief Get residue by legacy index (1-based)
     * 
     * Finds the residue that would be at the given legacy index when counting
     * in legacy order (PDB file order).
     * 
     * @param legacy_idx Legacy residue index (1-based)
     * @return Pointer to residue, or nullptr if not found
     */
    const Residue* get_residue_by_legacy_idx(int legacy_idx) const;

    /**
     * @brief Get legacy index for a residue
     * 
     * Returns the legacy index (1-based) for the given residue when counting
     * in legacy order (PDB file order).
     * 
     * @param residue The residue to find the index for
     * @return Legacy residue index (1-based), or 0 if not found
     */
    int get_legacy_idx_for_residue(const Residue* residue) const;

    /**
     * @brief Find chain by ID
     * @param chain_id Chain identifier
     * @return Optional chain if found
     */
    std::optional<Chain> find_chain(char chain_id) const {
        for (const auto& chain : chains_) {
            if (chain.chain_id() == chain_id) {
                return chain;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Convert to legacy JSON format (pdb_atoms record)
     */
    nlohmann::json to_json_legacy() const {
        nlohmann::json j;
        j["pdb_id"] = pdb_id_;
        j["num_atoms"] = static_cast<long>(num_atoms());
        j["num_residues"] = static_cast<long>(num_residues());
        j["num_chains"] = static_cast<long>(num_chains());

        // Collect all atoms from all chains
        j["atoms"] = nlohmann::json::array();
        long atom_idx = 1;
        for (const auto& chain : chains_) {
            for (const auto& residue : chain.residues()) {
                for (const auto& atom : residue.atoms()) {
                    auto atom_json = atom.to_json_legacy();
                    atom_json["atom_idx"] = atom_idx++;
                    j["atoms"].push_back(atom_json);
                }
            }
        }

        return j;
    }

    /**
     * @brief Create Structure from legacy JSON format (pdb_atoms record)
     */
    static Structure from_json_legacy(const nlohmann::json& j) {
        std::string pdb_id = j.value("pdb_id", "");
        Structure structure(pdb_id);

        if (!j.contains("atoms") || !j["atoms"].is_array()) {
            return structure;
        }

        // Group atoms by chain and residue
        std::map<std::pair<char, int>, std::vector<Atom>> residue_atoms;

        for (const auto& atom_json : j["atoms"]) {
            Atom atom = Atom::from_json_legacy(atom_json);
            char chain_id = atom.chain_id();
            int residue_seq = atom.residue_seq();
            std::string residue_name = atom.residue_name();

            std::pair<char, int> key = {chain_id, residue_seq};
            residue_atoms[key].push_back(atom);
        }

        // Create chains and residues
        std::map<char, Chain> chains;
        for (const auto& [key, atoms] : residue_atoms) {
            char chain_id = key.first;
            int residue_seq = key.second;

            if (atoms.empty()) {
                continue;
            }

            std::string residue_name = atoms[0].residue_name();
            Residue residue(residue_name, residue_seq, chain_id);

            for (const auto& atom : atoms) {
                residue.add_atom(atom);
            }

            if (chains.find(chain_id) == chains.end()) {
                chains[chain_id] = Chain(chain_id);
            }
            chains[chain_id].add_residue(residue);
        }

        // Add chains to structure
        for (auto& [chain_id, chain] : chains) {
            structure.add_chain(chain);
        }

        return structure;
    }

    /**
     * @brief Write atoms JSON to file (pdb_atoms format)
     * @param output_dir Output directory where pdb_atoms/<pdb_id>.json will be written
     */
    void write_atoms_json(const std::filesystem::path& output_dir) const {
        nlohmann::json record;
        record["num_atoms"] = num_atoms();
        record["atoms"] = nlohmann::json::array();
        
        for (const auto& chain : chains_) {
            for (const auto& residue : chain.residues()) {
                for (const auto& atom : residue.atoms()) {
                    record["atoms"].push_back(atom.to_json());  // Atom writes itself
                }
            }
        }
        
        // Write file
        std::filesystem::path file = output_dir / "pdb_atoms" / (pdb_id_ + ".json");
        std::filesystem::create_directories(file.parent_path());
        std::ofstream out(file);
        out << record.dump(2);
    }

    /**
     * @brief Convert to modern JSON format
     */
    nlohmann::json to_json() const {
        nlohmann::json j;
        j["pdb_id"] = pdb_id_;
        j["chains"] = nlohmann::json::array();
        for (const auto& chain : chains_) {
            j["chains"].push_back(chain.to_json());
        }
        return j;
    }

    /**
     * @brief Create Structure from modern JSON format
     */
    static Structure from_json(const nlohmann::json& j) {
        std::string pdb_id = j.value("pdb_id", "");
        Structure structure(pdb_id);

        if (j.contains("chains") && j["chains"].is_array()) {
            for (const auto& chain_json : j["chains"]) {
                structure.add_chain(Chain::from_json(chain_json));
            }
        }

        return structure;
    }

private:
    std::string pdb_id_;        // PDB identifier
    std::vector<Chain> chains_; // Chains in this structure
    // Note: Base pairs are NOT part of Structure - they are derived/calculated data
    // and should be managed separately (e.g., in algorithms/protocols)
};

} // namespace core
} // namespace x3dna
