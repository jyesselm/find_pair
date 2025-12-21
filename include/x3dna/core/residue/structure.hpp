/**
 * @file structure.hpp
 * @brief Structure class using polymorphic residue types
 */

#pragma once

#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <x3dna/core/residue/chain.hpp>

namespace x3dna {
namespace core {
namespace poly {


/**
 * @class Structure
 * @brief Represents a complete PDB structure with polymorphic chains and residues
 */
class Structure {
public:
    Structure() = default;

    explicit Structure(const std::string& pdb_id) : pdb_id_(pdb_id) {}

    // Move-only (chains contain unique_ptr)
    Structure(const Structure&) = delete;
    Structure& operator=(const Structure&) = delete;
    Structure(Structure&&) = default;
    Structure& operator=(Structure&&) = default;

    // Deep copy via clone()
    [[nodiscard]] Structure clone() const {
        Structure copy(pdb_id_);
        copy.chains_.reserve(chains_.size());
        for (const auto& chain : chains_) {
            copy.chains_.push_back(chain.clone());
        }
        copy.residue_record_types_ = residue_record_types_;
        return copy;
    }

    // === Identity ===
    [[nodiscard]] const std::string& pdb_id() const { return pdb_id_; }
    void set_pdb_id(const std::string& pdb_id) { pdb_id_ = pdb_id; }

    // === Chain access ===
    [[nodiscard]] size_t num_chains() const { return chains_.size(); }
    [[nodiscard]] size_t size() const { return chains_.size(); }
    [[nodiscard]] bool empty() const { return chains_.empty(); }

    [[nodiscard]] const Chain& operator[](size_t idx) const { return chains_[idx]; }
    [[nodiscard]] Chain& operator[](size_t idx) { return chains_[idx]; }

    [[nodiscard]] const Chain& at(size_t idx) const { return chains_.at(idx); }
    [[nodiscard]] Chain& at(size_t idx) { return chains_.at(idx); }

    // === Chain ownership ===
    void add_chain(Chain chain) {
        chains_.push_back(std::move(chain));
    }

    // === Iteration ===
    [[nodiscard]] auto begin() { return chains_.begin(); }
    [[nodiscard]] auto end() { return chains_.end(); }
    [[nodiscard]] auto begin() const { return chains_.begin(); }
    [[nodiscard]] auto end() const { return chains_.end(); }

    // === Counts ===
    [[nodiscard]] size_t num_residues() const {
        size_t count = 0;
        for (const auto& chain : chains_) {
            count += chain.num_residues();
        }
        return count;
    }

    [[nodiscard]] size_t num_atoms() const {
        size_t count = 0;
        for (const auto& chain : chains_) {
            count += chain.num_atoms();
        }
        return count;
    }

    // === Residue access ===
    [[nodiscard]] std::vector<IResidue*> all_residues() {
        std::vector<IResidue*> residues;
        for (auto& chain : chains_) {
            for (size_t i = 0; i < chain.size(); ++i) {
                residues.push_back(&chain[i]);
            }
        }
        return residues;
    }

    [[nodiscard]] std::vector<const IResidue*> all_residues() const {
        std::vector<const IResidue*> residues;
        for (const auto& chain : chains_) {
            for (size_t i = 0; i < chain.size(); ++i) {
                residues.push_back(&chain[i]);
            }
        }
        return residues;
    }

    // === Nucleotide access ===
    [[nodiscard]] std::vector<INucleotide*> nucleotides() {
        std::vector<INucleotide*> nts;
        for (auto& chain : chains_) {
            for (auto* nuc : chain.nucleotides()) {
                nts.push_back(nuc);
            }
        }
        return nts;
    }

    [[nodiscard]] std::vector<const INucleotide*> nucleotides() const {
        std::vector<const INucleotide*> nts;
        for (const auto& chain : chains_) {
            for (auto* nuc : chain.nucleotides()) {
                nts.push_back(nuc);
            }
        }
        return nts;
    }

    // === Find chain ===
    [[nodiscard]] Chain* find_chain(const std::string& chain_id) {
        for (auto& chain : chains_) {
            if (chain.chain_id() == chain_id) {
                return &chain;
            }
        }
        return nullptr;
    }

    [[nodiscard]] const Chain* find_chain(const std::string& chain_id) const {
        for (const auto& chain : chains_) {
            if (chain.chain_id() == chain_id) {
                return &chain;
            }
        }
        return nullptr;
    }

    // === Legacy index support ===
    void set_legacy_indices(
        const std::map<std::tuple<std::string, int, std::string, std::string>, int>& atom_idx_map,
        const std::map<std::tuple<std::string, int, std::string>, int>& residue_idx_map) {

        for (auto& chain : chains_) {
            for (size_t i = 0; i < chain.size(); ++i) {
                auto& residue = chain[i];
                const std::string& chain_id = residue.chain_id();
                int residue_seq = residue.seq_num();
                const std::string& insertion = residue.insertion();

                // Set legacy residue index
                auto residue_key = std::make_tuple(chain_id, residue_seq, insertion);
                auto residue_it = residue_idx_map.find(residue_key);
                if (residue_it != residue_idx_map.end() && residue_it->second > 0) {
                    residue.set_legacy_residue_idx(residue_it->second);
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

    // === Get residue by legacy index ===
    [[nodiscard]] IResidue* get_residue_by_legacy_idx(int legacy_idx) {
        for (auto& chain : chains_) {
            for (size_t i = 0; i < chain.size(); ++i) {
                if (chain[i].legacy_residue_idx() == legacy_idx) {
                    return &chain[i];
                }
            }
        }
        return nullptr;
    }

    [[nodiscard]] const IResidue* get_residue_by_legacy_idx(int legacy_idx) const {
        for (const auto& chain : chains_) {
            for (size_t i = 0; i < chain.size(); ++i) {
                if (chain[i].legacy_residue_idx() == legacy_idx) {
                    return &chain[i];
                }
            }
        }
        return nullptr;
    }

    // === Record type support ===
    void set_residue_record_type(const std::string& chain_id, int seq_num,
                                  const std::string& insertion, char record_type) {
        residue_record_types_[std::make_tuple(chain_id, seq_num, insertion)] = record_type;
    }

    [[nodiscard]] char get_residue_record_type(const std::string& chain_id, int seq_num,
                                                const std::string& insertion) const {
        auto key = std::make_tuple(chain_id, seq_num, insertion);
        auto it = residue_record_types_.find(key);
        return (it != residue_record_types_.end()) ? it->second : 'A';
    }

private:
    std::string pdb_id_;
    std::vector<Chain> chains_;
    std::map<std::tuple<std::string, int, std::string>, char> residue_record_types_;
};

} // namespace poly
} // namespace core
} // namespace x3dna
