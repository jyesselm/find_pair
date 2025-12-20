/**
 * @file chain.hpp
 * @brief Chain class representing a chain of residues in a PDB structure
 */

#pragma once

#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <x3dna/core/residue.hpp>

namespace x3dna {
namespace core {

/**
 * @class Chain
 * @brief Represents a chain of residues (typically a single polymer chain)
 */
class Chain {
public:
    /**
     * @brief Default constructor
     */
    Chain() = default;

    /**
     * @brief Constructor with chain ID
     * @param id Chain identifier (e.g., "A", "B", "AA" for CIF files)
     */
    explicit Chain(const std::string& id) : chain_id_(id) {}

    // Getters
    [[nodiscard]] const std::string& chain_id() const {
        return chain_id_;
    }
    [[nodiscard]] const std::vector<Residue>& residues() const {
        return residues_;
    }
    std::vector<Residue>& residues() {
        return residues_;
    }
    [[nodiscard]] size_t num_residues() const {
        return residues_.size();
    }

    // Container-like interface
    [[nodiscard]] auto begin() const {
        return residues_.begin();
    }
    [[nodiscard]] auto end() const {
        return residues_.end();
    }
    [[nodiscard]] auto begin() {
        return residues_.begin();
    }
    [[nodiscard]] auto end() {
        return residues_.end();
    }
    [[nodiscard]] size_t size() const {
        return residues_.size();
    }
    [[nodiscard]] bool empty() const {
        return residues_.empty();
    }
    [[nodiscard]] const Residue& operator[](size_t idx) const {
        return residues_[idx];
    }
    [[nodiscard]] Residue& operator[](size_t idx) {
        return residues_[idx];
    }

    /**
     * @brief Get total number of atoms in this chain
     */
    [[nodiscard]] size_t num_atoms() const {
        size_t count = 0;
        for (const auto& residue : residues_) {
            count += residue.num_atoms();
        }
        return count;
    }

    // Setters
    void set_chain_id(const std::string& id) {
        chain_id_ = id;
    }

    /**
     * @brief Add a residue to this chain
     */
    void add_residue(const Residue& residue) {
        residues_.push_back(residue);
    }

    /**
     * @brief Get sequence as one-letter code string
     * @return Sequence string (e.g., "ACGT")
     */
    [[nodiscard]] std::string sequence() const {
        std::string seq;
        for (const auto& residue : residues_) {
            char code = residue.one_letter_code();
            if (code != '?') {
                seq += code;
            }
        }
        return seq;
    }

    /**
     * @brief Get all nucleotides in this chain
     * @return Vector of nucleotide residues
     */
    [[nodiscard]] std::vector<Residue> nucleotides() const {
        std::vector<Residue> nts;
        for (const auto& residue : residues_) {
            if (residue.is_nucleotide()) {
                nts.push_back(residue);
            }
        }
        return nts;
    }

    /**
     * @brief Get residue by sequence number
     * @param seq_num Sequence number
     * @return Optional residue if found
     */
    [[nodiscard]] std::optional<Residue> find_residue(int seq_num) const {
        for (const auto& residue : residues_) {
            if (residue.seq_num() == seq_num) {
                return residue;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Convert to legacy JSON format
     */
    [[nodiscard]] nlohmann::json to_json_legacy() const {
        nlohmann::json j;
        j["chain_id"] = chain_id_;
        j["num_residues"] = static_cast<long>(residues_.size());
        j["residues"] = nlohmann::json::array();
        for (const auto& residue : residues_) {
            j["residues"].push_back(residue.to_json_legacy());
        }
        return j;
    }

    /**
     * @brief Create Chain from legacy JSON format
     */
    [[nodiscard]] static Chain from_json_legacy(const nlohmann::json& j) {
        std::string chain_id = j.value("chain_id", "");

        Chain chain(chain_id);

        if (j.contains("residues") && j["residues"].is_array()) {
            for (const auto& residue_json : j["residues"]) {
                chain.add_residue(Residue::from_json_legacy(residue_json));
            }
        }

        return chain;
    }

    /**
     * @brief Convert to modern JSON format
     */
    [[nodiscard]] nlohmann::json to_json() const {
        nlohmann::json j;
        j["chain_id"] = chain_id_;
        j["residues"] = nlohmann::json::array();
        for (const auto& residue : residues_) {
            j["residues"].push_back(residue.to_json());
        }
        return j;
    }

    /**
     * @brief Create Chain from modern JSON format
     */
    [[nodiscard]] static Chain from_json(const nlohmann::json& j) {
        std::string chain_id = j.value("chain_id", "");

        Chain chain(chain_id);

        if (j.contains("residues") && j["residues"].is_array()) {
            for (const auto& residue_json : j["residues"]) {
                chain.add_residue(Residue::from_json(residue_json));
            }
        }

        return chain;
    }

private:
    std::string chain_id_;          // Chain identifier (string for CIF compatibility)
    std::vector<Residue> residues_; // Residues in this chain
};

} // namespace core
} // namespace x3dna
