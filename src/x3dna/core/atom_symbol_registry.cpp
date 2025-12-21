#include <x3dna/core/atom_symbol_registry.hpp>
#include <x3dna/config/resource_locator.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <mutex>

using json = nlohmann::json;

namespace x3dna {
namespace core {

// Symbol to index mapping (matches legacy asym_idx values)
const std::map<std::string, int> AtomSymbolRegistry::SYMBOL_TO_IDX = {
    {"C", 1}, {"O", 2}, {"H", 3}, {"N", 4}, {"S", 5}, {"P", 6}
};

// Thread-safe lazy-loaded pattern registry
static std::map<std::string, std::string>& get_pattern_registry() {
    static std::map<std::string, std::string> patterns;
    static std::once_flag init_flag;

    std::call_once(init_flag, []() {
        std::filesystem::path config_file;
        try {
            config_file = config::ResourceLocator::config_file("atomlist.json");
        } catch (const std::runtime_error&) {
            // ResourceLocator not initialized or file not found - use empty registry
            // Fallback logic in get_symbol() will handle unknown atoms
            return;
        }

        std::ifstream file(config_file);
        if (!file.is_open()) {
            // File not found - use fallback logic only
            return;
        }

        try {
            json j = json::parse(file);

            for (const auto& [pattern, symbol] : j["patterns"].items()) {
                patterns[pattern] = symbol.get<std::string>();
            }
        } catch (const json::exception& e) {
            throw std::runtime_error("AtomSymbolRegistry: Error parsing atomlist.json: " +
                                     std::string(e.what()));
        }
    });

    return patterns;
}

const std::map<std::string, std::string>& AtomSymbolRegistry::get_patterns() {
    return get_pattern_registry();
}

std::string AtomSymbolRegistry::pad_atom_name(const std::string& atom_name) {
    std::string name = atom_name;
    if (name.length() < 4) {
        // Most nucleotide atoms are single-letter elements with format " XNN"
        // Prepend space for single-letter element names (length 1-3)
        if (name.length() <= 3 && !name.empty() && std::isupper(static_cast<unsigned char>(name[0]))) {
            name = " " + name;
        }
        while (name.length() < 4) {
            name += " ";
        }
    }
    return name.substr(0, 4);
}

std::string AtomSymbolRegistry::atom_name_to_pattern(const std::string& atom_name) {
    std::string padded = pad_atom_name(atom_name);
    std::string pattern = padded;

    // Convert to uppercase and replace non-alpha with '.'
    for (size_t i = 0; i < 4 && i < pattern.length(); i++) {
        char c = pattern[i];
        if (std::isalpha(static_cast<unsigned char>(c))) {
            pattern[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        } else {
            pattern[i] = '.';
        }
    }

    return pattern;
}

std::string AtomSymbolRegistry::get_symbol(const std::string& atom_name) {
    const auto& patterns = get_patterns();
    std::string pattern = atom_name_to_pattern(atom_name);

    // Try exact pattern match first
    auto it = patterns.find(pattern);
    if (it != patterns.end()) {
        return it->second;
    }

    // Fallback logic (matches legacy aname2asym behavior)
    bool unknown = (pattern == ".UNK");

    // Fallback case 1: Two-letter element at start (e.g., "FE..")
    if (pattern.length() >= 4 && pattern[0] != '.' && pattern[1] != '.' &&
        pattern[2] == '.' && pattern[3] == '.') {
        return std::string(1, pattern[1]);
    }

    // Fallback case 2: Single-letter element after dot (e.g., ".N..")
    if (pattern.length() >= 2 && pattern[0] == '.' && pattern[1] != '.' && !unknown) {
        return std::string(1, pattern[1]);
    }

    // Fallback case 3: Starts with H
    if (pattern.length() >= 1 && pattern[0] == 'H') {
        return "H";
    }

    // Unknown atom
    return "XX";
}

int AtomSymbolRegistry::get_atom_idx(const std::string& atom_name) {
    std::string symbol = get_symbol(atom_name);

    auto it = SYMBOL_TO_IDX.find(symbol);
    if (it != SYMBOL_TO_IDX.end()) {
        return it->second;
    }
    return 0;
}

bool AtomSymbolRegistry::contains_pattern(const std::string& pattern) {
    const auto& patterns = get_patterns();
    return patterns.find(pattern) != patterns.end();
}

} // namespace core
} // namespace x3dna
