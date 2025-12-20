/**
 * @file string_utils.hpp
 * @brief Centralized string utility functions for consistent string handling
 */

#pragma once

#include <string>
#include <algorithm>
#include <cctype>

namespace x3dna {
namespace core {

/**
 * @brief String utility functions for consistent trimming and comparison
 */
class StringUtils {
public:
    /**
     * @brief Trim whitespace from both ends of a string
     * @param str String to trim
     * @return Trimmed string
     */
    [[nodiscard]] static std::string trim(const std::string& str) {
        if (str.empty()) {
            return str;
        }

        size_t start = 0;
        size_t end = str.length();

        while (start < end && std::isspace(static_cast<unsigned char>(str[start]))) {
            ++start;
        }

        while (end > start && std::isspace(static_cast<unsigned char>(str[end - 1]))) {
            --end;
        }

        return str.substr(start, end - start);
    }

    /**
     * @brief Trim whitespace from the left side of a string
     * @param str String to trim
     * @return Left-trimmed string
     */
    [[nodiscard]] static std::string trim_left(const std::string& str) {
        if (str.empty()) {
            return str;
        }

        size_t start = 0;
        while (start < str.length() && std::isspace(static_cast<unsigned char>(str[start]))) {
            ++start;
        }

        return str.substr(start);
    }

    /**
     * @brief Trim whitespace from the right side of a string
     * @param str String to trim
     * @return Right-trimmed string
     */
    [[nodiscard]] static std::string trim_right(const std::string& str) {
        if (str.empty()) {
            return str;
        }

        size_t end = str.length();
        while (end > 0 && std::isspace(static_cast<unsigned char>(str[end - 1]))) {
            --end;
        }

        return str.substr(0, end);
    }

    /**
     * @brief Check if two strings are equal after trimming
     * @param a First string
     * @param b Second string
     * @return true if trimmed strings are equal
     */
    [[nodiscard]] static bool equals_trimmed(const std::string& a, const std::string& b) {
        return trim(a) == trim(b);
    }

    /**
     * @brief Convert string to uppercase
     * @param str String to convert
     * @return Uppercase string
     */
    [[nodiscard]] static std::string to_upper(const std::string& str) {
        std::string result = str;
        std::transform(result.begin(), result.end(), result.begin(),
                       [](unsigned char c) { return std::toupper(c); });
        return result;
    }

    /**
     * @brief Convert string to lowercase
     * @param str String to convert
     * @return Lowercase string
     */
    [[nodiscard]] static std::string to_lower(const std::string& str) {
        std::string result = str;
        std::transform(result.begin(), result.end(), result.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        return result;
    }
};

// Convenience free functions for common operations
[[nodiscard]] inline std::string trim(const std::string& str) {
    return StringUtils::trim(str);
}

[[nodiscard]] inline std::string trim_left(const std::string& str) {
    return StringUtils::trim_left(str);
}

[[nodiscard]] inline std::string trim_right(const std::string& str) {
    return StringUtils::trim_right(str);
}

} // namespace core
} // namespace x3dna
