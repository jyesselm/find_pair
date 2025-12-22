/**
 * @file protocol_base.hpp
 * @brief Base class for all protocols
 */

#pragma once

#include <x3dna/core/structure.hpp>
#include <x3dna/config/config_manager.hpp>
#include <memory>

namespace x3dna {
namespace protocols {

/**
 * @class ProtocolBase
 * @brief Abstract base class for all protocols
 */
class ProtocolBase {
public:
    virtual ~ProtocolBase() = default;

    /**
     * @brief Execute the protocol on a structure
     * @param structure The structure to process
     */
    virtual void execute(core::Structure& structure) = 0;

    /**
     * @brief Set the configuration manager
     * @param config Reference to ConfigManager
     */
    void set_config_manager(config::ConfigManager& config) {
        config_ = &config;
    }

    /**
     * @brief Get the configuration manager
     * @return Reference to ConfigManager
     */
    config::ConfigManager& config() const {
        if (config_) {
            return *config_;
        }
        // Return singleton instance if not set
        return config::ConfigManager::instance();
    }

protected:
    /**
     * @brief Constructor
     */
    ProtocolBase() = default;

    /**
     * @brief Copy constructor (deleted)
     */
    ProtocolBase(const ProtocolBase&) = default;

    /**
     * @brief Move constructor (deleted)
     */
    ProtocolBase(ProtocolBase&&) = default;

    /**
     * @brief Copy assignment (deleted)
     */
    ProtocolBase& operator=(const ProtocolBase&) = default;

    /**
     * @brief Move assignment (deleted)
     */
    ProtocolBase& operator=(ProtocolBase&&) = default;

    /**
     * @brief Configuration manager (may be null, falls back to singleton)
     */
    config::ConfigManager* config_ = nullptr;
};

} // namespace protocols
} // namespace x3dna
