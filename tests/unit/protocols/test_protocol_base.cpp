/**
 * @file test_protocol_base.cpp
 * @brief Unit tests for ProtocolBase
 */

#include <gtest/gtest.h>
#include <x3dna/protocols/protocol_base.hpp>
#include <x3dna/config/config_manager.hpp>
#include <x3dna/core/structure.hpp>

using namespace x3dna::protocols;
using namespace x3dna::config;
using namespace x3dna::core;

// Mock protocol implementation for testing
class MockProtocol : public ProtocolBase {
public:
    bool executed_ = false;
    Structure* last_structure_ = nullptr;
    
    void execute(Structure& structure) override {
        executed_ = true;
        last_structure_ = &structure;
    }
    
    bool has_config() const {
        return config_ != nullptr;
    }
};

class ProtocolBaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        protocol_ = std::make_unique<MockProtocol>();
        config_ = &ConfigManager::instance();
        config_->set_defaults();
    }
    
    std::unique_ptr<MockProtocol> protocol_;
    ConfigManager* config_;
    Structure structure_;
};

// Basic execution test
TEST_F(ProtocolBaseTest, Execute) {
    EXPECT_FALSE(protocol_->executed_);
    
    protocol_->execute(structure_);
    
    EXPECT_TRUE(protocol_->executed_);
    EXPECT_EQ(protocol_->last_structure_, &structure_);
}

// Configuration management tests
TEST_F(ProtocolBaseTest, SetConfigManager) {
    EXPECT_FALSE(protocol_->has_config());
    
    protocol_->set_config_manager(*config_);
    
    EXPECT_TRUE(protocol_->has_config());
    EXPECT_EQ(&protocol_->config(), config_);
}

TEST_F(ProtocolBaseTest, GetConfigWithoutSet) {
    // Should return singleton instance if not set
    auto& default_config = protocol_->config();
    EXPECT_EQ(&default_config, &ConfigManager::instance());
}

TEST_F(ProtocolBaseTest, GetConfigWithSet) {
    protocol_->set_config_manager(*config_);
    
    auto& retrieved_config = protocol_->config();
    EXPECT_EQ(&retrieved_config, config_);
}

TEST_F(ProtocolBaseTest, ConfigModification) {
    protocol_->set_config_manager(*config_);
    
    // Modify config through protocol
    protocol_->config().set_include_hetatm(true);
    protocol_->config().thresholds().max_dorg = 20.0;
    
    // Verify changes
    EXPECT_TRUE(config_->include_hetatm());
    EXPECT_DOUBLE_EQ(config_->thresholds().max_dorg, 20.0);
}

// Virtual destructor test (compile-time check)
TEST_F(ProtocolBaseTest, VirtualDestructor) {
    // This test verifies that ProtocolBase can be deleted through base pointer
    ProtocolBase* base_ptr = protocol_.release();
    EXPECT_NO_THROW(delete base_ptr);
}

// Multiple protocols with same config
TEST_F(ProtocolBaseTest, MultipleProtocolsSameConfig) {
    MockProtocol protocol1;
    MockProtocol protocol2;
    
    protocol1.set_config_manager(*config_);
    protocol2.set_config_manager(*config_);
    
    // Both should reference same config
    EXPECT_EQ(&protocol1.config(), &protocol2.config());
    EXPECT_EQ(&protocol1.config(), config_);
}

// Protocol without config (uses singleton)
TEST_F(ProtocolBaseTest, ProtocolWithoutConfig) {
    MockProtocol protocol;
    
    // Should use singleton
    auto& config = protocol.config();
    EXPECT_EQ(&config, &ConfigManager::instance());
    
    // Modifications should affect singleton
    config.set_legacy_mode(true);
    EXPECT_TRUE(ConfigManager::instance().legacy_mode());
}

