/**
 * @file test_standard_base_templates.cpp
 * @brief Unit tests for StandardBaseTemplates class
 */

#include <gtest/gtest.h>
#include <x3dna/algorithms/standard_base_templates.hpp>
#include <x3dna/core/structure.hpp>
#include <filesystem>

using namespace x3dna::algorithms;
using namespace x3dna::core;

class StandardBaseTemplatesTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Skip tests if templates directory doesn't exist
        if (!std::filesystem::exists("data/templates")) {
            GTEST_SKIP() << "Templates directory not found: data/templates";
        }
        templates_ = std::make_unique<StandardBaseTemplates>("data/templates");
    }

    std::unique_ptr<StandardBaseTemplates> templates_;
};

// Constructor tests
TEST_F(StandardBaseTemplatesTest, ConstructorWithPath) {
    StandardBaseTemplates t("data/templates");
    EXPECT_EQ(t.template_path(), std::filesystem::path("data/templates"));
}

// Template path tests
TEST_F(StandardBaseTemplatesTest, GetTemplatePath) {
    std::filesystem::path path_a = templates_->get_template_path(ResidueType::ADENINE);
    EXPECT_EQ(path_a.filename(), "Atomic_A.pdb");

    std::filesystem::path path_c = templates_->get_template_path(ResidueType::CYTOSINE);
    EXPECT_EQ(path_c.filename(), "Atomic_C.pdb");
}

// Template existence tests
TEST_F(StandardBaseTemplatesTest, TemplateExists) {
    // May or may not exist depending on whether files were copied
    // Just verify the method works without throwing
    EXPECT_NO_THROW({
        bool exists_a = templates_->template_exists(ResidueType::ADENINE);
        bool exists_c = templates_->template_exists(ResidueType::CYTOSINE);
        // Variables are used implicitly in the EXPECT_NO_THROW check
        (void)exists_a;
        (void)exists_c;
    });
}

// Template loading tests
TEST_F(StandardBaseTemplatesTest, LoadTemplate) {
    // Only test if template exists
    if (templates_->template_exists(ResidueType::ADENINE)) {
        Structure template_structure = templates_->load_template(ResidueType::ADENINE);
        EXPECT_GT(template_structure.num_atoms(), 0);

        // Load again - should use cache
        Structure template_structure2 = templates_->load_template(ResidueType::ADENINE);
        EXPECT_EQ(template_structure.num_atoms(), template_structure2.num_atoms());
    }
}

// Cache tests
TEST_F(StandardBaseTemplatesTest, ClearCache) {
    if (templates_->template_exists(ResidueType::ADENINE)) {
        (void)templates_->load_template(ResidueType::ADENINE);
        templates_->clear_cache();
        // Should still be able to load after clearing
        EXPECT_NO_THROW({
            Structure s = templates_->load_template(ResidueType::ADENINE);
            EXPECT_GT(s.num_atoms(), 0);
        });
    }
}

// Error handling tests
TEST_F(StandardBaseTemplatesTest, InvalidResidueType) {
    EXPECT_THROW({ (void)templates_->load_template(ResidueType::AMINO_ACID); }, std::invalid_argument);
}
