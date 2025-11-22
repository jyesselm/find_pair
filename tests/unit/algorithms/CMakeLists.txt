# Algorithms unit tests

add_executable(test_standard_base_templates
    test_standard_base_templates.cpp
)

target_link_libraries(test_standard_base_templates
    PRIVATE
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_standard_base_templates)

add_executable(test_ring_atom_matcher
    test_ring_atom_matcher.cpp
)

target_link_libraries(test_ring_atom_matcher
    PRIVATE
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_ring_atom_matcher)

add_executable(test_base_frame_calculator
    test_base_frame_calculator.cpp
)

target_link_libraries(test_base_frame_calculator
    PRIVATE
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_base_frame_calculator)

