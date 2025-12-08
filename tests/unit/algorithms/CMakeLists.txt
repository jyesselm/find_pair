# Algorithms unit tests

add_executable(test_standard_base_templates
    test_standard_base_templates.cpp
)

target_link_libraries(test_standard_base_templates
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_standard_base_templates)

add_executable(test_ring_atom_matcher
    test_ring_atom_matcher.cpp
)

target_link_libraries(test_ring_atom_matcher
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_ring_atom_matcher)

add_executable(test_base_frame_calculator
    test_base_frame_calculator.cpp
)

target_link_libraries(test_base_frame_calculator
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_base_frame_calculator)

# Base pair validator tests
add_executable(test_base_pair_validator
    test_base_pair_validator.cpp
)

target_include_directories(test_base_pair_validator PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_base_pair_validator
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_base_pair_validator)

# Base pair finder tests
add_executable(test_base_pair_finder
    test_base_pair_finder.cpp
)

target_include_directories(test_base_pair_finder PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_base_pair_finder
    PRIVATE
    x3dna
    gtest_main
)

# Link filesystem library if needed
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_base_pair_validator PRIVATE stdc++fs)
    target_link_libraries(test_base_pair_finder PRIVATE stdc++fs)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(test_base_pair_validator PRIVATE c++fs)
    target_link_libraries(test_base_pair_finder PRIVATE c++fs)
endif()

gtest_discover_tests(test_base_pair_finder)

# Residue index map tests
add_executable(test_residue_index_map
    test_residue_index_map.cpp
)

target_link_libraries(test_residue_index_map
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_residue_index_map)
