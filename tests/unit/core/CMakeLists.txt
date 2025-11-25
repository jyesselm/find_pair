# Core object unit tests

add_executable(test_atom
    test_atom.cpp
)

target_link_libraries(test_atom
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_atom)

add_executable(test_reference_frame
    test_reference_frame.cpp
)

target_link_libraries(test_reference_frame
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_reference_frame)

add_executable(test_residue
    test_residue.cpp
)

target_link_libraries(test_residue
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_residue)

add_executable(test_chain
    test_chain.cpp
)

target_link_libraries(test_chain
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_chain)

add_executable(test_structure
    test_structure.cpp
)

target_link_libraries(test_structure
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_structure)

add_executable(test_base_pair
    test_base_pair.cpp
)

target_link_libraries(test_base_pair
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_base_pair)

add_executable(test_parameters
    test_parameters.cpp
)

target_link_libraries(test_parameters
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_parameters)

add_executable(test_structure_legacy_order
    test_structure_legacy_order.cpp
)

target_link_libraries(test_structure_legacy_order
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_structure_legacy_order)

