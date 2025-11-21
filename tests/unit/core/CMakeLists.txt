# Core object unit tests

add_executable(test_atom
    test_atom.cpp
)

target_link_libraries(test_atom
    PRIVATE
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_atom)

add_executable(test_reference_frame
    test_reference_frame.cpp
)

target_link_libraries(test_reference_frame
    PRIVATE
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_reference_frame)

