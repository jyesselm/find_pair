# I/O unit tests

add_executable(test_pdb_parser
    test_pdb_parser.cpp
)

target_include_directories(test_pdb_parser PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(test_pdb_parser
    PRIVATE
    x3dna
    gtest
    gtest_main
)

gtest_discover_tests(test_pdb_parser)

