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
    gtest_main
)

gtest_discover_tests(test_pdb_parser)

add_executable(test_json_writer
    test_json_writer.cpp
)

target_link_libraries(test_json_writer
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_json_writer)

add_executable(test_json_reader
    test_json_reader.cpp
)

target_link_libraries(test_json_reader
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_json_reader)

add_executable(test_pdb_writer
    test_pdb_writer.cpp
)

target_link_libraries(test_pdb_writer
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_pdb_writer)

add_executable(test_input_file_parser
    test_input_file_parser.cpp
)

target_link_libraries(test_input_file_parser
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_input_file_parser)

add_executable(test_cif_parser
    test_cif_parser.cpp
)

target_link_libraries(test_cif_parser
    PRIVATE
    x3dna
    gtest_main
)

gtest_discover_tests(test_cif_parser)

