# Protocol unit tests

add_executable(test_protocol_base
    test_protocol_base.cpp
)

target_link_libraries(test_protocol_base
    x3dna
    gtest_main
)

gtest_discover_tests(test_protocol_base)

add_executable(test_find_pair_protocol
    test_find_pair_protocol.cpp
)

target_link_libraries(test_find_pair_protocol
    x3dna
    gtest_main
)

gtest_discover_tests(test_find_pair_protocol)

