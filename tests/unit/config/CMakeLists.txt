# Config unit tests

add_executable(test_config_manager
    test_config_manager.cpp
)

target_link_libraries(test_config_manager
    x3dna
    gtest_main
)

gtest_discover_tests(test_config_manager)

add_executable(test_resource_locator
    test_resource_locator.cpp
)

target_link_libraries(test_resource_locator
    x3dna
    gtest_main
)

gtest_discover_tests(test_resource_locator)

