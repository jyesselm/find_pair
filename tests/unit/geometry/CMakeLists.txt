# Geometry unit tests

add_executable(test_vector3d
    test_vector3d.cpp
)

target_link_libraries(test_vector3d
    x3dna
    gtest_main
)

gtest_discover_tests(test_vector3d)

add_executable(test_matrix3d
    test_matrix3d.cpp
)

target_link_libraries(test_matrix3d
    x3dna
    gtest_main
)

gtest_discover_tests(test_matrix3d)

add_executable(test_least_squares_fitter
    test_least_squares_fitter.cpp
)

target_link_libraries(test_least_squares_fitter
    x3dna
    gtest_main
)

gtest_discover_tests(test_least_squares_fitter)

