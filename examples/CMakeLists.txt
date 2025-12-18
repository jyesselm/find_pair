# Examples CMakeLists.txt

# Simple analysis example
add_executable(simple_analysis simple_analysis.cpp)
target_link_libraries(simple_analysis PRIVATE x3dna)

# You can add more examples here:
# add_executable(example_name source.cpp)
# target_link_libraries(example_name PRIVATE x3dna)
