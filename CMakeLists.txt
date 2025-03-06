cmake_minimum_required(VERSION 3.10)
project(FiniteElementMesh)

set(CMAKE_CXX_STANDARD 17)

# Enable testing
enable_testing()

# Add Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# Prevent GoogleTest from overriding our compiler/linker options
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Add main source files
add_library(mesh2d src/mesh2d.cpp)
add_executable(main src/main.cpp)
target_link_libraries(main mesh2d)

# Add tests
add_executable(tests tests/test_mesh2d.cpp)
target_link_libraries(tests mesh2d gtest gtest_main pthread)

# Register tests
add_test(NAME MeshTests COMMAND tests)