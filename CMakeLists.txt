cmake_minimum_required(VERSION 3.10)
project(FiniteElementMesh)

set(CMAKE_CXX_STANDARD 17)

add_executable(FiniteElementMesh src/main.cpp src/mesh2d.cpp)
