cmake_minimum_required(VERSION 3.10)
project(ParallelLinearProgrammingSolver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(MPI REQUIRED)

add_subdirectory(src)
add_subdirectory(tests)
